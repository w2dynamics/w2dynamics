#!/usr/bin/env python
"""Interface for the MaxEnt Solver"""
from __future__ import absolute_import, print_function, unicode_literals
import numpy as np
import h5py as hdf5
import sys
import optparse
import copy

from w2dyn.maxent import Solver
from w2dyn.auxiliaries import quantities as qttys

__version__ = (1,0)

def generateData(filename):   
   """
   Helper function to parse file data
   """
   with open(filename) as fdata:
      for line in fdata:
         data = np.fromstring(line, dtype=float, sep=" ")
         yield data

def readDat2Solver(fname,solver,offset):
   """
   Helper function to import data to solver object,
   the specification of solver.kernelMode determines
   if gtau/chi (3 columns) or giw (4 columns)
   is processed.
        
   the input parameters follow as:
   fname: input file name with data
   solver: maxent solver object with prefedined kernelMode,beta and Ncorr
   offset: number of comment lines to be skipped before data begins
   """
   
   line_iter = iter(generateData(fname))
   #exhausting generator to offset
   for i in range(offset):
      next(line_iter)
      
   for i, data in enumerate(line_iter, start=0):
      if(solver.kernelMode != 1):
         solver.corrGrid[i] = data[0]
         solver.corrReal[i] = data[1]
         solver.corrCov[i]  = data[2]
      else:
         solver.corrGrid[i] = data[0]
         solver.corrReal[i] = data[1]
         solver.corrImag[i] = data[2]
         solver.corrCov[i]  = data[3]

      solver.corr[i] = solver.corrReal[i] + solver.corrImag[i]*1j

      #initialize inverse diagonal covariance matrix
      solver.corrCov[i]=1./solver.corrCov[i]**2

      #define upper-limit for error, lower limit for covariance-matrix elements
      if (abs(solver.corrCov[i]) < 1e-5):
         solver.corrCov[i]=1e-5

def readHdf2Dat(tag,fin,outprefix,efunc,paramag=False,niter=-1,iiw=-1):
   """  Function to convert w2dynamics hdf5-file into data files
        accepted by the maxent program
        the (historic) file format consists of two comment lines,
        which are ignored, followed by data lines
         
        the input parameters follow as:
         tag: gtau or giw, determines what quantity is to be read out
         fin: hdf5 input file name
         outprefix: prefix for output files
         efunc: adjust the error manually (expects python expression)
         paramag: paramagnetic symmetrization
         niter: w2dyn iteration, if -1 then "dmft-last" is selected
         iiw: limit the number of iw frequencies on giw, symmetrically
              around 0, if iww=-1 all frequencies are used
   """
   hf = hdf5.File(fin,"r")
   
   file_version = tuple(hf.attrs["outfile-version"])
   beta = qttys.MetaQttyContainer("config", hf, file_version) \
            .select(qttys.SelectorPool(qttys.Selector("*beta"))) \
            .values()
   
   #we refer to tau or iw array as ximag in the following
   if tag == "gtau":
      ximag = hf[("axes", ".axes")[file_version[0]-1]]["taubin"][()]
   
   if tag == "giw-meas":
      ximag = hf[("axes", ".axes")[file_version[0]-1]]["iw"][()]
   
   #selecting the iteration
   if(niter==-1):
      diter = "dmft-last"
      try: 
         hf[diter]
      except KeyError:
         diter = "stat-last"
   else:
      diter = "dmft-" + "%03d" % niter
      try:
         hf[diter]
      except KeyError:
         diter = "stat-" + "%03d" % niter
   
   outlist = []
   
   for iineq, (all_gtaumean, all_gtauerr) in \
      enumerate(qttys.ineq_quantity(hf[diter], tag)):
         
         if paramag:
            all_gtaumean = all_gtaumean.mean(1)[:,None]
            all_gtauerr = np.sqrt((all_gtauerr**2).sum(1)[:,None])
            
         for iband, ispin in np.ndindex(all_gtaumean.shape[:2]):
            suffix = diter + "_" + str(iineq) + "_" + str(iband) + "_" + str(ispin) + ".dat"
            
            outlist.append(suffix)
            
            fdat = open(outprefix+suffix, "w")
            print(" ", file=fdat)
            print("# beta = %s, iter = %s, ineq = %s, band = %s, spin = %s" %
                  (beta[0], niter, iineq, iband, ispin), file=fdat)
            
            gtaumean = all_gtaumean[iband,ispin]

            #cutting frequency box symmetrically using iiw in case the number of freuqencies
            #in giw-meas is too large
            if tag == "giw-meas" and iiw != None:
               gtaumean = gtaumean[gtaumean.shape[0]//2-iiw//2:gtaumean.shape[0]//2+iiw//2]
               ximag = ximag[ximag.shape[0]//2-iiw//2:ximag.shape[0]//2+iiw//2]
                             
            gtauerr = all_gtauerr[iband, ispin]
            gtauerr[...] = efunc(ximag, gtaumean, gtauerr)
            
            
            for i in range(gtaumean.shape[0]):
               if tag == "gtau":
                  print(ximag[i], gtaumean[i], gtauerr[i], file=fdat)
               if tag == "giw-meas":
                  print(ximag[i],
                        -gtaumean[i].real, -gtaumean[i].imag, gtauerr[i],
                        file=fdat)
   return outlist

if __name__ == "__main__":
   """
   the following provides a possible interface to the maxent solver in the Solver.py
   class, default parameter values are currently set to gtau -> a(w) continuations.
   In case a continuation chi(iw)->chi(w) is desired, internal default parameters of
   Solver.py are used.
   """
   
   parser = optparse.OptionParser(usage = "%prog [OPTIONS] <hdf file> ...",
                               version = ".".join(map(str,__version__)))
   parser.add_option("--wmin", type="float", metavar="W", dest="wmin", default=-30.,
                  help="beginning of frequency window to use (default -30)")
   parser.add_option("--wmax", type="float", metavar="W", dest="wmax", default=30.,
                  help="end of frequency window to use (default 30)")
   parser.add_option("--wcount", type="int", metavar="N", dest="wcount", default=601,
                  help="number of frequencies to use (default 601)."
                        "For Greens functions/self-energies use 6*n+1, where n is an integer."
                        "The real frequency grid is divided into 3 intervals, symmetric around zero."
                        "The frequency w=0 is included, therefore wcount has to be odd.")
   parser.add_option("--nsmooth", type="int", metavar="N", dest="nsmooth", default=3, 
                  help="number of smoothening steps (default 3)")
   parser.add_option("--niter", type="int", metavar="N", dest="niter", default=-1,
                  help="hdf5-file dmft iteration number"
                       "if dmft-iteration does not exist, assume stat-iteration (default=-1)")
   parser.add_option("--gonly", action="store_true", dest="gonly", default=False, 
                  help="only write the Greens function to a file and exit")
   parser.add_option("--gerror", type="string", metavar="f(x, v, e)", default="e",
                  help="adjust the error manually (expects python expression)")
   parser.add_option("--kernelmode", type="int", metavar="N", dest="kmode", default=0,
                     help="kernelmode: 0-> gtau->A(w); 1-> giw->A(w);\
                     10-> chi(w)->A(w)=chi/omega; 11 -> chi(w)->A(w)=chi/(1-exp(...)) (default 0)")
   parser.add_option("--iiw", type="int", metavar="N", dest="iiw", default=None,
                     help="number of matsubaras from hdf5 to use for kernelMode 1")
   parser.add_option("--filetype", type="string", metavar="string", dest="ftype", default="hdf5",
                     help='"hdf5" to parse binary file, "dat" to parse ascii file')
   parser.add_option("--offset", type="int", metavar="N", dest="offset", default=2,
                     help="number of comment lines in ascii file before data")
   parser.add_option("--beta", type="float", metavar="N", dest="beta", default=None,
                     help="inverse temperature when reading from dat file")
   parser.add_option("--mom0", type="float", metavar="N",dest="mom0", default=0,
                     help="Zeroth moment (Hartree Term) of Self Energy (requires kernelmode=1 and filetype=dat)")
   parser.add_option("--mom1", type="float", metavar="N",dest="mom1", default=1,
                     help="First moment of Self Energy (requires kernelmode=1 and filetype=dat)")
   parser.add_option("--hfilling", type="int", metavar="N", dest="hfilling", default=0,
                     help="hfilling: 1-> enforces particle-hole symmetry in Green's functions")
   (options, args) = parser.parse_args()
   
   if not args:
      parser.error("expecting one or more HDF5/dat filenames")

   #loop over all hdf5 input filenames
   for filename in args:
      print("Processing file:", filename)
      
      if options.ftype == "hdf5":
         try:
            hf = hdf5.File(filename,"r")
         except IOError as e:
            parser.error("invalid HDF5 output file: %s" % filename)
            
         file_version = tuple(hf.attrs["outfile-version"])
         beta, paramag, magnetism = qttys.MetaQttyContainer("config", hf, file_version) \
                  .select(qttys.SelectorPool(qttys.Selector("*beta", "ParaMag", "general.magnetism"))) \
                  .values()
         print("beta = %s, paramag = %s, magnetism = %s" % (beta, paramag, magnetism))
         
         paramag = paramag or magnetism == 'para'  
         efunc = eval("lambda x, v, e: (%s)" % options.gerror)
      
         
         #creating files from hdf5, determining what iteration and quantity we are interested in
         if options.kmode==0:
            tag="gtau"
         elif options.kmode==1:
            tag="giw-meas"
         else:
            sys.exit("other kernel modes not yet parsed from hdf5, only read in from dat possible")
         
         outprefix = tag + "_"
         
         outlist = readHdf2Dat(tag,filename,outprefix,efunc,paramag,niter=options.niter,iiw=options.iiw)
         print(outlist)
         
         #only print gtau/giw, no maxent
         if options.gonly: continue
      
      elif options.ftype == "dat":
         outprefix = ""
         outlist = []
         outlist.append(filename)
         
         if options.beta == None:
            parser.error("no inverse temperature supplied")
         else:
            beta = options.beta
            
      for suffix in outlist:
               
         #counting number of lines
         with open(outprefix+suffix) as fdata:
            lines = sum(1 for line in fdata)
         
         #specifing comment lines and such
         offset = options.offset
         
         #creating maxent solver object with default window, logfiles will be named similar to file output
         if options.kmode == 0 or options.kmode == 1:
            
            #we contstruct a default equispaced omega grid
            wrange = options.wmax - options.wmin
            w1 = options.wmin + wrange/3
            w2 = options.wmax - wrange/3
            if options.wcount % 6 != 1:
              print('Warning: wcount should be a 6*n+1, e.g. 7,13,19')
              wcount_old=options.wcount
              options.wcount=1+6*(options.wcount//6 + 1) 
              print('Changing wcount from', wcount_old, 'to', options.wcount)

            solver = Solver.Green2Aw(kernelMode=options.kmode,beta=beta,Ncorr=lines-offset,\
                                    NGrid1=options.wcount//3,NGrid2=options.wcount//3+1,NGrid3=options.wcount//3,\
                                    wmin=options.wmin,w1=w1,w2=w2,wmax=options.wmax,NSmooth=options.nsmooth,\
                                    logprefix=outprefix+suffix)
         elif options.kmode == 10 or options.kmode == 11:    
            #FIXME: currently we use the internal default values, since they differ from the parser default values
            #at some point it will be necessary to play around with these values a bit
            print("ignoring input grid/frequency parameters", file=sys.stderr)
            solver = Solver.Chi2Aw(kernelMode=options.kmode,beta=beta,Ncorr=lines-offset,logprefix=outprefix+suffix)
         else:
            raise NotImplementedError("specified kernelMode unknown")
         
         #reading data from fname to object
         readDat2Solver(outprefix+suffix,solver,offset)
         
         #transfroming self energy into greens function
         if options.kmode==1:
             solver.corr = solver.corr-options.mom0
             solver.corr = solver.corr/options.mom1
             solver.corrCov = solver.corrCov*(options.mom1**2)
         
         #particle hole symmetry
         if options.hfilling==1:
            if options.kmode==0:
               solver.symmetrize()
            elif options.kmode==1:
               solver.corr = 1j*solver.corr.imag + 0.
            else:
               sys.exit("Half-Filling option only implemented for kernel 0 and 1")
             
         #computing the w-grid and the flat model
         solver.computeGrid()
         solver.computeModel()
         
         #passing the object members to the fortran maxent solver
         solver.computeSpec()
         
         #transfroming greens function into self energy
         if options.kmode==1:
             solver.spec = solver.spec*options.mom1
             
         #printing data
         solver.printSpec("aw_"+suffix)
         
         #converting aw to G(w) or Sigma(w)
         solver.computeGws()
         solver.gws = solver.gws + options.mom0
         solver.printGws("gws_"+suffix)
