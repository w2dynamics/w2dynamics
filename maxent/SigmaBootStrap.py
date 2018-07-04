 #!/opt/python/bin python
"""Helper script to generate a bootstrap error for the quantities (siw,gtau,gtau_worm ...)
   of the w2dynamics hdf5 output file. Generates plain text files which can be processed
   by MaxEnt.py afterwards.
   This helper script requires the scikits bootstrap library written by Constantine Evans
   
   Note that the current script is still in an experimental stage and may require some manual tuning.

   Author: Patrik Gunacker
"""

import numpy as np
import h5py as hdf5
import sys
import optparse, copy

import auxiliaries.quantities as qttys
import scikits.bootstrap as boot

__version__ = (1,0)

parser = optparse.OptionParser(usage = "%prog [OPTIONS] <hdf file> ...",
                               version = ".".join(map(str,__version__)))
parser.add_option("--quantity", type="string", metavar="string", dest="tag", default="siw",
                     help='QMC Quantity which will be bootstrapped')

(options, args) = parser.parse_args()
   
if not args:
   parser.error("expecting one or more HDF5/dat filenames")

for filename in args:      
   try:
      hf = hdf5.File(filename,"r")
   except IOError, e:
      parser.error("invalid HDF5 output file: %s" % filename)
               
   file_version = tuple(hf.attrs["outfile-version"])
   StatSteps = qttys.MetaQttyContainer("config", hf, file_version) \
                              .select(qttys.SelectorPool(qttys.Selector("*statisticsteps"))) \
                              .values()
                           
   print >> sys.stderr, "Number of Samples to be taken into account (StatSteps) = %s" % StatSteps[0]
   
   if StatSteps[0] == 0:
        print >> sys.stderr, "No Bins(StatSteps) given in %s" % filename
        continue
   
   NIneqs = hf["start/nneq/value"].value
   
   #we refer to tau or iw array as ximag in the following
   if options.tag == "gtau" or options.tag == "gtau_worm":
      freqval = hf[("axes", ".axes")[file_version[0]-1]]["taubin"].value
   
   elif options.tag == "siw" or options.tag == "giw":
      freqval = hf[("axes", ".axes")[file_version[0]-1]]["iw"].value
   else:
      sys.exit("Please specify x-axis manually in script")
      
   #loop over all ineqs
   for iineq in xrange(0,NIneqs):
      
      #full sample container
      niter = "stat-" + "001"
      #ineq = "ineq-" + "%03d" % iineq
      
      #only acces one inneq at a time, we can directly evaluate the content of the generator
      (val,err) = qttys.ineq_quantity(hf[niter], options.tag,ineq=iineq).next()
      
      #make sure the err array exits by overwriting it
      err = np.zeros_like(val,dtype=float)
      
      freqs = val.shape[-1]
      samples=np.zeros(shape=(StatSteps[0],freqs),dtype=complex)
      
      #loop over band and spins
      for iband, ispin in np.ndindex(val.shape[:2]):
            samples[0,:]=val[iband,ispin,:]
            
            #prepare filename
            fname = options.tag + "_" + str(iineq) + "_" + str(iband) + "_" + str(ispin) + ".dat"
            print >> sys.stderr, "Writing to %s" % fname
            
            fdat = file(fname, "w")
            print >> fdat, " "
            print >> fdat, "#ineq = %s, band = %s, spin = %s" % (iineq, iband, ispin)
            
            #add all other StatStep samples to array
            #FIXME: we should read out a single band spin at a time
            for isample in xrange(2,StatSteps[0]+1):
               niter = "stat-" + "%03d" % isample
               (val_add,err_add) = qttys.ineq_quantity(hf[niter], options.tag,ineq=iineq).next()
               samples[isample-1,:]=val_add[iband,ispin,:]
                        
            #calculate the mean
            val[iband,ispin,:]=np.mean(samples,axis=0)
            err_std=np.std(samples,axis=0)
            
            #calculate the bootstrap error bar for each freq
            for ifreq in xrange(0,freqs):
               #alpha sets the confidence interval to 1 sigma
               #for now we can only use 'bca' for either real or imaginary part
               ci_i = boot.ci( samples[:,ifreq], statfunction = lambda x: np.average(np.abs(x)), alpha=(1.-0.6827), n_samples=1000, method='bca', output='errorbar' )
               
               #add performance counter here
               
               #bootstrap gives us a lower and upper errorbar
               #we assume them to be almost equal such that a simple average is justified
               if options.tag == "gtau" or options.tag == "gtau_worm":
                  print >> fdat, freqval[ifreq], val[iband,ispin,ifreq].real, np.average(ci_i)
               elif options.tag == "siw" or options.tag == "giw":
                  print >> fdat, freqval[ifreq], val[iband,ispin,ifreq].real, val[iband,ispin,ifreq].imag, np.average(ci_i)
               
