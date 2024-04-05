#!/usr/bin/env python3
""" taylorfit_selfenergy.py [--help] [--iteration=i] <file.hdf5> 

Computes the Taylor fit to the real-axis self energy and writes it
back to the HDF5 file as new quantity "sfreal".
"""
from __future__ import print_function
import sys
from sys import exit, stderr as err
import optparse, re
##pw
#import matplotlib.pyplot as plt
from scipy.optimize import fmin
from scipy.optimize import fmin_bfgs
from scipy.optimize import leastsq
#import matplotlib.pyplot as plt
from configobj import ConfigObj

#import scipy as sp

#from scipy.optimize import minimize



def pol_expand(iw,param):
        #polynomial expansion in iw up to 3rd order
        return param[0]+param[1]*iw+param[2]*iw**2+param[3]*iw**3

def analyt_cont(w,rp,ip):
        #analytical continuation of polynomial on the real axis
        return rp[0]+ip[1]*w-rp[2]*w**2-ip[3]*w**3 \
              +1j*(ip[0]-rp[1]*w-ip[2]*w**2+rp[3]*w**3)

#def difffunc(param,targetval,iw):
        #return sum(abs(targetval-pol_expand(iw,param)))
def difffunc(param,targetval,iw):
        #difference function for minimization
        return (targetval-pol_expand(iw,param))**2
try:
    import numpy as np
    import h5py as hdf5
    h5ustrs = hdf5.special_dtype(vlen=(str
                                       if sys.version_info >= (3, )
                                       else unicode))

except ImportError:
    print("Error: script requires the h5py package to run", file=err)
    exit(3)
    
usage, _dummy, desc = __doc__.split("\n", 2)
parser = optparse.OptionParser(usage=usage, description=desc)
parser.add_option("-i", "--iteration", type="int", dest="iter", default=-1,
                  help="""iteration to work with (default last)""", metavar="itno")
## Add more options here
(options, args) = parser.parse_args()

try:
    filename = args.pop(0)
    hf = hdf5.File(filename, "r+")
except IndexError:
    parser.print_usage()
    print("Error: no file specified", file=err)
    exit(2)
except IOError:
    print("Error: unable to open HDF5 output file `%s'." % filename, file=err)
    exit(2)


#get configuration data
cfg = ConfigObj(list(hf["configfile"].value))

# get desired iteration
iterpat = re.compile(r"^(?:dmft|stat)-\d+$")
iterations = sorted([k for k in hf.keys() if iterpat.match(k)])
try:
    iter = hf[iterations[options.iter]]
except IndexError:
    print("Error: iteration not present (or no iterations at all)", file=err)
    exit(2)

# get desired quantities 
siw = iter["siw/value"].value
iwin = hf['axes/iw'].value
dc = hf['start/dc/value'].value

beta = float(cfg["beta"])
nat = int(cfg["NAt"])
norb = int(cfg["Nd"])

nf = 4  # the number of datapoints of the self energy for fitting
wmax = float(8)/beta #8k_B*T default w-range on real axis
dw = 0.001 #1 meV default increment on real axis
w=np.arange(-wmax,wmax,dw)
#print "w",w
Sfreal = np.zeros((nat,norb,2,w.size),dtype=complex)
iwpos = iwin.compress((iwin>0).flat) #positive frequency axis
tiw = iwpos[0:nf]
#idx1 =         (iwin>0).nonzero()

#plt.figure()

for iat in range(0,nat):

    for iorb in range(0,norb):
      
        for isp in range(0,1):
            #print iat
            locdc = dc[iat,iorb,isp]
            #print siw[iat,iorb,isp,1]
            #print iwin[1]
            #print locdc
            siwpos =  siw[iat,iorb,isp,:].compress((iwin>0).flat)+locdc
            #fit the real part of the Matsubara self energy
            tval = siwpos[0:nf].real
            param = [tval[0], 0, 0, 0]
            poptr = leastsq(difffunc, param, args=(tval, tiw))
            diffr = difffunc(poptr[0],tval,tiw)
            #fit the imaginary part of the Matsubara self energy
            tval = siwpos[0:nf].imag
            param = [0, tval[0]/tiw[0], 0, 0]
            popti = leastsq(difffunc, param, args=(tval, tiw))
            #print popti[0]
            diffi = difffunc(popti[0],tval,tiw)
        #print sum(diffr),sum(diffi)
            if sum(diffr)+sum(diffi) > 10**(-8):
                print("Warning: leastsq-fit for Taylor expansion returned with larger error than threshold")
                print("Error[real fit,imag fit]: ",sum(diffr),sum(diffi))
            Sfreal[iat,iorb,isp,:] = analyt_cont(w,poptr[0],popti[0])
            #this can be used to plot the fit
            #riw = tiw
            #riw = np.insert(riw, 0, 0)
            #fval = pol_expand(riw,poptr[0])
            #plt.plot(tiw,tval,'x',riw,fval,'--')
            #this can be used for plotting the real self energy
            #if iat == 0 and isp==0:
                  #plt.plot(w,Sfreal[iat,iorb,isp,:].real,'-',w,Sfreal[iat,iorb,isp,:].imag,'--')

#plt.legend(['Orb1 - Re Sf','Orb1 - Im Sf','Orb2 - Re Sf','Orb2 - Im Sf','Orb3 - Re Sf','Orb3 - Im Sf',\
            #'Orb4 - Re Sf','Orb4 - Im Sf','Orb5 - Re Sf','Orb5 - Im Sf'])
#plt.show()

# re-introduce into hdf5 file
iter.require_group("sfreal")
iter["sfreal"].attrs["desc"] = "Taylor-fitted self energy on real axis" 
iter["sfreal"].attrs.create("axes", ["ineq", "band", "spin", "wreal"],
                            dtype=h5ustrs)

hf.require_dataset("axes/wreal", shape=w.shape,dtype=float,data=w)
iter.require_dataset("sfreal/value", shape=Sfreal.shape,dtype=complex,data=Sfreal)
# iter.create_dataset("gtau-relerr/error", data=???)

## cleanup
hf.close()
for j in range(0,w.size):
        print('%12f %12f %12f %12f %12f %12f %12f %12f %12f %12f %12f %12f %12f %12f %12f %12f %12f %12f %12f %12f %12f' %
              (w[j],Sfreal[0,0,0,j].real,Sfreal[0,0,0,j].imag,
               Sfreal[0,1,0,j].real,Sfreal[0,1,0,j].imag,
               Sfreal[0,2,0,j].real,Sfreal[0,2,0,j].imag,
               Sfreal[0,3,0,j].real,Sfreal[0,3,0,j].imag,
               Sfreal[0,4,0,j].real,Sfreal[0,4,0,j].imag,
               Sfreal[0,0,0,j].real,Sfreal[0,0,0,j].imag,
               Sfreal[0,1,0,j].real,Sfreal[0,1,0,j].imag,
               Sfreal[0,2,0,j].real,Sfreal[0,2,0,j].imag,
               Sfreal[0,3,0,j].real,Sfreal[0,3,0,j].imag,
               Sfreal[0,4,0,j].real,Sfreal[0,4,0,j].imag  ))
      #print '%2f %3f %4d' % (x, x*x, x*x*x)



