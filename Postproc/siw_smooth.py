#!/usr/bin/env python
""" %prog [--help] [-i iter] [-u] [-t threshold] [-l maxorder] <file.hdf5> 

Computes siw-smooth and giw-smooth and adds it to output file
"""
from __future__ import print_function
import sys
from sys import exit, stderr as err
import optparse, re
try:
    import numpy as np
    import h5py as hdf5
    h5ustrs = hdf5.special_dtype(vlen=(str
                                       if sys.version_info >= (3, )
                                       else unicode))
except ImportError:
    print("error: script requires the h5py package to run", file=err)
    exit(3)

try:
    import w2dyn.auxiliaries.transform as tf
except ImportError:
    print("error: script cannot discover w2dynamics libraries", file=err)
    exit(3)
    

__version__ = "1.1"

usage, _dummy, desc = __doc__.split("\n", 2)
parser = optparse.OptionParser(usage=usage, description=desc, version="%prog "+__version__)

parser.add_option("-i", "--iteration", type="int", dest="iter", default=-1,
                  help="iteration to work with (default last)", metavar="n")
parser.add_option("-u", "--update", action="store_true", dest="update", default=False,
                  help="update old run if present")
parser.add_option("-t", "--threshold", type="float", dest="threshold",
                  help="use signal-to-noise ratio cutoff (in units of sigma)", metavar="r")
parser.add_option("-l", "--maxorder", type="int", dest="maxorder",
                  help="use legendre order cutoff (specify maximum order)", metavar="l")

(options, args) = parser.parse_args()

try:
    filename = args.pop(0)
    hf = hdf5.File(filename, "r+")
except IndexError:
    parser.error("no file specified")
except IOError:
    parser.error("unable to open HDF5 output file `%s'." % filename)

try:
    # Get desired iteration
    iterpat = re.compile(r"^(?:dmft|stat)-\d+$")
    iterations = sorted([k for k in hf.keys() if iterpat.match(k)])
    try:
        iter = hf[iterations[options.iter]]
    except IndexError:
        parser.error("iteration %d not present" % options.iter)
    
    # Checks for old run
    
    if "siw-smooth" in iter:   # replace with your quantity name
        if options.update:
            del iter["giw-smooth"]
            del iter["siw-smooth"]
        else:
            print("error: quantity already present (-u to override, -h for help)", file=err)
            exit(1)

    gleg = iter["gleg/value"].value
    gleg_error = iter["gleg/error"].value
    print("Initial number of coefficients:", gleg.shape[-1], file=err)

    # Use hard legendre order cut-off    
    if options.maxorder:
        print("Hard cut of %d coefficients" % (gleg.shape[-1]-options.maxorder-1), file=err)
        gleg = gleg[...,:options.maxorder+1]
        gleg_error = gleg_error[...,:options.maxorder+1]
    
    # Use soft sigma cut-off
    if options.threshold:
        noisy = np.abs(gleg)/gleg_error < options.threshold
        print("Number of cut noisy coefficients:", end=' ', file=err) 
        print(", ".join(map(str, np.sum(noisy, -1).flat)), file=err)
        gleg[noisy] = 0.

    iw = 1j * hf["axes/iw"].value
        
    gliw, gliw_error = tf.transform(-tf.leg2mat(gleg.shape[-1], iw.size), 
                                    gleg, gleg_error)
    
    grp = iter.create_group("giw-smooth")
    grp.attrs["desc"] = "Smoothened Green's function from Legendre" 
    grp.attrs.create("axes", ["ineq", "band", "spin", "iw"],
                     dtype=h5ustrs)
    grp.create_dataset("value", data=gliw)
    grp.create_dataset("error", data=gliw_error)
    
    # Now, compute G0(iw) in order to find the self energy
    fiw = iter["fiw/value"].value 
    
    if "mu-imp" in iter:
        g0iw_inv = iw - iter["mu-imp/value"][...,np.newaxis] - fiw[::-1]
    else:
        g0iw_inv = iw + iter["mu/value"] - fiw[::-1]
        
    sliw = g0iw_inv - 1/gliw
    
    grp = iter.create_group("siw-smooth")
    grp.attrs["desc"] = "Smoothened self energy from Legendre" 
    grp.attrs.create("axes", ["ineq", "band", "spin", "iw"],
                     dtype=h5ustrs)
    grp.create_dataset("value", data=sliw)

    print("success.", file=err)
finally:
    hf.close()  # cleanup
