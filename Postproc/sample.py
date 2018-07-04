#!/usr/bin/env python
""" %prog [--help] <file.hdf5> 

Add description of your script here. PLEASE DO THAT!
"""
from sys import exit, stderr as err
import optparse, re
try:
    import numpy as np
    import h5py as hdf5
except ImportError:
    print >> err, "error: script requires the h5py package to run"; exit(3)
    
__version__ = "1.1"

usage, _dummy, desc = __doc__.split("\n", 2)
parser = optparse.OptionParser(usage=usage, description=desc, version="%prog "+__version__)

parser.add_option("-i", "--iteration", type="int", dest="iter", default=-1,
                  help="iteration to work with (default last)", metavar="itno")
parser.add_option("-u", "--update", action="store_true", dest="update", default=False,
                  help="update old run if present")
# add more options here (see optparse documentation)

(options, args) = parser.parse_args()

try:
    filename = args.pop(0)
    hf = hdf5.File(filename, "r+")
except IndexError:
    parser.error("no file specified")
except IOError, e:
    parser.error("unable to open HDF5 output file `%s'." % filename)

try:
    # Get desired iteration
    iterpat = re.compile(r"^(?:dmft|stat)-\d+$")
    iterations = sorted([k for k in hf.iterkeys() if iterpat.match(k)])
    try:
        iter = hf[iterations[options.iter]]
    except IndexError:
        parser.error("iteration %d not present" % options.iter)
    
    # Checks for old run
    if "gtau-relerr" in iter:   # replace with your quantity name
        if options.update:
            del iter["gtau-relerr"]
        else:
            print >> err, "error: quantity already present (-u to override, -h for help)"
            exit(1)

    # HERE STARTS THE REAL WORK
    # Get quantities from iteration: use hf.py to find out their names.
    # "q/value" stores the mean, "q/error" stores the error, if present.  
    # .value selects the whole quantity array, but you can also use
    # square brackets to select part of it.
    
    gtau_mean = iter["gtau/value"].value
    gtau_error = iter["gtau/error"].value
    
    # Do computation: See numpy documentation on what you can do.
    
    gtau_rel_error = gtau_error/gtau_mean
    
    # Re-introduce into HDF5 file: this creates a new quantity in the
    # HDF5 file with your result and stores some metadata, which is
    # needed by hf.py to correctly display your data

    iter.create_group("gtau-relerr")
    iter["gtau-relerr"].attrs["desc"] = "Relative error of G(tau)" 
    iter["gtau-relerr"].attrs["axes"] = ["ineq", "band", "spin", "tau"]
    
    iter["gtau-relerr"].create_dataset("value", data=gtau_rel_error)
    # iter.create_dataset("gtau-relerr/error", data=???)

    print >> err, "success."
finally:
    hf.close()  # cleanup
