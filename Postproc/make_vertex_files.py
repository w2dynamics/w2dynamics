#!/usr/bin/env python3
""" %prog [--help] <file.hdf5> 

Generates vert_chi, vert_chi_pp and gm_wim files from HDF5
"""
from __future__ import division, print_function
from itertools import product as iter_product
from sys import exit, stderr as err
import optparse, re
try:
    import numpy as np
    import h5py as hdf5
except ImportError:
    print("error: script requires the h5py package to run", file=err)
    exit(3)
try:
    import w2dyn.auxiliaries.postprocessing as pp
except ImportError:
    print("error: script requires w2dynamics", file=err)
    exit(3)

__version__ = "1.0"

usage, _dummy, desc = __doc__.split("\n", 2)
parser = optparse.OptionParser(usage=usage, description=desc,
                               version="%prog "+__version__)
parser.add_option("-i", "--iteration", type="int", dest="iter", default=-1,
                  help="iteration to work with (default last)",
                  metavar="itno")
parser.add_option("-w", "--write-freq", action="store_true", default=False,
                  help="write frequencies instead of indices")
(options, args) = parser.parse_args()

try:
    filename = args.pop(0)
    hf = hdf5.File(filename, "r")
except IndexError:
    parser.error("no file specified")
except IOError:
    parser.error("unable to open HDF5 output file `%s'." % filename)

try:
    # Get desired iteration
    outfile_version = tuple(hf.attrs["outfile-version"])
    iterpat = re.compile(r"^(?:dmft|stat)-\d+$")
    iterations = sorted([k for k in hf.keys() if iterpat.match(k)])
    try:
        iter = hf[iterations[options.iter]]
    except IndexError:
        parser.error("iteration %d not present" % options.iter)
    
    print("Fetching quantities...", file=err)
    if outfile_version >= (2, 0):
        # new-style
        print("New-style file", file=err)
        iter = iter["ineq-001"]
        hfaxes = hf[".axes"]
        uu_slice = 0, 0, 0, 0
        ud_slice = 0, 0, 0, 1
        u_slice = 0, 0
        beta = hf[".config"].attrs["general.beta"]
    else:
        # old-style (DEPRECATED)
        print("Old-style file", file=err)
        hfaxes = hf["axes"]
        uu_slice = 0, 0, 0, 0, 0
        ud_slice = 0, 0, 0, 0, 1
        u_slice = 0, 0, 0
        entry = next(v for v in hf["configfile"].value if v.startswith("beta"))
        beta = float(entry.split("=")[1])

    giw = iter["giw/value"].value
    g4iw = iter["g4iw/value"].value
    g4iw_pp = iter["g4iw-pp/value"].value

    iw_one = hfaxes["iw"].value
    iw4f = hfaxes["iwf-g4"].value
    iw4b = hfaxes["iwb-g4"].value

finally:
    hf.close()  # cleanup

print("Generating quantities...", file=err)

# subtract the GG terms
chi_ph = pp.get_chi_ph(giw, g4iw, None, None)*beta
chi_pp = pp.get_chi_pp(giw, g4iw_pp, None, None)*beta

# transposition to move the bosonic frequency to the front
order = (2,0,1)

# select the appropriate slices
chi_ph_uu = chi_ph[uu_slice].transpose(order)
chi_ph_ud = chi_ph[ud_slice].transpose(order)
chi_pp_uu = chi_pp[uu_slice].transpose(order)
chi_pp_ud = chi_pp[ud_slice].transpose(order)

if options.write_freq:
    print("Writing frequencies ...", file=err)
    chi_fmt = "%17.10f" * 7 + "\n"
else:
    iw4b = np.arange(g4iw.shape[-1]) - g4iw.shape[-1]//2
    iw4f = np.arange(g4iw.shape[-2]) - g4iw.shape[-2]//2
    chi_fmt = "%7i" * 3 + "%20.10f" * 4 + "\n"

print("Building index ...", file=err)
iw, iv, ivp = np.broadcast_arrays(iw4b[:,None,None], iw4f[None,:,None],
                                  iw4f[None,None,:])

print("Writing vert_chi...", file=err)
f = file("vert_chi","w")
for tup in zip(iw.flat, iv.flat, ivp.flat,
               chi_ph_uu.real.flat, chi_ph_uu.imag.flat,
               chi_ph_ud.real.flat, chi_ph_ud.imag.flat):
    f.write(chi_fmt % tup)
f.close()

print("Writing vert_chi_pp...", file=err)
f = file("vert_chi_pp","w")
for tup in zip(iw.flat, iv.flat, ivp.flat,
               chi_pp_uu.real.flat, chi_pp_uu.imag.flat,
               chi_pp_ud.real.flat, chi_pp_ud.imag.flat):
    f.write(chi_fmt % tup)
f.close()

iwpos = slice(iw_one.size//2, None)

print("Writing gm_wim...", file=err)
f = file("gm_wim","w")
giw_uu = giw[u_slice + (iwpos,)]
for tup in zip(iw_one[iwpos], giw_uu.real, giw_uu.imag):
    print("".join("%20.10f" % f for f in tup), file=f)
f.close()

print("Complete.", file=err)
