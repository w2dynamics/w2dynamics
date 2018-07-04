#!/usr/bin/env python
""" %prog [--help] <file.hdf5>

Adds gloctau: Fourier transform of gloc with p-band error estimate added
"""
from sys import exit, stderr as err
import optparse, re
try:
    import numpy as np
    import h5py as hdf5
except ImportError:
    print >> err, "error: script requires the h5py package to run"; exit(3)

import auxiliaries.transform as tf

__version__ = "1.1"

usage, _dummy, desc = __doc__.split("\n", 2)
parser = optparse.OptionParser(usage=usage, description=desc, version="%prog "+__version__)

parser.add_option("-i", "--iteration", type="int", dest="iter", default=-1,
                  help="iteration to work with (default last)", metavar="itno")

(options, args) = parser.parse_args()

try:
    filename = args.pop(0)
    hf = hdf5.File(filename, "r+")
except IndexError:
    parser.error("no file specified")
except IOError, e:
    parser.error("unable to open HDF5 output file `%s'." % filename)

try:
    # get config
    config = hf[".config"].attrs
    beta = float(config["general.beta"])

    # Get desired iteration
    iterpat = re.compile(r"^(?:dmft|stat)-\d+$")
    iterations = sorted([k for k in hf.iterkeys() if iterpat.match(k)])

    try:
        iter_node = hf[iterations[options.iter]]
    except IndexError:
        parser.error("iteration %d not present" % options.iter)

    metainfo = hf[".quantities"].require_group("gloctau")
    metainfo.attrs["desc"] = "G_loc(tau) with p-band error estimate"
    metainfo.attrs["axes"] = ["ineq", "band", "spin", "tau"]

    taus = hf[".axes/tau"].value
    iw = 1j * hf[".axes/iw"].value

    for ineq_name, ineq_node in iter_node.iteritems():
        # Checks for old run
        if not ineq_name.startswith("ineq-"):
            continue

        if "gloctau" in ineq_node:   # replace with your quantity name
            print "WARNING: overriding"
            del ineq_node["gloctau"]

        gtau_mean = ineq_node["gtau/value"].value
        gtau_error = ineq_node["gtau/error"].value
        glociw = ineq_node["gloc/value"].value

        mat2tau_matrix = tf.mat2tau(beta, 'fermi', glociw.shape[-1], taus)
        gloctau = -tf.transform(mat2tau_matrix, glociw - 1/iw) + .5

        num_dbands = gtau_mean.shape[0]
        d_part = slice(None, num_dbands)
        p_part = slice(num_dbands, None)
        if not np.allclose(gloctau.imag, 0):
            raise RuntimeError("Gloctau is not real????")
        if not np.allclose(gloctau[d_part], gtau_mean):
            print "WARNING: d-part of G(tau) does not match with impurity"

        gloctau = gloctau.real
        gloctau_err = np.empty_like(gloctau)
        gloctau_err[d_part] = gtau_error
        gloctau_err[p_part] = gloctau_err[d_part].mean(0)[None]

        gloctau_group = ineq_node.create_group("gloctau")
        gloctau_group.create_dataset("value", data=gloctau)
        gloctau_group.create_dataset("error", data=gloctau_err)

    print >> err, "success."
finally:
    hf.close()  # cleanup
