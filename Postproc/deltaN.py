#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""" %prog [--help] <file.hdf5> 

Compute ΔN(k).
"""
from __future__ import print_function
import sys
from configobj import ConfigObj
#import optparse -- deprecated in favor of
import argparse
import re
from warnings import warn

try:
    import numpy as np
    import h5py as hdf5
    h5ustrs = hdf5.special_dtype(vlen=(str
                                       if sys.version_info >= (3, )
                                       else unicode))

except ImportError:
    print("error: script requires the h5py package to run", file=sys.stderr)
    sys.exit(3)

import w2dyn.auxiliaries.wien2k as wien2k
import w2dyn.auxiliaries.config as config
import w2dyn.auxiliaries.input as _input
import w2dyn.dmft.orbspin as orbspin 

__version__ = "0.2"

myname = 'delta-n'

usage, _dummy, desc = __doc__.split("\n", 2)
parser = argparse.ArgumentParser(description=desc, version="%prog "+__version__)

parser.add_argument('hf', type=config.Hdf5Type('r+', notify=sys.stderr), metavar='hdf5-file',
                    help='saved run, must include ‘axes/k-points’ (default: ‘latest’)')
parser.add_argument("-i", "--iteration", type=int, default=-1, metavar="ITNO",
                    help="iteration to work with (default last)")
parser.add_argument("-u", "--update", action='store_true',
                    help="update old run if present")
parser.add_argument("-C", "--config-file", default="Parameters.in",
                    help='config file to read (for ‘KListFile’), default: ‘Parameters.in’')
parser.add_argument("--debug", action='store_true',
                    help="store some additional information")

args = parser.parse_args()

cfg = config.get_cfg(args.config_file)

del parser

try:
    # Get desired iteration
    iterpat = re.compile(r"^(?:dmft|stat)-\d+$")
    iterations = sorted([k for k in args.hf.keys() if iterpat.match(k)])

    try:
        iter = args.hf[iterations[args.iteration]]
    except IndexError:
        print("iteration %d not present" % args.iteration, file=sys.stderr)
        sys.exit(1)

    # Checks for old run
    if myname in iter and not args.update:
        print("error: quantity already present (-u to replace)", file=sys.stderr)
        sys.exit(1)

    # HERE STARTS THE REAL WORK
    # Get quantities from iteration: use hf.py to find out their names.
    # "q/value" stores the mean, "q/error" stores the error, if present.  
    # .value selects the whole quantity array, but you can also use
    # square brackets to select part of it.

    has_spin_orbit = cfg["General"]["DOS"] == 'ReadInSO'
    if has_spin_orbit:
       print("Spin-Orbit Hamiltonians not supported - adapt wien2k.py", file=sys.stderr)
       sys.exit()

    hkfile = file(cfg["General"]["HkFile"], "r")
    Hk, kpoints = _input.read_hamiltonian(hkfile, has_spin_orbit)
    beta        = args.hf['.config'].attrs['general.beta']
    w           = args.hf['.axes/iw']          .value
    try:       # FIXME: this is wrong for DMFT iterations before last!
        mu_dmft = args.hf['finish/mu/value']   .value
    except KeyError:
        warn("`finish/mu' not found, using iteration/mu instead")
        mu_dmft = iter   ['mu/value']          .value
    mu_lda      = args.hf['start/lda-mu/value'].value
    
    nat         = len(cfg["Atoms"])
    nneq       = len([ineq for ineq in args.hf['start'].keys() if ineq[0:4]==u'ineq'])

    nw = len(w)
    Siw = np.zeros((0, 2, nw))
    DC  = np.zeros((0, 2))
    istart = 0
    for i in range(1,nneq+1):
        ineq = "ineq-%03d"%i
        s = iter[ineq]['siw/value'].value
        nds = args.hf['.config'].attrs['atoms.'+str(i)+'.nd']
        nps = args.hf['.config'].attrs['atoms.'+str(i)+'.np']
        d = iter['dc-latt/value'].value[istart:nds+nps]
        d=orbspin.extract_diagonal(d)
        # We upcast the selfenergy to (nbands x 2 x nw) here.  This
        # could be optimized by passing Np(i) to delta_N and upcasting
        # from Nd(i) to Nbands for each w individually.
        Siw = np.vstack((Siw, s, np.zeros((nps, 2, nw))))
        DC  = np.vstack((DC,  d))
        istart = nds+nps

    wien2k.init(cfg, kpoints, sys.stdout)

    if args.debug:
        (dnk, g1, g2) = wien2k.delta_N_debug(beta, w, mu_lda, mu_dmft, Hk, Siw, DC, nat, nneq)
    else:
        dnk = wien2k.delta_N(beta, w, mu_lda, mu_dmft, Hk, Siw, DC, nat, nneq)

    # Re-introduce into HDF5 file: this creates a new quantity in the
    # HDF5 file with your result and stores some metadata, which is
    # needed by hf.py to correctly display your data

    if myname in iter and args.update:
        del iter[myname]
    else:
        g = args.hf['.quantities'].create_group(myname)
        g.attrs['desc'] = "DMFT-induced change of occupancy matrix"
        g.attrs.create('axes', ["lda-band", "lda-band", "spin", "kpoint"],
                       dtype=h5ustrs)

    iter.create_group(myname).create_dataset('value', data=dnk)

    if args.debug:
        iter[myname].create_dataset('g',   data=g1)
        iter[myname].create_dataset('gks', data=g2)

finally:
    args.hf.close()  # cleanup
