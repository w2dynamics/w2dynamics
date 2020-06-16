#!/usr/bin/env python
"""Helper script to generate a bootstrap error for the quantities (siw, gtau,
gtau_worm ...) of the w2dynamics hdf5 output file. Generates plain text files
which can be processed by MaxEnt.py afterwards.

This helper script requires the scikits bootstrap library written by
Constantine Evans.

Note that the current script is still in an experimental stage and may require
some manual tuning.

Author: Patrik Gunacker
"""
from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
import numpy as np
import h5py as hdf5
import sys
import optparse
import multiprocessing
from multiprocessing import Pool

import auxiliaries.quantities as qttys
import scikits.bootstrap as boot


def ci_eval(samples):
    # alpha sets the confidence interval to 1 sigma
    # bootstrap gives us a lower and upper errorbar
    # we assume them to be almost equal such that
    # a simple average is justified
    return np.average(boot.ci(samples,
                              statfunction=(lambda x:
                                            np.average(np.abs(x))),
                              alpha=(1.-0.6827),
                              n_samples=5000,
                              method='bca',
                              output='errorbar'))


class Enumerator:
    def __init__(self):
        pass

    def __getitem__(self, i):
        return int(i)


if __name__ == '__main__':
    __version__ = (1, 0)

    parser = optparse.OptionParser(usage="%prog [OPTIONS] <hdf file> ...",
                                   version=".".join(map(str, __version__)))
    parser.add_option("--quantity", type="string", metavar="string",
                      dest="tag", default="siw",
                      help='QMC Quantity which will be bootstrapped')
    parser.add_option("--dmft-steps-begin", type="int",
                      dest="dmft_steps_begin", default=-1)
    parser.add_option("--dmft-steps-end", type="int",
                      dest="dmft_steps_end", default=-1)
    parser.add_option("--stat-steps-begin", type="int",
                      dest="stat_steps_begin", default=1)
    parser.add_option("--stat-steps-end", type="int",
                      dest="stat_steps_end", default=-1)
    parser.add_option("--filter", dest="filterfile")

    (options, args) = parser.parse_args()

    if not args:
        parser.error("expecting one or more HDF5/dat filenames")

    pool = Pool(processes=multiprocessing.cpu_count())

    for filename in args:
        try:
            hf = hdf5.File(filename, "r")
        except IOError:
            parser.error("invalid HDF5 output file: %s" % filename)

        file_version = tuple(hf.attrs["outfile-version"])
        StatSteps = qttys.MetaQttyContainer("config", hf, file_version) \
                         .select(qttys.
                                 SelectorPool(qttys.
                                              Selector("*statisticsteps"))) \
                         .values()[0]
        DMFTsteps = qttys.MetaQttyContainer("config", hf, file_version) \
                         .select(qttys.
                                 SelectorPool(qttys.
                                              Selector("*dmftsteps"))) \
                         .values()[0]

        if (options.stat_steps_begin <= 0
           or options.stat_steps_begin > StatSteps):
            options.stat_steps_begin = StatSteps + 1
        if options.stat_steps_end <= 0 or options.stat_steps_end > StatSteps:
            options.stat_steps_end = StatSteps + 1
        statsteps = options.stat_steps_end - options.stat_steps_begin
        if statsteps < 0:
            statsteps = 0

        if (options.dmft_steps_begin <= 0
           or options.dmft_steps_begin > DMFTsteps):
            options.dmft_steps_begin = DMFTsteps + 1
        if options.dmft_steps_end <= 0 or options.dmft_steps_end > DMFTsteps:
            options.dmft_steps_end = DMFTsteps + 1
        dmftsteps = options.dmft_steps_end - options.dmft_steps_begin
        if dmftsteps < 0:
            dmftsteps = 0

        filteriter = []
        filterstat = 0
        filterdmft = 0
        if options.filterfile is not None:
            with open(options.filterfile, "r") as filterfile:
                for line in filterfile:
                    iterstr = line.strip()
                    filteriter.append(iterstr)
                    if iterstr.startswith("stat-"):
                        statsteps -= 1
                        filterstat += 1
                    elif iterstr.startswith("dmft-"):
                        dmftsteps -= 1
                        filterdmft += 1
                    else:
                        raise RuntimeWarning("Iteration to be filtered not "
                                             "recognized: {}".format(iterstr))

        print("Number of Samples to be taken into account = {}"
              .format(statsteps + dmftsteps), file=sys.stderr)

        if statsteps + dmftsteps == 0:
            print("No Bins given in {}".format(filename), file=sys.stderr)
            continue

        NIneqs = hf[".config"].attrs["general.nat"]

        # we refer to tau or iw array as ximag in the following
        if options.tag == "gtau" or options.tag == "gtau_worm":
            freqval = hf[("axes", ".axes")[file_version[0]-1]]["taubin"].value
        elif options.tag == "siw" or options.tag == "giw":
            freqval = hf[("axes", ".axes")[file_version[0]-1]]["iw"].value
        elif options.tag == "ntau-n0" or options.tag == "sztau-sz0":
            freqval = hf[("axes", ".axes")[file_version[0]-1]]["tausus"].value
        else:
            print("Warning: ascending integers will be used as x-axis")
            freqval = Enumerator()

        # loop over all ineqs
        for iineq in range(0, NIneqs):

            # full sample container
            # niter = "stat-" + "001"
            # ineq = "ineq-" + "%03d" % iineq

            # only acces one ineq at a time, we can directly evaluate the
            # content of the generator
            try:
                niter = "stat-" + "%03d" % options.stat_steps_begin
                val = next(qttys.ineq_quantity(hf[niter],
                                               options.tag,
                                               ineq=iineq))[0]
            except KeyError:
                niter = "dmft-" + "%03d" % options.dmft_steps_begin
                val = next(qttys.ineq_quantity(hf[niter],
                                               options.tag,
                                               ineq=iineq))[0]

            freqs = val.shape[-1]
            samples = np.zeros(shape=(statsteps + dmftsteps, freqs),
                               dtype=complex)

            # loop over band and spins
            for indices in np.ndindex(val.shape[:-1]):
                # prepare filename
                fname = options.tag + "_" + str(iineq) + "_" +\
                    "_".join(str(x) for x in indices) + ".dat"
                print("Writing to {}".format(fname), file=sys.stderr)

                fdat = open(fname, "w")
                print("#ineq = {}, indices = {}".format(iineq, indices),
                      file=fdat)

                # add all other samples to array
                # FIXME: we should read out a single band spin at a time
                offset = 0
                for isample in range(0, statsteps + dmftsteps
                                     + filterstat + filterdmft):
                    if isample < statsteps + filterstat:
                        niter = "stat-" + "%03d" % (options.stat_steps_begin
                                                    + isample)
                    else:
                        niter = "dmft-" + "%03d" % (options.dmft_steps_begin
                                                    + isample - statsteps)
                    if niter in filteriter:
                        offset += 1
                        continue
                    val_add = (next(qttys.ineq_quantity(hf[niter],
                                                        options.tag,
                                                        ineq=iineq))[0]
                               [indices][:])
                    samples[isample - offset, :] = val_add

                # calculate the mean
                val[indices][:] = np.mean(samples, axis=0)

                results = []
                # calculate the bootstrap error bar for each freq
                for ifreq in range(0, freqs):
                    results.append(pool.apply_async(ci_eval,
                                                    (samples[:, ifreq], )))

                for ifreq in range(0, freqs):
                    if options.tag == "siw" or options.tag == "giw":
                        print(freqval[ifreq],
                              val[indices + (ifreq, )].real,
                              val[indices + (ifreq, )].imag,
                              results[ifreq].get(),
                              file=fdat)
                    else:
                        print(freqval[ifreq],
                              val[indices + (ifreq, )].real,
                              results[ifreq].get(),
                              file=fdat)
