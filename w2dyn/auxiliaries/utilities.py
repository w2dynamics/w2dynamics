"""Various supporting functionality and common functionality of the
driver scripts cthyb and DMFT.py"""
from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
import numpy as np


def decomplexify(arr):
    """Takes a complex array-like arr and returns it as a real array with
    an axis for real/imag. part appended.
    """
    return np.concatenate((np.real(arr)[..., np.newaxis],
                           np.imag(arr)[..., np.newaxis]),
                          axis=-1)


def diagonal_covariance(dsample):
    """Takes a complex DistributedSample with data dimensions (iw, band,
    spin, band, spin) and returns the band/spin-diagonal real
    covariance of the mean with dimensions (band, spin, iw1, part1,
    iw2, part2)
    """
    diagsample = (dsample
                  # to (iw, spin, spin, band)
                  .apply(lambda qty:
                         np.diagonal(qty,
                                     axis1=1,
                                     axis2=3))
                  # to (iw, band, spin)
                  .apply(lambda qty:
                         np.diagonal(qty,
                                     axis1=1,
                                     axis2=2))
                  # positive iw only
                  .apply(lambda qty: qty[qty.shape[0]//2:]))
    return np.array([[diagsample
                      .apply(lambda qty, band=b, spin=s: qty[:, band, spin])
                      .apply(decomplexify)
                      .cov_of_mean()
                      for s in range(diagsample.local.shape[-1])]
                     for b in range(diagsample.local.shape[-2])])
