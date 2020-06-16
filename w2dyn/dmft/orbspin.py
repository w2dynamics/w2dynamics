"""Auxiliary functions for working with spins and orbitals"""
from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from warnings import warn
import numpy as np

from . import _compat as _linalg

STD_TOL = 1e-8

def allzeros(arr, tolerance=STD_TOL):
    return (np.abs(arr) <= tolerance).all()

def is_spindiagonal(arr, tolerance=STD_TOL):
    chk = arr.copy()
    spin = np.arange(arr.shape[-1])
    chk[:, :, spin, :, spin] = 0
    return allzeros(chk, tolerance)

def is_paramag(arr, tolerance=STD_TOL):
    spdiag = np.diagonal(arr, 0, 2, 4)
    return allzeros(spdiag[..., 1:] - spdiag[..., :1], tolerance)

def is_diagonal(arr, tolerance=STD_TOL):
    """Check if array is band/spin diagonal"""
    nflavours = arr.shape[-2] * arr.shape[-1]
    flavour = np.arange(nflavours)
    arr = arr.reshape(-1, nflavours, nflavours)
    offdiag = arr.copy()
    offdiag[:, flavour, flavour] = 0  # operate on copy
    return allzeros(offdiag, tolerance)

def set_diagonal(arr, diag_values):
    """Sets the band/spin diagonal of a (..., band, spin, band, spin) array"""
    nflavours = arr.shape[-2] * arr.shape[-1]
    flavour_idx = np.arange(nflavours)
    view = arr.view()
    view.shape = -1, nflavours, nflavours
    view[:, flavour_idx, flavour_idx] = diag_values

def promote_diagonal(diag):
    """Construct (..., band, spin, band, spin) array from band/spin-diagonal"""
    orig_shape = diag.shape
    nflavours = orig_shape[-2] * orig_shape[-1]
    diag = diag.reshape(-1, nflavours)
    ziw_full = np.zeros(diag.shape + (nflavours,), diag.dtype)
    flavour = np.arange(nflavours)
    ziw_full[:, flavour, flavour] = diag
    return ziw_full.reshape(orig_shape + orig_shape[-2:])

def extract_diagonal(arr):
    """Extract band/spin-diagonal from (..., band, spin, band, spin) array"""
    orig_shape = arr.shape
    nflavours = orig_shape[-2] * orig_shape[-1]
    arr = arr.reshape(-1, nflavours, nflavours)
    # note that the diagonal will be the last axis, so we need to transpose
    ziw_diag = arr.diagonal(0, -2, -1).copy()
    return ziw_diag.reshape(orig_shape[:-2])

def trace(arr):
    return arr.trace(0, -3, -1).trace(0, -2, -1)

def warn_offdiagonal(arr, tolerance=STD_TOL):
    """Checks for off-diagonal terms, then prints a warning if there are any"""
    if is_diagonal(arr, tolerance): return
    # Now we have some problem, provide a detailled report
    # Analyse diagonal entries
    diag_magn = np.abs(extract_diagonal(arr))
    diag_mean = diag_magn.mean()
    diag_max = diag_magn.max()
    # Analyse off-diagonal entries
    nflavours = arr.shape[-2] * arr.shape[-1]
    noff = np.prod(arr.shape[:-4]) * nflavours * (nflavours - 1)
    off_magn = np.abs(arr)
    set_diagonal(off_magn, 0)
    nflagged = (~(off_magn <= tolerance)).sum()
    where_max = np.unravel_index(off_magn.argmax(), off_magn.shape)
    off_max = off_magn[where_max]
    off_mean = off_magn.mean()
    # Print warning
    warn("\n\tDetected %d (%.2f%%) non-zero offdiagonal elements (tol = %g):"
         "\n\tOffdiagonal magnitudes   mean: %12g   max: %12g"
         "\n\tDiagonal magnitudes      mean: %12g   max: %12g"
         "\n\tIndices of maximum offdiagonal: %s" %
         (nflagged, 100*nflagged/noff, tolerance, off_mean, off_max, diag_mean,
          diag_max, where_max),
         UserWarning, 2)

def invert(arr):
    """Performs band/spin inverse of (..., band, spin, band, spin) array"""
    res_shape = arr.shape
    nflavours = res_shape[1] * res_shape[2]
    arr = arr.reshape(-1, nflavours, nflavours)
    return _linalg.inv(arr).reshape(res_shape)

def diag_invert(arr):
    """Inversion for an array that is band/spin-diagonal."""
    res_shape = arr.shape
    nflavours = res_shape[1] * res_shape[2]
    arr = arr.reshape(-1, nflavours, nflavours)
    flavour = np.arange(nflavours)
    inv = np.zeros_like(arr)
    inv[:, flavour, flavour] = 1/arr.diagonal(0, 1, 2)
    return inv.reshape(res_shape)

def symm_spins(arr):
    """Symmetrise Hermitean array over the spin dimensions"""
    if arr.shape[-1] != 2:
        raise NotImplementedError("Only works for 2 spins for now")

    new_arr = np.zeros_like(arr)
    new_arr[...,0,:,0] = arr[...,0,:,0] + arr[...,1,:,1]
    new_arr[...,1,:,1] = arr[...,0,:,0] + arr[...,1,:,1]
    new_arr[...,0,:,1] = arr[...,0,:,1] + arr[...,1,:,0].conj()
    new_arr[...,1,:,0] = arr[...,1,:,0] + arr[...,0,:,1].conj()
    new_arr /= 2.
    return new_arr

def multiply(arr1,arr2):
    """
    multiply two arrays of shape (nfreq, b, s, b, s) in orbital space
    """
    nfreq = arr1.shape[0]
    nb = arr1.shape[1]
    ns = arr1.shape[2]
    arr1 = arr1.reshape((nfreq, nb*ns, nb*ns))
    arr2 = arr2.reshape((nfreq, nb*ns, nb*ns))

    arr3 = np.einsum('vij,vjk->vik',arr1,arr2)
    return arr3.reshape((nfreq, nb, ns, nb, ns))

