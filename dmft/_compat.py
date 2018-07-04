"""Compatiblity module for numpy"""
from warnings import warn
import numpy as np

numpy_version = tuple(map(int, np.__version__.split(".")[:3]))

def diagonal(arr, offset=0, axis1=0, axis2=1):
    """Works around diagonal breaking result shape for zero-size arrays"""
    if arr.size:
        return arr.diagonal(offset, axis1, axis2)
    # zero-size array
    shape = (arr.shape[:axis1] + arr.shape[axis1+1:axis2] +
             arr.shape[axis2+1:] + (arr.shape[axis1],))
    return np.empty(shape, arr.dtype)

def fromiter(iterator, dtype, shape):
    """Generalises `numpy.fromiter()` to multi-dimesional arrays.

    Instead of the number of elements, the parameter `shape` has to be given,
    which contains the shape of the output array. The first dimension may be
    `-1`, in which case it is inferred from the iterator.

    Note that in NumPy (< 1.8), errors in the iterator are not forwarded, but
    instead reported as "itertor too short" since fromiter does not only catch
    StopIteration, but any exception class (SciPy issue #1999).
    """
    res_shape = shape[1:]
    if not res_shape:
        # Fallback to the "normal" fromiter in case of 1-D arrays
        return np.fromiter(iterator, dtype, shape[0])

    # This wrapping of the iterator is necessary because when used with the
    # field trick, np.fromiter does not enforce consistency of the shapes
    # returned with the '_' field and silently cuts additional elements.
    def shape_checker(iterator, res_shape):
        for value in iterator:
            if value.shape != res_shape:
                raise ValueError("shape of returned object %s does not match"
                                 " given shape %s" % (value.shape, res_shape))
            yield value,

    return np.fromiter(shape_checker(iterator, res_shape),
                       [("_", dtype, res_shape)], shape[0])["_"]

if numpy_version >= (1, 8):
    from numpy.linalg import inv as inv
    from numpy.linalg import eigvals as eigvals
    from numpy.linalg import eigh as eigh

else:
    warn("Emulating numpy.linalg functions for old version (< 1.8).\n"
         "This will result in lower performance.", UserWarning, 2)

    def inv(a):
        """`numpy.linalg.inv` over the last two dimensions"""
        a = np.asarray(a)
        if a.ndim == 2: return np.linalg.inv(a)
        return np.reshape([np.linalg.inv(p) for p
                           in a.reshape(-1, a.shape[-2], a.shape[-1])],
                          a.shape)

    def eigvals(a):
        """`numpy.linalg.eigvals` over the last two dimensions"""
        a = np.asarray(a)
        if a.ndim == 2: return np.linalg.eigvals(a)
        return np.reshape([np.linalg.eigvals(p) for p
                        in a.reshape(-1, a.shape[-2], a.shape[-1])],
                        a.shape[:-1])


    def eigh(a):
        """`numpy.linalg.eigh` over the last two dimensions"""
        a = np.asarray(a)
        if a.ndim == 2: return np.linalg.eigh(a)
        ew, ev = zip(*(np.linalg.eigh(p) for p in
                       a.reshape(-1, a.shape[-2], a.shape[-1])))
        return np.reshape(ew, a.shape[:-1]), np.reshape(ev, a.shape)

