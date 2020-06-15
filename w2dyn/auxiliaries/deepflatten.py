"""
This module flattens nested python lists with numbers and numpy arrays as
entries. `types` and `shapes` extracts information from the input while
`flatten` destroyes it. With the output of `types` and `shapes` `restore` can be
called with a modified array to restore the original shape.
"""

import numbers
import numpy as np

def types(x):
    """Returns the argument's type or, if it is a potentially nested non-ndarray
    iterable, its constituents' types in the form of a potentially nested
    list."""
    # numpy arrays can give us the type directly
    if type(x) is np.ndarray: return x.dtype
    # ordinary numbers can give us the type directly
    if isinstance(x, numbers.Number): return type(x)
    # we assume x is iterable and call `types` recursively
    result = []
    for i in x: result.append(types(i))
    return result

def shapes(x):
    """Returns the argument's shape or, if it is a potentially nested non-ndarray
    iterable, its constituents' shapes in the form of a potentially nested
    list."""
    # numpy arrays can give us the shape directly
    if type(x) is np.ndarray: return x.shape
    # we use `1` as a symbol for an ordinary number
    if isinstance(x, numbers.Number): return 1
    # we assume x is iterable and call `shapes` recursively
    result = []
    for i in x: result.append(shapes(i))
    return result

def flatten(x):
    """Flattens an ndarray or a potentially nested list containing ndarrays and
    numbers into one one-dimensional ndarray output."""
    # numpy arrays can be flattened directly
    if type(x) is np.ndarray: return x.flatten()
    # ordinary numbers are flat
    if isinstance(x, numbers.Number): return np.array([x])
    # we assume x is iterable and call `flatten` recursively to append
    # everything to the result
    result = np.array([])
    for i in x: result = np.append(result, flatten(i))
    return result

def restore(x, shapes, types):
    """Restores a potentially nested list of ndarrays and numbers from its
    flattened form, its shapes, and types."""
    def restore_(x, shapes, types):
        # tuples represent the shape of numpy arrays
        if isinstance(shapes, tuple):
            size = np.prod(shapes)
            return x[:size].reshape(shapes).astype(types), size
        # we used `1` as a symbol for an ordinary number
        if shapes == 1:
            return types(x[0]), 1
        # we assume x is iterable and call `restore` recursively
        result = []
        offset = 0
        for s, t in zip(shapes, types):
            tmp = restore_(x[offset:], s, t)
            result.append(tmp[0])
            offset += tmp[1]
        return result, offset
    return restore_(x, shapes, types)[0]

