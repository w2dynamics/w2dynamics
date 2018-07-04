"""Ensures compatibility of the code with python versions >= 2.4"""

def alt_iter_product(*args, **kwds):
    "Alternative fast implementation of product for python < 2.6"
    def cycle(sequence, uplevel):
        while True:
            vals = next(uplevel)   # will forward uplevel's StopIteration 
            it = iter(sequence)
            try:
                while True: yield vals + (next(it),)
            except StopIteration:
                pass
    step = iter(((),))
    for pool in map(tuple, args)*kwds.get('repeat', 1):
        step = cycle(pool, step)
    return step

try:
    from itertools import product as iter_product
except ImportError:
    iter_product = alt_itertools_product
