"""@package onlinestatistics"""
import numpy as np
import itertools

def nan_like(a):
    x = np.empty_like(a)
    x[...] = np.nan
    return x

class Sample(object):
    """ A statistical sample with first and second moment """
    def __init__(self, mean, sumqdiff=0, n=1):
        self.mean = mean
        self.sumqdiff = sumqdiff
        self.n = n
        
    def get_variance(self, ddof=0):
        """ Returns the current variance """
        assert ddof >= 0, "DDoF must be non-negative"
        if self.n <= ddof:
            return nan_like(self.sumqdiff)
        return self.sumqdiff/float(self.n - ddof)

    variance = property(lambda self: self.get_variance(),
                        doc="Biased variance of the sample")
    
    variance_unbiased = property(lambda self: self.get_variance(1),
                                 doc="Unbiased variance of the sample")
    
    stddev = property(lambda self: np.sqrt(self.get_variance()), 
                      doc="Biased standard deviation of the sample")
    
    stddev_unbiased = property(lambda self: np.sqrt(self.get_variance(1)), 
                               doc="Unbiased standard deviation of the sample")


class Aggregate(Sample):
    """ An on-line statistical aggregator for mean and variance.

    Implementation of an on-line algorithm for statistical variables, where mean
    and variance are updated continuously as new values are added, without the
    need to provide the whole sample for the computation. 

    This is intended for situations where the sample is too big for the memory or
    some continuous monitoring of the quantities is required. This class occupies
    at most four times the memory of a single data value at any time.

    The algorithm is stable and based on:
      [1] D.E. Knuth 1999, The Art of Computer Programming, Vol. 2
      [2] B.P. Welford 1962, Technometrics, Vol. 4(3), pp. 419-420
    """
    def __init__(self):
        """ Initialises an empty set of values """
        Sample.__init__(self, None, None, 0)

    def reset(self):
        """ Empties the set of values and clears the aggregates """
        self.__init__()

    def add(self, p):
        """ Adds an element to the sample, updating the aggregates """
        if not isinstance(p, Sample):
            p = Sample(p, None, 1)
        
        if self.n == 0:
            # First entry is handled specially
            self.n = p.n
            self.mean = p.mean
            if p.sumqdiff is not None:
                self.sumqdiff = p.sumqdiff
            else:
                self.sumqdiff = np.abs(np.zeros_like(p.mean)) # ensure realness
        else:
            # Generalising Knuth's method [1] to complex variables, but replacing
            # the variance addition with Welford's original formula [2] to avoid
            # multiplication with potentially big array and complex cancellation.
            oldcount = self.n
            delta = p.mean - self.mean
            self.n += p.n
            adjust = float(p.n)/self.n
            
            self.mean += delta * adjust
            self.sumqdiff += np.square(np.abs(delta)) * (oldcount*adjust)
            if p.sumqdiff is not None:
                self.sumqdiff += p.sumqdiff

    def add_all(self, ps):
        """ Adds all elements of the iterators """
        for p in ps:
            self.add(p)

    error = Sample.stddev_unbiased