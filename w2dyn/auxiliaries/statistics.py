import numpy as np
from mpi4py import MPI as mpi


def join(*samples):
    """Joins a set of DistributedSample instances to one big sample"""
    if not samples:
        raise ValueError("Must pass list of at least one sample")

    mpi_comm = samples[0].mpi_comm
    if any(sample.mpi_comm != mpi_comm for sample in samples):
        raise ValueError("MPI communicator must agree for all samples")

    local = np.concatenate(tuple(sample.local for sample in samples), axis=0)
    ntotal = sum(sample.ntotal for sample in samples)
    return DistributedSample(local, mpi_comm, ntotal)


class DistributedSample:
    """Sample distributed over MPI ranks.

    A sample is basically an array, where the zeroth dimension
    corresponds to the observations, and (optional) other dimensions
    correspond to the components of a single observation. This class
    models a situation where the samples are distributed over MPI ranks,
    e.g.::

        array index:     | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 |
                            -------------   ---------   ---------
        stored in:            rank 0         rank 1      rank 2

    Functions for calculating various statistical quantities across the
    whole distributed sample are supported.
    """

    def __init__(self, local, mpi_comm, ntotal=None):
        """Create a new distributed sample.

        :arg local:
            (n,...)-array corresponding to n locally stored samples
        :arg mpi_comm:
            MPI communicator instance
        :arg ntotal:
            total number of samples (inferred by default)
        """

        self.mpi_comm = mpi_comm
        self.local = np.asarray(local)
        if self.local.ndim == 0:
            self.local = self.local.ravel()
        self.nlocal = self.local.shape[0]

        if ntotal is not None:
            if ntotal < self.nlocal:
                raise ValueError("ntotal must at least cover local")
            if self.mpi_comm.Get_size() == 1 and ntotal != self.nlocal:
                raise ValueError("ntotal inconsistent for single-core")

            self.ntotal = ntotal
        else:
            # Figure out total number of bins via MPI
            self.ntotal = np.array(self.nlocal, np.int64)
            self._reduce_inplace(self.ntotal)

    def apply(self, func, inplace=False):
        """Transform all samples by applying function `func` to the data.

        Note that in comparison to the input object `self` the
        returned DistributedSample instance will simply have all
        samples `x` replaced with the transformed samples `func(x)`,
        which will bias statistical estimates such as the mean if
        `func` is not linear. Consider using something like
        `DistributedJackknife.transform` instead if this is
        undesirable.

        :arg func:
            Function taking and returning an array. The dimensions of
            the returned array should be the same for any possibly
            contained sample.
        :arg inplace:
            Whether to transform the samples contained in self.local or
            return a new DistributedSample containing the transformed
            values
        """
        if inplace:
            self.local = np.array([func(self.local[i, ...])
                                   for i in range(self.local.shape[0])])
            return self
        else:
            return DistributedSample(np.array([func(self.local[i, ...])
                                               for i in range(self.local.shape[0])]),
                                     self.mpi_comm,
                                     self.ntotal)

    def _reduce_inplace(self, arr):
        # Perform MPI in-place reduce-to-all.
        self.mpi_comm.Allreduce(mpi.IN_PLACE, arr, op=mpi.SUM)

    def sum(self):
        """Return sum of sample values across bins."""
        if self.local is None:
            raise RuntimeError("No data in sample")
        result = np.asarray(self.local.sum(0))
        self._reduce_inplace(result)
        return result

    def _sqdiff(self, mean, dtype=float):
        # Return sum of squared differences to a given mean.
        sqdiff = np.asarray(
            sum(np.abs(bin_i - mean, dtype=dtype)**2 for bin_i in self.local))
        self._reduce_inplace(sqdiff)
        return sqdiff

    def _bidiff(self, mean, dtype=float):
        """Return outer product of element differences to a given mean
        (bilinear in differences).
        """
        mean = np.asarray(mean, dtype=dtype)
        assert self.local.shape[1:] == mean.shape
        flatmean = mean.ravel()
        flatdata = np.reshape(self.local, (self.local.shape[0], -1))
        diff = np.einsum('ij,ik->jk',
                         (flatdata - flatmean),
                         np.conj(flatdata - flatmean))
        diff = np.reshape(diff,
                          mean.shape * 2)
        self._reduce_inplace(diff)
        return diff

    def mean(self):
        """Calculate and return the sample mean."""
        return self.sum() / self.ntotal

    def var(self, mean=None, ddof=1):
        """Calculate and return the sample variance :math:`s^2`.

            .. math::
                s^2 = \\frac{\\sum_i (y_i-\\bar{y})^2}{n-\\textrm{ddof}}

        :arg mean:
            If ``mean`` is not given, it is calculated.
        :arg ddof:
            delta degrees of freedom
        """
        if mean is None:
            result = self._sqdiff(self.mean())
        else:
            result = self._sqdiff(mean)
        result /= (self.ntotal - ddof)
        return result

    def stddev(self, var=None, ddof=1):
        """Calculate and return the sample standard deviation :math:`s`.

            .. math::
                s = \\sqrt{\\frac{\\sum_i (y_i-\\bar{y})^2}{n-
                    \\textrm{ddof}}}

        If the sample variance is given,

            .. math::
                s = \\sqrt{s^2}.

        :arg var:
            sample variance :math:`s^2`
        :arg ddof:
            delta degrees of freedom

        .. note::
            While the sample variance :math:`s^2` with ``ddof=1`` is an
            unbiased estimator, the sample standard deviation :math:`s`
            is not. It generally underestimates the standard deviation.
        """
        if var is None:
            var = self.var(ddof=ddof)
            return np.sqrt(var, out=var)
        else:
            return np.sqrt(var)

    def stderr(self, stddev=None, ddof=1):
        """Calculate and return the sample standard error of the mean
        :math:`\\textrm{SEM}`.

            .. math::
                \\textrm{SEM} = \\sqrt{\\frac{\\sum_i (y_i-\\bar{y})^2}
                {n(n-\\textrm{ddof})}}

        if the standard deviation is given,

            .. math::
                \\textrm{SEM} = \\frac{s}{\\sqrt{n}}.

        :arg stddev:
            sample standard deviation :math:`s`
        :arg ddof:
            delta degrees of freedom

        .. note::
            While the sample variance :math:`s^2` with ``ddof=1`` is an
            unbiased estimator, the sample standard error of the mean
            :math:`\\textrm{SEM}` is not. It generally underestimates
            the standard error of the mean.
        """
        if stddev is None:
            result = self._sqdiff(self.mean())
            result /= self.ntotal * (self.ntotal - ddof)
            return np.sqrt(result, out=result)
        else:
            return stddev / np.sqrt(self.ntotal)

    def cov(self, mean=None, ddof=1):
        """Calculate and return the sample covariance.

            .. math::
                K_{YY,jk} = \\frac{\\sum_i (y_{j,i}-\\bar{y_j})
                                           {(y_{k,i}-\\bar{y_k})}^*}
                                  {n-\\textrm{ddof}}

        :arg mean:
            If ``mean`` is not given, it is calculated.
        :arg ddof:
            delta degrees of freedom
        """
        result = self._bidiff(mean if mean is not None else self.mean())
        result /= self.ntotal - ddof
        return result

    def cov_of_mean(self, cov=None):
        """Calculate and return the covariance estimate of the sample mean.

        Return full covariance matrix of the sample mean. Its diagonal
        entries are equal to the squared standard error of the mean
        `stderr()`, whereas by comparison the diagonal entries of the
        sample covariance returned by `cov()` are equal to the sample
        variance, i.e. the square of the sample standard deviation.

        Note that similar to `stddev()` and `stderr()`, `cov()` and
        `cov_of_mean()` serve distinct purposes: `cov()` computes the
        naive estimate of the (ensemble) covariance matrix from the
        sample, whereas this function is the covariance estimate of
        the estimator that yields the sample mean.

        :arg cov:
            sample covariance :math:`K_{YY}`, calculated if not given
        """
        if cov is None:
            cov = self.cov()
        return cov / self.ntotal


class DistributedJackknife:
    """A class calculating the pseudovalues of the jackknife estimator.

    The jackknife is a resampling method that allows to estimate the
    variance and bias of a parameter. This can further be used to obtain
    the bias-corrected jackknife estimator. The input sample is assumed
    to be a :class:`DistributedSample`.

    Given :math:`n` input samples :math:`x_i` and a transformation
    function :math:`f()`, the goal is to estimate :math:`y=f(x)`, with
    :math:`x` and :math:`y` being the true values. A simple estimate
    would just be the transformed input mean :math:`f(\\bar{x})`, with
    :math:`\\bar{x}` being the input sample mean. In general this
    carries a bias of order :math:`1/n`. Using jackknife resampling this
    leading :math:`1/n` term can be estimated and removed, yielding the
    bias-corrected jackknife estimator.

    First, the so-called leave-one-out statistics :math:`x'_i` are
    calculated by

        .. math::
            x'_i = (n\\bar{x} - x_i) / (n-1).

    These are then transformed into bias-corrected pseudovalues
    :math:`y'_i` according to the following formular

        .. math::
            y'_i = n f(\\bar{x}) - (n-1) f(x'_i),

    The jackknife estimator :math:`\\hat{y}_{\\textrm{jk}}` can now be
    obtained simply by calculating the sample mean of the :math:`y'_i`.
    Other statistical quantities like the variance for example can also
    be calculated from the pseudovalues.
    """
    def __init__(self, sample, inplace=False):
        """Calculate the leave-one-out statistics :math:`x'_i`.

        :arg sample:
            the distributed input sample :math:`x_i`
        :arg inplace:
            If ``True`` the leave-one-out statistics are calculated in
            place of the input sample, otherwise a copy is created.
        """

        self.mpi_comm = sample.mpi_comm
        self.nlocal = sample.nlocal
        self.ntotal = sample.ntotal

        ssum = sample.sum()

        if self.ntotal < 2:
            # Handle the case n<2 specially, for which the leave-out statistics
            # do not exist.
            self.leaveoneout = None
        else:
            if inplace:
                self.leaveoneout = sample.local
                sample.local = None
            else:
                self.leaveoneout = sample.local.copy()

            # Compute leaveout statistics
            self.leaveoneout *= -1
            self.leaveoneout += ssum
            self.leaveoneout /= self.ntotal - 1

        self._mean = ssum
        self._mean /= self.ntotal

    def mean(self): return self._mean

    def transform(self, f):
        """Calculate and return the pseudovalues :math:`y'_i`.

        :arg f:
            transformation function :math:`f()`
        """
        if self.ntotal < 2:
            return DistributedSample([f(self._mean)], self.mpi_comm, self.ntotal)

        ymean = f(self._mean)
        yi = np.empty((self.nlocal,) + ymean.shape, ymean.dtype)
        yi[...] = self.ntotal * ymean
        del ymean

        for i in range(self.nlocal):
            yi[i] -= (self.ntotal - 1) * f(self.leaveoneout[i])

        return DistributedSample(yi, self.mpi_comm, self.ntotal)
