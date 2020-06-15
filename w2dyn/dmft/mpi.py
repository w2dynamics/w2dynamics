"""Collection of MPI-related things"""
from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
import sys
import numpy as np

from mpi4py import MPI as mpi
# NOTE: For a comprehensive API documentation of mpi4py, refer to any one
#       parallel universe where they actually bothered to write more than just
#       repeating the name of the respective function.


DEBUG = True
MPI_COMM_WORLD = mpi.COMM_WORLD

class MPIProxy(object):
    """Object that transparently wraps around another object, providing a
    starting point for MPI-parallelised versions of other objects.

    This is basically a trivial implementation of the proxy pattern [1]: it
    wraps an inner object `_inner` in a transparent way, i.e., it redirects all
    methods and attributes to this instance, such that the inner and the proxy
    behave identically, so the calling code should not use rank-based forks.

    This is in itself not very useful. However, derived classes can substitute
    some or all (single-core) methods with a MPI-parallelised versions that
    expose the same interface, thus providing a drop-in replacement for the
    non-parallelised version.  Note that every proxy should be used identically
    on all cores.

    MPIProxy also provides a `_controller`, which is an instance of
    `MPIProxy.Controller` and keeps track of the MPI topology. Derived proxies
    may inherit from this controller to provide attributes and methods that aid
    MPI parallelisation.

    [1]: E Gamma et al. - Design Patterns. ISBN 0-201-63361-2
    """
    __slots__ = "_inner", "_controller"

    class Controller(object):
        def __init__(self, mpi_comm, mpi_root):
            global MPI_COMM_WORLD
            if mpi_comm is None: mpi_comm = MPI_COMM_WORLD
            self.comm = mpi_comm
            self.root = mpi_root
            self.rank = mpi_comm.Get_rank()
            self.size = mpi_comm.Get_size()
            self.is_root = self.rank == mpi_root

        def rdebug(self, fmt, *params):
            """Write debugging message to STDERR, but only on root"""
            if DEBUG and self.is_root:
                sys.stderr.write("MPI: %s\n" % (str(fmt) % params))

        def on_root(self, func):
            """Execute something on root only, but broadcast the result"""
            if self.is_root: result = func()
            else: result = None
            return self.comm.bcast(result)

    class Aggregate(object):
        def __init__(self, name):
            self.name = name

        def __get__(self, obj, objtype=None):
            if obj is None: return self   # descriptor protocol
            # Aggregate data from nodes
            comm = obj._controller.comm
            rankres = getattr(obj._inner, self.name)
            return comm.allreduce(rankres)

        def __set__(self, obj, value): raise NotImplementedError()

    def __init__(self, inner, mpi_comm=None, mpi_root=0, controller=None):
        if controller is None:
            controller = self.Controller(mpi_comm, mpi_root)
        MPIProxy._inner.__set__(self, inner)
        MPIProxy._controller.__set__(self, controller)

    def __getattr__(self, attr):
        """Forward all other attributes to the inner class"""
        return getattr(self._inner, attr)

    def __setattr__(self, attr, value):
        """Forward all other attributes to the inner class"""
        return setattr(self._inner, attr, value)

class MPIFrequencyParallel(MPIProxy):
    """Object which distributes work on frequencies via MPI.

    The instance wraps around another instance (`inner`) and attempts
    to load-balance the work of functions over multiple cores using
    MPI (mpi4py). This is achieved by splitting the Matsubara axis into chunks.
    Example: assume 100 frequencies are distributed over 12 cores, then then
    the first 8 ranks receive 8 frequencies each, and the last 4 ranks receive
    9 frequencies each.

    All methods of this class should be called simultaneously on all cores.
    The total result (full Matsubara axis) will be returned on every core, not
    just the root.  Attributes and methods that are not in the class are
    forwarded to the encapsulated instance `inner`.
    """
    __slots__ = ()

    class Controller(MPIProxy.Controller):
        def __init__(self, mpi_comm, mpi_root):
            MPIProxy.Controller.__init__(self, mpi_comm, mpi_root)
            self.niw = None

        def set_niw(self, niw):
            if self.niw == niw: return
            rank = self.comm.Get_rank()
            nranks = self.comm.Get_size()
            iw_per_rank = niw // nranks
            iw_excess = niw - iw_per_rank * nranks

            self.niw = niw
            self.sizes = iw_per_rank * np.ones(nranks, int)
            if iw_excess:
                self.sizes[-iw_excess:] += 1

            slice_ends = self.sizes.cumsum()
            self.slices = list(map(slice, slice_ends - self.sizes, slice_ends))

            self.myniw = self.sizes[rank]
            self.myslice = self.slices[rank]
            self.myshare = float(self.myniw) / self.niw
            self.rdebug("distributing %d frequencies (%d per rank)",
                        niw, iw_per_rank)

        def allgather(self, rank_result):
            if self.myniw != rank_result.shape[0]:
                raise ValueError("Illegal frequency gathering operation")
            tot_shape = (self.niw,) + rank_result.shape[1:]
            tot_result = np.empty(tot_shape, rank_result.dtype)
            tot_result[...] = np.nan
            other_dims = np.prod(rank_result.shape[1:])

            # The sizes argument needs the total number of elements rather than
            # just the first axis. The type argument is inferred.
            self.comm.Allgatherv(rank_result,
                                 [tot_result, self.sizes*other_dims])
            return tot_result

        def parts(self, rank_result):
            for rank, curr_slice in enumerate(list(self.slices)):
                curr_result = self.comm.bcast(rank_result, rank)
                yield curr_slice, curr_result


class MPILattice(MPIFrequencyParallel):
    def gloc(self, iw, mu, siw):
        """Perform the k-space integral, load-balanced over cores."""
        self._controller.set_niw(iw.size)
        myslice = self._controller.myslice
        rank_result = self._inner.gloc(iw[myslice], mu, siw[myslice])
        return self._controller.allgather(rank_result)

    # -- gmodel() does not need to be parallelised --

    def trace(self, iw, siw, smom, strategy='local'):
        """For every core, return an object that allows to compute a series of
        traces of the Green's function for different chemical potentials"""
        return self.TraceFactory(self, iw, siw, smom, strategy)

    class TraceFactory:
        """A class matching the requirements of `tr_gloc()`"""
        def __init__(self, lattice, iw, siw, smom, strategy):
            self._controller = lattice._controller
            self._controller.set_niw(iw.size)
            self._inner = lattice._inner.trace(iw[self._controller.myslice],
                                               siw[self._controller.myslice],
                                               smom, strategy)

        def trgloc(self, mu):
            rank_result = self._inner.trgloc(mu)
            return self._controller.allgather(rank_result)

        def trgmodel(self, mu):
            rank_result = self._inner.trgmodel(mu)
            return self._controller.allgather(rank_result)


