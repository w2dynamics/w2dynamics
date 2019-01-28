"""Package abstracting details of the lattice model."""

from __future__ import division
import numpy as np
import scipy.ndimage
import scipy.signal
import scipy.optimize
import scipy.integrate
from warnings import warn

from ..auxiliaries import transform as tr

import _compat as _linalg
import orbspin as orbspin
import h5py
import copy

def fermi(x, beta, warn_complex=True):
    """Fermi function"""
    arg = np.asarray(float(beta) * x)
    result = np.empty_like(arg)

    if warn_complex and not (np.abs(arg.imag) < 1e-6).all():
        warn("Fermi function called with complex number.\n"
             "This will give the right result, but may indicate a problem\n"
             "with the calling code.", UserWarning, 2)

    # Create a piece-wise function in the argument: in the asymptotic region,
    # where Re(beta*x) > 50, the Fermi function is approximated by a step.
    # This is to quiet the warning for invalid divisions. np.where does not
    # do the trick because it does not short-circuit. In principle, one could
    # rewrite the function in terms of the tanh, however one runs into numpy
    # bug #5518 for complex numbers.
    inner_region = np.asarray(np.abs(arg.real) < 50)
    result[...] = arg < 0
    result[inner_region] = 1/(1. + np.exp(arg[inner_region]))
    return result

def dfermi(x, beta):
    """Derivative of the Fermi function `df(x)/dx`"""
    return -.25 * beta / np.cosh(-.5 * beta * x)**2

def rmap(values, mapping):
    """Performs an approximate reverse mapping of values to indices"""
    if (mapping.argsort() != np.arange(mapping.size)).any():
        raise ValueError("illegal mapping")
    atol = np.diff(mapping).min()/10
    mapping = np.repeat(mapping[:,None], 2, -1)
    mapping += -atol, +atol
    indices = mapping.reshape(-1).searchsorted(values, 'left')
    if not (indices % 2).all():
        problem = (indices % 2 == 0).nonzero()
        print values[problem], "\n", mapping
        raise ValueError("Inverse mapping failed: illegal values")
    return indices // 2


# ----- lattice base class -----

class Lattice:
    """Lattice which is aware of the way it performs the Dyson equation.

    A lattice class provides an abstraction of the way the self-energy is
    plugged into the Dyson equation to give either the local Green's function
    or the band/spin trace over it (useful for computing densities).

    Additionally, define the following read-only attributes or properties:
      - `hloc`:      local Hamiltonian as `(band, spin, band, spin)` array
      - `hmom2`:     second moment of Hamiltonian (same shape as `hloc`)
    """
    def __init__(self, beta, norbitals, nspins):
        if not (norbitals > 0 and nspins > 0):
            raise ValueError("norbitals and nspins must be positive")
        if beta <= 0:
            raise ValueError("beta must be positive")

        # nspins is thought to be a starting point for extension to Nambu space
        self.norbitals = norbitals
        self.nspins = nspins
        self.nflavours = norbitals * nspins
        self.beta = beta
        self.eye = np.eye(self.nflavours).reshape(norbitals, nspins,
                                                  norbitals, nspins)

    def compute_dos(self, fixmu, mu=None, totdens=None):
        """Compute non-interacting fillings and density of states.

        This function fills the following fields in the lattice class:
          - `w`:         frequency/energy axis for density of states
          - `dos`:       density-of-states as `(band, spin, w)` array
          - `mu`:        non-interacting chemical potential
          - `densities`: fillings at `mu` as `(band, spin)` array
        """
        raise NotImplementedError()

    def gloc(self, iw, mu, siw, siw_gw=None):
        """Compute the local Green's function from the quantitiy `A(iw)`.

        In principle, for this lattice, perform the sum:

                G(iw) = \sum_k (iw + \mu - H(k) - S(iw))^{-1},  where

        S(iw) is supplied in `siw` as 5-dimensional NumPy array, where the axes
        correspond to frequency, band, spin, band, spin.  The local Green's
        function G(iw) returned has the same dimensions as S(iw).
        """
        raise NotImplementedError()

    def gmodel(self, iw, mu, smom, smom_gw=None, strategy='local'):
        """Compute a model Green's function and its analytical sum.

        For this lattice, perform the (approximate) sum:

                G_m(iw) = \sum_k (iw + \mu - H_k - S)^{-1},  where
                      S = \Sigma_0 + \Sigma_1/iw + ...

        S is supplied in `smom` as 5-dimensional NumPy array, where the axes
        correspond to moment, band, spin, band, spin. Also compute the density
        matrix `D` analytically, i.e. over the whole Matsubara frequency range:

                      D = \sum_{iw=-\infty}^\infty G_m(iw)
        """
        if strategy == 'local':
            return self._model_integrate(iw, mu, self.hloc, smom, smom_gw)
        else:
            raise NotImplemented("unknown strategy")

    def trace(self, iw, siw, smom, strategy='local', use_gw=False, siw_gw=None, smom_gw=None):
        r"""Return factory for computing spin/band-traces for different mu.

        This class returns a factory (generating object), which can be
        called with different shifts of the chemical potential to yield the
        trace of the local Green's function:

              N(iw, \Delta\mu) = \sum_k \tr (iw + mu - H_k - S(iw))^{-1}

        This is useful in computing total filling, which is then just given by
        the frequency sum.  While this can in principle be computed by repeated
        calls to `gloc()`, for many Hamiltonians there exists a more efficient
        solution for a set of traces.
        """
        return Lattice.TraceFactory(self, iw, siw, smom, strategy, use_gw, siw_gw, smom_gw)

    class TraceFactory:
        """Functor for computing traces of Gloc(iw) for different mu.

        This functor satisfies the requirements to `Lattice.trace()`.  It
        provides a naive and slow implementation of the trace by just wrapping
        over the corresponding `Lattice.gloc()`.
        """
        def __init__(self, lattice, iw, siw, smom, strategy, use_gw=False, siw_gw=None, smom_gw=None):
            self.lattice = lattice
            self.iw = iw
            self.siw = siw
            self.smom = smom
            self.strategy = strategy
            self.use_gw = use_gw
            self.siw_gw = siw_gw
            self.smom_gw = smom_gw
            
        def trgloc(self, mu):
            """Return spin/band trace of the local Green's function"""
            gloc = self.lattice.gloc(self.iw, mu, self.siw, self.siw_gw)
            return gloc.trace(0,1,3).trace(0,1,2)

        def trgmodel(self, mu):
            """Return spin/band trace of the model Green's function and the
               corresponding density matrix"""
            gmodel, densmatrix = self.lattice.gmodel(self.iw, mu, self.smom, self.smom_gw)
            return gmodel.trace(0,1,3).trace(0,1,2), \
                   densmatrix.trace(0,0,2).trace(0,0,1)

    def _model_integrate(self, iw, mu, h, smom, smom_GW=None):
        """Helper function for computing model and densities"""
        if smom.shape[0] != 1:
            warn("only taking into account zeroth moment yet", UserWarning, 2)

        # ==> GW START
        if smom_GW is not None:
          if np.sum(smom) != 0:
            smom_GW_loc = smom_GW.mean(0)   # Must be ZERO for k-Average
            gmom0 = h + smom[0] + smom_GW_loc - mu*self.eye
          else:
            gmom0 = h + smom[0] - mu*self.eye
        else:
          gmom0 = h + smom[0] - mu*self.eye
        # ==> GW END

        gmom0 = gmom0.reshape(self.nflavours, self.nflavours)
        gvals, gbasis = np.linalg.eig(gmom0)
        iw = iw[:, np.newaxis]
        ginvbasis = np.linalg.inv(gbasis)

        giw = np.einsum('ij,wj,jl->wil', gbasis, 1/(1j*iw - gvals), ginvbasis)
        dens = np.einsum('ij,j,jl->il',
                         gbasis, fermi(gvals, self.beta), ginvbasis)

        tshape = -1, self.norbitals, self.nspins, self.norbitals, self.nspins
        return giw.reshape(tshape), dens.reshape(tshape[1:])

# ---- concrete lattices ----

class KspaceHamiltonian(Lattice):
    """Simple lattice given by a k-space Hamiltonian `H(k)`.

    Does not consider any symmetries and just performs a "naive" k summation
    instead of tetrahedronal fancitude.
    """
    def __init__(self, beta, hk, check_herm=False):
        if hk.ndim != 5:
            raise ValueError("H(k) must be (k,band,spin,band,spin)-array")
        if hk.shape[1] != hk.shape[3] or hk.shape[2] != hk.shape[4]:
            raise ValueError("H(k) must be a square matrix")
        if check_herm and not np.allclose(hk.conj().transpose(0,3,4,1,2), hk):
            warn("Hamiltonian is not Hermitean", RuntimeWarning, 2)

        # If the spin-dimension of the Hamiltonian is 2, then we have a full
        # spin-orbit calculation.
        self.nkpoints, norbitals = hk.shape[:2]
        self.spinorbit = (hk.shape[2] == 2)

        # nspins is thought to be a starting point for extension to Nambu space
        # H(k) will however not be depending on that, therefore we hardcode it
        # to 2 above.
        Lattice.__init__(self, beta, norbitals, nspins=2)

        if not self.spinorbit:
            self.hk = np.zeros((self.nkpoints, self.norbitals, 2,
                                self.norbitals, 2), hk.dtype)
            # in case of spin-independent Hamiltonians, we copy twice to the
            # diagonal
            self.hk[:, :, :1, :, :1] = hk
            self.hk[:, :, 1:2, :, 1:2] = hk
        else:
            self.hk = hk
            
        # set moments of H(k). Note that we need the second moment, rather
        # than the second central moment (see Markus' thesis, Eq. (4.63b)
        self.hloc = self.hk.mean(0)
        self.hmom2 = np.einsum("kasbt,kbtcu->ascu", self.hk,
                               self.hk)/self.nkpoints

    def compute_dos(self, fixmu, mu=None, totdens=None):
        """Compute non-interacting density of states by electron counting"""
        # Because it's the sanest choice for double-counting reasons but then
        # the user has to be aware that totdens is indeed important.
        fixmu = True

        nkpoints, nbands = self.hk.shape[:2]
        nflavours = 2 * nbands
        hk = self.hk.reshape(nkpoints, nflavours, nflavours)

        ee, ev = _linalg.eigh(hk)
        ev = ev.transpose(0,2,1)  # turn eigenbasis to eigenvectors

         # Flatten k dimension, collect all eigenvalues in one array
        ee = ee.reshape(nkpoints*nflavours)
        ev = ev.reshape(nkpoints*nflavours, nflavours)

        # Sort eigenvectors according to eigenvalues (ascending)
        index = ee.argsort()
        ee = ee[index]
        ev = ev[index]

        if not fixmu:
            # calculate at fixed chemical potential
            if mu is None: raise ValueError("must supply mu")
            mu_pos = np.searchsorted(ee, mu)
            self.mu = mu
        else:
            # find chemical potential by counting electrons
            if totdens is None: raise ValueError("must supply totdens")
            mu_pos = int(round(nkpoints * totdens))
            if mu_pos >= len(ee):
                warn("chemical potential outside LDA DOS")
                self.mu = ee[-1]
            elif mu_pos == 0:
                self.mu = ee[0]
            else:
                # try to model a possible gap
                self.mu = (ee[mu_pos] + ee[mu_pos-1])/2

        # Integrate eigenvectors on the energy axis up to the chemical potential
        self.densities = np.sum(np.absolute(ev[:mu_pos,:])**2/float(nkpoints),
                                axis=0).reshape(nbands, 2)

        # Try to estimate "natural" discretisation of the energy space
        discr = np.diff(ee)
        discr = discr[discr > 1e-10].mean()
        sigma = .05 #* discr  # on average, take smoothen over 2 energy levels

        # boundaries of the non-interacting density of states
        wborder = 3*sigma
        wstep = .001 # put into small bins, and let the filter work its magic
        wmin = np.around(np.amin(ee) - wborder, 3)
        wmax = np.around(np.amax(ee) + wborder, 3) + 1e-4
        w = np.arange(wmin, wmax, wstep)

        #wstep=10.*abs(wmax-wmin)/float(self.kpoints.shape[0])
        # Binning ee for dos
        # Use historgram in a loop because it only supports 1-D weights
        dos = [np.histogram(ee, w, weights=np.abs(ev[:,iflv])**2, density=True)[0]
               for iflv in range(nflavours)]
        dos = np.transpose(dos)  # -> w, flavour
        dos = scipy.ndimage.filters.gaussian_filter1d(dos, sigma/wstep, axis=0)

        self.w = w[:-1]
        self.dos = dos.reshape(self.w.size, nbands, 2)   # spin-dependency

    def gloc(self, iw, mu, siw, siw_gw=None):
        """Compute the local Green's function from the quantitiy `S(iw)`."""
        if not siw.size:   # some numpy does not like zero-sized arrays
            return siw.copy()

        ziw = (1j*iw + mu)[:, None, None, None, None] * self.eye - siw

        if self.spinorbit or not orbspin.is_spindiagonal(siw):
            # the Hamiltonian is a full spin-orbit Hk
            hk = self.hk.reshape(self.nkpoints, self.nflavours, self.nflavours)
            ziw = ziw.reshape(-1, self.nflavours, self.nflavours)

            glociw = sum(_linalg.inv(ziw - h_kfix) for h_kfix in hk)
            glociw /= self.nkpoints
            glociw = glociw.reshape(-1, self.norbitals, self.nspins,
                                        self.norbitals, self.nspins)
            if siw_gw is not None:
              raise NotImplemented("GW not implemented for Spin-Orbit Coupling!")
              
        else:
            # The Hamiltonian is not spin-dependent and the ziw is diagonal in
            # spin, so the inversion in the spin space is trivial.
            # self.hk is always "big" therefore, we extract the upup block here
            hk_upup = self.hk[:, :, 0, :, 0]
            # Here we take the main spin-diagonal --> iw, band, band, spin
            ziw_diag = ziw.diagonal(0, 2, 4).transpose(3, 0, 1, 2) # ->s,iw,b,b

            # In the paramagnetic case, we restrict ourselves to up-up block
            paramag = orbspin.is_paramag(siw)
            if siw_gw is not None:
              
              if np.sum(siw) != 0:
                siw_gw = siw_gw.diagonal(0,3,5).transpose(4,1,0,2,3)
                if paramag:
                  siw_gw = siw_gw[:1]
                hk_upup = hk_upup + siw_gw
                # Move the k-component into the first attribute of the array: [k,1,w,b1,b2]
                hk_upup = hk_upup.transpose(2,0,1,3,4)
                # At this point hk_upup and siw_diag are identical in their last 4 arguments: [1,w,b1,b2]


            if paramag:
                ziw_diag = ziw_diag[:1]

            gloc_diag = sum(_linalg.inv(ziw_diag - h_kfix) for h_kfix
                            in hk_upup)
            gloc_diag /= self.nkpoints

            glociw = np.zeros_like(ziw)
            for i in range(self.nspins):
                if paramag:
                    glociw[:,:,i,:,i] = gloc_diag[0]
                else:
                    glociw[:,:,i,:,i] = gloc_diag[i]

        return glociw

    def trace(self, iw, siw, smom, strategy='local', use_gw=False, siw_gw=None, smom_gw=None):
        """Functor for computing traces of Gloc(iw) for different mu."""

        if use_gw == 1:
          return self.TraceFactory(self, iw, siw, smom, strategy, use_gw, siw_gw, smom_gw)
        else:
          try:
              return self.EVTraceFactory(self, iw, siw, smom, strategy)
          except MemoryError:
              warn("Unable to use improved functor due to memory constraints, \n"
                  "falling back to naive trace.  Expect lower performance of \n"
                  "mu search.", RuntimeWarning, 2)
              return self.TraceFactory(self, iw, siw, smom, strategy)


    class EVTraceFactory(Lattice.TraceFactory):
        """Computing a set trace by changing to the eigenbasis of A - H

        This functor exploits two facts: (1) the trace of a matrix `tr(S)` is
        invariant under a change of basis and (2) the eigenbasis of `S` is the
        same as the eigenbasis of `S + c1`. Therefore:

            N(iw, mu) = \sum_k \tr (z(iw) - H(k) + \mu)^{-1}
                      = \sum_k \tr (S(iw,k) + \mu)^{-1}
                      = \sum_k \tr B(iw,k) (E(iw,k) + \mu)^{-1} B^{-1}(iw,k)
                      = \sum_k \sum_n 1/(e(iw,k,n) + \mu)

        where `e` just marks the eigenvalues of `ziw` and the `hk` combined.
        This has the same complexity `O(N*k**3)` as a "naive" computation, but
        scales as `O(N*k)` for subsequent calls.
        """
        def get_eigvals(self, ziw, hk):
            # This returns a (niw * nkpoints * nbands * nspins) array, which is
            # huge, so it is preallocated to ensure that it fails immediately
            # when there is insufficient memory and filled with an iterator to
            # avoid copying.  The argument to eigvals is bigger still by a
            # factor of (nbands * nspins), so iteration must be used.  We
            # pre-initialise the eigv array to ensure that the memory is used
            # and the code fails with OOM as early as possible.
            niw = ziw.shape[0]
            result_shape = np.broadcast(ziw[:,None], hk).shape[:-1]
            eigv = np.empty(result_shape, ziw.dtype)
            eigv[...] = np.nan
            for iiw, z_fixiw in enumerate(ziw):
                eigv[iiw] = _linalg.eigvals(z_fixiw - hk)
            return eigv.reshape(niw, -1)

        def __init__(self, lattice, iw, siw, smom, strategy):
            Lattice.TraceFactory.__init__(self, lattice, iw, siw, smom, strategy)
            nflavours = lattice.nflavours

            if not siw.size:   # some numpy does not like zero-sized arrays
                self.eigvals = np.zeros((0, 1))
                self.multiplicity = 1

            ziw = 1j * iw[:, None, None, None, None] * lattice.eye - siw
            if lattice.spinorbit or (not orbspin.is_spindiagonal(siw)):
                # the Hamiltonian is a full spin-orbit Hk
                hk = lattice.hk.reshape(-1, nflavours, nflavours)
                ziw = ziw.reshape(-1, nflavours, nflavours)
                self.eigvals = self.get_eigvals(ziw, hk)
                self.multiplicity = 1

            else:
                # The Hamiltonian is not spin-dependent and the ziw is diagonal
                # in spin, so the inversion in the spin space is trivial.
                # self.hk is always "big" therefore, we extract the upup block
                # --> we arrive at a kpoint, (spin), band, band array
                hk_upup = lattice.hk[:, np.newaxis, :, 0, :, 0]

                # Here we take the main spin-diagonal --> iw, spin, band, band
                ziw_diag = ziw.diagonal(0, 2, 4).transpose(0, 3, 1, 2)
                # In the paramagnetic case, we restrict ourselves to up-up
                paramag = orbspin.is_paramag(siw)
                if paramag:
                    ziw_diag = ziw_diag[:, :1]
                self.eigvals = self.get_eigvals(ziw_diag, hk_upup)
                self.multiplicity = lattice.nspins**paramag

            self.multiplicity /= float(lattice.nkpoints)

        def trgloc(self, mu):
            # We will typically operate at a low memory setting here because
            # of the huge eigenvalues array, so we split the calculation by
            # fermionic frequency. This is typically not a problem
            # performance-wise since Niw ~ O(1000), which is fast to loop.
            return np.array([(1./(mu + eigvals_iw)).sum() for eigvals_iw
                             in self.eigvals]) * self.multiplicity



class Bethe(Lattice):
    """Lattice with a semi-elliptic density of states.

    Note that this class does not support off-diagonal self-energies.
    """
    def _get_dos(self, w, integrated=False):
        eps = (w - self.crystalfield)/self.d
        eps = np.clip(eps, -1, 1)
        if integrated:
            # The integrated DOS
            return .5 + (eps * np.sqrt(1 - eps**2) + np.arcsin(eps))/np.pi
        else:
            # Semi-elliptic DOS
            return 2/(np.pi*self.d) * np.sqrt(1 - eps**2)

    def __init__(self, beta, half_bw, crystalfield=None):
        # crystalfield marks the centers of the DOS (hmean)
        if crystalfield is None:
            crystalfield = np.zeros_like(half_bw)

        self.crystalfield = np.asarray(crystalfield)
        self.d = np.asarray(half_bw)

        norbitals, nspins = self.d.shape
        Lattice.__init__(self, beta, norbitals, nspins)

        # add dummy iw dimension to satisfy requirements for promote_diagonal
        tstar2 = self.d**2 / 4.
        self.hloc = orbspin.promote_diagonal(self.crystalfield)
        self.hmom2 = orbspin.promote_diagonal(tstar2) + self.hloc**2

    def compute_dos(self, fixmu, mu=None, totdens=None):
        # find boundaries of the non-interacting problem
        wstep = self.d.mean()/20
        dmax = self.d.max() + wstep
        wmax = self.crystalfield.max() + dmax + wstep
        wmin = self.crystalfield.min() - dmax

        # calculate density of states
        self.w = np.arange(wmin, wmax, wstep)
        self.dos = self._get_dos(self.w[:,None,None])

        if fixmu:
            # Find mu by numerically inverting the integrated density of states
            self.mu = scipy.optimize.brentq(
                            lambda w: self._get_dos(w, True).sum() - totdens,
                            wmin, wmax)
        else:
            self.mu = mu

        self.densities = self._get_dos(self.mu, True)

    def gloc(self, iw, mu, siw, siw_gw=None):
        ziw = (1j*iw + mu)[:, None, None, None, None] * self.eye - siw
        orbspin.warn_offdiagonal(siw)
        zeta = orbspin.extract_diagonal(ziw - self.hloc)
        d2 = self.d[np.newaxis]**2
        glociw = 2*zeta/d2 * (1 - np.sqrt(1 - d2/zeta**2))
        glociw = orbspin.promote_diagonal(glociw)
        return glociw

class DensityOfStates(Lattice):
    """Lattice represented as density of states.

    Note that this class does not support off-diagonal self-energies.
    """
    @staticmethod
    def get_dos_bethe(w, d=2., w0=0., z=None):
        """Computes density of states for the Bethe lattice of coordination Z.

        Given a half bandwitdth `d` (this implies proper scaling of the
        hopping) and a band centre of `w0`, returns the DoS of the Bethe
        lattice for a given coordination number `z`. If z is omitted, this
        gives the well-known semi-elliptic DoS.
        """
        if np.any(d <= 0):
            raise ValueError("half bandwidth d is negative")
        eps = np.clip((w - w0)/d, -1, 1)
        dos = 2/(np.pi * d) * np.sqrt(1 - eps**2)

        if z is not None:   # correct for finite dimension
            if np.any(z < 2):
                raise ValueError("coordination number must be at least 2")
            dos /= z/(z - 1) - eps**2/(4 * z)

        return dos

    def __init__(self, beta, w, dos):
        try:
            nw, norbitals, nspins = dos.shape
        except ValueError:
            raise ValueError("dos must be a nw, norbital, nspins array")
        if w.shape != (nw,):
            raise ValueError("w and dos are not consistent in shape")

        Lattice.__init__(self, beta, norbitals, nspins)

        integrate = lambda fw: scipy.integrate.trapz(fw, w)
        cum_int = lambda fw: scipy.integrate.cumtrapz(fw, w, initial=0)

        if (dos < 0).any() or dos.imag.any():
            warn("dos is not a positive real function.", UserWarning, 2)

        dos_w_last = dos.transpose(1,2,0)
        norm = integrate(dos_w_last)
        if not np.allclose(norm, 1):
            warn("dos seems not to be properly normalised.\n Norms: %s" % norm,
                 UserWarning, 2)

        self.w = w
        self.dos = dos
        self.int_dos = cum_int(dos_w_last).transpose(2,0,1)

        self.hloc = integrate(dos_w_last * w)
        self.hmom2 = integrate(dos_w_last * w**2)

        self.hloc = orbspin.promote_diagonal(self.hloc)
        self.hmom2 = orbspin.promote_diagonal(self.hmom2)

    def compute_dos(self, fixmu, mu=None, totdens=None):
        # We don't have to compute the DOS since it is the input of this thing
        if fixmu:
            # Find mu by numerically inverting the integrated density of states
            nw = self.w.size
            int_dos_tot = self.int_dos.reshape(nw, -1).sum(-1)
            mu_pos = int_dos_tot.searchsorted(totdens, 'left')

            if mu_pos == nw:
                warn("ill-defined LDA chemical potential - outside LDA DOS")
                self.mu = self.w[-1]
            elif mu_pos == 0:
                warn("ill-defined LDA chemical potential - outside LDA DOS")
                self.mu = self.w[0]
            else:
                # try to model a possible gap
                self.mu = self.w[mu_pos-1:mu_pos+1].mean()
        else:
            self.mu = mu
            mu_pos = self.w.searchsorted(mu, 'left')

        # compute also densities
        if mu_pos == nw:
            self.densities = self.int_dos[-1]
        elif mu_pos == 0:
            self.densities = self.int_dos[0]
        else:
            self.densities = self.int_dos[mu_pos-1:mu_pos+1].mean(0)

    def gloc(self, iw, mu, siw, siw_gw=None):
        orbspin.warn_offdiagonal(siw)
        siw = orbspin.extract_diagonal(siw)
        zeta = (1j*iw + mu)[:, None, None] - siw

        # compute  integral de N(e)/(zeta - e)
        glociw = self.integrator(self.dos.transpose(1, 2, 0)[None,:,:,:]/
                                 (zeta[:,:,:,None] - self.w[None,None,None,:]))
        glociw = orbspin.promote_diagonal(glociw)
        return glociw


class NanoLattice(KspaceHamiltonian):

#GS: attempt to define a set_leadsiw method inside the KspaceHamiltonian class
    def get_leadsiw(self, leadsw, w_hyb, nleads, beta, niw):
        norb = self.norbitals
        nsp  = self.nspins
        if nleads:
            leadsiw = np.zeros((niw, nleads, norb, nsp, norb, nsp), dtype=complex)
            #wre2mat wants the total # of Matsubara, positive + negative
            # we pass the hybridization spectral density to the transformation
            leadtrf = tr.transform(tr.wre2mat(beta, w_hyb, 'fermi', niw), -leadsw.imag/np.pi)
            leadtrf = leadtrf.transpose(5, 0, 1, 2, 3, 4)   # iw_n, lead, band, spin, band, spin
            if not self.spinorbit:
#GS:            self.leadsiw = np.zeros_like(leadsiw)
                # in case of spin-independent Hamiltonians, we copy twice to the diagonal
                leadsiw[:, :, :, :1, :, :1] = leadtrf
                leadsiw[:, :, :, 1:2, :, 1:2] = leadtrf
            else:
                leadsiw=leadtrf
#           leadsiw.real += self.hk[0,...].real   #adding a (band, spin, band, spin) array to a (w, lead, band, spin, band, spin) one
#GS:!!!!!!  leadsiw.real += self.hk.mean(0).real  #here I am less sure about the right indices, but it seems to work. Ask Markus!
            leadsmom = tr.transform(tr.wre_int(beta, w_hyb), -leadsw.imag/np.pi)
            print "leadsmom = ", leadsmom.shape
            leadsmom = leadsmom.transpose(5, 0, 1, 2, 3, 4)
#           print "leadsmom = ", leadsmom
        else:
            leadsiw = np.zeros((niw, 1, norb, nsp, norb, nsp), dtype=complex)
            leadsmom = np.zeros((1, norb, nsp, norb, nsp))
        return leadsiw, leadsmom


    def __init__(self, hk, beta, niw, leadsw, w_hyb, nleads, check_herm=False,
                 deltino=1e-1):
#GS:def __init__(self, hk, leadsw, w_hyb, nleads, beta, check_herm=False):


#GS:    set up of Hk exactly as in KspaceHamiltonian.
#GS:	Moments of H(k) need not to be redifined, right (Angelo)?? Check Eq. B.18 of Angelo's PhD Thesis
#GS: 	In the moments of G0 there is however also a contribution from V (Eq. B.21), which will have to be taken into account at some point
        KspaceHamiltonian.__init__(self, beta, hk, check_herm=True)

        print "norbitals, nspins, nflavours = ", self.norbitals , self.nspins , self.nflavours
        print "total # of Matsubara (positive + negative) = ", niw

        # store the Delta(w) for the computation of the DOS
        self.leadsw = leadsw.transpose(5, 0, 1, 2, 3, 4)
        self.w = w_hyb
        self.deltino = deltino

#GS:	add also a check that # of k-points is really 1 (ask Angelo)

#GS:    we now call read_ImHyb as well as the transformation function to get leadsiw. We should then define moments from leadsiw
#GSAV   if nleads:
#GSAV       print "check: after read_ImHyb", leadsw.shape, nleads
#GS:        add a check that the dimension of the lead matrix is compatible with nflavours

#GS:        initialize self.leadsiw --> this is not optimal, as if readleads=False, then you would have no attribut leadsiw for the nano-class
        self.leadsiw, leadsmom = self.get_leadsiw(leadsw, w_hyb, nleads, beta, niw)
        self.iwleads = tr.matfreq(beta, 'fermi', niw)

#GS:        now transform each of these to matsubara and then fill self.leadsiw in
#GS:        leadsiw = tr.transform(tr.wre2mat(beta, w_hyb-hk[0,0,0,0,0].real, 'fermi', 800), leadsw)  # <-- this is needed if ImDelta(w) is symmetric
#GS!!!!!!   leadsiw = tr.transform(tr.wre2mat(beta, w_hyb, 'fermi', niw), leadsw)  #wre2mat wants the total # of Matsubara, positive + negative
#GS:        leadsiw = tr.wre2mat_quad(beta, w_hyb, leadsw, 'fermi', 2) # replace 10 with a parameter and don't forget to add hk (or hk.mean?)

#GS: In order to use it as argument for the function transform, the frequency must be the last argument:
#GS!!!!!!   leadsiw = leadsiw.transpose(5, 0, 1, 2, 3, 4)   # w, lead, band, spin, band, spin
        print "leadsiw -> ", self.leadsiw.shape

#GS:        add hk[0:...] to the real part of leadsiw
#           leadsiw.real += hk[0,...].real   #adding a (band, spin, band, spin) array to a (w, lead, band, spin, band, spin) one
#GS:!!!!!   leadsiw.real += hk.mean(0).real  #here I am less sure about the right indices, but it seems to work. Ask Markus!

#GS:        add a check that the imaginary part of hk is zero.....

#GS:    leadsiw needs to be put in the right form, if spin_orbit=false
#GS:!!!!!   if not self.spinorbit:
#GS:            self.leadsiw = np.zeros_like(leadsiw)
                # in case of spin-independent Hamiltonians, we copy twice to the diagonal
#GS:!!!!!       self.leadsiw[:, :, :, :1, :, :1] = leadsiw
#GS:!!!!!       self.leadsiw[:, :, :, 1:2, :, 1:2] = leadsiw
#GS:!!!!!   else:
#GS:!!!!!       self.leadsiw=leadsiw
#GS:!!!!else:
#GS:!!!!!   self.leadsiw = 0

        # set "non-interacting moments": those from H(k) *and* from leadsiw
        self.hloc = self.hk.mean(0)
        self.hmom2 = np.einsum("kasbt,kbtcu->ascu", self.hk,
                               self.hk)/self.nkpoints - self.hloc**2 + leadsmom.sum(1)
        
        print "LMOM"
        print leadsmom.sum(1)

    def compute_dos(self, fixmu, mu=None, totdens=None):
#AV:    as the gloc method is written for the Matsubara representation,
#       we pass an ad-hoc frequency to get the retarded Green's function
        fake_iw = -1j * self.w + self.deltino
#       the effective leads' self-energy is given by the sum over nleads
#       we also need to subtract the crystal fields, already included in Delta
#       at the moment we subtract H(k) or should we rather subtract just the diagonal elements epsilon_d?
#       this depends on the definition of the Delta obtained from the DFT calculation
        sw_leads = self.leadsw.sum(1) - self.hk.mean(0)
#       get the retarded Green's function and extract the dos from the diagonal elements
        gloc_w = KspaceHamiltonian.gloc(self, fake_iw, 0, sw_leads)
        self.dos = -1/np.pi * orbspin.extract_diagonal(gloc_w.imag)

        siw_leads = self.leadsiw.sum(1)
        def get_totdens_for_mu(mu):
            gloc_iw = KspaceHamiltonian.gloc(self, self.iwleads, mu, siw_leads)
            #np.savetxt("hk_test", self.hk[0,:,0,:,0])
            gloc_model = self.eye * 1/(1j*self.iwleads + mu)[:,None,None,None,None]
            dens_model = fermi(-mu, self.beta) * self.eye
            dens_matrix = (gloc_iw - gloc_model).sum(0)/self.beta + dens_model
            self.densities = orbspin.extract_diagonal(dens_matrix)
            return self.densities.sum()

        if fixmu:
            self.mu = scipy.optimize.brentq(
                            lambda mu: get_totdens_for_mu(mu) - totdens,
                            self.w.min(), self.w.max())
        else:
            self.mu = mu
            get_totdens_for_mu(mu)


    def gloc(self, iw, mu, siw, siw_gw=None):
        arg = siw + self._get_leadsiw(iw)
        return KspaceHamiltonian.gloc(self, iw, mu, arg)


    def trace(self, iw, siw, smom, strategy='local'):
        # FIXME (@Giorgio) does the moment also need to be corrected?!
        # Giorgio: see Angelo's PhD thesis Eq. B23: the V^2-like moment 
        # of leadsiw enters at the same level of the 1/iw_n moment of Sigma. 
        # We should include it if we go one step "higher" in the moments.
        # Actually, why don't we do it? ;-)
        arg = siw + self._get_leadsiw(iw)
        return KspaceHamiltonian.trace(self, iw, arg, smom, strategy)


    def _get_leadsiw(self, iw):
        return self.leadsiw.sum(1)[rmap(iw, self.iwleads)]
