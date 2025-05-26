"""Package abstracting impurity problems and the way of solving it.

Provides three basic classes for the dealing with the impurity picture of the
DMFT self-consistency loop. Their relation can be visualised as follows:

                           +----------------+
        ------------------>| ImpuritySolver |----------------->
          ImpurityProblem  +----------------+  ImpurityResult
"""
from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
import time
from warnings import warn
import numpy as np

from w2dyn.auxiliaries.statistics import (DistributedSample,
                                          DistributedJackknife)
import w2dyn.dmft.orbspin as orbspin


def lattice_convention(qtty):
    """Function to be used for tranposing and reshaping three-dimensional
    band/spin-diagonal arrays with (band, spin) as first two
    dimensions as obtained from the solver for some quantities into
    full five-dimensional band/spin-matrices with (band, spin, band,
    spin) as last four dimensions as used in the DMFT code.
    """
    return orbspin.promote_diagonal(qtty.transpose(2, 0, 1))


class DummyMpiCommunicator:
    def __init__(self): pass
    def Get_size(self): return 1
    def Get_rank(self): return 0
    def allreduce(self, rank_result): return rank_result
    def allgather(self, rank_result): return rank_result


class ImpurityProblem:
    """Specification of a single-impurity Anderson problem.

    Single-impurity Anderson problem, specified by the following attributes:
      - `beta`:        Inverse temperature
      - `g0inviw`:     Inverse of the Weiss field `G_0^{-1}(iw)` in Matsubara
      - `fiw`:         Hybridisation function in Matsubara
      - `ftau`:        Hybridisation function in imaginary time
      - `muimp`:       single-particle impurity Hamiltonian as correction to mu
      - `interaction`: d-d electron interaction
      - `paramag`:     true if system should be forced para-magnetic
      - `symmetries`:  symmetries of the system
      - `screening`:   screened interaction
    """
    def __init__(self, beta, g0inviw, fiw, fmom, ftau, muimp, interaction,
                 phonon_g=None, phonon_omega0=None, symmetry_moves=(),
                 paramag=False, screening=None):
        """Creates new impurity problem from all parameters"""
        # assign to problem
        self.beta = float(beta)
        self.g0inviw = np.asarray(g0inviw)
        self.fiw = np.asarray(fiw)
        self.fmom = np.asarray(fmom)
        self.ftau = np.asarray(ftau)
        self.muimp = np.asarray(muimp)
        self.interaction = interaction
        self.paramag = bool(paramag)
        self.symmetry_moves = symmetry_moves

        # extract geometries and discretisations
        self.niw, self.norbitals, self.nspins = fiw.shape[:3]
        self.nflavours = self.norbitals * self.nspins
        self.nftau = ftau.shape[0]

        # phonon stuff
        self.use_phonons = phonon_g is not None or phonon_omega0 is not None
        if not self.use_phonons:
            phonon_g = phonon_omega0 = ()

        self.phonon_g = np.asarray(phonon_g)
        self.phonon_omega0 = np.asarray(phonon_omega0)

        # screened interaction
        self.use_screening = screening is not None
        if not self.use_screening:
            screening = np.zeros((self.nftau, self.norbitals, self.nspins,
                                  self.norbitals, self.nspins), np.cdouble)
        self.screening = screening
        self.verify()

    def verify(self):
        """Verifies consistency of the impurity parameters"""
        def assert_dims(qtty_name, *req_shape):
            shape = getattr(self, qtty_name).shape
            if shape != req_shape: raise RuntimeError(
                "quantity %s of impurity problem has wrong dimensions.\n"
                "expected: %s, actual: %s" % (qtty_name, shape, req_shape))
        # shorthands
        norb, nsp = self.norbitals, self.nspins
        # hard checks
        if self.beta <= 0: raise RuntimeError("beta is smaller than 0")
        assert_dims("g0inviw", self.niw, norb, nsp, norb, nsp)
        assert_dims("fiw", self.niw, norb, nsp, norb, nsp)
        assert_dims("ftau", self.nftau, norb, nsp, norb, nsp)
        assert_dims("muimp", norb, nsp, norb, nsp)
        assert_dims("screening", self.nftau, norb, nsp, norb, nsp)
        # TODO: "soft" checks like if diag(ftau) >= 0 etc.
        pass

class ImpurityResult:
    """Result of an Anderson impurity problem"""

    def __init__(self, problem, giw, gsigmaiw=None, smom=None,
                 g_diagonal_only=False, **other):

        #if g_diagonal_only:
        #    if not orbspin.is_diagonal(giw, 0):
        #        raise ValueError("Giw contains offdiagonals even though"
        #                         "the result is marked as diagonal")
        #    if not (gsigmaiw is None or orbspin.is_diagonal(gsigmaiw, 0)):
        #        raise ValueError("GSigmaiw contains offdiagonals even though"
        #                         "the result is marked as diagonal")

        self.problem = problem
        self.giw = giw
        self.gsigmaiw = gsigmaiw
        self.smom = smom
        self.g_diagonal_only = g_diagonal_only
        self.other = other

    def postprocessing(self, siw_method="dyson", smom_method="estimate"):
        """Set self-energy `self.siw` and its moment `self.smom`.

        Parameters
        ----------
        siw_method : {'dyson', 'improved_worm', 'symmetric_improved_worm'}, optional
            whether to calculate self-energy from Dyson equation,
            from improved estimators: `self.gsigmaiw/self.giw`,
            or from symmetric improved estimators: `(\Sigma_H + \theta) / (1 + \mathcal{G}(\Sigma_H + \theta))`
        smom_method : {'extract', 'estimate'}, optional
            whether to calculate the moment from self-energy data ('extract')
            or from calculated densities ('estimate').

        """
        # This function is disjoint from the initialization, because it is a
        # "derived" quantity that should be computed after averaging over bins
        self.siw = self.calculate_siw(siw_method)
        self.smom = self.calculate_smom(smom_method)

    def calculate_siw(self, method):
        """Calculate the self-energy from impurity data.

        Parameters
        ----------
        method : {'dyson', 'improved_worm', 'symmetric_improved_worm'}, optional
            whether to calculate self-energy from Dyson equation or
            from improved estimators: `self.gsigmaiw/self.giw`,
            or from symmetric improved estimators: `(\Sigma_H + \theta) / (1 + \mathcal{G}(\Sigma_H + \theta))`

        Returns
        -------
        siw : complex ndarray
            The self-energy.

        Raises
        ------
        NotImplementedError
            If a different method is given.

        """
        method = method.lower()  # ignore case
        if method == "dyson":
            if self.g_diagonal_only:
                # If the solver is only capable of supplying the spin/orbital-
                # diagonal terms of the Green's function, we need to make sure
                # that also only the diagonal terms of G_0 are taken into
                # account in the Dyson equation, otherwise Sigma erroneously
                # "counteracts" these terms.  As an added performance benefit,
                # we can use 1/giw rather than the inversion in this case.
                get_siw = lambda giw: orbspin.promote_diagonal(
                            orbspin.extract_diagonal(self.problem.g0inviw)
                            - 1/orbspin.extract_diagonal(giw))
            else:
                get_siw = lambda giw: self.problem.g0inviw - orbspin.invert(giw)

        elif method == "improved":
            if self.gsigmaiw is None:
                warn("Cannot compute improved estimators - GSigmaiw missing\n"
                     "Falling back to Dyson equation", UserWarning, 2)
                return self.calculate_siw("dyson")
            raise NotImplementedError()  # FIXME

        elif method == 'improved_worm':
            if self.gsigmaiw is None:
                warn("Cannot compute improved estimators - GSigmaiw missing\n"
                     "Falling back to Dyson equation", UserWarning, 2)
                return self.calculate_siw("dyson")

            if self.g_diagonal_only:
                get_siw = lambda giw: orbspin.promote_diagonal(
                    orbspin.extract_diagonal(self.problem.g0inviw)
                    - 1./orbspin.extract_diagonal(giw))
            else:
                self.siw = self.problem.g0inviw - orbspin.invert(self.giw)
                raise NotImplementedError("Offdiagonal worm improved \n"
                                          "estimators not implemented")

        elif method == 'symmetric_improved_worm':

            if self.g_diagonal_only:
                #pseudo_se = orbspin.promote_diagonal(self.qqiw) + self.smom.value
                #se = pseudo_se / (1.
                #        + orbspin.promote_diagonal(orbspin.invert(self.problem.g0inviw)) * pseudo_se)
                #return se
                get_siw = lambda giw: orbspin.promote_diagonal(
                    orbspin.extract_diagonal(self.problem.g0inviw)
                    - 1./orbspin.extract_diagonal(giw))
            else:
                self.siw = self.problem.g0inviw - orbspin.invert(self.giw)
                raise NotImplementedError("Offdiagonal worm improved \n"
                                          "estimators not implemented")
        else:
            raise ValueError("unknown siw_method: %s" % method)

        giw_jk_input = DistributedJackknife(self.giw)
        return giw_jk_input.transform(get_siw)

    def calculate_smom(self, method):
        """Calculate the moment of the self-energy.

        Parameters
        ----------
        method : {'extract', 'estimate'}, optional
            whether to calculate the moment from self-energy data ('extract')
            or from calculated densities ('estimate').

        Returns
        -------
        smom : float ndarray
            The moment of the self-energy.

        Raises
        ------
        NotImplementedError
            If a different method is given.

        """
        method = method.lower()
        if method == "extract":
            warn("Extracting moment of self-energy from the data",
                 UserWarning, 2)
            return DistributedSample(
                self.siw.local[:, :2].real.copy(),
                self.siw.mpi_comm
            )
        if method == "estimate":
            if self.smom is None:
                warn("Cannot compute smom estimate - rho1 and rho2 missing\n"
                     "Falling back to extraction", UserWarning, 2)
                return self.calculate_smom("extract")
            return self.smom
        raise ValueError("unknown smom_method: %s" % method)


class StatisticalImpurityResult(ImpurityResult):
    pass


class ImpuritySolver(object):
    """Solver that takes an ImpurityProblem and yields an ImpurityResult."""
    def __init__(self, config, **params):
        """Initialises and configures the solver"""
        self.config = config
        self.abort = False

    def set_problem(self, problem, iter_type):
        """Sets the problem the impurity solver should solve"""
        self.problem = problem

    def solve(self, mpi_rank=0):
        """Instructs the solver to solve the impurity problem set and return the
           result as an instance of ImpurityResult."""
        raise NotImplementedError()
