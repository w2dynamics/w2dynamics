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
import scipy.optimize as opt
import sys

from w2dyn.auxiliaries.CTQMC import mctqmc as ctqmc
import w2dyn.auxiliaries.transform as tf
import w2dyn.auxiliaries.postprocessing as postproc
import w2dyn.auxiliaries.compound_index as ci
import w2dyn.auxiliaries.statistics as statistics
from w2dyn.auxiliaries.statistics import (DistributedSample,
                                          DistributedJackknife)
import w2dyn.dmft.dynamicalU as dynamicalU
# import copy

import w2dyn.dmft.orbspin as orbspin

# deprecated function time.clock removed in python >= 3.8
if not hasattr(time, 'clock'):
    time.clock = time.perf_counter


def lattice_convention(qtty):
    """Function to be used for tranposing and reshaping three-dimensional
    band/spin-diagonal arrays with (band, spin) as first two
    dimensions as obtained from the solver for some quantities into
    full five-dimensional band/spin-matrices with (band, spin, band,
    spin) as last four dimensions as used in the DMFT code.
    """
    return orbspin.promote_diagonal(qtty.transpose(2, 0, 1))


class CtHybConfig:
    def __init__(self, taus, orbs, spins, cas, hashybs,
                 outer_sst, outer_state):
        self.taus = taus
        self.orbs = orbs
        self.spins = spins
        self.cas = cas
        self.hashybs = hashybs
        self.outer_sst = outer_sst
        self.outer_state = outer_state

    def __str__(self):
        return ("Monte Carlo configuration:\n" +
                "outer sst: {}\n".format(self.outer_sst) +
                "outer state: {}\n".format(self.outer_state) +
                "operators: {}".format(np.vstack((self.taus,
                                                  self.orbs,
                                                  self.spins,
                                                  self.cas,
                                                  self.hashybs))))

    @classmethod
    def get_from_ctqmc(cls, solver):
        noper, outer_sst, outer_state = solver.get_mc_config_scalars()
        if noper > 0:
            return cls(*solver.get_mc_config(noper))
        else:
            return cls(np.array((), dtype=np.double),
                       np.array((), dtype=np.int),
                       np.array((), dtype=np.int),
                       np.array((), dtype=np.int),
                       np.array((), dtype=np.int),
                       outer_sst, outer_state)

    def set_to_ctqmc(self, solver):
        if self.taus.size > 0:
            solver.set_mc_config(self.taus.size,
                                 self.taus, self.orbs, self.spins, self.cas,
                                 self.hashybs,
                                 self.outer_sst, self.outer_state)
        else:
            solver.set_empty_mc_config(self.outer_sst, self.outer_state)

    @classmethod
    def load_from_file(cls, file):
        with np.load(file) as f:
            return cls(f["taus"],
                       f["orbs"],
                       f["spins"],
                       f["cas"],
                       f["hashybs"],
                       f["outer_sst"][0],
                       f["outer_state"][0])

    def save_to_file(self, file):
        np.savez(file,
                 taus=self.taus,
                 orbs=self.orbs,
                 spins=self.spins,
                 cas=self.cas,
                 hashybs=self.hashybs,
                 outer_sst=np.array([self.outer_sst]),
                 outer_state=np.array([self.outer_state]))


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
                                  self.norbitals, self.nspins), np.complex)
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
            return statistics.DistributedSample(
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


class CtHybSolver(ImpuritySolver):
    """Wrapper for Fortran CT-HYB(M-V) solver for diagonal hybridisation.

    This class wraps around the Fortran module `auxiliaries.CTQMC.mctqmc`. This
    module exposes a CT-QMC solver in the hybridisation expansion, implemented
    both using the Krylov subspace and matrix-vector technique.

    Note that by being a Fortran module, it effectively makes this class a
    singleton, so you cannot run multiple instances of CtHybSolver on a single
    core in parallel.
    """
    def __init__(self, config, seed, Uw, Uw_Mat, epsn, interactive=True, mpi_comm=None):
        super(CtHybSolver, self).__init__(config)
        self.mpi_comm = mpi_comm
        self.mpi_rank = 0
        if mpi_comm is not None:
            self.mpi_rank = mpi_comm.Get_rank()

        self.seed = seed + self.mpi_rank
        self.param_string = self.config_to_pstring(config["QMC"])
        self.g_diagonal_only = (config["QMC"]["offdiag"] == 0)
        self.Uw = Uw
        self.Uw_Mat = Uw_Mat
        self.epsn = epsn
        ctqmc.init_rnd(seed)
        ctqmc.fancy_prog = interactive

        # START DYNAMICAL U
        if self.Uw == 1:
          Nd = config["Atoms"]["1"]["Nd"]
          beta = config["General"]["beta"]
          Nftau = config["QMC"]["Nftau"]
          umatrixdummy = np.zeros((Nd*2,Nd*2,Nd*2,Nd*2))
          self._Retarded2Shifts = dynamicalU.Replace_Retarded_Interaction_by_a_Shift_of_instantaneous_Potentials(beta, umatrixdummy, self.Uw_Mat)
          screeningdummy = np.zeros((Nd,2,Nd,2,Nftau))
          self.screening = self._Retarded2Shifts.create_tau_dependent_UV_screening_function(Nftau, beta, screeningdummy)
          print("         ****************************")
          print("         ******* U(w) Message *******")
          print("         ****************************")
          print(" ")

    def set_problem(self, problem, compute_fourpoint=0):
        # problem
        self.problem = problem
        problem_params = self.config_to_pstring({
            "Hamiltonian": problem.interaction.name,
            "beta": problem.beta,
            "Nd": problem.norbitals,
            "ParaMag": int(problem.paramag),
            "QuantumNumbers": " ".join(problem.interaction.quantum_numbers),
            "Phonon": int(problem.use_phonons),
            "g_phonon": " ".join(problem.phonon_g),
            "omega0": " ".join(problem.phonon_omega0),
            "Screening": int(problem.use_screening),
            "Uw": self.Uw,
            })
        symm_params = self.symm_to_pstring(problem.symmetry_moves)
        all_params = self.param_string + problem_params + symm_params
        ctqmc.init_paras(all_params)
        ctqmc.fourpnt = compute_fourpoint

        # prepare parameters
        if self.config["QMC"]["offdiag"] == 0:
           self.ftau=orbspin.promote_diagonal(orbspin.extract_diagonal(self.problem.ftau))

        self.ftau = self.problem.ftau.transpose(1,2,3,4,0)

        # CT-HYB uses a different convention for muimp
        #self.muimp = -self.prepare_input(self.problem.muimp)
        self.muimp = -self.problem.muimp
        self.umatrix = self.problem.interaction.u_matrix.reshape(
                             self.problem.nflavours, self.problem.nflavours,
                             self.problem.nflavours, self.problem.nflavours)
        if self.Uw == 0:
            self.screening = self.problem.screening.transpose(1, 2, 3, 4, 0).real

        # START DYNAMICAL U
        if self.Uw == 1:
            #print " "
            #print "         ****************************"
            #print "         ****** U(w) is active ******"
            #print "         ****************************"
            ### ===> UPDATE This is now done in the init above.
            ### ===> UPDATE _Retarded2Shifts = dynamicalU.Replace_Retarded_Interaction_by_a_Shift_of_instantaneous_Potentials(problem.beta, self.umatrix, self.Uw_Mat)

            ### ===> UPDATE The umatrix shifting is now performed in config.py.
            ### ===> UPDATE The umatrix shifting is now performed in config.py.
            ### ===> # This should be executed for all atoms in the first iteration
            ### ===> try:
            ### ===>    self.iteration_ID    # Not defined for first DMFT iteration
            ### ===> except:
            ### ===>    self.umatrix_original = copy.copy(self.umatrix)
            ### ===>    self.iteration_ID = 1
            ### ===> if (self.umatrix==self.umatrix_original).all():
            ### ===>    self.umatrix = _Retarded2Shifts.shift_instantaneous_densitydensity_potentials(self.umatrix)
            ### ===> UPDATE The umatrix shifting is now performed in config.py.
            ### ===> UPDATE The umatrix shifting is now performed in config.py.

            ### ===> UPDATE This is now done in the init above.
            ### ===> UPDATE self.screening = _Retarded2Shifts.create_tau_dependent_UV_screening_function(problem.nftau, problem.beta, self.screening)

            # Perform chemical potential shift in each DMFT iteration
            if problem.norbitals == 1 and self.epsn == 0:
                self.muimp = self._Retarded2Shifts.shift_chemical_potential_for_single_band_case(self.muimp, self.umatrix)
            #else:
                #print "  ==> U(w): No shift of mu for mu-search"
                #if problem.norbitals > 1 and self.epsn==0:
                #print "===>>> U(w): mu-Search (epsn > 0) is recommended for the chosen problem."
        # END DYNAMICAL U


    def solve(self, iter_no, mc_config_inout=None):
        sys.stdout.flush()
        sys.stderr.flush()
        if self.mpi_comm is not None:
            self.mpi_comm.Barrier()

        # Call CT-QMC solver
        ctqmc.simid = self.mpi_rank
        time_qmc = time.clock()
        ctqmc.init_solver(self.umatrix, self.ftau, self.muimp, self.screening, iter_no+1)

        firstrun = True
        if type(mc_config_inout) is list and len(mc_config_inout) > 0:
            firstrun = False
            mc_config_inout.pop().set_to_ctqmc(ctqmc)

        if self.config["QMC"]["segment"] == 0:
            ctqmc.set_taudiffmax(self.config["QMC"]["TaudiffMax"]
                                 if self.config["QMC"]["TaudiffMax"] > 0.0
                                 else self.calibrate_taudiff_max())

        ctqmc.ctqmc_warmup(self.config["QMC"]["Nwarmups"]
                           if firstrun or self.config["QMC"]["Nwarmups2Plus"] < 0
                           else self.config["QMC"]["Nwarmups2Plus"])
        ctqmc.ctqmc_measure(1,1)
        time_qmc = time.clock() - time_qmc
        if type(mc_config_inout) is list:
            mc_config_inout.append(CtHybConfig.get_from_ctqmc(ctqmc))

        if ctqmc.aborted == 2:  # no data, kill run
            raise KeyboardInterrupt("CTQMC was interrupted by signal")
        elif ctqmc.aborted == 1:
            self.abort = True
            warn("CT-QMC was aborted before all measurements were conducted.\n"
                 "Expect larger errorbars.", UserWarning, 2)

        # Get result set from the module variables.
        result = self.get_resultset(self.config["QMC"])
        result["time-qmc"] = time_qmc

        # Extract Green's function to use for Dyson equation and transform it
        # from (band, spin, iw) to (iw, band, spin, band, spin) array for use
        # in DMFT self-consistency cycle.
        giw = self.get_giw(self.config, self.problem, result)
        if self.config["QMC"]["offdiag"] == 0:
            giw = lattice_convention(giw)
        else:
            giw = np.ascontiguousarray(giw.transpose(4, 0, 1, 2, 3))
        giw = DistributedSample([giw], self.mpi_comm)


        # Extract improved estimators
        gsigmaiw = self.get_gsigmaiw(self.config, self.problem, result)
        if gsigmaiw is not None:
            gsigmaiw = lattice_convention(gsigmaiw)
            gsigmaiw = DistributedSample([gsigmaiw], self.mpi_comm)

        # Extract moments of the self-energy from 1- and 2-particle reduced
        # density matrix
        try:
            rho1 = result["rho1"]
            rho2 = result["rho2"]
        except KeyError:
            smom = None
        else:
            umatrix_ctr = self.problem.interaction.u_matrix.reshape(
                                                [self.problem.nflavours] * 4)
            smom = postproc.get_siw_mom(umatrix_ctr, rho1, rho2)
            smom = lattice_convention(smom)
            smom = DistributedSample([smom], self.mpi_comm)

        result = {key: DistributedSample([value], self.mpi_comm)
                  for (key, value) in result.items()}

        # Construct result object from result set and collect it from all nodes
        result = StatisticalImpurityResult(self.problem, giw, gsigmaiw, smom,
                                           self.g_diagonal_only, **result)

        # Cleanup
        ctqmc.dest_ctqmc()
        ctqmc.dest_paras()

        return result


    def solve_component(self,iter_no,isector,icomponent,mc_config_inout=[]):
        """ This solver returns the worm estimator for a given component. 
           This results in several advantages:

              1. the balancing between worm space and partition function space
                 (i.e. eta-search) is done for specified component, and not in
                 an integrated way

              2. memory requirements are lowered for 2P quantities
                 (averaging can be done in python)

              3. (following 2.) no need for task-local files, postprocessing
                 formula to calculate H2->chi_conn can be implemented directly

           Convention of sampling spaces (Sector):

           Sector=1: Z (partition function, always sampled)
           Sector=2: G1 (one particle Green's function)
           Sector=3: GSigma (one particle improved estimator)
           Sector=4: G2 (two particle Green's function)
           Sector=5: H2 (two particle improved estimator)
           Sector=6: P2PH (single-freqeuncy 2P-GF in PH channel)
           Sector=7: P2PP (single-freqeuncy 2P-GF in PP channel)
           Sector=8: P3PH (two-freqeuncy 2P-GF in PH channel)
           Sector=9: P3PP (two-freqeuncy 2P-GF in PP channel)
           Sector=10: QQ (one-particle symmetric IE)
           Sector=11: QQQQ (two-particle symmetric IE)
           Sector=12: NQQdag
           Sector=13: QQdd
           Sector=14: Ucaca (P2 with U matrices)
           Sector=15: Uccaa (P2pp with U matrices)
           Sector=16: QUDdag
        """
        sys.stdout.flush()
        sys.stderr.flush()
        if self.mpi_comm is not None:
            self.mpi_comm.Barrier()

        # Call CT-QMC solver
        ctqmc.simid = self.mpi_rank
        time_qmc = time.clock()
        ctqmc.init_solver(self.umatrix, self.ftau, self.muimp, self.screening, iter_no+1)

        if type(mc_config_inout) is list and len(mc_config_inout) > 0:
            mc_config_inout[0].set_to_ctqmc(ctqmc)

        ctqmc.set_taudiffmax(self.config["QMC"]["TaudiffMax"]
                             if self.config["QMC"]["TaudiffMax"] > 0.0
                             else self.calibrate_taudiff_max())

        # first warm-up
        if not mc_config_inout:
           ctqmc.ctqmc_warmup(self.config["QMC"]["Nwarmups"])
           mc_config_inout.append(CtHybConfig.get_from_ctqmc(ctqmc))

        if self.config["QMC"]["Nwarmups2Plus"] < 0:
           warn("Nwarmups2Plus not set; assuming 0.1*Nwarmups for worm")
           WormWarmups=0.1*self.config["QMC"]["Nwarmups"]
        else:
           WormWarmups=self.config["QMC"]["Nwarmups2Plus"]

        if int(self.config['QMC']['WormSearchEta']) == 1:
            # Patrik's eta search: bisection (?)
            ctqmc.wormeta[:] = np.float(self.config['QMC']['WormEta'])
            max_iter = 50
            for i_try in range(max_iter):
                steps_worm = np.array(0, dtype=np.int)
                steps_z = np.array(0, dtype=np.int)
                ctqmc.count_steps(isector, icomponent, WormWarmups, steps_worm, steps_z)
                if not (steps_worm + steps_z == WormWarmups):
                    raise ValueError('{} steps in worm, {} steps in Z, '
                                     'some steps missing!'.format(steps_worm, steps_z))
                res = abs(float(steps_worm - steps_z)) / float(steps_worm + steps_z)
                print('At eta={}: {} steps in worm, {} steps in Z -> {}'.format(ctqmc.wormeta[0], steps_worm, steps_z, res))
                if res > 0.1:
                    if i_try == max_iter - 1 or steps_worm == 0:
                        print('No suitable value of eta found, setting it to 0.')
                        ctqmc.wormeta[:] = 0.
                        break
                    else:
                        ctqmc.wormeta[:] = ctqmc.wormeta[:] * float(steps_z) / float(steps_worm)
                else:
                    break
            eta_stat = DistributedSample([ctqmc.wormeta[0]], self.mpi_comm)
            eta_mean, eta_err = eta_stat.mean(), eta_stat.stderr()
            print('final eta on core: {}'.format(ctqmc.wormeta[0]))
            print('final eta={}, stdev={}'.format(eta_mean, eta_err))

        elif int(self.config['QMC']['WormSearchEta']) == 2:
            # Josef's eta search: brentq root finding
            def eta_root(eta):
                '''
                Eta is determined by the root of this function. 
                The ctqmc module function count_steps is called 
                to get the number of steps in Z and worm space
                for a given eta. 
                '''
                steps_worm = np.array(0, dtype=np.int)
                steps_z = np.array(0, dtype=np.int)
                ctqmc.wormeta[:] = eta # TODO: is it impossible to set a single array element?
                ctqmc.count_steps(isector, icomponent, WormWarmups, steps_worm, steps_z)
                if not (steps_worm + steps_z == WormWarmups):
                    raise ValueError('{} steps in worm, {} steps in Z, '
                                     'some steps missing!'.format(steps_worm, steps_z))
                res = float(steps_worm - steps_z) / float(steps_worm + steps_z)
                print('At eta={}: {} steps in worm, {} steps in Z -> {}'.format(eta, steps_worm, steps_z, res))
                return res

            bracket_min = 1e-4
            bracket_max = 1e2
            try:
                eta = opt.brentq(eta_root, bracket_min, bracket_max, xtol=0.001, rtol=0.001, maxiter=50, disp=True)
            except ValueError:
                try:  # for practical purposes, enlarging once should be sufficient.
                    print('Trying larger bracket for eta-search.')
                    eta = opt.brentq(eta_root, 0.01 * bracket_min, 100. * bracket_max, xtol=0.001, rtol=0.001, maxiter=50, disp=True)
                except ValueError:
                    print('No suitable value of eta found, setting it to 0.')
                    eta = 0.

            eta_stat = DistributedSample([eta], self.mpi_comm)
            eta_mean, eta_err = eta_stat.mean(), eta_stat.stderr()
            print('final eta on core: {}'.format(eta))
            print('final eta={}, stdev={}'.format(eta_mean, eta_err))

        else:
            eta_mean, eta_err = float(self.config['QMC']['WormEta']), 0.

        ctqmc.wormeta[:] = eta_mean # TODO: is it impossible to set a single array element?

        ctqmc.ctqmc_worm_warmup(WormWarmups,isector,icomponent)
        ctqmc.ctqmc_measure(isector,icomponent)

        time_qmc = time.clock() - time_qmc

        if ctqmc.aborted == 2:  # no data, kill run
           raise KeyboardInterrupt("CTQMC was interrupted by signal")
        elif ctqmc.aborted == 1:
           self.abort = True
           warn("CT-QMC was aborted before all measurements were conducted.\n"
           "Expect larger errorbars.", UserWarning, 2)


        # set dummy Green's function
        giw = np.zeros_like(self.problem.g0inviw)

        # Get result for the component specified - normalize to component
        result = {}

        if isector == 4:
           if self.config['QMC']['G4ph'] == 1:
               result["g4iw-worm/" + str(icomponent).zfill(5)] = \
                    ctqmc.g4iw_worm.copy()/(self.problem.nflavours**4)
           if self.config['QMC']['G4pp'] == 1:
               result["g4iwpp-worm/" + str(icomponent).zfill(5)] = \
                    ctqmc.g4iwpp_worm.copy()/(self.problem.nflavours**4)

        if isector == 5:
           result["h4iw-worm/" + str(icomponent).zfill(5)] = \
                ctqmc.h4iw_worm.copy()/(self.problem.nflavours**4)

        if isector == 6 and self.config["QMC"]["WormMeasP2iwPH"] == 1:
           result["p2iw-worm/" + str(icomponent).zfill(5)] = \
                ctqmc.p2iw_worm.copy()/(self.problem.nflavours**4)

        if isector == 6 and self.config["QMC"]["WormMeasP2tauPH"] == 1:
           result["p2tau-worm/" + str(icomponent).zfill(5)] = \
                ctqmc.p2tau_worm.copy()/(self.problem.nflavours**4)

        if isector == 7 and self.config["QMC"]["WormMeasP2iwPP"] == 1:
           result["p2iwpp-worm/" + str(icomponent).zfill(5)] = \
                ctqmc.p2iwpp_worm.copy()/(self.problem.nflavours**4)

        if isector == 7 and self.config["QMC"]["WormMeasP2tauPP"] == 1:
           result["p2taupp-worm/" + str(icomponent).zfill(5)] = \
                ctqmc.p2taupp_worm.copy()/(self.problem.nflavours**4)

        if isector == 8:
           result["p3iw-worm/" + str(icomponent).zfill(5)] = \
                ctqmc.p3iw_worm.copy()/(self.problem.nflavours**4)

        if isector == 9:
           result["p3iwpp-worm/" + str(icomponent).zfill(5)] = \
                ctqmc.p3iwpp_worm.copy()/(self.problem.nflavours**4)

        if isector == 11:
           result["qqqq-worm/" + str(icomponent).zfill(5)] = \
                ctqmc.qqqq_worm.copy()/(self.problem.nflavours**4)

        if isector == 12:
            result["nqqdag-worm/{:05}".format(icomponent)] = \
                ctqmc.nqqdag_worm.copy()/(self.problem.nflavours**4)

        if isector == 13:
            result["qqdd-worm/{:05}".format(icomponent)] = \
                ctqmc.qqdd_worm.copy()/(self.problem.nflavours**4)

        if isector == 14 and self.config['QMC']['WormMeasUcaca'] == 1:
            result["ucaca-worm/{:05}".format(icomponent)] = \
                ctqmc.ucaca_worm.copy()/(self.problem.nflavours**4)

        if isector == 14 and self.config['QMC']['WormMeasUcacatau'] == 1:
            result["ucacatau-worm/{:05}".format(icomponent)] = \
                ctqmc.ucacatau_worm.copy()/(self.problem.nflavours**4)

        if isector == 15 and self.config['QMC']['WormMeasUccaa'] == 1:
            result["uccaa-worm/{:05}".format(icomponent)] = \
                ctqmc.uccaa_worm.copy()/(self.problem.nflavours**4)

        if isector == 15 and self.config['QMC']['WormMeasUccaatau'] == 1:
            result["uccaatau-worm/{:05}".format(icomponent)] = \
                ctqmc.uccaatau_worm.copy()/(self.problem.nflavours**4)

        if isector == 16 and self.config['QMC']['WormMeasQUDdag'] == 1:
            result['quddag-worm/{:05}'.format(icomponent)] = \
                    ctqmc.quddag_worm.copy()/(self.problem.nflavours**4)

        # Get set of result quantities NOT grouped by component
        result_gen = self.get_resultset(self.config["QMC"])
        result_gen = StatisticalImpurityResult(self.problem, giw, None, None,
                                               self.g_diagonal_only,
                                               **result_gen)

        # Get set of result quantities grouped by component and update
        # with previously extracted component of sector-defining
        # quantity
        result_comp = self.get_resultset_worm(str(icomponent).zfill(5),
                                              self.config['QMC'])
        result_comp.update(result)
        result_comp = StatisticalImpurityResult(self.problem, giw, None, None,
                                                self.g_diagonal_only,
                                                **result_comp)

        # Cleanup
        ctqmc.dest_ctqmc()
        ctqmc.dest_paras()
        return result_gen, result_comp

    def solve_comp_stats(self, iter_no, isector, icomponent,
                         mc_config_inout=[]):
        result_gen, result_comp = self.solve_component(iter_no,
                                                       isector,
                                                       icomponent,
                                                       mc_config_inout)
        result_gen.other = {key: DistributedSample([value], self.mpi_comm)
                            for key, value in result_gen.other.items()}
        result_comp.other = {key: DistributedSample([value], self.mpi_comm)
                             for key, value in result_comp.other.items()}
        return result_gen, result_comp

    def solve_worm(self, iter_no, log_function=None):
        """
        Solve all non-vanishing components of G or GSigma or QQ by component sampling.
        Then convert the component list to standard format (matrix in orbitals).
        From this, compute the impurity Green's function and put it in a StatisticalImpurityResult.
        Furthermore, put all worm quantities to a separate StatisticalImpurityResult.
        Both of them are returned.
        :return: (StatisticalImpurityResult, StatisticalImpurityResult)
        """
        # get the components to be sampled
        conf_comp = self.config["QMC"]["WormComponents"]
        if conf_comp is not None: # use the ones specified in the parameters file
            components = map(int, conf_comp)
        if self.config['QMC']['offdiag'] == 0: # generate all diagonal components
            components = [ci.GFComponent(bandspin=(iflav, iflav), n_ops=2, n_bands=self.problem.nflavours//2).index
                          for iflav in range(self.problem.nflavours)]
        else: # generate the full list
            components = 1 + np.arange(self.problem.nflavours**2) # one-based

        # introduce booleans for the measurement mode
        meas_g = self.config["QMC"]["WormMeasGiw"] != 0
        meas_gs = self.config["QMC"]["WormMeasGSigmaiw"] != 0
        meas_qq = self.config["QMC"]["WormMeasQQ"] != 0

        if meas_g and not meas_gs and not meas_qq:
            sector = 2
            sector_name = 'giw-worm'
        elif not meas_g and meas_gs and not meas_qq:
            sector = 3
            sector_name = 'gsigmaiw-worm'
        elif not meas_g and not meas_gs and meas_qq:
            sector = 10
            sector_name = 'qqiw-worm'
        else:
            raise ValueError("Only one of (WormMeasGiw, "
                             "WormMeasGSigmaiw, WormMeasQQ) can be chosen")

        # perform a solve_component worm sampling run for each component.
        # in this run, the specified worm estimator is measured,
        # and also the densities in Z space, since we need them for moments etc.
        result_list_z = []
        result_list_worm = []
        conf_container = []
        for component in components:
            log_function("sampling component {}".format(component))
            self.set_problem(self.problem, self.config["QMC"]["FourPnt"])
            res_z, res_worm = self.solve_component(iter_no, sector, component, conf_container)
            result_list_z.append((component, res_z))
            result_list_worm.append((component, res_worm))


        # merge worm results from all component runs
        # (they are distinguishable, since their names always contain the component number,
        # e.g. qqiw-worm/00001)
        worm_dict_complete = {}
        for res in result_list_worm:
            worm_dict_complete.update(res[1].other)
        worm_dict_complete = {key: DistributedSample([value], self.mpi_comm)
                              for key, value in worm_dict_complete.items()}

        # Z-sampling quantities, such as occ, are sampled together with each component.
        # To exploit all measurements, we take the average.
        def component_avg(qtty_name, result_list, mpi_comm):
            val_list = []
            for res in result_list:
                qtty = res[1].other[qtty_name]
                val_list.append(DistributedSample([qtty], mpi_comm))
            return statistics.join(*val_list)

        # average Z-sampled quantities from all worm component runs
        result_z_other = {key: component_avg(key, result_list_z, self.mpi_comm)
                          for key in result_list_z[0][1].other.keys()}

        # Extract moments of the self-energy from 1- and 2-particle reduced
        # density matrix
        umatrix_ctr = self.problem.interaction.u_matrix.reshape(
                                            [self.problem.nflavours] * 4)
        smom = postproc.get_siw_mom(np.asarray(umatrix_ctr),
                                    np.asarray(result_z_other['rho1'].mean()),
                                    np.asarray(result_z_other['rho2'].mean()))
        smom = lattice_convention(smom)
        smom = DistributedSample([smom], self.mpi_comm)

        if meas_g:
            # convert component list of g to matrix
            giw = self.get_giw(self.config, self.problem, result_list_worm)
            gsigmaiw = None # not needed here
        if meas_gs:
            # convert component list of gsigma to matrix,
            # then apply IE equation to calculate giw
            gsigmaiw = self.get_gsigmaiw_worm(self.config, self.problem, result_list_worm)
            giw = self.giw_from_gsigmaiw_worm(self.config, self.problem, gsigmaiw)
        if meas_qq:
            # convert component list of qq to matrix,
            # then apply SIE equation to calculate giw
            gsigmaiw = None # not needed here
            qqiw = self.get_qqiw_worm(self.config, self.problem, result_list_worm)
            giw = self.giw_from_qqiw_worm(self.config, self.problem, qqiw, smom.local[0][0])

        giw = DistributedSample([giw], self.mpi_comm)

        return (StatisticalImpurityResult(self.problem, giw, gsigmaiw, smom, self.g_diagonal_only,
                                          **result_z_other),
                StatisticalImpurityResult(self.problem, giw, gsigmaiw, smom, self.g_diagonal_only,
                                          **worm_dict_complete))


        # -------- auxiliary helper methods ------------

    @classmethod
    def config_to_pstring(cls, config):
        return "".join("%s = %s #\n" % (key, value)
                       for key, value in sorted(config.items()))

    @classmethod
    def symm_to_pstring(cls, symmetries):
        s = "NSymMove = %d #\n" % len(symmetries)
        s +=  "".join("SymMove%02d = %s #\n" % (i+1, " ".join(map(str,vals+1)))
                      for i, vals in enumerate(symmetries))
        return s

    @classmethod
    def prepare_input_offdiag(cls, qtty):
        orbspin.warn_offdiagonal(qtty)
        #qtty = orbspin.extract_diagonal(qtty)
        if (np.abs(qtty.imag) > 1e-10).any():
            warn("Quantity has non-vanishing imaginary part", UserWarning, 2)
        qtty = qtty.real
        if qtty.ndim == 5:
            axes_list = [-4, -3, -2, -1] + list(range(qtty.ndim - 4))
            qtty = qtty.transpose(axes_list)
        return qtty

    @classmethod
    def prepare_input(cls, qtty):
        orbspin.warn_offdiagonal(qtty)
        qtty = orbspin.extract_diagonal(qtty)
        if (np.abs(qtty.imag) > 1e-10).any():
            warn("Quantity has non-vanishing imaginary part", UserWarning, 2)
        qtty = qtty.real
        if qtty.ndim > 2:
            axes_list = [-2, -1] + list(range(qtty.ndim - 2))
            qtty = qtty.transpose(axes_list)
        return qtty

    @classmethod
    def get_giw(cls, config, problem, result):
        #typ = "legendre"  # FIXME: QMC config does not contain FTtype
        typ = config["General"]["FTType"]
        if config["QMC"]["offdiag"] != 0:
            if typ == "legendre" or typ == "legendre_full":
                typ = "legendre_full"
            elif typ == "none_worm":
                pass
            else:
                warn("Chosen FTType unsupported for non-diagonal hybridization, "
                     "using legendre_full instead")
                typ = "legendre_full"
        if typ == "none":
            giw = -result["giw-meas"]
        elif typ == "none_worm":
            giw = np.zeros_like(problem.g0inviw).transpose((1,2,3,4,0))
            for icomp, res in result:
                comp_ind = ci.GFComponent(index=icomp, n_ops=2,
                    n_bands=problem.nflavours//problem.nspins).bsbs()
                giw[comp_ind] = res.other['giw-worm/{:05}'.format(icomp)]
            giw = giw.transpose((4,0,1,2,3))
        elif typ == "plain":
            ntau = config["QMC"]["Ntau"]
            tf_matrix = tf.tau2mat(problem.beta, ntau, 'fermi',
                                   problem.niw)
            giw = tf.transform(tf_matrix, result["gtau"], None)
        elif typ == "legendre":
            nleg = config["QMC"]["NLegOrder"]
            tf_matrix = -tf.leg2mat(nleg, problem.niw)
            giw = tf.transform(tf_matrix, result["gleg"], None,
                               onmismatch=tf.Truncate.tofirst, warn=False)
        elif typ == "legendre_full":
            nleg = config["QMC"]["NLegOrder"]
            tf_matrix = -tf.leg2mat(nleg, problem.niw)
            nbands = result["gleg-full"].shape[0]
            giw = np.zeros(shape=(nbands, 2, nbands, 2, problem.niw),
                           dtype=complex)

            for b in range(0,nbands):
               for s in range(0,2):
                  giw[b,s,:,:,:] = tf.transform(tf_matrix, result["gleg-full"][b,s,:,:,:], None,
                                     onmismatch=tf.Truncate.tofirst, warn=False)

            giw += giw.transpose(2, 3, 0, 1, 4)
            giw /= 2.0

        else:
            raise RuntimeError("Invalid transform type: `%s'" % typ)

        return giw   # copy not necessary here.

    @classmethod
    def get_gsigmaiw_worm(cls, config, problem, result):
        gsigma = np.zeros_like(problem.g0inviw).transpose((1,2,3,4,0))
        for icomp, res in result:
            comp_ind = ci.GFComponent(index=icomp, n_ops=2,
                n_bands=problem.nflavours//problem.nspins).bsbs()
            gsigma[comp_ind] = res.other['gsigmaiw-worm/{:05}'.format(icomp)]
        return gsigma.transpose((4,0,1,2,3))

    @classmethod
    def get_qqiw_worm(cls, config, problem, result):
        qq = np.zeros_like(problem.g0inviw).transpose((1, 2, 3, 4, 0))
        for icomp, res in result:
            comp_ind = ci.GFComponent(index=icomp, n_ops=2,
                    n_bands=problem.nflavours // problem.nspins).bsbs()
            qq[comp_ind] = res.other['qqiw-worm/{:05}'.format(icomp)]
        return qq.transpose((4, 0, 1, 2, 3))

    @classmethod
    def get_gsigmaiw(cls, config, problem, result):
        gsigmaiw = None
        # re-enable this when/if Z-space GSigma measurement is fixed
        # if config["QMC"]["MeasGSigmaiw"] == 1:
        #     gsigmaiw = result["gsigmaiw"]
        #     del result["gsigmaiw"] # have to delete it here, since gsigmaiw is passed explicitly to StatisticalImpurityResult
        return gsigmaiw

    @classmethod
    def giw_from_gsigmaiw_worm(cls, config, problem, gsigma):
        """
        Calculate impurity Green's function by improved estimator formula,
        Eq. (9) in PRB 100, 075119 (2019).
        """
        g0iw = orbspin.invert(problem.g0inviw)
        giw = g0iw + orbspin.multiply(g0iw, gsigma)
        return giw

    @classmethod
    def giw_from_qqiw_worm(cls, config, problem, qq, mom0):
        """
        Calculate impurity Green's function by symmetric improved estimator formula,
        Eq. (13) in PRB 100, 075119 (2019).
        """
        g0iw = orbspin.invert(problem.g0inviw)
        giw = g0iw + orbspin.multiply(g0iw, orbspin.multiply(qq + mom0, g0iw))
        return giw

    def calibrate_taudiff_max(self):
        phase1pairnum = 1000
        phase2taugrid = 100
        phase2pairnum = 1000
        nprogress = 100000
        ctqmc.init_counters()
        ctqmc.ctqmc_calibrate(True,
                              phase1pairnum, phase2taugrid, phase2pairnum,
                              nprogress)
        acctaulist = ctqmc.accpairtau
        acctaulist.sort()
        taudiffmax = DistributedSample([np.array((acctaulist[int(0.995 * phase1pairnum)],))],
                                       self.mpi_comm).mean()
        ctqmc.set_taudiffmax(taudiffmax)
        ctqmc.init_counters()
        ctqmc.ctqmc_calibrate(False,
                              phase1pairnum, phase2taugrid, phase2pairnum,
                              nprogress)
        return self.estimate_taudiff_max(self.config,
                                         self.mpi_comm,
                                         ctqmc.accpair,
                                         taudiffmax,
                                         phase2taugrid)

    @staticmethod
    def estimate_taudiff_max(config, mpi_comm, accept_pair, taudiffmax, phase2taugrid):
        accept_pair = DistributedSample([accept_pair],
                                        mpi_comm).mean()

        def expdecay(p, taus, acc):
            return acc - p[0] * np.exp(p[1] * taus)

        estimates = np.empty((accept_pair.shape[0],
                              accept_pair.shape[1]),
                             dtype=np.double)
        for b in range(accept_pair.shape[0]):
            for s in range(accept_pair.shape[1]):
                fit = opt.leastsq(expdecay,
                                  (accept_pair[b, s, 0], -1),
                                  (tf.tau_bins(taudiffmax,
                                               phase2taugrid,
                                               "centre"),
                                   accept_pair[b, s]))
                estimates[b, s] = np.log(1 - 0.99) / fit[0][1]

        return np.amax(estimates)

    @classmethod
    def get_resultset(cls, qmc_config):
        # Build up result array from QMC data
        norbitals = ctqmc.gtau.shape[0]
        result = {
            "gtau": ctqmc.gtau,
            "gtau-full": ctqmc.gtau_full,
            "rhist": ctqmc.rhisto,
            "lhist": ctqmc.lhisto,
            "hist": ctqmc.histo,
            "ssts-states": ctqmc.statessuperstates,
            "energies-eigenstates": ctqmc.eigenstatesenergies,
            "occbasis-mapping": ctqmc.occbasismapping,
            "hist-sst": ctqmc.outersuperstatehisto,
            "hist-state": ctqmc.outerstatehisto,
            "sign-sst": ctqmc.signhistosuperstates,
            "sign-state": ctqmc.signhistostates,
            "contrib-sst": ctqmc.tracecontribsuperstates,
            "contrib-state": ctqmc.tracecontribstates,
            "occ": ctqmc.occ,
            "gleg": ctqmc.gleg,
            "gleg-full": ctqmc.gleg_full,
            "accept-ins": ctqmc.accadd,
            "accept-ins4": ctqmc.accadd4,
            "accept-rem": ctqmc.accrem,
            "accept-rem4": ctqmc.accrem4,
            "accept-pair-tau": ctqmc.accpair,
            "accept-glob": ctqmc.accglob,
            "accept-shift": ctqmc.accshift,
            "accept-flavourchange": ctqmc.accflavc,
            "sign": ctqmc.mean_sign,
            "time-warmup": ctqmc.time_warmup,
            "time-simulation": ctqmc.time_sim,
            "time-sampling": ctqmc.time_sim_steps,
        }
        if qmc_config["Gtau_mean_step"] != 0:
            result["gtau-mean-step"] = ctqmc.gtau_mean_step
        if qmc_config["Gtau_mid_step"] != 0:
            result["gtau-mid-step"] = ctqmc.gtau_mid_step
        if qmc_config["sign_step"] != 0:
            result["sign-step"] = ctqmc.sign_step
        if qmc_config["MeasGiw"] != 0 and qmc_config["offdiag"] == 0:
            result["giw-meas"] = ctqmc.giw
        # if qmc_config["MeasGSigmaiw"] != 0:  # enable if GSigma Z-meas fixed
        #     result["gsigmaiw"] = ctqmc.gsigmaiw
        if qmc_config["MeasG2iw"] != 0:
            result["g2iw"] = ctqmc.g2iw
        if qmc_config["MeasDensityMatrix"] != 0 and qmc_config["segment"] == 0:
            result["densitymatrix"] = ctqmc.densitymatrix
            result["rho2"] = ctqmc.rho2.reshape((norbitals, 2)*4)
            result["rho1"] = ctqmc.rho1.reshape((norbitals, 2)*2)
        if qmc_config["MeasExpResDensityMatrix"] != 0:
            result["expresdensitymatrix"] = ctqmc.expresdensitymatrix
        if qmc_config["FourPnt"] & 1:
            result["g4tau"] = ctqmc.g4tau
        if qmc_config["FourPnt"] & 2:
            result["g4leg"] = ctqmc.g4leg
        if qmc_config["FourPnt"] & 4:
            result["g4iw"] = ctqmc.g4iw
        if qmc_config["FourPnt"] & 4 and qmc_config["MeasG4iwPP"] != 0:
            result["g4iw-pp"] = ctqmc.g4iw_pp
        if qmc_config["MeasurementTiming"] >= 0:
            result["time-giw"] = ctqmc.timings_giw
            result["time-g4iw-add"] = ctqmc.timings_g4iw_add
            result["time-g4iw-ft"] = ctqmc.timings_g4iw_ft
        if qmc_config["segment"] != 0:
            result["ntau-n0"] = ctqmc.ntau_n0
            result["hist-seg"]= ctqmc.histo_seg
        else:
            if qmc_config["MeasSusz"] != 0:
                result["ntau-n0"] = ctqmc.ntau_n0
            result["hist"]= ctqmc.histo

        # Attempt to read blocks only if chosen accumulator provides blocks
        if qmc_config["AccumGtau"] == 2:
            result["gtau-blocks"] = ctqmc.gtau_blocks

        # Note that for safety, we need to copy the data since the ctqmc
        # cleanup routine makes those quantities invalid.
        result = dict((k, result[k].copy()) for k in result)

        return result


    @classmethod
    def get_resultset_worm(cls, grp_str, qmc_config):
        # Build up result array from QMC data
        result = {
            "hist/"+grp_str: ctqmc.histo,
            "accept-ins/"+grp_str: ctqmc.accadd,
            "accept-ins4/"+grp_str: ctqmc.accadd4,
            "accept-rem/"+grp_str: ctqmc.accrem,
            "accept-rem4/"+grp_str: ctqmc.accrem4,
            "accept-pair-tau/"+grp_str: ctqmc.accpair,
            "accept-glob/"+grp_str: ctqmc.accglob,
            "accept-shift/"+grp_str: ctqmc.accshift,
            "accept-flavourchange/"+grp_str: ctqmc.accflavc,
            "accept-worm-ins/"+grp_str: ctqmc.accwormadd,
            "accept-worm-rem/"+grp_str: ctqmc.accwormrem,
            "accept-worm-rep/"+grp_str: ctqmc.accwormrep,
            "sign/"+grp_str: ctqmc.mean_sign,
            "time-warmup/"+grp_str: ctqmc.time_warmup,
            "time-simulation/"+grp_str: ctqmc.time_sim,
            "time-sampling/"+grp_str: ctqmc.time_sim_steps,
            "steps-worm-partition/"+grp_str: np.array(ctqmc.cntsampling, np.double),
            "worm-eta/"+grp_str: ctqmc.wormeta,
            # The conversion to double for large integers is necessary because
            # otherwise squaring it for the calculation of errors produces
            # negative values and nan errorbars.
        }
        if qmc_config["WormMeasGiw"] != 0:
            result.update({"giw-worm" + '/' + grp_str: ctqmc.giw_worm})
        if qmc_config["WormMeasGtau"] != 0:
            result.update({"gtau-worm" + '/' + grp_str: ctqmc.gtau_worm})
        if qmc_config["WormMeasGSigmaiw"] != 0:
            result.update({"gsigmaiw-worm" + '/' + grp_str: ctqmc.gsigmaiw_worm})
        if qmc_config["WormMeasQQ"] != 0:
            result.update({"qqiw-worm" + '/' + grp_str: ctqmc.qq_worm})
        # FIXME: all the following quantities are also extracted in
        # solve_component directly, but with additional normalization!
        # can we remove the following lines, and what about the
        # preceding ones? why are components of these four quantities
        # extracted here and not in solve_component, or vice versa?
        if qmc_config["WormMeasQQQQ"] != 0:
            result.update({"qqqq-worm" + '/' + grp_str: ctqmc.qqqq_worm})
        if qmc_config["WormMeasP2iwPH"] != 0:
            result.update({"p2iw-worm" + '/' + grp_str: ctqmc.p2iw_worm})
        if qmc_config["WormMeasP2iwPP"] != 0:
            result.update({"p2iwpp-worm" + '/' + grp_str: ctqmc.p2iwpp_worm})
        if qmc_config["WormMeasP2tauPH"] != 0:
            result.update({"p2tau-worm" + '/' + grp_str: ctqmc.p2tau_worm})
        if qmc_config["WormMeasP2tauPP"] != 0:
            result.update({"p2taupp-worm" + '/' + grp_str: ctqmc.p2taupp_worm})
        if qmc_config["WormMeasP3iwPH"] != 0:
            result.update({"p3iw-worm" + '/' + grp_str: ctqmc.p3iw_worm})
        if qmc_config["WormMeasP3iwPP"] != 0:
            result.update({"p3iwpp-worm" + '/' + grp_str: ctqmc.p3iwpp_worm})
        # Note that for safety, we need to copy the data since the ctqmc
        # cleanup routine makes those quantities invalid.
        result = dict((k, result[k].copy()) for k in result)

        return result
