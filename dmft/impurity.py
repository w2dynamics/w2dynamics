"""Package abstracting impurity problems and the way of solving it.

Provides three basic classes for the dealing with the impurity picture of the
DMFT self-consistency loop. Their relation can be visualised as follows:

                           +----------------+
        ------------------>| ImpuritySolver |----------------->
          ImpurityProblem  +----------------+  ImpurityResult
"""
import time
from warnings import warn
import numpy as np
import scipy.optimize as opt

from ..auxiliaries.CTQMC import mctqmc as ctqmc
from ..auxiliaries import transform as tf
from ..auxiliaries import postprocessing as postproc
import dynamicalU as dynamicalU
import copy

import orbspin as orbspin


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

def reduce_stderr(mpi_comm, qtty_rank):
    mpi_size = mpi_comm.Get_size()
    
    if mpi_size == 1:
        # no need for such stuff
        qtty_err = np.empty_like(np.abs(qtty_rank))
        qtty_err[...] = np.nan
        return qtty_rank, qtty_err

    # make sure we have a numpy array
    qtty_rank = np.asarray(qtty_rank)
    qtty_mean = np.zeros_like(qtty_rank)
    
    # sum of the quantity
    mpi_comm.Allreduce(qtty_rank, qtty_mean)

    # sum abs-squared of the quantity
    qtty_rank = np.abs(qtty_rank)**2
    qttysq_mean = np.zeros_like(qtty_rank)
    mpi_comm.Allreduce(qtty_rank, qttysq_mean)

    # compute mean and standard deviation
    qtty_mean /= mpi_size
    qttysq_mean /= mpi_size
    qttysq_mean = np.sqrt(qttysq_mean - np.abs(qtty_mean)**2)

    # Get the standard error rather than the sample variance. This assumes that
    # the samples are uncorrelated, which is true in practise since they come
    # from completely different MC runs.
    qttysq_mean /= np.sqrt(mpi_size - 1.)
    return qtty_mean, qttysq_mean


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
    def postprocessing(self, siw_method="improved", smom_method="estimate"):
        # This function is disjoint from the initialisation, because it is a
        # "derived" quantity that should be computed after averaging over bins
        fallback = False

        if siw_method == "dyson":
            if self.g_diagonal_only:
                # If the solver is only capable of supplying the spin/orbital-
                # diagonal terms of the Green's function, we need to make sure
                # that also only the diagonal terms of G_0 are taken into
                # account in the Dyson equation, otherwise Sigma erroneously
                # "counteracts" these terms.  As an added performance benefit,
                # we can use 1/giw rather than the inversion in this case.
                self.siw = orbspin.promote_diagonal(
                                orbspin.extract_diagonal(self.problem.g0inviw) -
                                1/orbspin.extract_diagonal(self.giw))
            else:
                self.siw = self.problem.g0inviw - orbspin.invert(self.giw)

        elif siw_method == "improved":
            if self.gsigmaiw is None:
                warn("Cannot compute improved estimators - GSigmaiw missing\n"
                     "Falling back to Dyson equation", UserWarning, 2)
                siw_method = "dyson"
                fallback = True

            raise NotImplementedError()  # FIXME

        elif siw_method == 'improved_worm':
            if self.gsigmaiw is None:
                warn("Cannot compute improved estimators - GSigmaiw missing\n"
                     "Falling back to Dyson equation", UserWarning, 2)
                siw_method = "dyson"
                fallback = True
            
            if self.g_diagonal_only: 
                #TODO: check gsigma sign
                self.siw = -orbspin.promote_diagonal(
                              orbspin.extract_diagonal(self.gsigmaiw) / 
                              orbspin.extract_diagonal(self.giw))
            else:
                raise NotImplementedError("Offdiagonal worm improved \n"
                                          "estimators not implemented")

        else:
            raise ValueError("unknown siw_method: %s" % siw_method)

        if smom_method == "extract":
            warn("Extracting moment of self-energy from the data",
                 UserWarning, 2)
            self.smom = self.siw[:1].real.copy()
        elif smom_method == "estimate":
            if self.smom is None:
                warn("Cannot compute smom estimate - rho1 and rho2 missing\n"
                     "Falling back to extraction", UserWarning, 2)
                smom_method = "extract"
                fallback = True
        else:
            raise ValueError("unknown smom_method: %s" % smom_method)

        if fallback:
            self.postprocessing(siw_method, smom_method)

    def __init__(self, problem, giw, gsigmaiw=None, smom=None,
                 g_diagonal_only=False, **other):

        if g_diagonal_only:
            if not orbspin.is_diagonal(giw, 0):
                raise ValueError("Giw contains offdiagonals even though"
                                 "the result is marked as diagonal")
            if not (gsigmaiw is None or orbspin.is_diagonal(gsigmaiw, 0)):
                raise ValueError("GSigmaiw contains offdiagonals even though"
                                 "the result is marked as diagonal")

        self.problem = problem
        self.giw = np.asarray(giw)
        self.gsigmaiw = gsigmaiw
        self.smom = smom
        self.g_diagonal_only = g_diagonal_only
        self.other = other

class StatisticalImpurityResult(ImpurityResult):
    def collect(self, mpi_comm):
        # Provide an estimate for the error by averaging over bins
        self.giw, self.giw_err = reduce_stderr(mpi_comm, self.giw)
        if self.gsigmaiw is not None:
            self.gsigmaiw, self.gsigmaiw_err = reduce_stderr(mpi_comm, self.gsigmaiw)
        if self.smom is not None:
            self.smom, self.smom_err = reduce_stderr(mpi_comm, self.smom)

        # the reduce the others:
        other = {}
        for qtty_name, qtty_value in self.other.iteritems():
            qtty_mean, qtty_error = reduce_stderr(mpi_comm, qtty_value)
            qtty_value = np.rec.fromarrays((qtty_mean, qtty_error),
                                           names=("value", "error"))
            other[qtty_name] = qtty_value
        self.other = other

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
          print "         ****************************"
          print "         ******* U(w) Message *******"
          print "         ****************************"
          print " "

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
        def lattice_convention(qtty):
            return orbspin.promote_diagonal(qtty.transpose(2,0,1))

        # Call CT-QMC solver
        ctqmc.simid = self.mpi_rank
        time_qmc = time.clock()
        ctqmc.init_solver(self.umatrix, self.ftau, self.muimp, self.screening, iter_no+1)

        firstrun = True
        if type(mc_config_inout) is list and len(mc_config_inout) > 0:
            firstrun = False
            mc_config_inout.pop().set_to_ctqmc(ctqmc)

        if self.config["QMC"]["segment"] == 0:
            if self.config["QMC"]["TaudiffMax"] <= 0:
                ctqmc.init_counters()
                while ctqmc.apt_index <= 1000:
                    ctqmc.ctqmc_calibrate(10000, True)
                acctaulist = ctqmc.accpairtau
                acctaulist.sort()
                taudiffmax = reduce_stderr(self.mpi_comm,
                                           np.array((acctaulist[995],)))[0]
                ctqmc.set_taudiffmax(taudiffmax)
                ctqmc.init_counters()
                ctqmc.ctqmc_calibrate(100000, False)
                ctqmc.set_taudiffmax(self.estimate_taudiff_max(self.config,
                                                               self.mpi_comm,
                                                               ctqmc.accpair,
                                                               taudiffmax))
            else:
                ctqmc.set_taudiffmax(self.config["QMC"]["TaudiffMax"])
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
            ### this has to be done later in offiagonal case!
            #giw=giw.transpose(4,0,1,2,3)
            pass

        # Extract improved estimators
        gsigmaiw = self.get_gsigmaiw(self.config, self.problem, result)
        if gsigmaiw is not None:
            gsigmaiw = lattice_convention(gsigmaiw)

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


        # Construct result object from result set and collect it from all nodes
        result = StatisticalImpurityResult(self.problem, giw, gsigmaiw, smom,
                                           self.g_diagonal_only, **result)

        if self.mpi_comm is not None:
            result.collect(self.mpi_comm)


        # Cleanup
        ctqmc.dest_ctqmc()
        ctqmc.dest_paras()

        ### I have to reshape the giw at the very end, since the shitty mpi will not 
        ### process the giw array in its full form?!
        if self.config["QMC"]["offdiag"] == 1:
            result.giw=result.giw.reshape(self.problem.norbitals, 2,
                                          self.problem.norbitals, 2, self.problem.niw)
            result.giw=result.giw.transpose(4,0,1,2,3)

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
        """
        
        def lattice_convention(qtty):
            return orbspin.promote_diagonal(qtty.transpose(2,0,1))

        # Call CT-QMC solver
        ctqmc.simid = self.mpi_rank
        time_qmc = time.clock()
        ctqmc.init_solver(self.umatrix, self.ftau, self.muimp, self.screening, iter_no+1)
        
        if type(mc_config_inout) is list and len(mc_config_inout) > 0:
            mc_config_inout[0].set_to_ctqmc(ctqmc)

        if self.config["QMC"]["TaudiffMax"] <= 0:
            ctqmc.init_counters()
            while ctqmc.apt_index <= 1000:
                ctqmc.ctqmc_calibrate(10000, True)
            acctaulist = ctqmc.accpairtau
            acctaulist.sort()
            taudiffmax = reduce_stderr(self.mpi_comm,
                                       np.array((acctaulist[995],)))[0]
            ctqmc.set_taudiffmax(taudiffmax)
            ctqmc.init_counters()
            ctqmc.ctqmc_calibrate(100000, False)
            ctqmc.set_taudiffmax(self.estimate_taudiff_max(self.config,
                                                           self.mpi_comm,
                                                           ctqmc.accpair,
                                                           taudiffmax))
        else:
            ctqmc.set_taudiffmax(self.config["QMC"]["TaudiffMax"])
       
        #first warm-up 
        if not mc_config_inout:
           ctqmc.ctqmc_warmup(self.config["QMC"]["Nwarmups"])
           mc_config_inout.append(CtHybConfig.get_from_ctqmc(ctqmc))

        if self.config["QMC"]["Nwarmups2Plus"] < 0:
           warn("Nwarmups2Plus not set; assuming 0.1*Nwarmups for worm")
           WormWarmups=0.1*self.config["QMC"]["Nwarmups"]
        else:
           WormWarmups=self.config["QMC"]["Nwarmups2Plus"]

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
        if isector == 2 and self.config["QMC"]["WormMeasGiw"] == 1:
           result["giw-worm/" + str(icomponent).zfill(5)] = \
                ctqmc.giw_worm.copy()/(self.problem.nflavours**2)

        if isector == 2 and self.config["QMC"]["WormMeasGtau"] == 1:
           result["gtau-worm/" + str(icomponent).zfill(5)] = \
                ctqmc.gtau_worm.copy()/(self.problem.nflavours**2)

        if isector == 3:
           result["gsigmaiw-worm/" + str(icomponent).zfill(5)] = \
                ctqmc.gsigmaiw_worm.copy()/(self.problem.nflavours**2)

        if isector == 4:
           result["g4iw-worm/" + str(icomponent).zfill(5)] = \
                ctqmc.g4iw_worm.copy()/(self.problem.nflavours**4)

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
         
        # Construct result object from result set and collect it from all nodes        
        result = StatisticalImpurityResult(self.problem, giw, None, None,
                                           self.g_diagonal_only, **result)

        if self.mpi_comm is not None:
            result.collect(self.mpi_comm)
        
        #auxiliary data
        result_aux = self.get_resultset_worm(str(icomponent).zfill(5))
        result_aux = StatisticalImpurityResult(self.problem, giw, None, None,
                                               self.g_diagonal_only, **result_aux)

        if self.mpi_comm is not None:
            result_aux.collect(self.mpi_comm)

        # Cleanup
        ctqmc.dest_ctqmc()
        ctqmc.dest_paras()
        return result, result_aux

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
            axes_list = [-4, -3, -2, -1] + range(qtty.ndim - 4)
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
            axes_list = [-2, -1] + range(qtty.ndim - 2)
            qtty = qtty.transpose(axes_list)
        return qtty

    @classmethod
    def get_giw(cls, config, problem, result):
        #typ = "legendre"  # FIXME: QMC config does not contain FTtype
        typ = config["General"]["FTType"]
        if config["QMC"]["offdiag"]==1 and typ=="legendre":
            typ = "legendre_full"
        if typ == "none":
            giw = -result["giw-meas"]
        elif typ == "none_worm":
            giw = result["giw-worm"]
            giw = orbspin.extract_diagonal(giw.transpose(4,0,1,2,3)) \
                         .transpose(1,2,0)
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
            nbands=result["gleg-full"].shape[0]
            giw=np.zeros(shape=(nbands,2,nbands,2,problem.niw),dtype=complex)

            for b in range(0,nbands):
               for s in range(0,2):
                  giw[b,s,:,:,:] = tf.transform(tf_matrix, result["gleg-full"][b,s,:,:,:], None,
                                     onmismatch=tf.Truncate.tofirst, warn=False)

                  #### the old one
                  #tmp = tf.transform(tf_matrix, result["gleg-full"][b,s,:,:,:], None,
                                     #onmismatch=tf.Truncate.tofirst, warn=False)
                  #for b2 in range(0,nbands):
                     #for s2 in range(0,2):
                        #for w in range(0,problem.niw):
                           #print "w", w
                           #giw[b,s,b2,s2,w] = tmp[b2,s2,w]

            giw=giw.reshape(nbands*2,nbands*2,problem.niw)

            giw_new=np.zeros_like(giw,dtype=complex)
            for i in range(0,problem.niw):
               tmp=giw[:,:,i]
               giw_new[:,:,i]=0.5*(tmp.transpose(1,0)+tmp)

            giw=giw_new

        else:
            raise RuntimeError("Invalid transform type: `%s'" % typ)

        return giw   # copy not necessary here.

    @classmethod
    def get_gsigmaiw(cls, config, problem, result):
        gsigmaiw = None
        if config["QMC"]["MeasGSigmaiw"] == 1:
            gsigmaiw = result["gsigmaiw"]
        if config["QMC"]["WormMeasGSigmaiw"] == 1:
            if gsigmaiw is not None:
                warn("Worm GSigmaiw overwriting Z GSigmaiw", UserWarning, 2)

            gsigmaiw = result["gsigmaiw-worm"]
            gsigmaiw = orbspin.extract_diagonal(gsigmaiw.transpose(4,0,1,2,3)) \
                                    .transpose(1,2,0)

        return gsigmaiw

    @classmethod
    def estimate_taudiff_max(cls, config, mpi_comm, accept_pair, taudiffmax):
        accept_pair = reduce_stderr(mpi_comm, accept_pair)[0]

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
                                               1000,
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
            "rho2": ctqmc.rho2.reshape((norbitals, 2)*4),
            "rho1": ctqmc.rho1.reshape((norbitals, 2)*2),
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
            "accept-worm-ins": ctqmc.accwormadd,
            "accept-worm-rem": ctqmc.accwormrem,
            "accept-worm-rep": ctqmc.accwormrep,
            "sign": ctqmc.mean_sign,
            "time-warmup": ctqmc.time_warmup,
            "time-simulation": ctqmc.time_sim,
            "time-sampling": ctqmc.time_sim_steps,

            # The conversion to double for large integers is necessary because
            # otherwise squaring it for the calculation of errors produces
            # negative values and nan errorbars.
            "steps-worm-partition": np.array(ctqmc.cntsampling, np.double),
            "worm-eta": ctqmc.wormeta,
        }
        if qmc_config["Gtau_mean_step"] != 0:
            result["gtau-mean-step"] = ctqmc.gtau_mean_step
        if qmc_config["Gtau_mid_step"] != 0:
            result["gtau-mid-step"] = ctqmc.gtau_mid_step
        if qmc_config["sign_step"] != 0:
            result["sign-step"] = ctqmc.sign_step
        if qmc_config["MeasGiw"] != 0:
            result["giw-meas"] = ctqmc.giw
        if qmc_config["WormMeasGiw"] != 0:
            result["giw-worm"] = ctqmc.giw_worm
        if qmc_config["WormMeasGtau"] != 0:
            result["gtau-worm"] = ctqmc.gtau_worm
        if qmc_config["MeasG2iw"] != 0:
            result["g2iw"] = ctqmc.g2iw
        if qmc_config["MeasDensityMatrix"] != 0 or qmc_config["Eigenbasis"] != 0:
            result["densitymatrix"] = ctqmc.densitymatrix
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
        if qmc_config["WormMeasGSigmaiw"] != 0:
            result["gsigmaiw-worm"] = ctqmc.gsigmaiw_worm
        if qmc_config["WormMeasP2iwPH"] != 0:
            result["p2iw-worm"] = ctqmc.p2iw_worm
        if qmc_config["WormMeasP2iwPP"] != 0:
            result["p2iwpp-worm"] = ctqmc.p2iwpp_worm
        if qmc_config["WormMeasP2tauPH"] != 0:
            result["p2tau-worm"] = ctqmc.p2tau_worm
        if qmc_config["WormMeasP2tauPP"] != 0:
            result["p2taupp-worm"] = ctqmc.p2taupp_worm
        if qmc_config["WormMeasP3iwPH"] != 0:
            result["p3iw-worm"] = ctqmc.p3iw_worm
        if qmc_config["WormMeasP3iwPP"] != 0:
            result["p3iwpp-worm"] = ctqmc.p3iwpp_worm
        if qmc_config["segment"] != 0:
            result["ntau-n0"] = ctqmc.ntau_n0
            result["hist-seg"]= ctqmc.histo_seg
        else:
            if qmc_config["MeasSusz"] != 0:
                result["ntau-n0"] = ctqmc.ntau_n0
            result["hist"]= ctqmc.histo

        # Note that for safety, we need to copy the data since the ctqmc
        # cleanup routine makes those quantities invalid.
        result = dict((k, v.copy()) for k,v in result.iteritems())

        return result


    @classmethod
    def get_resultset_worm(cls,grp_str):
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
        }
        # Note that for safety, we need to copy the data since the ctqmc
        # cleanup routine makes those quantities invalid.
        result = dict((k, v.copy()) for k,v in result.iteritems())

        return result
