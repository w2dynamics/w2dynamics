from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
import time
from warnings import warn
from textwrap import dedent
from sys import stdout, stderr
from pathlib import Path
from os import chdir
from tempfile import mkstemp

import numpy as np
import scipy
import scipy.optimize as opt

import w2dyn.auxiliaries.transform as tf
import w2dyn.auxiliaries.postprocessing as postproc
from w2dyn.auxiliaries.input import write_u_matrix
from w2dyn.auxiliaries.statistics import DistributedSample
from w2dyn.dmft.selfcons import iw_to_tau_fast
from w2dyn.dmft.impurity import (ImpuritySolver,
                                 StatisticalImpurityResult)

import w2dyn.dmft.orbspin as orbspin


def lattice_convention(qtty):
    """Function to be used for tranposing and reshaping three-dimensional
    band/spin-diagonal arrays with (band, spin) as first two
    dimensions as obtained from the solver for some quantities into
    full five-dimensional band/spin-matrices with (band, spin, band,
    spin) as last four dimensions as used in the DMFT code.
    """
    return orbspin.promote_diagonal(qtty.transpose(2, 0, 1))


class EDIpySolver(ImpuritySolver):
    """EDIpy solver"""
    def __init__(self, config, seed=0, Uw=0, Uw_Mat=0, epsn=0, interactive=False, mpi_comm=None):
        super(EDIpySolver, self).__init__(config)

        from edipy2 import global_env
        self.ed = global_env

        self.mpi_comm = mpi_comm
        self.mpi_rank = 0
        if mpi_comm is not None:
            self.mpi_rank = mpi_comm.Get_rank()

        self.g_diagonal_only = (config["QMC"]["offdiag"] == 0)

    def config_to_edipack(self):
        c = self.config

        def l(boolvar):
            return "T" if boolvar else "F"

        def f(floatvar):
            return f"{floatvar:e}".replace('e', 'd')

        umatrixname = None
        if self.mpi_rank == 0:
            umatrixfd, umatrixname = mkstemp(".restart", "umat.", Path.cwd(), True)
            with open(umatrixfd, "x") as file:
                write_u_matrix(file,
                               np.flip(self.problem.interaction.u_matrix, axis=(1, 3, 5, 7)),
                               force_spin=True)
            # paths containing '/' must not be passed to read_input
            # and .restart is appended automatically
            umatrixname = Path(umatrixname).stem
        umatrixname = self.mpi_comm.bcast(umatrixname, root=0)

        params = dedent(f"""\
        NORB={self.problem.norbitals}                 !Number of impurity orbitals (max 5).
        NBATH={c["CI"]["npoles"]}                     !Number of bath sites:(normal=>Nbath per orb)(hybrid=>Nbath total)(replica/general=>Nbath=Nreplica/Ngeneral)
        NSPIN=2                                       !Number of spin degeneracy (max 2)
        NPH=0                                         !Max number of phonons allowed (cut off)
        BATH_TYPE={c["EDIPACK"]["bath_type"]}         !flag to set bath type: normal (1bath/imp), hybrid(1bath), replica(1replica/imp), general(replica++)
        ULOC=2.000000000,0.d0,0.d0,0.d0,0.d0          !Values of the local interaction per orbital (max 5)
        UST=0.d0                                      !Value of the inter-orbital interaction term
        JH=0.d0                                       !Hunds coupling
        JX=0.d0                                       !S-E coupling
        JP=0.d0                                       !P-H coupling
        NLOOP=100                                     !Max number of DMFT iterations.
        NSUCCESS=1                                    !Number of successive iterations below threshold for convergence
        DMFT_ERROR=1.000000000E-05                    !Error threshold for DMFT convergence
        SB_FIELD=1.000000000E-01                      !Value of a symmetry breaking field for magnetic solutions.
        DELTASC=2.000000000E-02                       !Value of the SC symmetry breaking term.
        BETA={f(self.problem.beta)}                      !Inverse temperature, at T=0 is used as a IR cut-off.
        XMU=0.d0                                      !Chemical potential. If HFMODE=T, xmu=0 indicates half-filling condition.
        G_PH={",".join([f(0)]*self.problem.norbitals)} !Electron-phonon coupling density constant
        W0_PH=0.d0                                    !Phonon frequency
        A_PH=0.d0                                     !Forcing field coupled to phonons displacement operator
        GPHFILE=NONE                                  !File of Phonon couplings. Put NONE to use only density couplings.
        SPIN_FIELD_X={",".join([f(0)]*self.problem.norbitals)} !magnetic field per orbital coupling to X-spin component
        SPIN_FIELD_Y={",".join([f(0)]*self.problem.norbitals)} !magnetic field per orbital coupling to Y-spin component
        SPIN_FIELD_Z={",".join([f(0)]*self.problem.norbitals)} !magnetic field per orbital coupling to Z-spin component
        EXC_FIELD=0.d0,0.d0,0.d0,0.d0                 !external field coupling to exciton order parameters
        PAIR_FIELD={",".join([f(0)]*self.problem.norbitals)} !pair field per orbital coupling to s-wave order parameter component
        CHISPIN_FLAG=F                                !Flag to activate spin susceptibility calculation.
        CHIDENS_FLAG=F                                !Flag to activate density susceptibility calculation.
        CHIPAIR_FLAG=F                                !Flag to activate pair susceptibility calculation.
        CHIEXCT_FLAG=F                                !Flag to activate excitonis susceptibility calculation.
        ED_MODE=normal                                !Flag to set ED type: normal=normal, superc=superconductive, nonsu2=broken SU(2)
        ED_FINITE_TEMP={l(c["EDIPACK"]["finite_temp"])} !flag to select finite temperature method. note that if T then lanc_nstates_total must be > 1
        ED_SECTORS=F                                  !flag to reduce sector scan for the spectrum to specific sectors +/- ed_sectors_shift.
        ED_SECTORS_SHIFT=1                            !shift to ed_sectors
        ED_SPARSE_H=T                                 !flag to select  storage of sparse matrix H (mem--, cpu++) if TRUE, or direct on-the-fly H*v product (mem++, cpu--) if FALSE
        ED_TOTAL_UD=T                                 !flag to select which type of quantum numbers have to be considered: T (default) total Nup-Ndw, F orbital based Nup-Ndw
        ED_TWIN=F                                     !flag to reduce (T) or not (F,default) the number of visited sector using twin symmetry.
        ED_READ_UMATRIX=T                             !flag to read (T) or not (F,default) the two-body operators from an external file.
        ED_OBS_ALL=F                                  !flag to print observables for every loop.
        ED_SOLVE_OFFDIAG_GF={l(not self.g_diagonal_only)} !flag to select the calculation of the off-diagonal impurity GF. this is T by default if bath_type/=normal
        ED_PRINT_SIGMA=F                              !flag to print impurity Self-energies
        ED_PRINT_G=F                                  !flag to print impurity Greens function
        ED_PRINT_G0=F                                 !flag to print non-interacting impurity Greens function
        ED_PRINT_CHISPIN=T                            !flag to print impurity spin susceptibility
        ED_PRINT_CHIDENS=T                            !flag to print impurity dens susceptibility
        ED_PRINT_CHIPAIR=T                            !flag to print impurity pair susceptibility
        ED_PRINT_CHIEXCT=T                            !flag to print impurity exct susceptibility
        ED_ALL_G=T                                    !flag to evaluate all the components of the impurity Green`s functions irrespective of the symmetries
        ED_VERBOSE={c["EDIPACK"]["verbosity"]}        !Verbosity level: 0=almost nothing --> 5:all. Really: all
        ED_HW_BATH=2.000000000                        !half-bandwidth for the bath initialization: flat in -ed_hw_bath:ed_hw_bath
        ED_OFFSET_BATH=1.000000000E-01                !offset for the initialization of diagonal terms in replica/general bath: -offset:offset
        LMATS={c["QMC"]["Niw"]}                       !Number of Matsubara frequencies.
        LREAL={c["CI"]["Nrealw"]}                     !Number of real-axis frequencies.
        LTAU={c["QMC"]["Ntau"]}                       !Number of imaginary time points.
        LFIT={c["CI"]["niw_fit"]}                     !Number of Matsubara frequencies used in the \Chi2 fit.
        LPOS=100                                      !Number of points for the lattice PDF.
        NREAD=0.d0                                    !Objective density for fixed density calculations.
        NERR=1.000000000E-04                          !Error threshold for fixed density calculations.
        NDELTA=1.000000000E-01                        !Initial step for fixed density calculations.
        NCOEFF=1.000000000                            !multiplier for the initial ndelta read from a file (ndelta-->ndelta*ncoeff).
        WINI={f(c["CI"]["minrealw"])}                    !Smallest real-axis frequency
        WFIN={f(c["CI"]["minrealw"])}                    !Largest real-axis frequency
        XMIN=-3.000000000                             !Smallest position for the lattice PDF
        XMAX=3.000000000                              !Largest position for the lattice PDF
        RDM_FLAG=T                                    !Flag to activate RDM calculation.
        CHISPIN_FLAG=F                                !Flag to activate spin susceptibility calculation.
        CHISPIN_FLAG=F                                !Flag to activate spin susceptibility calculation.
        CHIDENS_FLAG=F                                !Flag to activate density susceptibility calculation.
        CHIPAIR_FLAG=F                                !Flag to activate pair susceptibility calculation.
        CHIEXCT_FLAG=F                                !Flag to activate excitonis susceptibility calculation.
        HFMODE=F                                      !Flag to set the Hartree form of the interaction (n-1/2). see xmu.
        EPS={f(c["CI"]["imag_eta"])}                     !Broadening on the real-axis.
        CUTOFF={f(c["EDIPACK"]["cutoff"])}               !Spectrum cut-off, used to determine the number states to be retained.
        GS_THRESHOLD=1.000000000E-09                  !Energy threshold for ground state degeneracy loop up
        LANC_METHOD=arpack                            !select the lanczos method to be used in the determination of the spectrum. ARPACK (default), LANCZOS (T=0 only), DVDSON (no MPI)
        LANC_NSTATES_SECTOR=2                         !Initial number of states per sector to be determined.
        LANC_NSTATES_TOTAL=1                          !Initial number of total states to be determined.
        LANC_NSTATES_STEP=2                           !Number of states added to the spectrum at each step.
        LANC_NCV_FACTOR=10                            !Set the size of the block used in Lanczos-Arpack by multiplying the required Neigen (Ncv=lanc_ncv_factor*Neigen+lanc_ncv_add)
        LANC_NCV_ADD=0                                !Adds up to the size of the block to prevent it to become too small (Ncv=lanc_ncv_factor*Neigen+lanc_ncv_add)
        LANC_NITER=512                                !Number of Lanczos iteration in spectrum determination.
        LANC_NGFITER=200                              !Number of Lanczos iteration in GF determination. Number of momenta.
        LANC_TOLERANCE=1.000000000E-18                !Tolerance for the Lanczos iterations as used in Arpack and plain lanczos.
        LANC_DIM_THRESHOLD=1024                       !Min dimension threshold to use Lanczos determination of the spectrum rather than Lapack based exact diagonalization.
        CG_METHOD=0                                   !Conjugate-Gradient method: 0=NumericalRecipes, 1=minimize.
        CG_GRAD=0                                     !Gradient evaluation method: 0=analytic (default), 1=numeric.
        CG_FTOL=1.000000000E-05                       !Conjugate-Gradient tolerance.
        CG_STOP=0                                     !Conjugate-Gradient stopping condition: 0-2, 0=C1.AND.C2, 1=C1, 2=C2 with C1=|F_n-1 -F_n|<tol*(1+F_n), C2=||x_n-1 -x_n||<tol*(1+||x_n||).
        CG_NITER=500                                  !Max. number of Conjugate-Gradient iterations.
        CG_WEIGHT=1                                   !Conjugate-Gradient weight form: 1=1.0, 2=1/n , 3=1/w_n.
        CG_SCHEME=delta                               !Conjugate-Gradient fit scheme: delta or weiss.
        CG_NORM=elemental                             !Conjugate-Gradient norm definition: elemental (default) or frobenius.
        CG_POW=2                                      !Fit power for the calculation of the generalized distance as |G0 - G0and|**cg_pow
        CG_MINIMIZE_VER=F                             !Flag to pick old/.false. (Krauth) or new/.true. (Lichtenstein) version of the minimize CG routine
        CG_MINIMIZE_HH=1.000000000E-04                !Unknown parameter used in the CG minimize procedure.
        JZ_BASIS=F                                    !Flag to enable the Jz basis
        JZ_MAX=F                                      !Whether to cutoff Jz
        JZ_MAX_VALUE=1000.000000000                   !Maximum Jz
        SECTORFILE=sectors                            !File where to retrieve/store the sectors contributing to the spectrum.
        HFILE=hamiltonian                             !File where to retrieve/store the bath parameters.
        HLOCFILE=inputHLOC.in                         !File read the input local H.
        UMATRIX_FILE={umatrixname}                    !File read the two-body operator list from.
        PRINT_INPUT_VARS=F                            !Flag to toggle console printing of input variables list
        LOGFILE=6                                     !LOG unit."""
                          )

        paramname = None
        if self.mpi_rank == 0:
            paramfd, paramname = mkstemp(".conf", "inED.", Path.cwd(), True)
            with open(paramfd, "x") as file:
                file.write(params)
            paramname = Path(paramname).name  # paths containing dirs must not be passed to read_input
        paramname = self.mpi_comm.bcast(paramname, root=0)
        self.ed.read_input(paramname)
        self.mpi_comm.Barrier()

    def set_problem(self, problem, compute_fourpoint=0):
        # problem
        self.problem = problem

        # remove hybridization and one-particle Hamiltonian off-diagonals in diagonal calculation
        if self.g_diagonal_only:
            self.ftau = orbspin.promote_diagonal(orbspin.extract_diagonal(self.problem.ftau))
            self.fiw = orbspin.promote_diagonal(orbspin.extract_diagonal(self.problem.fiw))
            self.muimp = orbspin.promote_diagonal(orbspin.extract_diagonal(self.problem.muimp))
        else:
            self.ftau = self.problem.ftau
            self.fiw = self.problem.fiw
            self.muimp = self.problem.muimp

        # move tau-axis to last
        self.ftau = self.ftau.transpose(1,2,3,4,0)
        self.fiw = self.fiw.transpose(1,2,3,4,0)

        # solver uses a different convention for muimp
        self.muimp = -self.muimp
        self.umatrix = self.problem.interaction.u_matrix.reshape(
                             self.problem.nflavours, self.problem.nflavours,
                             self.problem.nflavours, self.problem.nflavours)

    def solve(self, iter_no, step_cache=None, prefixdir=None):
        stdout.flush()
        stderr.flush()

        def log(*args, **kwargs):
            print(time.strftime("%y-%m-%d %H:%M:%S"), *args, **kwargs)

        time_qmc = time.perf_counter()

        # Solver-specific cache
        firstrun = True
        if type(step_cache) is list and len(step_cache) > 0:
            firstrun = False
            oldcache = step_cache.pop(0)
        else:
            oldcache = {}

        newcache = {}
        if type(step_cache) is list:
            step_cache.append(newcache)

        # transposition and flip by convention
        fiw = np.transpose(self.fiw, (1, 3, 0, 2, 4))[:, :, :, :, ::-1]
        # w2dynamics has spins in (down, up) order
        fiw = np.flip(fiw, axis=(0, 1))
        omega = 2.0 * np.pi * (np.arange(-fiw.shape[-1]//2, fiw.shape[-1]//2) + 0.5) / self.problem.beta

        fiw_pos = fiw[..., fiw.shape[-1]//2:]
        omega_pos = omega[omega.size//2:]

        if prefixdir is not None:
            w2d_wdir = Path.cwd()
            ed_wdir = Path(prefixdir)
            ed_wdir.mkdir(parents=True, exist_ok=True)
            chdir(ed_wdir)

        self.config_to_edipack()
        self.ed.set_hloc(np.flip(self.muimp.transpose(1, 3, 0, 2), axis=(0, 1)).astype(complex))

        bath = self.ed.init_solver()
        if "bath" in oldcache and oldcache["bath"].shape == bath.shape:
            bath = oldcache["bath"]

        if iter_no == 0 and self.config["CI"]["initial_bath"] is not None:
            bath = np.array([float(x) for x in self.config["CI"]["initial_bath"]])
        else:
            bath = self.ed.chi2_fitgf(fiw_pos.astype(complex), bath, ispin=0)
            bath = self.ed.chi2_fitgf(fiw_pos.astype(complex), bath, ispin=1)
        log(f"BATH: {bath}")

        # WRITE FIT RESULTS TO RESULT DICT
        result = {}

        fiw_fit = self.ed.get_delta(1.0j * omega, bath, 5)

        result["fiw-fit-error-rss"] = np.sqrt(np.sum(np.abs(fiw - fiw_fit)**2))
        result["fiw-fit-error-max"] = np.amax(np.abs(fiw - fiw_fit))

        log(f"Bath fit: whole range tot. mismatch {result['fiw-fit-error-rss']}")
        log(f"Bath fit: whole range max mismatch {result['fiw-fit-error-max']}", flush=True)

        fiw_fit = np.flip(fiw_fit, axis=(0, 1))
        fiw_fit = np.transpose(fiw_fit, (2, 0, 3, 1, 4))
        result["fiw-fit"] = fiw_fit

        self.problem.g0inviw = (1.0j * (omega[:, None, None, None, None]
                                         * np.eye(self.problem.norbitals)[None, :, None, :, None]
                                         * np.eye(2)[None, None, :, None, :])
                                - self.muimp[None, :, :, :, :]
                                - np.transpose(fiw_fit, (4, 0, 1, 2, 3)))

        result["g0inviw-imp"] = np.transpose(self.problem.g0inviw, (1, 2, 3, 4, 0))
        result["g0iw-imp"] = np.transpose(orbspin.invert(self.problem.g0inviw), (1, 2, 3, 4, 0))
        # FIXME: adapt to edipack if possible
        # result["fiw-bath-energies"] = np.array([e for e, _ in bathparams])
        # result["fiw-bath-hybvecs"] = np.stack([np.reshape(v, (self.problem.norbitals, 2))
        #                                        for _, v in bathparams], axis=0)

        self.ed.solve(bath)
        newcache["bath"] = bath

        # FIXME: adapt to edipack if possible
        # totalocc = sum(qci.Operator.op_number(index, nf=nf)
        #                for index in range(nf))
        # totalocc_res = totalocc.braket(wfgs, wfgs)
        # log(f"Found ground state with total expected filling {totalocc_res:.5f}")
        # result["aim-total-occ"] = totalocc_res

        # totalsz = sum((-1 + 2 * (index % 2)) * qci.Operator.op_number(index, nf=nf) for index in range(nf))
        # totalsz_res = totalsz.braket(wfgs, wfgs)
        # log(f"Found ground state with total expected spin {totalsz_res:.5f}")
        # result["aim-total-sz"] = totalsz_res

        # energy = H.braket(wfgs, wfgs)
        # log(f"Ground state total energy: {energy}")
        # result["aim-total-gs-energy"] = energy

        time_qmc = time.perf_counter() - time_qmc
        result["time-qmc"] = time_qmc

        # Extract output quantities
        occ = np.full((self.problem.norbitals, 2) * 2, np.nan, dtype=np.float64)

        soccs = np.zeros((self.problem.norbitals, 2), dtype=np.float64)
        soccs[...] = self.ed.get_dens()[:, np.newaxis] / 2.0
        soccs += self.ed.get_mag("z")[:, np.newaxis] * np.array([-0.5, 0.5])[np.newaxis, :]

        for o in range(self.problem.norbitals):
            for s in (0, 1):
                occ[o, s, o, s] = soccs[o, s]

        doccs = self.ed.get_docc()

        for o in range(self.problem.norbitals):
            for s in (0, 1):
                occ[o, s, o, 1 - s] = doccs[o]

        result["occ"] = occ
        # result["rho1"] = np.reshape(rho1, (self.problem.norbitals, 2) * 2)
        # result["rho2"] = np.reshape(rho2, (self.problem.norbitals, 2) * 4)

        giw = np.transpose(
            np.flip(self.ed.get_gimp(ishape=5, axis='m'), axis=(0, 1)),
            (2, 0, 3, 1, 4)
        )
        giw = np.concatenate((np.conj(giw[..., ::-1]), giw), axis=4)
        result["gomega"] = np.transpose(
            np.flip(self.ed.get_gimp(ishape=5, axis='r'), axis=(0, 1)),
            (2, 0, 3, 1, 4)
        )
        result["somega"] = np.transpose(
            np.flip(self.ed.get_sigma(ishape=5, axis='r'), axis=(0, 1)),
            (2, 0, 3, 1, 4)
        )

        gtau_ft = iw_to_tau_fast(giw, self.config["QMC"]["Ntau"], self.problem.beta, axis=-1)
        result["gtau-ft"] = gtau_ft

        giw = np.ascontiguousarray(giw.transpose(4, 0, 1, 2, 3))
        giw = DistributedSample([giw],
                                self.mpi_comm)

        self.ed.finalize_solver()

        if prefixdir is not None:
            chdir(w2d_wdir)

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
            smom = DistributedSample([smom],
                                     self.mpi_comm)

        result = {key: DistributedSample([value],
                                         self.mpi_comm)
                  for (key, value) in result.items()}

        # Construct result object from result set and collect it from all nodes
        result = StatisticalImpurityResult(self.problem, giw, None, smom,
                                           self.g_diagonal_only, **result)

        return result


    def solve_component(self,iter_no,isector,icomponent,mc_config_inout=[]):
        raise NotImplementedError()

    def solve_worm(self, iter_no, log_function=None):
        raise NotImplementedError()
