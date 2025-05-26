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

        # FIXME: auto-determine
        self.g_diagonal_only = False

    def config_to_edipack(self, ed_mode):
        c = self.config

        def l(boolvar):
            return "T" if boolvar else "F"

        def f(floatvar):
            return f"{floatvar:e}".replace('e', 'd')

        params = dedent(f"""\
        NORB={self.problem.norbitals}
        NBATH={c["EDIPACK"]["NBATH"]}
        NSPIN=2
        NPH=0
        BATH_TYPE={c["EDIPACK"]["BATH_TYPE"]}
        ULOC=0.000000000,0.d0,0.d0,0.d0,0.d0
        UST=0.d0
        JH=0.d0
        JX=0.d0
        JP=0.d0
        NLOOP=100
        NSUCCESS=1
        DMFT_ERROR=1.000000000E-05
        SB_FIELD=0.0
        DELTASC=0.0
        BETA={f(self.problem.beta)}
        XMU=0.d0
        G_PH={",".join([f(0)]*self.problem.norbitals)}
        W0_PH=0.d0
        A_PH=0.d0
        GPHFILE=NONE
        SPIN_FIELD_X={",".join([f(0)]*self.problem.norbitals)}
        SPIN_FIELD_Y={",".join([f(0)]*self.problem.norbitals)}
        SPIN_FIELD_Z={",".join([f(0)]*self.problem.norbitals)}
        EXC_FIELD=0.d0,0.d0,0.d0,0.d0
        PAIR_FIELD={",".join([f(0)]*self.problem.norbitals)}
        CHISPIN_FLAG=F
        CHIDENS_FLAG=F
        CHIPAIR_FLAG=F
        CHIEXCT_FLAG=F
        ED_MODE={ed_mode}
        ED_FINITE_TEMP={l(c["EDIPACK"]["ED_FINITE_TEMP"])}
        ED_SECTORS=F
        ED_SECTORS_SHIFT=1
        ED_SPARSE_H=T
        ED_TOTAL_UD=T
        ED_TWIN={l(c["EDIPACK"]["ED_TWIN"])}
        ED_READ_UMATRIX=F
        ED_USE_KANAMORI=F
        ED_OBS_ALL=F
        ED_SOLVE_OFFDIAG_GF={l(not self.g_diagonal_only)}
        ED_PRINT_SIGMA=F
        ED_PRINT_G=F
        ED_PRINT_G0=F
        ED_PRINT_CHISPIN=T
        ED_PRINT_CHIDENS=T
        ED_PRINT_CHIPAIR=T
        ED_PRINT_CHIEXCT=T
        ED_ALL_G=T
        ED_VERBOSE={c["EDIPACK"]["ED_VERBOSE"]}
        ED_HW_BATH=2.000000000
        ED_OFFSET_BATH=1.000000000E-01
        LMATS={c["QMC"]["Niw"]}
        LREAL={c["EDIPACK"]["LREAL"]}
        LTAU={c["QMC"]["Ntau"]}
        LFIT={c["EDIPACK"]["LFIT"]}
        LPOS=100
        NREAD=0.d0
        NERR=1.000000000E-04
        NDELTA=1.000000000E-01
        NCOEFF=1.000000000
        WINI={f(c["EDIPACK"]["WINI"])}
        WFIN={f(c["EDIPACK"]["WFIN"])}
        XMIN=-3.000000000
        XMAX=3.000000000
        RDM_FLAG=T
        CHISPIN_FLAG=F
        CHISPIN_FLAG=F
        CHIDENS_FLAG=F
        CHIPAIR_FLAG=F
        CHIEXCT_FLAG=F
        HFMODE=F
        EPS={f(c["EDIPACK"]["EPS"])}
        CUTOFF={f(c["EDIPACK"]["CUTOFF"])}
        GS_THRESHOLD={f(c["EDIPACK"]["GS_THRESHOLD"])}
        LANC_METHOD={c["EDIPACK"]["LANC_METHOD"]}
        LANC_NSTATES_SECTOR={c["EDIPACK"]["LANC_NSTATES_SECTOR"]}
        LANC_NSTATES_TOTAL={c["EDIPACK"]["LANC_NSTATES_TOTAL"]}
        LANC_NSTATES_STEP={c["EDIPACK"]["LANC_NSTATES_STEP"]}
        LANC_NCV_FACTOR={c["EDIPACK"]["LANC_NCV_FACTOR"]}
        LANC_NCV_ADD={c["EDIPACK"]["LANC_NCV_ADD"]}
        LANC_NITER={c["EDIPACK"]["LANC_NITER"]}
        LANC_NGFITER={c["EDIPACK"]["LANC_NGFITER"]}
        LANC_TOLERANCE={f(c["EDIPACK"]["LANC_TOLERANCE"])}
        LANC_DIM_THRESHOLD={c["EDIPACK"]["LANC_DIM_THRESHOLD"]}
        CG_METHOD={c["EDIPACK"]["CG_METHOD"]}
        CG_GRAD={c["EDIPACK"]["CG_GRAD"]}
        CG_FTOL={f(c["EDIPACK"]["CG_FTOL"])}
        CG_STOP={c["EDIPACK"]["CG_STOP"]}
        CG_NITER={c["EDIPACK"]["CG_NITER"]}
        CG_WEIGHT={c["EDIPACK"]["CG_WEIGHT"]}
        CG_SCHEME={c["EDIPACK"]["CG_SCHEME"]}
        CG_NORM={c["EDIPACK"]["CG_NORM"]}
        CG_POW={c["EDIPACK"]["CG_POW"]}
        CG_MINIMIZE_VER={l(c["EDIPACK"]["CG_MINIMIZE_VER"])}
        CG_MINIMIZE_HH={f(c["EDIPACK"]["CG_MINIMIZE_HH"])}
        JZ_BASIS={l(c["EDIPACK"]["JZ_BASIS"])}
        JZ_MAX={l(c["EDIPACK"]["JZ_MAX"])}
        JZ_MAX_VALUE={f(c["EDIPACK"]["JZ_MAX_VALUE"])}
        SECTORFILE=sectors
        HFILE=hamiltonian
        HLOCFILE=inputHLOC.in
        UMATRIX_FILE=umatrix
        PRINT_INPUT_VARS=F
        LOGFILE=6
        """)

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

        # hermitize explicitly (check whether this is really needed)
        self.problem.interaction.u_matrix = 0.5 * (
            self.problem.interaction.u_matrix
            + np.conj(np.transpose(self.problem.interaction.u_matrix,
                                   (4, 5, 6, 7, 0, 1, 2, 3)))
        )

    def solve(self, iter_no, step_cache=None, prefixdir=None):
        stdout.flush()
        stderr.flush()

        def log(*args, **kwargs):
            if self.mpi_rank == 0:
                print(time.strftime("%y-%m-%d %H:%M:%S"), *args, **kwargs)

        time_solver = time.perf_counter()

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

        # F[o, s, o, s, w] -> Delta[s, s, o, o, w]
        fiw = np.conj(np.transpose(self.fiw, (1, 3, 0, 2, 4)))
        g0iw = np.transpose(orbspin.invert(self.problem.g0inviw), (2, 4, 1, 3, 0))
        # w2dynamics has spins in (down, up) order
        # flip w2d to ed spin order
        fiw = np.flip(fiw, axis=(0, 1))
        g0iw = np.flip(g0iw, axis=(0, 1))
        omega = 2.0 * np.pi * (np.arange(-fiw.shape[-1]//2, fiw.shape[-1]//2) + 0.5) / self.problem.beta

        fiw_pos = fiw[..., fiw.shape[-1]//2:]
        g0iw_pos = g0iw[..., g0iw.shape[-1]//2:]
        omega_pos = omega[omega.size//2:]

        if prefixdir is not None:
            w2d_wdir = Path.cwd()
            ed_wdir = Path(prefixdir)
            ed_wdir.mkdir(parents=True, exist_ok=True)
            chdir(ed_wdir)

        # check for spin-offdiagonals
        if (all(np.allclose(self.muimp[:, s1, :, 1 - s1], 0)
                for s1 in (0, 1))
            and all(np.allclose(self.fiw[:, s1, :, 1 - s1, :], 0)
                    for s1 in (0, 1))
            and np.allclose(np.imag(self.muimp), 0)):
            ed_mode = "normal"
        else:
            ed_mode = "nonsu2"

        self.config_to_edipack(ed_mode)
        # flip w2d to ed spin order
        self.ed.set_hloc(np.flip(self.muimp.transpose(1, 3, 0, 2), axis=(0, 1)).astype(complex))

        bath = self.ed.init_solver()
        if "bath" in oldcache and oldcache["bath"].shape == bath.shape:
            bath = oldcache["bath"]

        self.ed.reset_umatrix()
        # flip w2d to ed spin order
        ed_umat = np.flip(self.problem.interaction.u_matrix, axis=(1, 3, 5, 7))
        nonzero_umat_indices = np.transpose(
            np.nonzero(ed_umat)
        )
        for i in nonzero_umat_indices:
            value = ed_umat[tuple(i)]
            for icol, val in enumerate(i):
                if icol % 2 == 1: # spin: 0 to down, 1 to up
                    i[icol] = ord('d') + (ord('u') - ord('d')) * val
            self.ed.add_twobody_operator(i[0], chr(i[1]),
                                         i[2], chr(i[3]),
                                         i[4], chr(i[5]),
                                         i[6], chr(i[7]),
                                         value.real)  # accepts only real

        if iter_no == 0 and self.config["CI"]["initial_bath"] is not None:
            bath = np.array([float(x) for x in self.config["CI"]["initial_bath"]])
        elif self.config["EDIPACK"]["CG_SCHEME"] == "delta":
            bath = self.ed.chi2_fitgf(fiw_pos.astype(complex), bath, ispin=0)
            bath = self.ed.chi2_fitgf(fiw_pos.astype(complex), bath, ispin=1)
        elif self.config["EDIPACK"]["CG_SCHEME"] == "weiss":
            bath = self.ed.chi2_fitgf(g0iw_pos.astype(complex), bath, ispin=0)
            bath = self.ed.chi2_fitgf(g0iw_pos.astype(complex), bath, ispin=1)
        log(f"BATH: {bath}")

        # WRITE FIT RESULTS TO RESULT DICT
        result = {}

        fiw_fit = self.ed.get_delta(1.0j * omega, bath, 5)

        result["fiw-fit-error-rss"] = np.sqrt(np.sum(np.abs(fiw - fiw_fit)**2))
        result["fiw-fit-error-max"] = np.amax(np.abs(fiw - fiw_fit))

        log(f"Bath fit: whole range tot. mismatch {result['fiw-fit-error-rss']}")
        log(f"Bath fit: whole range max mismatch {result['fiw-fit-error-max']}", flush=True)

        # flip ed to w2d spin order
        fiw_fit = np.flip(fiw_fit, axis=(0, 1))
        fiw_fit = np.conj(np.transpose(fiw_fit, (2, 0, 3, 1, 4)))
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

        time_solver = time.perf_counter() - time_solver
        result["time-qmc"] = time_solver

        def op_c(orb, sp, norb=self.problem.norbitals):
            # reverse only orb, spin compensated by ed -> w2d
            orb = 1 - orb
            mat = np.zeros(2**(2 * norb) - 1, dtype=np.float64)
            mat[0::2] = 1
            mat = np.diagflat(mat, k=1)
            mat = np.reshape(mat, (2,) * (2 * norb * 2))
            mat = np.moveaxis(mat,
                              (norb * 2 - 1, 2 * norb * 2 - 1),
                              (orb + sp * norb, norb * 2 + orb + sp * norb))
            return np.reshape(mat, (2**(2 * norb),
                                    2**(2 * norb)))

        def op_cdag(orb, sp, norb=self.problem.norbitals):
            # reverse only orb, spin compensated by ed -> w2d
            orb = 1 - orb
            mat = np.zeros(2**(2 * norb) - 1, dtype=np.float64)
            mat[0::2] = 1
            mat = np.diagflat(mat, k=-1)
            mat = np.reshape(mat, (2,) * (2 * norb * 2))
            mat = np.moveaxis(mat,
                              (norb * 2 - 1, 2 * norb * 2 - 1),
                              (orb + sp * norb, norb * 2 + orb + sp * norb))
            return np.reshape(mat, (2**(2 * norb),
                                    2**(2 * norb)))

        try:
            rdm = np.loadtxt("reduced_density_matrix.ed")
            # could be indexed with occupation numbers if reshaped like this:
            # rdm = np.reshape(rdm, (2,) * (2 * self.problem.norbitals * 2))  # factors of 2: row/col and spin

            occ = np.zeros((self.problem.norbitals, 2) * 2, dtype=np.float64)
            rho1 = np.zeros((self.problem.norbitals, 2) * 2, dtype=np.complex128)
            rho2 = np.zeros((self.problem.norbitals, 2) * 4, dtype=np.complex128)
            for o1 in range(self.problem.norbitals):
                for s1 in (0, 1):
                    for o2 in range(self.problem.norbitals):
                        for s2 in (0, 1):
                            rho1[o1, s1, o2, s2] = np.trace(rdm @ op_cdag(o1, s1) @ op_c(o2, s2))
                            if o1 == o2 and s1 == s2:
                                occ[o1, s1, o2, s2] = (
                                    rho1[o1, s1, o2, s2].real
                                )
                            for o3 in range(self.problem.norbitals):
                                for s3 in (0, 1):
                                    for o4 in range(self.problem.norbitals):
                                        for s4 in (0, 1):
                                            rho2[o1, s1, o2, s2, o3, s3, o4, s4] = np.trace(
                                                rdm
                                                @ op_cdag(o1, s1)
                                                @ op_cdag(o2, s2)
                                                @ op_c(o3, s3)
                                                @ op_c(o4, s4)
                                            )
                                            if (o1 == o4 and s1 == s4
                                                and o3 == o2 and s3 == s2
                                                and o1 != o2 and s1 != s2):
                                                occ[o1, s1, o2, s2] = (
                                                    rho2[o1, s1, o2, s2, o3, s3, o4, s4].real
                                                )

            result["occ"] = occ
            result["rho1"] = rho1
            result["rho2"] = rho2

            # occbasis state index, orb, sp -> occupation
            obm = np.zeros((2**(2 * self.problem.norbitals), self.problem.norbitals, 2), dtype=int)
            divisor = 2
            for s in (0, 1):  # descending in ed order = ascending in w2d order
                for o in range(self.problem.norbitals - 1, -1, -1):
                    obm[1::divisor, o, s] = 1
                    divisor *= 2

            result["occbasis-mapping"] = obm
            result["densitymatrix"] = rdm
        except Exception as e:
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

        # flip ed to w2d spin order
        giw = np.transpose(
            np.flip(self.ed.get_gimp(ishape=5, axis='m'), axis=(0, 1)),
            (2, 0, 3, 1, 4)
        )
        giw = np.concatenate(
            (np.conj(np.transpose(giw, (2, 3, 0, 1, 4))[..., ::-1]),
             giw),
            axis=4
        )
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
