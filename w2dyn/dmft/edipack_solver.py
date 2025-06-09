from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
import time
from textwrap import dedent
from sys import stdout, stderr
from pathlib import Path
from os import chdir
from tempfile import mkstemp
from threading import Lock

import numpy as np

import w2dyn.auxiliaries.postprocessing as postproc
import w2dyn.dmft.orbspin as orbspin
from w2dyn.auxiliaries.statistics import DistributedSample
from w2dyn.dmft.selfcons import iw_to_tau_fast
from w2dyn.dmft.impurity import (ImpuritySolver,
                                 StatisticalImpurityResult)


replicaerr = ValueError(
    "EDIPACK.BATH_TYPE is replica "
    "but no bath basis file has been provided.\n"
    "Set EDIPACK.bathbasisfile to the path of a file in "
    "EDIpack bath Hamiltonian file format, i.e.\n"
    "  i\n"
    "  V_0 lambda_00 ... lambda_0i\n"
    "  ... ...       ... ...\n"
    "  V_r lambda_r0 ... lambda_ri\n"
    "\n"
    "  mat_0_00_re mat_0_00_im mat_0_01_re ...\n"
    "  ...\n"
    "  mat_0_m0_re ...\n"
    "\n"
    "  mat_1_00_re ...\n"
    "  ...\n"
    "\n"
    "  ...\n"
    "\n"
    "  mat_i_00_re ...\n"
    "with i the number or basis matrices, r the number of "
    "replicas, m the number of orbitals (SU(2)) or twice "
    "the number of orbitals (with explicit spin) and V ignored.\n"
    "In case of explicit spin, set EDIPACK.basisconv_w2d to "
    "True if the matrices are given in "
    "w2dynamics spin convention rather than EDIpack spin "
    "convention."
)


generalerr = ValueError(
    "EDIPACK.BATH_TYPE is general "
    "but no bath basis file has been provided.\n"
    "Set EDIPACK.bathbasisfile to the path of a file in EDIpack "
    "bath Hamiltonian file format, i.e.\n"
    "  i\n"
    "  V_00 ... V_0m lambda_00 ... lambda_0i\n"
    "  ... ...       ... ...\n"
    "  V_r0 ... V_rm lambda_r0 ... lambda_ri\n"
    "\n"
    "  mat_0_00_re mat_0_00_im mat_0_01_re ...\n"
    "  ...\n"
    "  mat_0_m0_re ...\n"
    "\n"
    "  mat_1_00_re ...\n"
    "  ...\n"
    "\n"
    "  ...\n"
    "\n"
    "  mat_i_00_re ...\n"
    "with i the number or basis matrices, r the number of "
    "replicas, m the number of orbitals (SU(2)) or twice "
    "the number of orbitals (with explicit spin) and V ignored.\n"
    "In case of explicit spin, set EDIPACK.basisconv_w2d to "
    "True if the matrices are given "
    "in w2dynamics spin convention rather than EDIpack spin "
    "convention."
)


def lattice_convention(qtty):
    """Function to be used for tranposing and reshaping three-dimensional
    band/spin-diagonal arrays with (band, spin) as first two
    dimensions as obtained from the solver for some quantities into
    full five-dimensional band/spin-matrices with (band, spin, band,
    spin) as last four dimensions as used in the DMFT code.
    """
    return orbspin.promote_diagonal(qtty.transpose(2, 0, 1))


class EDIpackSolver(ImpuritySolver):
    """EDIpack ED solver"""
    # Lock to prevent concurrent use of EDIpack through multiple calls
    # to EDIpackSolver.solve in one process
    solver_lock = Lock()

    def __init__(self, config, seed=0, Uw=0, Uw_Mat=0, epsn=0,
                 interactive=False, mpi_comm=None):
        super(EDIpackSolver, self).__init__(config)

        from edipack2py import global_env
        self.ed = global_env

        self.mpi_comm = mpi_comm
        self.mpi_rank = 0
        if mpi_comm is not None:
            self.mpi_rank = mpi_comm.Get_rank()

        # FIXME: auto-determine
        self.g_diagonal_only = False

    def config_to_edipack(self, ed_mode):
        c = self.config
        edc = self.config["EDIPACK"]

        def l(boolvar):
            return "T" if boolvar else "F"

        def f(floatvar):
            return f"{floatvar:e}".replace('e', 'd')

        params = dedent(f"""\
        NORB={self.problem.norbitals}
        NBATH={edc["NBATH"]}
        NSPIN=2
        NPH={edc["NPH"]}
        BATH_TYPE={edc["BATH_TYPE"]}
        UST=0.d0
        JH=0.d0
        JX=0.d0
        JP=0.d0
        NLOOP=100
        NSUCCESS=1
        DMFT_ERROR=1.000000000E-05
        SB_FIELD={f(edc["SB_FIELD"])}
        DELTASC=0.0
        BETA={f(self.problem.beta)}
        XMU=0.d0
        G_PH={",".join(f(x)
              for x
              in (edc["G_PH"]
                  if edc["G_PH"] is not None
                  else [0] * self.problem.norbitals))}
        W0_PH={f(edc["W0_PH"])}
        A_PH={f(edc["A_PH"])}
        GPHFILE={edc["GPHFILE"]}
        SPIN_FIELD_X={",".join([f(0)] * self.problem.norbitals)}
        SPIN_FIELD_Y={",".join([f(0)] * self.problem.norbitals)}
        SPIN_FIELD_Z={",".join([f(0)] * self.problem.norbitals)}
        EXC_FIELD=0.d0,0.d0,0.d0,0.d0
        PAIR_FIELD={",".join([f(0)] * self.problem.norbitals)}
        CHISPIN_FLAG={l(edc["CHISPIN_FLAG"])}
        CHIDENS_FLAG={l(edc["CHIDENS_FLAG"])}
        CHIPAIR_FLAG={l(edc["CHIPAIR_FLAG"])}
        CHIEXCT_FLAG={l(edc["CHIEXCT_FLAG"])}
        ED_MODE={ed_mode}
        ED_FINITE_TEMP={l(edc["ED_FINITE_TEMP"])}
        ED_SECTORS={l(edc["ED_SECTORS"])}
        ED_SECTORS_SHIFT={edc["ED_SECTORS_SHIFT"]}
        ED_SPARSE_H={l(edc["ED_SPARSE_H"])}
        ED_TOTAL_UD={l(edc["ED_TOTAL_UD"])}
        ED_TWIN={l(edc["ED_TWIN"])}
        ED_READ_UMATRIX=F
        ED_USE_KANAMORI=F
        ED_OBS_ALL=F
        ED_SOLVE_OFFDIAG_GF={l(not self.g_diagonal_only)}
        ED_PRINT_SIGMA=F
        ED_PRINT_G=F
        ED_PRINT_G0=F
        ED_PRINT_CHISPIN={l(edc["ED_PRINT_CHISPIN"])}
        ED_PRINT_CHIDENS={l(edc["ED_PRINT_CHIDENS"])}
        ED_PRINT_CHIPAIR={l(edc["ED_PRINT_CHIPAIR"])}
        ED_PRINT_CHIEXCT={l(edc["ED_PRINT_CHIEXCT"])}
        ED_ALL_G=T
        ED_VERBOSE={edc["ED_VERBOSE"]}
        ED_HW_BATH=2.000000000
        ED_OFFSET_BATH=1.000000000E-01
        LMATS={c["QMC"]["Niw"]}
        LREAL={edc["LREAL"]}
        LTAU={c["QMC"]["Ntau"]}
        LFIT={edc["LFIT"]}
        LPOS=100
        NREAD=0.d0
        NERR=1.000000000E-04
        NDELTA=1.000000000E-01
        NCOEFF=1.000000000
        WINI={f(edc["WINI"])}
        WFIN={f(edc["WFIN"])}
        XMIN=-3.000000000
        XMAX=3.000000000
        RDM_FLAG=T
        CHISPIN_FLAG=F
        CHISPIN_FLAG=F
        CHIDENS_FLAG=F
        CHIPAIR_FLAG=F
        CHIEXCT_FLAG=F
        HFMODE=F
        EPS={f(edc["EPS"])}
        CUTOFF={f(edc["CUTOFF"])}
        GS_THRESHOLD={f(edc["GS_THRESHOLD"])}
        LANC_METHOD={edc["LANC_METHOD"]}
        LANC_NSTATES_SECTOR={edc["LANC_NSTATES_SECTOR"]}
        LANC_NSTATES_TOTAL={edc["LANC_NSTATES_TOTAL"]}
        LANC_NSTATES_STEP={edc["LANC_NSTATES_STEP"]}
        LANC_NCV_FACTOR={edc["LANC_NCV_FACTOR"]}
        LANC_NCV_ADD={edc["LANC_NCV_ADD"]}
        LANC_NITER={edc["LANC_NITER"]}
        LANC_NGFITER={edc["LANC_NGFITER"]}
        LANC_TOLERANCE={f(edc["LANC_TOLERANCE"])}
        LANC_DIM_THRESHOLD={edc["LANC_DIM_THRESHOLD"]}
        CG_METHOD={edc["CG_METHOD"]}
        CG_GRAD={edc["CG_GRAD"]}
        CG_FTOL={f(edc["CG_FTOL"])}
        CG_STOP={edc["CG_STOP"]}
        CG_NITER={edc["CG_NITER"]}
        CG_WEIGHT={edc["CG_WEIGHT"]}
        CG_SCHEME={edc["CG_SCHEME"]}
        CG_NORM={edc["CG_NORM"]}
        CG_POW={edc["CG_POW"]}
        CG_MINIMIZE_VER={l(edc["CG_MINIMIZE_VER"])}
        CG_MINIMIZE_HH={f(edc["CG_MINIMIZE_HH"])}
        JZ_BASIS={l(edc["JZ_BASIS"])}
        JZ_MAX={l(edc["JZ_MAX"])}
        JZ_MAX_VALUE={f(edc["JZ_MAX_VALUE"])}
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
            with open(paramfd, "w") as file:
                file.write(params)
            # paths containing dirs must not be passed to read_input
            paramname = Path(paramname).name
        paramname = self.mpi_comm.bcast(paramname, root=0)
        self.ed.read_input(paramname)
        self.mpi_comm.Barrier()

    def set_problem(self, problem, compute_fourpoint=0):
        # problem
        self.problem = problem

        # remove hybridization and one-particle Hamiltonian
        # off-diagonals in diagonal calculation
        if self.g_diagonal_only:
            self.ftau = orbspin.promote_diagonal(
                orbspin.extract_diagonal(self.problem.ftau)
            )
            self.fiw = orbspin.promote_diagonal(
                orbspin.extract_diagonal(self.problem.fiw)
            )
            self.muimp = orbspin.promote_diagonal(
                orbspin.extract_diagonal(self.problem.muimp)
            )
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
                print(time.strftime("%y-%m-%d %H:%M:%S"),
                      *args,
                      **{"flush": True, **kwargs})

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
        g0iw = np.transpose(
            orbspin.invert(self.problem.g0inviw),
            (2, 4, 1, 3, 0)
        )
        # w2dynamics has spins in (down, up) order
        # flip w2d to ed spin order
        fiw = np.flip(fiw, axis=(0, 1))
        g0iw = np.flip(g0iw, axis=(0, 1))
        omega = 2.0 * np.pi * (np.arange(-fiw.shape[-1]//2,
                                         fiw.shape[-1]//2)
                               + 0.5) / self.problem.beta

        fiw_pos = fiw[..., fiw.shape[-1]//2:]
        g0iw_pos = g0iw[..., g0iw.shape[-1]//2:]

        # check for spin-offdiagonals
        if (all(np.allclose(self.muimp[:, s1, :, 1 - s1], 0)
                for s1 in (0, 1))
            and all(np.allclose(self.fiw[:, s1, :, 1 - s1, :], 0)
                    for s1 in (0, 1))
            and np.allclose(np.imag(self.muimp), 0)):
            ed_mode = "normal"
            log("Setting ED_MODE=normal after Hamiltonian "
                "and hybridization check")
        else:
            ed_mode = "nonsu2"
            log("Setting ED_MODE=nonsu2 after Hamiltonian "
                "and hybridization check")

        def check_bathbasis_nonsu2(bm):
            if (not np.allclose(bm[0, 1], 0)
                or not np.allclose(bm[1, 0], 0)
                or not np.allclose(bm[0, 0], bm[1, 1])):
                return True
            return False

        if (bath_type := self.config["EDIPACK"]["BATH_TYPE"]) in (
                "replica",
                "general"
        ):
            if f"{bath_type}bath" in oldcache:
                basis, lambdas = oldcache[f"{bath_type}bath"]
            else:
                if self.config["EDIPACK"]["bathbasisfile"] is None:
                    raise replicaerr if bath_type == "replica" else generalerr
                log(f"Reading basis for BATH_TYPE={bath_type} from "
                    f"file {self.config['EDIPACK']['bathbasisfile']}")
                with open(self.config["EDIPACK"]["bathbasisfile"], "r") as f:
                    basis, lambdas, _ = read_bath_basis_file(
                        f,
                        general=(False if bath_type == "replica" else True),
                        w2d_spin_convention=self.config["EDIPACK"]["basisconv_w2d"],
                        norb=self.problem.norbitals,
                        spin_autodetect=True
                    )
            if check_bathbasis_nonsu2(basis):
                ed_mode = "nonsu2"
                log("Setting ED_MODE=nonsu2 after bath basis check")
            newcache[f"{bath_type}bath"] = (basis, lambdas)
            if lambdas.shape[0] != self.config["EDIPACK"]["NBATH"]:
                raise ValueError("Bath basis file inconsistent with NBATH")

        actual_gphfile_param = self.config["EDIPACK"]["GPHFILE"]
        if self.config["EDIPACK"]["GPHFILE"] != "NONE":
            log(f"Reading GPHFILE {self.config['EDIPACK']['GPHFILE']}")
            with open(self.config["EDIPACK"]["GPHFILE"], "rb") as f:
                gphcontent = f.read()
        else:
            gphcontent = None

        if prefixdir is not None:
            w2d_wdir = Path.cwd()
            log(f"Changing to ED runtime directory {prefixdir}")
            ed_wdir = Path(prefixdir)
            ed_wdir.mkdir(parents=True, exist_ok=True)
            chdir(ed_wdir)

        EDIpackSolver.solver_lock.acquire()

        if gphcontent is not None:
            if self.mpi_rank == 0:
                gphfd, gphpath = mkstemp(prefix="gphfile", dir=Path.cwd())
                with open(gphfd, "wb") as f:
                    f.write(gphcontent)
            else:
                gphpath = None
            self.mpi_comm.bcast(gphpath, root=0)
            self.config["EDIPACK"]["GPHFILE"] = Path(gphpath).name

        self.config_to_edipack(ed_mode)
        self.config["EDIPACK"]["GPHFILE"] = actual_gphfile_param
        # flip w2d to ed spin order
        self.ed.set_hloc(np.flip(
            self.muimp.transpose(1, 3, 0, 2),
            axis=(0, 1)
        ).astype(complex))

        if self.config["EDIPACK"]["BATH_TYPE"] == "replica":
            log("Setting basis for replica bath")
            self.ed.set_hreplica(basis, lambdas)
        elif self.config["EDIPACK"]["BATH_TYPE"] == "general":
            log("Setting basis for general bath")
            self.ed.set_hgeneral(basis, lambdas)

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
            bath = np.array([float(x)
                             for x in self.config["CI"]["initial_bath"]])
        elif self.config["EDIPACK"]["CG_SCHEME"] == "delta":
            log("Fitting bath to hybridization function Delta")
            bath = self.ed.chi2_fitgf(fiw_pos.astype(complex), bath, ispin=0)
            bath = self.ed.chi2_fitgf(fiw_pos.astype(complex), bath, ispin=1)
        elif self.config["EDIPACK"]["CG_SCHEME"] == "weiss":
            log("Fitting bath to Weiss field G0")
            bath = self.ed.chi2_fitgf(g0iw_pos.astype(complex), bath, ispin=0)
            bath = self.ed.chi2_fitgf(g0iw_pos.astype(complex), bath, ispin=1)
        log(f"Bath array: {bath}")

        # WRITE FIT RESULTS TO RESULT DICT
        result = {}

        fiw_fit = self.ed.get_delta(1.0j * omega, bath, 5)

        result["fiw-fit-error-rss"] = np.sqrt(np.sum(np.abs(fiw - fiw_fit)**2))
        result["fiw-fit-error-max"] = np.amax(np.abs(fiw - fiw_fit))

        log(f"Bath root of summed squares mismatch over the whole frequency range: "
            f"{result['fiw-fit-error-rss']}")
        log(f"Bath maximum mismatch over the whole frequency range: "
            f"{result['fiw-fit-error-max']}", flush=True)

        # flip ed to w2d spin order
        fiw_fit = np.flip(fiw_fit, axis=(0, 1))
        fiw_fit = np.conj(np.transpose(fiw_fit, (2, 0, 3, 1, 4)))
        result["fiw-fit"] = fiw_fit

        edipack_g0 = np.flip(
            self.ed.get_g0and(1.0j * omega, bath, ishape=5),
            axis=(0, 1)
        ).transpose(4, 2, 0, 3, 1)
        edipack_g0inviw = orbspin.invert(edipack_g0)

        self.problem.g0inviw = edipack_g0inviw

        # FIXME: figure out problem here
        # self.problem.g0inviw = (
        #     1.0j * (omega[:, None, None, None, None]
        #             * np.eye(self.problem.norbitals)[None, :, None, :, None]
        #             * np.eye(2)[None, None, :, None, :])
        #     - self.muimp[None, :, :, :, :]
        #     - np.transpose(fiw_fit, (4, 0, 1, 2, 3))
        # )

        result["g0inviw-imp"] = np.transpose(self.problem.g0inviw,
                                             (1, 2, 3, 4, 0))
        # result["g0iw-imp"] = np.transpose(
        #     orbspin.invert(self.problem.g0inviw),
        #     (1, 2, 3, 4, 0)
        # )
        result["g0iw-imp"] = np.transpose(edipack_g0, (1, 2, 3, 4, 0))

        # FIXME: adapt to edipack if possible
        # result["fiw-bath-energies"] = np.array([e for e, _ in bathparams])
        # result["fiw-bath-hybvecs"] = np.stack(
        #     [np.reshape(v, (self.problem.norbitals, 2))
        #      for _, v in bathparams],
        #     axis=0
        # )

        log("Performing exact diagonalization")
        self.ed.solve(bath)
        newcache["bath"] = bath

        # FIXME: adapt to edipack if possible
        # totalocc = sum(qci.Operator.op_number(index, nf=nf)
        #                for index in range(nf))
        # totalocc_res = totalocc.braket(wfgs, wfgs)
        # log(f"Found ground state with total expected filling "
        #     f"{totalocc_res:.5f}")
        # result["aim-total-occ"] = totalocc_res

        # totalsz = sum((-1 + 2 * (index % 2))
        #               * qci.Operator.op_number(index, nf=nf)
        #               for index in range(nf))
        # totalsz_res = totalsz.braket(wfgs, wfgs)
        # log(f"Found ground state with total expected spin {totalsz_res:.5f}")
        # result["aim-total-sz"] = totalsz_res

        # energy = H.braket(wfgs, wfgs)
        # log(f"Ground state total energy: {energy}")
        # result["aim-total-gs-energy"] = energy

        time_solver = time.perf_counter() - time_solver
        result["time-qmc"] = time_solver

        self.mpi_comm.Barrier()  # prevent undeterministic occ problems

        def op_c(orb, sp, norb=self.problem.norbitals):
            # reverse only orb, spin compensated by ed -> w2d
            orb = norb - 1 - orb
            halfshape = (2,) * (norb * 2)  # occupation number basis state dims
            mat = np.zeros(halfshape + halfshape, dtype=np.float64)
            for i in np.ndindex(halfshape):
                # iterate over column / ingoing
                if i[orb + sp * norb] == 0:
                    # skip if arg orb/sp unoccupied
                    continue
                # make index with arg orb/sp unoccupied in j
                j = list(i)
                j[orb + sp * norb] = 0
                j = tuple(j)
                # sign: -1 if odd number of preceding orb/sp occupied
                mat[j + i] = 1 - 2 * (np.sum(i[:orb + sp * norb]) % 2)
            return np.reshape(mat, (2**(2 * norb),
                                    2**(2 * norb)))

        def op_cdag(orb, sp, norb=self.problem.norbitals):
            return np.transpose(op_c(orb, sp, norb))

        try:
            log("Processing reduced density matrix")

            rdm = np.loadtxt("reduced_density_matrix.ed")
            # could be indexed with occupation numbers if reshaped like this:
            # factors of 2: row/col and spin
            # rdm = np.reshape(rdm, (2,) * (2 * self.problem.norbitals * 2))

            occ = np.zeros((self.problem.norbitals, 2) * 2,
                           dtype=np.float64)
            rho1 = np.zeros((self.problem.norbitals, 2) * 2,
                            dtype=np.complex128)
            rho2 = np.zeros((self.problem.norbitals, 2) * 4,
                            dtype=np.complex128)
            for o1 in range(self.problem.norbitals):
                for s1 in (0, 1):
                    for o2 in range(self.problem.norbitals):
                        for s2 in (0, 1):
                            rho1[o1, s1, o2, s2] = np.trace(rdm
                                                            @ op_cdag(o1, s1)
                                                            @ op_c(o2, s2))
                            if o1 == o2 and s1 == s2:
                                occ[o1, s1, o2, s2] = (
                                    rho1[o1, s1, o2, s2].real
                                )
                            for o3 in range(self.problem.norbitals):
                                for s3 in (0, 1):
                                    for o4 in range(self.problem.norbitals):
                                        for s4 in (0, 1):
                                            rho2[o1, s1, o2, s2,
                                                 o3, s3, o4, s4] = np.trace(
                                                rdm
                                                @ op_cdag(o1, s1)
                                                @ op_cdag(o2, s2)
                                                @ op_c(o3, s3)
                                                @ op_c(o4, s4)
                                            )
                                            if (o1 == o4 and s1 == s4
                                                and o3 == o2 and s3 == s2
                                                and (o1 != o2 or s1 != s2)):
                                                occ[o1, s1, o2, s2] = (
                                                    rho2[o1, s1, o2, s2,
                                                         o3, s3, o4, s4].real
                                                )

            result["occ"] = occ
            result["rho1"] = rho1
            result["rho2"] = rho2

            # occbasis state index, orb, sp -> occupation
            obm = np.zeros((2**(2 * self.problem.norbitals),
                            self.problem.norbitals, 2), dtype=int)
            divisor = 2
            for s in (0, 1):  # descending in ed order = ascending in w2d order
                for o in range(self.problem.norbitals - 1, -1, -1):
                    obm[1::divisor, o, s] = 1
                    divisor *= 2

            result["occbasis-mapping"] = obm
            result["densitymatrix"] = rdm
        except Exception as e:
            log("Falling back to direct extraction of occupations")

            occ = np.full((self.problem.norbitals, 2) * 2, np.nan,
                          dtype=np.float64)

            soccs = np.zeros((self.problem.norbitals, 2), dtype=np.float64)
            soccs[...] = self.ed.get_dens()[:, np.newaxis] / 2.0
            soccs += (self.ed.get_mag("z")[:, np.newaxis]
                      * np.array([-0.5, 0.5])[np.newaxis, :])

            for o in range(self.problem.norbitals):
                for s in (0, 1):
                    occ[o, s, o, s] = soccs[o, s]

            doccs = self.ed.get_docc()

            for o in range(self.problem.norbitals):
                for s in (0, 1):
                    occ[o, s, o, 1 - s] = doccs[o]

            result["occ"] = occ

        log("Extracting Green's function and self-energy")
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

        gtau_ft = iw_to_tau_fast(giw, self.config["QMC"]["Ntau"],
                                 self.problem.beta, axis=-1)
        result["gtau-ft"] = gtau_ft

        giw = np.ascontiguousarray(giw.transpose(4, 0, 1, 2, 3))
        giw = DistributedSample([giw],
                                self.mpi_comm)

        self.ed.finalize_solver()
        EDIpackSolver.solver_lock.release()

        if prefixdir is not None:
            log("Leaving ED runtime directory for current iteration")
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


def read_bath_basis_file(f,
                         general=False,
                         spin=False,
                         w2d_spin_convention=False,
                         norb=None,
                         spin_autodetect=False):
    """Read the basis matrices and coefficients for a replica /
    general bath from a file object containing data in the EDIpack
    bath hamiltonian file format.

    Parameters
    ----------
    f: file object to read from
    general: general bath (otherwise replica)
    spin: spin dimensions present (in EDIpack-generated if non-SU(2))
    w2d_spin_convention: assume faster running spin dim in order down, up
    norb: number of orbitals (for size checking)
    spin_autodetect: autodetect spin (only allowed if norb provided)
    """
    while (nextline := f.readline().strip()).startswith("#"):
        pass  # ignore comment header(s)
    nbasis = int(nextline)  # size of basis

    # coefficient lines (possibly flavor dependent hybridization V
    # followed by basis coefficients lambda)
    coefflines = []
    while (nextline := f.readline().strip()) != '':
        coefflines.append(nextline)
    coefficients = np.genfromtxt(coefflines)

    matdim = None
    if norb is not None:
        nflav = 2 * norb
        if general and spin_autodetect:
            if coefficients.shape[1] == (nflav + nbasis):
                spin = True
            elif coefficients.shape[1] == (norb + nbasis):
                spin = False
            else:
                raise ValueError("bath basis file malformed: "
                                 "wrong number of columns: "
                                 f"expected {norb} (V) + {nbasis} (lambda) "
                                 f"or {nflav} (V) + {nbasis} (lambda), "
                                 f"found {coefficients.shape[1]}")
            matdim = norb if not spin else nflav
        elif not spin_autodetect:
            matdim = norb if not spin else nflav
        flavdim = 1 if not general else norb if not spin else nflav

        if coefficients.shape[1] != (flavdim + nbasis):
            raise ValueError("bath basis file malformed: "
                             "wrong number of columns: "
                             f"expected {flavdim} (V) + {nbasis} (lambda), "
                             f"found {coefficients.shape[1]}")
    else:
        if spin_autodetect:
            raise ValueError("spin_autodetect but norb is None")
        flavdim = coefficients.shape[1] - nbasis
        nflav = None if not general else flavdim * 2 if not spin else flavdim
        norb = None if not general else nflav // 2
        matdim = None if not general else norb if not spin else nflav

    vs = coefficients[:, :flavdim]
    lambdas = coefficients[:, flavdim:]
    nreplica = coefficients.shape[0]

    if general and spin and w2d_spin_convention:
        vs = np.transpose(np.flip(np.reshape(vs, (nreplica, norb, 2)), axis=2),
                          (0, 2, 1))
    elif general and spin:
        vs = np.reshape(vs, (nreplica, 2, norb))
    elif not general:
        vs = vs[:, 0]

    # basis matrices
    basis = []
    for i in range(nbasis):
        matrixlines = []
        while (nextline := f.readline().strip()) != '':
            # complex numbers in actual EDIpack output are "(real,
            # imag)", we want only numbers and whitespace
            matrixlines.append(nextline.translate(str.maketrans("(),", "   ")))

        basismatrix = np.genfromtxt(matrixlines)
        basismatrix = basismatrix[:, 0::2] + 1.0j * basismatrix[:, 1::2]

        if not general and spin_autodetect:
            matdim = basismatrix.shape[0]
            if matdim == nflav:
                spin = True
            elif matdim == norb:
                spin = False
            else:
                raise ValueError("bath basis file malformed: "
                                 f"matrix {i} has {matdim} rows, "
                                 f"expected {norb} or {nflav} rows")


        if basismatrix.shape[0] != basismatrix.shape[1]:
            raise ValueError("bath basis file malformed: "
                             f"matrix {i} not square")
        if matdim is None:
            matdim = basismatrix.shape[0]
            nflav = matdim if spin else 2 * matdim
            norb = nflav // 2
        if basismatrix.shape[0] != matdim:
            raise ValueError("bath basis file malformed: "
                             f"matrix {i} has {basismatrix.shape[0]} rows, "
                             f"expected {matdim} rows")

        # transform into (spin, spin, orb, orb) with spin order up, down
        if spin and w2d_spin_convention:
            basismatrix = np.transpose(
                np.flip(
                    np.reshape(basismatrix,
                               (norb, 2, norb, 2)),
                    axis=(1, 3)
                ),
                (1, 3, 0, 2)
            )
        elif spin:
            basismatrix = np.reshape(basismatrix, (2, norb, 2, norb))
        else:
            bm = np.zeros((2, norb, 2, norb), dtype=basismatrix.dtype)
            bm[0, :, 0, :] = basismatrix[...]
            bm[1, :, 1, :] = basismatrix[...]
            basismatrix = bm

        basis.append(basismatrix)

    basis = np.transpose(np.array(basis), (1, 3, 2, 4, 0))
    return basis, lambdas, vs
