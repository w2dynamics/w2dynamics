from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
import time
import subprocess
import pathlib
import tempfile
from warnings import warn
from numbers import Number
from threading import Lock
from collections import deque
import numpy as np
import scipy
import scipy.optimize as opt
import sys

import w2dyn
import w2dyn.auxiliaries.transform as tf
import w2dyn.auxiliaries.postprocessing as postproc
from w2dyn.auxiliaries.statistics import DistributedSample
from w2dyn.dmft.selfcons import iw_to_tau_fast
from w2dyn.dmft.impurity import (ImpuritySolver,
                                 StatisticalImpurityResult)

import w2dyn.dmft.orbspin as orbspin


class CISolver(ImpuritySolver):
    """CI solver"""
    def __init__(self, config, seed=0, Uw=0, Uw_Mat=0, epsn=0, interactive=False, mpi_comm=None):
        super(CISolver, self).__init__(config)
        self.mpi_comm = mpi_comm
        self.mpi_rank = 0
        if mpi_comm is not None:
            self.mpi_rank = mpi_comm.Get_rank()

        # self.seed = seed + self.mpi_rank
        self.g_diagonal_only = False

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

        # CT-HYB uses a different convention for muimp
        self.muimp = -self.muimp
        self.umatrix = self.problem.interaction.u_matrix.reshape(
                             self.problem.nflavours, self.problem.nflavours,
                             self.problem.nflavours, self.problem.nflavours)

    def solve(self, iter_no, step_cache=None, write_ham=True, only_write_ham=True):
        sys.stdout.flush()
        sys.stderr.flush()

        def log(*args, **kwargs):
            print(time.strftime("%y-%m-%d %H:%M:%S"), *args, **kwargs)

        time_qmc = time.perf_counter()

        # Read potentially passed-in values for initialization of
        # energies and hybridization coefficients
        firstrun = True
        if type(step_cache) is list and len(step_cache) > 0:
            firstrun = False
            oldcache = step_cache.pop(0)
            oldbath = oldcache["bath"]
            log(f"Passing saved old config on to fitting procedure:")
            log(f"{oldbath=}")
        else:
            oldcache = {}
            oldbath = None

        newcache = {}
        if type(step_cache) is list:
            step_cache.append(newcache)

        # transposition and flip by convention
        fiw = np.transpose(np.reshape(self.fiw, (self.problem.nflavours, self.problem.nflavours, -1)), (1, 0, 2))[:, :, ::-1]
        ftau = np.reshape(self.ftau, (self.problem.nflavours, self.problem.nflavours, -1))
        omega = 2.0 * np.pi * (np.arange(-fiw.shape[-1]//2, fiw.shape[-1]//2) + 0.5) / self.problem.beta

        fiw_pos = fiw[:, :, fiw.shape[-1]//2:]
        omega_pos = omega[omega.size//2:]

        def fcmp(energies, xmat, omega):
            return np.sum(xmat[:, :, :, None] / (1.0j * omega[None, None, None, :] - energies[:, None, None, None]), axis=0)

        npoles = self.config["CI"]["npoles"]
        fit_method = self.config["CI"]["bath_fit_method"]

        def spin_symmetrize(fiw):
            fiw = np.reshape(fiw, (self.problem.norbitals, 2,
                                   self.problem.norbitals, 2,
                                   -1))
            if not (np.allclose(fiw[:, 0, :, 1, :], 0.)
                    and np.allclose(fiw[:, 1, :, 0, :], 0.)):
                log("WARNING: off-diagonal entries in fiw ignored due "
                    "to assume_spin_symmetry")
            if not np.allclose(fiw[:, 0, :, 0, :], fiw[:, 1, :, 1, :]):
                log("WARNING: significant differences between spin diagonals "
                    "of fiw ignored due to assume_spin_symmetry")
            fiw = np.mean(fiw, axis=(1, 3))
            return fiw

        if self.config["CI"]["assume_spin_symmetry"]:
            fiw_pos = spin_symmetrize(fiw_pos)

        if iter_no == 0 and self.config["CI"]["initial_bath"] is not None:
            energies = np.array([float(x) for x in self.config["CI"]["initial_bath"]])
            energies, vs = np.split(energies, (energies.size // (self.problem.nflavours + 1),))
            vs = np.reshape(vs, (energies.size, self.problem.nflavours))
            bathparams, oldbath = CISolver.bathparams_from_energies_vs(energies, vs), (energies, vs)
            log(f"Using bath from configuration file: {bathparams=}")
        else:
            log(f"Fitting bath ({fit_method=})", flush=True)

            if fit_method == 'adapol':
                bathparams, oldbath = self.bath_fit_adapol(npoles, omega_pos, fiw_pos)

            elif fit_method == 'sdr_bfgs':
                bathparams, oldbath = self.bath_fit_sdr(npoles, omega_pos, fiw_pos, oldbath)
                log(f"After SDR optimization: whole range tot. mismatch {deltaobj_total(omega, fiw, bathparams)}")
                log(f"After SDR optimization: whole range max mismatch {deltaobj_maxnorm(omega, fiw, bathparams)}", flush=True)

                energies = np.array([e for e, _ in bathparams])
                vs = np.reshape([v for _, vs in bathparams for v in vs], (energies.size, -1))

                bathparams, _ = self.bath_fit_bfgs(energies.size, omega_pos, fiw_pos, (energies, vs), sort=False)
                log(f"After BFGS optimization: whole range tot. mismatch {deltaobj_total(omega, fiw, bathparams)}")
                log(f"After BFGS optimization: whole range max mismatch {deltaobj_maxnorm(omega, fiw, bathparams)}", flush=True)

                energies, xmats = energies_xmats_from_bathparams(bathparams)
                energies = np.reshape(energies, [-1, self.problem.nflavours])
                # assert np.allclose(np.mean(energies, axis=1)[:, None], energies)
                # FIXME: oldbath for non-diagonal case
                xdiags = np.sum(np.reshape(np.diagonal(xmats, 0, 1, 2), (-1, self.problem.nflavours, self.problem.nflavours)), axis=1)
                oldbath = (np.mean(energies, axis=1), xdiags)

            elif fit_method == 'sdr':
                bathparams, oldbath = self.bath_fit_sdr(npoles, omega_pos, fiw_pos, oldbath)

            elif fit_method == 'bfgs':
                bathparams, oldbath = self.bath_fit_bfgs(npoles, omega_pos, fiw_pos, oldbath)

        newcache["bath"] = oldbath

        if self.config["CI"]["assume_spin_symmetry"]:
            new_bps = []
            for energy, vs in bathparams:
                newvs = np.zeros((2 * vs.size,), dtype=vs.dtype)
                newvs[0::2] = vs
                new_bps.append((energy, newvs))
                newvs = np.zeros((2 * vs.size,), dtype=vs.dtype)
                newvs[1::2] = vs
                new_bps.append((energy, newvs))
            bathparams = new_bps


        # WRITE FIT RESULTS TO RESULT DICT
        result = {}
        result["fiw-fit-error-rss"] = deltaobj_total(omega, fiw, bathparams)
        result["fiw-fit-error-max"] = deltaobj_maxnorm(omega, fiw, bathparams)

        log(f"After optimization: whole range tot. mismatch {result['fiw-fit-error-rss']}")
        log(f"After optimization: whole range max mismatch {result['fiw-fit-error-max']}", flush=True)

        result["fiw-fit"], result["ftau-fit"] = [
            np.reshape(
                fit,
                (self.problem.norbitals, 2,
                 self.problem.norbitals, 2,
                 lastdim)
            )
            for fit, lastdim in
            zip(
                tf.hybr_from_sites(
                    [energy
                     for energy, _ in bathparams],
                    np.vstack([hybvec
                               for _, hybvec in bathparams]),
                    omega,
                    np.linspace(0, self.problem.beta, self.problem.nftau)
                ),
                (omega.size, self.problem.nftau)
            )
        ]

        self.problem.g0inviw = (1.0j * (omega[:, None, None, None, None]
                                         * np.eye(self.problem.norbitals)[None, :, None, :, None]
                                         * np.eye(2)[None, None, :, None, :])
                                - self.muimp[None, :, :, :, :]
                                - np.transpose(result["fiw-fit"], (4, 0, 1, 2, 3)))

        result["g0inviw-imp"] = np.transpose(self.problem.g0inviw, (1, 2, 3, 4, 0))
        result["g0iw-imp"] = np.transpose(orbspin.invert(self.problem.g0inviw), (1, 2, 3, 4, 0))
        result["fiw-bath-energies"] = np.array([e for e, _ in bathparams])
        result["fiw-bath-hybvecs"] = np.stack([np.reshape(v, (self.problem.norbitals, 2))
                                               for _, v in bathparams], axis=0)

        # MAKE HAMILTONIAN
        nf = self.problem.nflavours + len(bathparams)  # number of single-fermion states
        H = Operator(nf=nf)

        muimp_f = np.reshape(self.muimp, (self.muimp.shape[0] * self.muimp.shape[1],
                                          self.muimp.shape[2] * self.muimp.shape[3]))
        quadimpcount, quartimpcount, quadbathcount, hybcount = 0, 0, 0, 0
        for flav1 in range(self.problem.nflavours):
            for flav2 in range(self.problem.nflavours):
                if not np.allclose((coeff := muimp_f[flav1, flav2]), 0):
                    H += np.real_if_close(coeff) * Operator.op_cdag(flav1, nf=nf) @ Operator.op_c(flav2, nf=nf)
                    quadimpcount += 1
                for flav3 in range(self.problem.nflavours):
                    for flav4 in range(self.problem.nflavours):
                        if not np.allclose((coeff := self.umatrix[flav1, flav2, flav3, flav4]), 0):
                            H += 0.5 * np.real_if_close(coeff) * (
                                Operator.op_cdag(flav1, nf=nf) @ Operator.op_cdag(flav2, nf=nf)
                                @ Operator.op_c(flav4, nf=nf) @ Operator.op_c(flav3, nf=nf)
                            )
                            quartimpcount += 1

        bathoffset = self.problem.nflavours
        for energy, vvec in bathparams:
            if not np.allclose(energy, 0):
                H += np.real(energy) * Operator.op_cdag(bathoffset, nf=nf) @ Operator.op_c(bathoffset, nf=nf)
                quadbathcount += 1
            for i, v in enumerate(vvec):
                if not np.allclose(v, 0):
                    H += np.real_if_close(v) * Operator.op_cdag(bathoffset, nf=nf) @ Operator.op_c(i, nf=nf)
                    H += np.real_if_close(v) * Operator.op_cdag(i, nf=nf) @ Operator.op_c(bathoffset, nf=nf)
                    hybcount += 2
            bathoffset += 1

        log(f"Elements added to Hamiltonian: {quadimpcount} bilin. imp., {quartimpcount} quadrilin. imp, {quadbathcount} quadratic bath, {hybcount} bilin. hyb.")

        nbathst = nf - self.problem.nflavours
        log(f"Number of fermionic single-particle states: {nf}")
        log(f"Of which bath states: {nbathst}")

        if write_ham:
            log(f"Writing Hamiltonian to file: {write_ham}")
            H.write_quanty_file(".", write_ham)
            if only_write_ham:
                result = StatisticalImpurityResult(self.problem,
                                                   None, None, None,
                                                   self.g_diagonal_only,
                                                   **result)
                return result

        raise NotImplementedError("not available in this version")

    def bath_fit_sdr(self, npoles, omega, fiw, oldcfg=None):
        if self.g_diagonal_only:
            return self.bath_fit_sdr_diagonal(npoles, omega, fiw, oldcfg)
        else:
            return self.bath_fit_sdr_general(npoles, omega, fiw, oldcfg)


    def bath_fit_sdr_general(self, npoles, omega, fiw, oldcfg=None):
        import cvxpy as cvx

        # Fit bath using semi-definite relaxation algorithm (PRB 101, 035143)

        if (niw_fit := self.config["CI"]["niw_fit"]) is not None and niw_fit < omega.size:
            grid = np.floor(loggrid(omega.size, niw_fit) - 1).astype(np.intp)
            print(f"Using reduced frequency grid for fit, {grid=}")
            omega_fit = omega[grid]
            fiw_fit = fiw[:, :, grid]
        else:
            omega_fit = omega
            fiw_fit = fiw

        # Initialize energies and hybridization coefficients if no
        # values were passed in
        if oldcfg is None:
            energies = np.zeros(npoles, dtype=np.float64)
            if True:
                # Use AAA for energy initialization (cf. PRB 107, 075151)
                poles = aaa_poles(np.trace(fiw), 1.0j*omega, nspmax=(npoles + 1), rtol=0)
                print(f"Using real parts of {poles.size} poles from AAA algo.: {poles=}")
                energies[:poles.size] = np.real(poles[:])
            else:
                omegacenter = (np.amin(omega) + np.amax(omega))/2
                omegaquarterwidth = (np.amax(omega) - np.amin(omega))/4
                energies = np.linspace(omegacenter - omegaquarterwidth, omegacenter + omegaquarterwidth, npoles)
            xval = np.broadcast_to(np.eye(self.problem.nflavours)[None, :, :],
                                   (npoles, self.problem.nflavours, self.problem.nflavours))
        else:
            energies, xval = oldcfg

        print(f"{energies=}")
        print(f"{xval=}")

        def objfun(energies, xval):
            return np.sqrt(np.sum(
                np.abs(fiw_fit
                       - np.sum(xval[:, :, :, None]
                                / (1.0j * omega_fit[None, None, None, :]
                                   - energies[:, None, None, None]), axis=0))**2
            ))

        def objjac(energies, xval):
            iwn = 1.0j * omega_fit

            # delta_fit[iflavor, jflavor, iomega]
            delta_fit = (xval[:, :, :, None] / (iwn[:] - energies[:, None, None, None])).sum(0)
            chi = np.sqrt(np.sum(np.abs(fiw_fit - delta_fit)**2))

            # derivatives wrt epsk
            jac = - np.real(np.sum(np.conj(fiw_fit - delta_fit)[None, :, :, :] * xval[:, :, :, None]/(iwn[:] - energies[:, None, None, None])**2, axis=(1, 2, 3))) / chi

            return jac

        lastobjval = np.inf
        fulldiff_current = np.inf
        i = 0
        miniter = self.config["CI"]["sdr_min_iter"]
        maxiter = self.config["CI"]["sdr_max_iter"]
        converged = 0
        eps = self.config["CI"]["epsbatheigval"]

        # optimize w.r.t. X matrices (corr. to hyb. V) using cvxpy
        psdx = [cvx.Variable((self.problem.nflavours, self.problem.nflavours), PSD=True)
                for _ in range(npoles)]
        constraints = [X == X.T for X in psdx]

        while (converged < 5 and i <= maxiter) or i < miniter:
            i += 1

            # Enforce positive semidefinite inputs
            allxevals, allxevecs = np.linalg.eigh(xval)
            xval = np.einsum('abc,ac,adc->abd', allxevecs, np.maximum(allxevals, 0), np.conj(allxevecs))

            for lastx, xvar in zip(xval, psdx):
                xvar.value = lastx

            xobjective = cvx.Minimize(
                sum(cvx.sum(cvx.abs(fiw_fit[:, :, i] - sum(cvx.multiply(x, 1.0/(1.0j * w - l))
                                                           for x, l in zip(psdx, energies)))**2)
                    for i, w in enumerate(omega_fit))
            )
            xoptprob = cvx.Problem(xobjective, constraints)
            xoptprob.solve()
            print(f"{xoptprob.status=}, {xoptprob.value=}")
            xval = np.array([x.value for x in psdx])

            print(f"SDR iter {i} after X opt: mismatch {objfun(energies, xval)}")
            print(f"{np.sort(np.linalg.eigvalsh(xval).flatten())=}")


            loptres = opt.minimize(objfun,
                                   energies,
                                   args=xval,
                                   jac=objjac,
                                   method="BFGS")
            energies = np.array(loptres.x)
            print(f"{np.sort(energies)=}")
            if (fitchange := loptres.fun - lastobjval) > 0:
                warn("Discrete bath fit got worse in the last iteration",
                     UserWarning, 2)
            print(f"SDR iter {i} after E opt: mismatch {loptres.fun}, change {fitchange}")

            fulldiff_last = fulldiff_current
            fulldiff_current = deltaobj_total(
                omega, fiw,
                bathparams_from_xmats(xval, energies, 0, trim=False)
            )
            print(f"End of SDR iteration {i}: "
                  f"whole range tot. mismatch {fulldiff_current}, "
                  f"change {fulldiff_current - fulldiff_last}")
            maxdiff = deltaobj_maxnorm(omega, fiw, bathparams_from_xmats(xval, energies, eps))
            print(f"End of SDR iteration {i}: whole range max mismatch {maxdiff}")

            if abs(loptres.fun - lastobjval) < self.config["CI"]["epsbathfit"]:
                converged += 1
            else:
                converged = 0
            lastobjval = loptres.fun

        esort = np.argsort(energies)
        energies = energies[esort]
        xval = xval[esort, :, :]

        bathparams = bathparams_from_xmats(xval, energies, eps, log=print)

        return bathparams, (energies, xval)


    def bath_fit_sdr_diagonal(self, npoles, omega, fiw, oldcfg=None, recover_zeros=False):
        import cvxpy as cvx

        # Fit bath using semi-definite relaxation algorithm (PRB 101, 035143)

        if (niw_fit := self.config["CI"]["niw_fit"]) is not None and niw_fit < omega.size:
            grid = np.floor(loggrid(omega.size, niw_fit) - 1).astype(np.intp)
            print(f"Using reduced frequency grid for fit, {grid=}")
            omega_fit = omega[grid]
            fiw_fit = fiw[:, :, grid]
        else:
            omega_fit = omega
            fiw_fit = fiw

        # Initialize energies and hybridization coefficients if no
        # values were passed in
        if oldcfg is None:
            energies = np.zeros(npoles, dtype=np.float64)
            xdiag = np.ones((npoles, self.problem.nflavours), dtype=np.float64)
            if True:
                # Use AAA for energy initialization (cf. PRB 107, 075151)
                poles = aaa_poles(np.trace(fiw), 1.0j*omega, nspmax=(npoles + 1), rtol=0)
                print(f"Using real parts of {poles.size} poles from AAA algo.: {poles=}")
                energies[:poles.size] = np.real(poles[:])

                # Use AAA for energy initialization (cf. PRB 107, 075151)
                # poles = aaa_poles(np.trace(fiw), 1.0j*omega, nspmax=(npoles + 1), rtol=0)
                # xdiag = np.zeros((npoles, self.problem.nflavours), dtype=np.float64)
                # for i, (esub, vsub) in enumerate(zip(np.array_split(energies, self.problem.nflavours),
                #                                      np.array_split(xdiag, self.problem.nflavours))):
                #     poles = aaa_poles(fiw[i, i, :], 1.0j*omega, nspmax=(esub.size + 1), rtol=0)
                #     print(f"Using real parts of {poles.size} poles from AAA algo.: {poles=}")
                #     esub[:] = np.real(poles[:])
                #     vsub[:, :] = 1.0
            else:
                omegacenter = (np.amin(omega) + np.amax(omega))/2
                omegaquarterwidth = (np.amax(omega) - np.amin(omega))/4
                energies = np.linspace(omegacenter - omegaquarterwidth, omegacenter + omegaquarterwidth, npoles)
        else:
            energies, xdiag = oldcfg

        print(f"{energies=}")
        print(f"{xdiag=}")

        fiw_fit = np.transpose(np.diagonal(fiw_fit))

        def objfun(energies, xdiag):
            return np.sqrt(np.sum(
                np.abs(fiw_fit
                       - np.sum(xdiag[:, :, None]
                                / (1.0j * omega_fit[None, None, :]
                                   - energies[:, None, None]), axis=0))**2
            ))

        def objjac(energies, xdiag):
            iwn = 1.0j * omega_fit

            # delta_fit[iflavor, jflavor, iomega]
            delta_fit = (xdiag[:, :, None] / (iwn[:] - energies[:, None, None])).sum(0)
            chi = np.sqrt(np.sum(np.abs(fiw_fit - delta_fit)**2))

            # derivatives wrt epsk
            jac = - np.real(np.sum(np.conj(fiw_fit - delta_fit)[None, :, :] * xdiag[:, :, None]/(iwn[:] - energies[:, None, None])**2, axis=(1, 2))) / chi

            return jac

        def objfun_singlesite(e_xdiag, nsite, energies, xdiag):
            energies[nsite] = e_xdiag[0]
            xdiag[nsite, :] = e_xdiag[1:]
            return objfun(energies, xdiag)

        lastobjval = np.inf
        fulldiff_current = np.inf
        i = 0
        miniter = self.config["CI"]["sdr_min_iter"]
        maxiter = self.config["CI"]["sdr_max_iter"]
        converged = 0
        eps = self.config["CI"]["epsbatheigval"]

        # optimize w.r.t. X matrices (corr. to hyb. V) using cvxpy
        px = cvx.Variable((npoles, self.problem.nflavours), nonneg=True)

        while (converged < 5 and i <= maxiter) or i < miniter:
            i += 1

            # Enforce non-negative inputs
            px.value = np.maximum(0, xdiag)

            # FIXME: try to use cvx.Parameter instead of fiw, w, l
            xobjective = cvx.Minimize(
                cvx.norm(
                    cvx.abs(
                        np.transpose(fiw_fit[:, :])
                        - cvx.matmul(1.0/(1.0j * omega_fit[:, None] - energies[None, :]), px)
                    ),
                    "fro"
                )
            )
            xoptprob = cvx.Problem(xobjective)
            xoptprob.solve()
            print(f"{xoptprob.status=}, {xoptprob.value=}")
            xdiag = px.value

            print(f"SDR iter {i} after X opt: mismatch {objfun(energies, xdiag)}")
            print(f"{xdiag=}")

            # Try to recover from zeros
            if recover_zeros:
                for l in range(energies.size):
                    if np.allclose(0, xdiag[l]):
                        soptret = opt.minimize(objfun_singlesite,
                                               np.concatenate((energies[l, None],
                                                               xdiag[l, :]), axis=0),
                                               args=(l, energies, xdiag))
                        energies[l] = soptret.x[0]
                        xdiag[l] = soptret.x[1:]

            loptres = opt.minimize(objfun,
                                   energies,
                                   args=xdiag,
                                   jac=objjac,
                                   method="BFGS")
            energies = np.array(loptres.x)
            print(f"{np.sort(energies)=}")
            if (fitchange := loptres.fun - lastobjval) > 0:
                warn("Discrete bath fit got worse in the last iteration",
                     UserWarning, 2)
            print(f"SDR iter {i} after E opt: mismatch {loptres.fun}, change {fitchange}")

            fulldiff_last = fulldiff_current
            fulldiff_current = deltaobj_total(
                omega, fiw,
                bathparams_from_xdiags(xdiag, energies, 0, trim=False)
            )
            print(f"End of SDR iteration {i}: "
                  f"whole range tot. mismatch {fulldiff_current}, "
                  f"change {fulldiff_current - fulldiff_last}")
            maxdiff = deltaobj_maxnorm(omega, fiw, bathparams_from_xdiags(xdiag, energies, eps))
            print(f"End of SDR iteration {i}: whole range max mismatch {maxdiff}")

            if abs(loptres.fun - lastobjval) < self.config["CI"]["epsbathfit"]:
                converged += 1
            else:
                converged = 0
            lastobjval = loptres.fun

        esort = np.argsort(energies)
        energies = energies[esort]
        xdiag = xdiag[esort, :]

        bathparams = bathparams_from_xdiags(xdiag, energies, eps, trim=False, log=print)

        return bathparams, (energies, xdiag)


    def bath_fit_bfgs(self, npoles, omega, fiw, oldcfg=None, sort=True):
        # Fit bath using scipy-provided optimization algorithms

        if (niw_fit := self.config["CI"]["niw_fit"]) is not None and niw_fit < omega.size:
            grid = np.floor(loggrid(omega.size, niw_fit) - 1).astype(np.intp)
            print(f"Using reduced frequency grid for fit, {grid=}")
            omega_fit = omega[grid]
            fiw_fit = fiw[:, :, grid]
        else:
            omega_fit = omega
            fiw_fit = fiw

        # Initialize energies and hybridization coefficients if no
        # values were passed in
        if oldcfg is None:
            energies = np.zeros(npoles, dtype=np.float64)
            if True:
                # Use AAA for energy initialization (cf. PRB 107, 075151)
                # poles = aaa_poles(np.trace(fiw), 1.0j*omega, nspmax=(npoles + 1), rtol=0)
                vs = np.zeros((npoles, self.problem.nflavours), dtype=np.float64)
                for i, (esub, vsub) in enumerate(zip(np.array_split(energies, self.problem.nflavours),
                                                     np.array_split(vs, self.problem.nflavours))):
                    poles = aaa_poles(fiw[i, i, :], 1.0j*omega, nspmax=(esub.size + 1), rtol=0)
                    print(f"Using real parts of {poles.size} poles from AAA algo.: {poles=}")
                    esub[:] = np.real(poles[:])
                    vsub[:] = 1.0
            else:
                omegacenter = (np.amin(omega) + np.amax(omega))/2
                omegaquarterwidth = (np.amax(omega) - np.amin(omega))/4
                energies = np.linspace(omegacenter - omegaquarterwidth, omegacenter + omegaquarterwidth, npoles)
                vs = np.ones((npoles, self.problem.nflavours))
        else:
            energies, vs = oldcfg

        print(f"{energies=}")
        print(f"{vs=}")

        optres = opt.minimize(deltaobj_opt_realvs,
                              np.concatenate((energies, vs.flatten())),
                              args=(energies.size, omega_fit, fiw_fit),
                              jac=deltajac_opt_realvs,
                              method='BFGS',
                              options={'maxiter': 10000,
                                       'disp': True})

        energies, vs = np.split(np.array(optres.x), (energies.size,))
        vs = np.reshape(vs, (energies.size, self.problem.nflavours))

        if sort:
            bathparams = CISolver.bathparams_from_energies_vs(energies, vs)
        else:
            print(f"{energies=}")
            print(f"{vs=}")
            bathparams = [(e, v) for e, v in zip(energies, vs)]

        print(f"End of BFGS optimization: tot. mismatch {deltaobj_total(omega_fit, fiw_fit, bathparams)} = {optres.fun}")
        print(f"End of BFGS optimization: whole range tot. mismatch {deltaobj_total(omega, fiw, bathparams)}")
        print(f"End of BFGS optimization: whole range max mismatch {deltaobj_maxnorm(omega, fiw, bathparams)}")

        return bathparams, (energies, vs)

    @staticmethod
    def bathparams_from_energies_vs(energies, vs):
        esort = np.argsort(energies)
        energies = energies[esort]
        vs = vs[esort, :]
        print(f"{energies=}")
        print(f"{vs=}")

        bathparams_down = [(e, v) for e, v in zip(energies, vs)
                           if np.sum(np.abs(v[0::2])**2) >= np.sum(np.abs(v[1::2])**2)]
        bathparams_up = [(e, v) for e, v in zip(energies, vs)
                         if np.sum(np.abs(v[0::2])**2) < np.sum(np.abs(v[1::2])**2)]
        if len(bathparams_down) != len(bathparams_up):
            print(f"WARNING: Bath level spins seem imbalanced {len(bathparams_down)=} {len(bathparams_up)=}",
                  file=sys.stderr, flush=True)
        else:
            print("Bath level spins balanced")
        bathparams = []
        for i in range(max(len(bathparams_down), len(bathparams_up))):
            if i < len(bathparams_down):
                bathparams.append(bathparams_down[i])
            if i < len(bathparams_up):
                bathparams.append(bathparams_up[i])

        return bathparams

    def bath_fit_adapol(self, npoles, omega, fiw, oldcfg=None):
        import adapol
        energies, vs, error, func = adapol.hybfit(fiw.transpose(2, 0, 1),
                                                  1.0j * omega,
                                                  Np=npoles,
                                                  svdtol=1e-10,
                                                  verbose=True)

        esort = np.argsort(energies)
        energies = energies[esort]
        vs = vs[esort, :]
        print(f"{energies=}")
        print(f"{vs=}")

        bathparams_down = [(e, v) for e, v in zip(energies, vs)
                           if np.sum(np.abs(v[0::2])**2) >= np.sum(np.abs(v[1::2])**2)]
        bathparams_up = [(e, v) for e, v in zip(energies, vs)
                         if np.sum(np.abs(v[0::2])**2) < np.sum(np.abs(v[1::2])**2)]
        if len(bathparams_down) != len(bathparams_up):
            print(f"WARNING: Bath level spins seem imbalanced {len(bathparams_down)=} {len(bathparams_up)=}",
                  file=sys.stderr, flush=True)
        else:
            print("Bath level spins balanced")
        bathparams = []
        for i in range(max(len(bathparams_down), len(bathparams_up))):
            if i < len(bathparams_down):
                bathparams.append(bathparams_down[i])
            if i < len(bathparams_up):
                bathparams.append(bathparams_up[i])

        print(f"End of adapol optimization: reported error {error}")
        print(f"End of adapol optimization: whole range tot. mismatch {deltaobj_total(omega, fiw, bathparams)}")
        print(f"End of adapol optimization: whole range max mismatch {deltaobj_maxnorm(omega, fiw, bathparams)}")

        return bathparams, (energies, vs)


# Fit function fnvals given at (complex) sample points omega using a
# rational approximation of numerator and denominator degrees at most
# nspmax, less if relative tolerance rtol reached earlier.
#
# Returns:
# - Poles of the rational approximation
#
# Reference:
# - Y. Nakatsukasa, O. Sète, and L. N. Trefethen, The AAA Algorithm
#   for Rational Approximation, SIAM J. Sci. Comput. 40, A1494 (2018).
def aaa_poles(fnvals, omega, nspmax=100, rtol=1e-8):
    selectedsps = np.zeros(fnvals.shape, dtype=np.bool_)
    fnrat = np.full_like(fnvals, np.mean(fnvals))
    for nsp in range(nspmax):
        nextspindex = np.argmax(np.abs(fnrat - fnvals))
        if nsp > 0 and np.all(np.abs((fnrat[nextspindex]
                                      - fnvals[nextspindex])
                                     /fnvals[nextspindex]) < rtol):
            break
        selectedsps[nextspindex] = True
        C = 1.0/(omega[np.logical_not(selectedsps)][:, np.newaxis]
                 - omega[selectedsps][np.newaxis, :])
        w = np.linalg.svd(
            fnvals[np.logical_not(selectedsps)][:, np.newaxis] * C
            - C * fnvals[selectedsps][np.newaxis, :],
            full_matrices=False
        )[-1][-1, :]
        N = np.dot(C, w*fnvals[selectedsps])
        D = np.dot(C, w)
        fnrat[np.logical_not(selectedsps)] = N/D
        fnrat[selectedsps] = fnvals[selectedsps]
        # print(f"in AAA iter for {nsp=}")
    polemat = np.diagflat(np.concatenate(([0.0j], omega[selectedsps]), axis=0))
    polemat[0, 1:] = w
    polemat[1:, 0] = 1
    poles = scipy.linalg.eig(
        polemat,
        np.diagflat(np.concatenate(([0], np.ones((w.size))), axis=0))
    )[0]
    poles = poles[np.isfinite(poles)]
    # determine residuals "numerically" as done in the reference...
    # polecrosses = 1e-6*np.exp(0.5j*np.pi*np.array(range(4)))[np.newaxis, :] + poles[:, np.newaxis]
    # Cpc = 1.0/(np.reshape(polecrosses, (-1,))[:, np.newaxis] - omega[selectedsps][np.newaxis, :])
    # resid = np.sum(np.reshape(np.dot(Cpc, w*fnvals[selectedsps]) / np.dot(Cpc, w), (-1, 4)), axis=1)/4
    return poles


def loggrid(nfreq_tot, nfreq_red):
    def grid(x, t):
        return 1 - (1 - t**(x - 1))/(t - 1)

    x = opt.brentq(lambda x: (grid(nfreq_red, x) - nfreq_tot), 1.00001, 10)
    return grid(np.arange(nfreq_red) + 1, x)


def bathparams_from_xmats(xval, energies, eps, trim=True, log=lambda x, *a, **k: None):
    allxevals, allxevecs = np.linalg.eigh(xval)
    log(f"{allxevals=}")

    bathparams_down = []
    bathparams_up = []

    for xevals, xunit, energy in zip(allxevals, allxevecs, energies):
        # negative eigenvalue zeroing intentional
        xevals[xevals < eps] = 0
        log(f"{np.count_nonzero(xevals)} eigenvalues of {xevals.size} for energy {energy}")

        # pull diagonal matrix into unitary matrix such that x = umat @ umat.T
        umat = xunit @ np.diagflat(np.sqrt(xevals))
        if trim and np.allclose(umat, 0):
            log(f"No bath sites with energy {energy} added")
            continue

        # each column is one v vector with one entry per impurity flavor
        #for v in np.unstack(umat, axis=1)
        for i, v in enumerate(np.hsplit(umat, list(range(1, umat.shape[1])))):
            if not (trim and np.allclose(v, 0)):
                if np.allclose(np.sum(np.abs(v[0::2])**2), np.sum(np.abs(v[1::2])**2)):
                    # in order to still assign them equally if
                    # e.g. zeros are not trimmed or strange cases
                    if i % 2 == 0:
                        bathparams_down.append((energy, np.reshape(v, (-1,))))
                    else:
                        bathparams_up.append((energy, np.reshape(v, (-1,))))
                elif np.sum(np.abs(v[0::2])**2) >= np.sum(np.abs(v[1::2])**2):
                    bathparams_down.append((energy, np.reshape(v, (-1,))))
                else:
                    bathparams_up.append((energy, np.reshape(v, (-1,))))

    if len(bathparams_down) != len(bathparams_up):
        log(f"WARNING: Bath level spins seem imbalanced {len(bathparams_down)=} {len(bathparams_up)=}",
            file=sys.stderr, flush=True)
    else:
        log("Bath level spins balanced")

    bathparams = []
    for i in range(max(len(bathparams_down), len(bathparams_up))):
        if i < len(bathparams_down):
            bathparams.append(bathparams_down[i])
        if i < len(bathparams_up):
            bathparams.append(bathparams_up[i])

    log(f"{bathparams=}")
    return bathparams


def bathparams_from_xdiags(xdmat, energies, eps, trim=True, log=lambda x, *a, **k: None):
    bathparams_down = []
    bathparams_up = []

    for xdiags, energy in zip(xdmat, energies):
        # negative eigenvalue zeroing intentional
        xdiags[xdiags < eps] = 0
        log(f"{np.count_nonzero(xdiags)} diagonals of {xdiags.size} for energy {energy}")

        if trim and np.allclose(xdiags, 0):
            log(f"No bath sites with energy {energy} added")
            continue

        # pull diagonal matrix into unitary matrix such that x = umat @ umat.T
        xdiags = np.diagflat(np.sqrt(xdiags))

        # each column is one v vector with one entry per impurity flavor
        #for v in np.unstack(umat, axis=1)
        for i, v in enumerate(np.hsplit(xdiags, list(range(1, xdiags.shape[1])))):
            if not (trim and np.allclose(v, 0)):
                if np.allclose(np.sum(np.abs(v[0::2])**2), np.sum(np.abs(v[1::2])**2)):
                    # in order to still assign them equally if
                    # e.g. zeros are not trimmed or strange cases
                    if i % 2 == 0:
                        bathparams_down.append((energy, np.reshape(v, (-1,))))
                    else:
                        bathparams_up.append((energy, np.reshape(v, (-1,))))
                elif np.sum(np.abs(v[0::2])**2) >= np.sum(np.abs(v[1::2])**2):
                    bathparams_down.append((energy, np.reshape(v, (-1,))))
                else:
                    bathparams_up.append((energy, np.reshape(v, (-1,))))

    if len(bathparams_down) != len(bathparams_up):
        log(f"WARNING: Bath level spins seem imbalanced {len(bathparams_down)=} {len(bathparams_up)=}",
            file=sys.stderr, flush=True)
    else:
        log("Bath level spins balanced")

    bathparams = []
    for i in range(max(len(bathparams_down), len(bathparams_up))):
        if i < len(bathparams_down):
            bathparams.append(bathparams_down[i])
        if i < len(bathparams_up):
            bathparams.append(bathparams_up[i])

    log(f"{bathparams=}")
    return bathparams


def energies_xmats_from_bathparams(bathparams):
    energies = np.zeros(len(bathparams), dtype=np.float64)
    xmats = np.zeros((len(bathparams), bathparams[0][1].size, bathparams[0][1].size), dtype=np.float64)
    for i, (energy, v) in enumerate(bathparams):
        energies[i] = energy
        xmats[i] = np.outer(v, v)
    return energies, xmats


def epsk_vki_from_bathparams(bathparams):
    epsk = np.zeros(len(bathparams), dtype=np.float64)
    vki = np.zeros((len(bathparams), bathparams[0][1].size), dtype=np.float64)
    for k, (eps, vi) in enumerate(bathparams):
        epsk[k] = eps
        vki[k, :] = vi
    return epsk, vki


def deltajac_opt_realvs(energies_vs_array, n, omega, fiw):
    # energies_vs_array: 1-dim array of energies followed by vs
    # n_omega_fiw_tuple: tuple of (n, omega, fiw)
    # n: number of energies (and of vs)
    # omega: matsubara frequencies of fiw
    # fiw: hybridization function (nflavor, nflavor, nmatsubara)
    epsk = energies_vs_array[:n]
    vki = np.reshape(energies_vs_array[n:],
                    (n, -1))
    jac = np.zeros(energies_vs_array.size, dtype=np.float64)
    jacv = np.zeros_like(vki)
    iwn = 1.0j * omega
    vijk = vki.T.conj()[:,None,:] * vki.T[None,:,:]

    # delta_fit[iflavor, jflavor, iomega]
    delta_fit = (vijk[:,:,None,:] / (iwn[:,None] - epsk)).sum(-1)
    chi = np.sqrt(np.sum(np.abs(fiw - delta_fit)**2))

    # derivatives wrt epsk
    jac[:n] = - np.real(np.sum(np.conj(fiw - delta_fit)[:, :, :, None] * vijk[:, :, None, :]/(iwn[:, None] - epsk)**2, axis=(0, 1, 2))) / chi

    # wrt V
    for i in range(delta_fit.shape[0]):
        for j in range(delta_fit.shape[1]):
            if i == j:
                jacv[:, i] += -2 * np.real(np.sum(np.conj(fiw[i, i] - delta_fit[i, i])[:, None] * vki[None, :, i] / (iwn[:, None] - epsk), axis=0)) / chi
            else:
                jacv[:, i] += - np.real(np.sum(np.conj(fiw[i, j] - delta_fit[i, j])[:, None] * vki[None, :, j] / (iwn[:, None] - epsk), axis=0)) / chi
                jacv[:, j] += - np.real(np.sum(np.conj(fiw[i, j] - delta_fit[i, j])[:, None] * vki[None, :, i] / (iwn[:, None] - epsk), axis=0)) / chi
    jac[n:] = np.reshape(jacv, (-1))

    return jac


def deltaobj_opt_realvs(energies_vs_array, n, omega, fiw):
    # energies_vs_array: 1-dim array of energies followed by vs
    # n: number of energies (and of vs)
    # omega: matsubara frequencies of fiw
    # fiw: hybridization function (nflavor, nflavor, nmatsubara)
    return np.sqrt(np.sum(np.abs(
        fiw
        - tf.hybr_from_sites(
            energies_vs_array[:n],
            np.reshape(energies_vs_array[n:],
                       (n, -1)),
            omega,
            np.array((0,))  # for unused F(tau)
        )[0]
    )**2))


def deltaobj_total(omega, fiw, bathparams):
    # omega: 1-dim array of matsubara frequencies of fiw
    # fiw: 3-dim array of hybridization function, dims. (nflavor, nflavor, nmatsubara)
    # bathparams: seq. of (energy, hybvec)
    return deltaobj(*epsk_vki_from_bathparams(bathparams), omega, fiw)


def deltaobj(epsk, vki, omega, fiw):
    iwn = 1.0j * omega

    vijk = vki.T.conj()[:,None,:] * vki.T[None,:,:]
    return np.sqrt(np.sum(np.abs(fiw - (vijk[:,:,None,:] / (iwn[:,None] - epsk)).sum(-1))**2))


def deltaobj_maxnorm(omega, fiw, bathparams):
    # omega: 1-dim array of matsubara frequencies of fiw
    # fiw: 3-dim array of hybridization function, dims. (nflavor, nflavor, nmatsubara)
    # bathparams: seq. of (energy, hybvec)
    return np.amax(np.abs(
        fiw
        - tf.hybr_from_sites(
            [energy
             for energy, _ in bathparams],
            np.vstack([hybvec
                       for _, hybvec in bathparams]),
            omega,
            np.linspace(0, 1, 2)  # for unused F(tau)
        )[0]
    ))


### Quanty interface
_last_nf = None
def _auto_nf(nf=None):
    global _last_nf
    if nf is not None:
        assert nf == int(nf) > 0, "nf must be a positive integer"
        _last_nf = int(nf)
    return _last_nf


QUANTYREPL = r"""
function main()
   io.stdout:write("\nREPL READY\n")
   io.stdout:flush()
   io.stderr:flush()

   while true do

      local instr = ""
      while true do
	 local nextline = io.stdin:read("*L")
	 if nextline == nil then
	    do return end
	 elseif nextline == "EOF\n" then
	    break
	 else
	    instr = instr .. nextline
	 end
      end

      local fn = load(instr, instr, "t")
      instr = ""
      a, b = pcall(fn)
      if not a then
	 io.stderr:write(b)
      end
      io.stdout:write("\nEOF\n")
      io.stdout:flush()
      io.stderr:write("\nEOF\n")
      io.stderr:flush()
      collectgarbage()

   end
end

main()
"""


class _QuantyManager:
    _singleton = None

    @classmethod
    def get_singleton(cls):
        if cls._singleton is None:
            cls._singleton = _QuantyManager()
        return cls._singleton

    def __init__(self):
        self._qexe = pathlib.Path("Quanty")
        if not self._qexe.exists():
            self._qexe = pathlib.Path(w2dyn.__file__).parent.parent / "Quanty"
        self._qexe = self._qexe.absolute()
        with tempfile.NamedTemporaryFile(mode="w+",
                                         prefix="repl",
                                         suffix="lua") as f:
            f.write(QUANTYREPL)
            f.flush()
            self._process = subprocess.Popen(
                [
                    self._qexe,
                    f.name
                ],
                text=True,
                bufsize=1,
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE
            )

            self._quanty_preamble = ""
            line = self._process.stdout.readline()
            while line != "REPL READY\n":
                self._quanty_preamble += line
                line = self._process.stdout.readline()

        self._varcounter = 0

        self._process_lock = Lock()
        self._defer_destroy_deque_lock = Lock()
        self._defer_destroy_deque = deque()

    def new_var(self, obj):
        self._varcounter += 1
        varname = f"pyautovar{self._varcounter}"
        return varname

    def destroy_var(self, obj):
        if not self._process_lock.acquire(blocking=False):
            if not self._defer_destroy_deque_lock.acquire(blocking=False):
                raise RuntimeError("destroy_var can't acquire either lock")
            self._defer_destroy_deque.append(obj._qvar)
            self._defer_destroy_deque_lock.release()
        else:
            self._process.stdin.write(f"{obj._qvar} = nil\n")
            self._quanty_finalize_command()
            self._quanty_get_output()
            self._quanty_check_error()
            self._process_lock.release()

    def _destroy_deferred(self):
        while True:
            with self._defer_destroy_deque_lock:
                try:
                    nextvar = self._defer_destroy_deque.popleft()
                except IndexError:
                    break
            self.run_command(f"{nextvar} = nil\n")

    def _quanty_finalize_command(self):
        """Write input terminator."""
        self._process.stdin.write("\nEOF\n")
        self._process.stdin.flush()

    def _quanty_check_error(self):
        error = ''
        line = self._process.stderr.readline()
        while line not in ("EOF\n", ""):
            error += line
            line = self._process.stderr.readline()

        error = error.strip()
        if len(error) > 0:
            raise RuntimeError(
                f"Got error message from Quanty: {error}"
            )

    def _quanty_get_output(self):
        output = ''
        line = self._process.stdout.readline()
        while line not in ("EOF\n", ""):
            output += line
            line = self._process.stdout.readline()

        return output

    def var_from_quanty_file(self, file, target):
        with self._process_lock:
            with open(file, "r") as f:
                # write only our variable name instead of that in the
                # first line, then just pipe the rest to quanty in
                # ~page-sized chunks
                firstline = f.readline()
                self._process.stdin.write(f"{target._qvar}")
                self._process.stdin.write(firstline[firstline.index('='):])
                self._process.stdin.flush()
                while len(chunk := f.read(4096)) > 0:
                    self._process.stdin.write(chunk)
                    self._process.stdin.flush()
            self._quanty_finalize_command()
            output = self._quanty_get_output()
            self._quanty_check_error()
            return output
        self._destroy_deferred()

    def run_command(self, command, *args):
        with self._process_lock:
            self._process.stdin.write(command.format(*[arg._qvar for arg in args]))
            self._quanty_finalize_command()
            output = self._quanty_get_output()
            self._quanty_check_error()
        self._destroy_deferred()
        return output

    def stream_command_input(self, prelude_inputgen, *args):
        with self._process_lock:
            command, input_generator = prelude_inputgen
            self._process.stdin.write(command.format(*[arg._qvar for arg in args]))
            self._quanty_finalize_command()

            for inputstr in input_generator():
                self._process.stdin.write(inputstr)

            output = self._quanty_get_output()
            self._quanty_check_error()
        self._destroy_deferred()
        return output

    def parse_command_output(self, cmd_linefn_finfn, *args):
        with self._process_lock:
            command, process_line, finish = cmd_linefn_finfn
            self._process.stdin.write(command.format(*[arg._qvar for arg in args]))
            self._quanty_finalize_command()

            line = self._process.stdout.readline()
            while line != "EOF\n":
                process_line = process_line(line)
                line = self._process.stdout.readline()

            self._quanty_check_error()
        self._destroy_deferred()
        return finish()

    def _set_number(self, varname, number):
        if number.imag == 0:
            return self.run_command(
                f"{varname} = {number.real}\n"
            )
        else:
            return self.run_command(
                f"{varname} = {number.real} + I * {number.imag}\n"
            )

    def _get_number(self, varname):
        res = _qm().parse_command_output(
            (f"pyqmtempvar1 = {varname} + I * 0\n"
             "print(string.format('%25.17e', pyqmtempvar1.real))\n"
             "print(string.format('%25.17e', pyqmtempvar1.imag))\n"
             "pyqmtempvar1 = nil\n",
             *parse_nums()),
        )
        return (res[0] if res[1] == 0 else res[0] + 1.0j * res[1])

    def _execute_luafile(self, filename):
        return self.run_command(f"dofile('{filename}')\n")

    def __del__(self):
        if self.__class__._singleton is self:
            self.__class__._singleton = None
        try:
            self._process.stdin.write("os.exit(true)\nEOF\n")
            self._process.stdin.flush()
            try:
                self._process.communicate(timeout=5)
            except subprocess.TimeoutExpired:
                self._process.kill()
        except BrokenPipeError:
            self._process.kill()

def _qm():
    return _QuantyManager.get_singleton()


class _QuantyObject:
    def __init__(self, *args, **kwargs):
        self._qvar = _qm().new_var(self)
        super().__init__(*args, **kwargs)

    def clear(self):
        _qm().destroy_var(self)

    def write_quanty_file(self, dir, varname=None):
        if varname is None:
            return _qm().run_command(
                f"io.put('{dir}', '{{}}')\n",
                self
            )
        else:
            tempname = "pytempvar1" if varname.strip() != "pytempvar1" else "pytempvar2"
            return _qm().run_command(
                f"{tempname} = {varname}\n"
                f"{varname} = {{}}\n"
                f"io.put('{dir}', '{varname}')\n"
                f"{varname} = {tempname}\n"
                f"{tempname} = nil\n",
                self
            )

    def __del__(self):
        try:
            self.clear()
        except (TypeError, BrokenPipeError):
            pass


class Operator(_QuantyObject):
    def __init__(self, nf=None, *args, **kwargs):
        """
        nf: Number of fermionic single-particle states
        """
        nf = _auto_nf(nf)
        super().__init__(*args, **kwargs)
        self.nf = nf
        _qm().run_command(
            f"{{}} = NewOperator({self.nf}, 0, {{{{{{{{0}}}}}}}})\n",
            self
        )

    def get_nf(self):
        return self.nf

    def print_from_quanty(self):
        print(_qm().run_command("{}.Print()\n", self))

    @staticmethod
    def op_c(index, nf=None):
        nf = _auto_nf(nf)
        res = Operator(nf)
        _qm().run_command(
            f"{{}} = NewOperator('An', {nf}, {index})\n",
            res
        )
        return res

    @staticmethod
    def op_cdag(index, nf=None):
        nf = _auto_nf(nf)
        res = Operator(nf)
        _qm().run_command(
            f"{{}} = NewOperator('Cr', {nf}, {index})\n",
            res
        )
        return res

    @staticmethod
    def op_number(index, nf=None):
        nf = _auto_nf(nf)
        res = Operator(nf)
        _qm().run_command(
            f"{{}} = NewOperator('Number', {nf}, {index}, {index})\n",
            res
        )
        return res

    def __add__(self, other):
        if isinstance(other, Operator):
            res = self.__class__(self.nf)
            _qm().run_command("{} = {} + {}\n", res, self, other)
            return res
        elif isinstance(other, Number):
            res = self.__class__(self.nf)
            if other.imag == 0.0:
                _qm().run_command(f"{{}} = {{}} + {other.real}\n", res, self)
            else:
                _qm().run_command(f"{{}} = {{}} + ({other.real} +  I * {other.imag})\n",
                                  res, self)
            return res
        else:
            return NotImplemented

    def __radd__(self, other):
        return self.__add__(other)

    def __sub__(self, other):
        if isinstance(other, Operator):
            res = self.__class__(self.nf)
            _qm().run_command("{} = {} - {}\n", res, self, other)
            return res
        elif isinstance(other, Number):
            res = self.__class__(self.nf)
            if other.imag == 0.0:
                _qm().run_command(f"{{}} = {{}} - {other.real}\n", res, self)
            else:
                _qm().run_command(f"{{}} = {{}} - ({other.real} +  I * {other.imag})\n",
                                  res, self)
            return res
        else:
            return NotImplemented

    def __rsub__(self, other):
        if isinstance(other, Number):
            res = self.__class__(self.nf)
            if other.imag == 0.0:
                _qm().run_command(f"{{}} = {other.real} - {{}}\n", res, self)
            else:
                _qm().run_command(f"{{}} = ({other.real} +  I * {other.imag}) - {{}}\n",
                                  res, self)
            return res
        else:
            return NotImplemented

    def __mul__(self, other):
        if isinstance(other, Number):
            res = self.__class__(self.nf)
            if other.imag == 0.0:
                _qm().run_command(f"{{}} = {other.real} * {{}}\n", res, self)
            else:
                _qm().run_command(f"{{}} = ({other.real} + I * {other.imag}) * {{}}\n",
                                  res, self)
            return res
        else:
            return NotImplemented

    def __truediv__(self, other):
        if isinstance(other, Number):
            res = self.__class__(self.nf)
            if other.imag == 0.0:
                _qm().run_command(f"{{}} = {{}} / {other.real}\n", res, self)
            else:
                _qm().run_command(f"{{}} = {{}} / ({other.real} + I * {other.imag})\n",
                                  res, self)
            return res
        else:
            return NotImplemented

    def __rmul__(self, other):
        return self.__mul__(other)

    def __neg__(self):
        return -1.0 * self

    def __pos__(self):
        return self

    def __matmul__(self, other):
        if isinstance(other, Operator):
            res = self.__class__(self.nf)
            _qm().run_command("{} = {} * {}\n", res, self, other)
            return res
        else:
            return NotImplemented

    def __deepcopy__(self, memo=None):
        res = self.__class__(self.nf)
        _qm().run_command("{} = Clone({})\n", res, self)
        return res
