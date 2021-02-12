#!/usr/bin/env python
"""Program for the DMFT self-consistency loop"""
from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
import os
import os.path
import errno
import glob
import random
from subprocess import Popen, PIPE
import sys
import time
import optparse
import warnings

import numpy as np

import w2dyn.auxiliaries as aux
from w2dyn.auxiliaries import transform as tf
from w2dyn.auxiliaries import hdfout
from w2dyn.auxiliaries import config
from w2dyn.auxiliaries.utilities import diagonal_covariance

from w2dyn.dmft import impurity
from w2dyn.dmft import lattice
from w2dyn.dmft import atoms
from w2dyn.dmft import interaction
from w2dyn.dmft import selfcons
from w2dyn.dmft import orbspin
from w2dyn.dmft import worm
from w2dyn.dmft import mixing

def git_revision():
    try:
        return os.environ["DMFT_GIT_REVISION"]
    except KeyError:
        pass
    try:
        return Popen(["git", "rev-parse", "HEAD"], stdout=PIPE,
                     stderr=open(os.path.devnull, 'w'),
                     cwd=os.path.dirname(os.path.realpath(__file__))
                     ).communicate()[0]
    except Exception:
        return None

# MPI initialisation
use_mpi = True
if use_mpi:
    from w2dyn.dmft import mpi
    mpi_comm = mpi.MPI_COMM_WORLD
    mpi_rank = mpi_comm.Get_rank()
    mpi_size = mpi_comm.Get_size()

    def mpi_abort(type, value, traceback):
        sys.__excepthook__(type, value, traceback)
        sys.stderr.write("Error: Exception at top-level on rank %s "
                         "(see previous output for error message)\n" %
                         mpi_rank)
        sys.stderr.flush()
        mpi_comm.Abort(1)

    sys.excepthook = mpi_abort

    def mpi_on_root(func):
        """Executes a function only on the root node"""
        if mpi_rank == 0: result = func()
        else: result = None
        return mpi_comm.bcast(result)
else:
    mpi_comm = None
    mpi_rank = 0
    mpi_size = 1
    mpi_on_root = lambda func: func()


mpi_iamroot = mpi_rank == 0

def my_show_warning(message, category, filename, lineno, file=None, line=None):
    if not mpi_iamroot: return
    message = str(message).replace("\n", "\n\t")
    sys.stderr.write("\nWARNING: %s\n\t%s triggered at %s:%s\n\n" %
                     (message, category.__name__, filename, lineno))

warnings.showwarning = my_show_warning

if mpi_iamroot:
    log = lambda s, *a: sys.stderr.write(str(s) % a + "\n")
    rerr = sys.stderr
else:
    log = lambda s, *a: None
    rerr = open(os.devnull, "w")

# Print banners and general information
log(aux.BANNER, aux.CODE_VERSION_STRING, aux.CODE_DATE)
log("Running on %d core%s", mpi_size, " s"[mpi_size > 1])
log("Calculation started %s", time.strftime("%c"))

# Parse positional arguments
key_value_args, argv = config.parse_pairs(sys.argv[1:])
parser = optparse.OptionParser(usage="%prog [key=[value] ...] [FILE ...]",
                               description=__doc__,
                               version="%prog " + aux.CODE_VERSION_STRING)
prog_options, argv = parser.parse_args(argv)
if len(argv) == 0:
    log("No config file name given, using `Parameters.in' ...")
    cfg_file_name = 'Parameters.in'
elif len(argv) == 1:
    cfg_file_name = argv[0]
else:
    parser.error("Expecting exactly one filename")

cfg =  mpi_on_root(lambda: config.get_cfg(cfg_file_name, key_value_args,
                                          err=rerr))

del cfg_file_name, key_value_args, argv, parser
# Finished argument parsing, cfg now contains the configuration for the run

total_time1=time.time()

# create output file and write metadata
log("Writing basic data for the run ...")
output = hdfout.HdfOutput(cfg, git_revision(), mpi_comm=mpi_comm)

# construct lattice
log("Constructing lattice ...")
mylattice = mpi_on_root(lambda: config.lattice_from_cfg(cfg))

# compute non-interacting chemical potential
log("Computing non-interacting density of states ...")
natoms = cfg["General"]["NAt"]
epsn = cfg["General"]["EPSN"]
mu = cfg["General"]["mu"]
niw = 2*cfg["QMC"]["Niw"]
fix_mu = epsn > 0  # 'fix' like 'repair', not like 'frozen'
mu_method = cfg["General"]["mu_search"]
totdens = cfg["General"]["totdens"] * natoms
if fix_mu and mu_method == "kappa":
    last_mudiff = cfg["General"]["initial_last_mudiff"]
    kappa_forced_sign = cfg["General"]["kappa_force_sign"]
    kmuiter = selfcons.KappaMuExtrapolator(totdens,
                                           (last_mudiff
                                            if last_mudiff != 0.0
                                            else None),
                                           kappa_forced_sign)
mylattice.compute_dos(fix_mu, mu, totdens)
log("LDA chemical potential mu = %.6g, total electrons N = %.6g",
    mylattice.mu, mylattice.densities.sum().real)

# If we want to fix the chemical potential, then start from the LDA one
if fix_mu:
    mu = mylattice.mu

# generate atoms
log("Generating list of atoms ...")
latt_type = cfg["General"]["DOS"]
norbitals = mylattice.norbitals
nspins = mylattice.nspins
atom_list = config.atomlist_from_cfg(cfg, norbitals)

epseq = cfg["General"]["EPSEQ"]
equiv = atoms.check_equivalence(None, atom_list,
        lambda at1, at2: (at1.nd == at2.nd and at1.typ == at2.typ and
                          at1.se_shift == at2.se_shift and
                          interaction.similar(at1.dd_int, at2.dd_int, epseq))
        )
log("Equivalence before G0 check: %s",  equiv)

# Compute G0
beta = cfg["General"]["beta"]
iwf = tf.matfreq(beta, 'fermi', niw)
paramag = cfg["General"]["magnetism"] == "para"

eqmaxfreq = cfg["General"]["eqmaxfreq"]
eqslice = slice(iwf.searchsorted(0), iwf.searchsorted(eqmaxfreq)+1)
siw_zero = np.zeros((iwf.size, mylattice.norbitals, mylattice.nspins,
                     mylattice.norbitals, mylattice.nspins), complex)

log("Checking equivalence of d-d blocks (checking %d frequencies) ...",
    eqslice.stop - eqslice.start)
chk_lattice = mylattice
if use_mpi:
    chk_lattice = mpi.MPILattice(chk_lattice)
glociw_eq = chk_lattice.gloc(iwf[eqslice], mylattice.mu, siw_zero[eqslice])
equiv = atoms.check_equivalence(equiv, atom_list,
        lambda at1, at2: np.allclose(at1.d_downfold(glociw_eq),
                                     at2.d_downfold(glociw_eq), atol=epseq)
        )
log("Equivalence after G0 check: %s", equiv)

if cfg["General"]["equiv"] is not None:
    equiv_new = np.asarray(cfg["General"]["equiv"],dtype=int)

    ### some checks of sanity of equiv_new
    if equiv_new.size != len(atom_list):
        parser.error("length of equiv array must match number of atoms!")
    #if (equiv_new >= len(atom_list)).any():
        #parser.error("invalid atom in array equiv!")
    if (equiv_new < 0).any():
        parser.error("invalid atom in array equiv!")
    if (equiv_new > np.arange(len(atom_list))).any():
        parser.error("wrong order of array equiv; they have to be in ascending order!")

    log("Forcing equivalence to: %s", equiv_new)
    equiv = equiv_new[:]


log("Checking if d-d blocks are diagonal in spin and orbital ...")
epsoffdiag = cfg["General"]["EPSOFFDIAG"]
for iatom, atom in enumerate(atom_list):
    g0iw_dd = atom.d_downfold(glociw_eq).copy()
    orbspin.warn_offdiagonal(g0iw_dd, tolerance=epsoffdiag)

log("Generating inequivalent atoms ...")
ineq_indices = np.unique(equiv)
ineq_list = []
for ineq_index in ineq_indices:
    clones = [at for iat, at in enumerate(atom_list)
              if equiv[iat] == ineq_index]
    curr_ineq = atoms.InequivalentAtom(clones)
    ineq_list.append(curr_ineq)
    log("  inequivalent atom %d: %d clones, starting: %s",
        ineq_index, len(clones)-1, ", ".join(str(c.start) for c in clones))

output.ineq_list = ineq_list
output.write_quantity("lda-mu", mylattice.mu)
output.write_quantity("lda-dens", mylattice.densities)
output.write_quantity("lda-dos", mylattice.dos.transpose(1,2,0))
output.write_quantity("h-mean", mylattice.hloc)
output.write_quantity("h-mom2", mylattice.hmom2)
output.write_quantity("hk-mean", mylattice.hloc)
# FIXME: hack
if mpi_iamroot:
    output.file[".axes"].create_dataset("w-dos", data=mylattice.w)
if isinstance(mylattice, lattice.NanoLattice):
    output.write_quantity("leadsiw-full", mylattice.leadsiw)

# generate inter-atom U
log("Constructing inter-atom U matrix ...")
udd_full, udp_full, upp_full = atoms.construct_ufull(atom_list,
                                                     cfg["General"]["Uw"],
                                                     cfg["General"]["Uw_Mat"])

log("Constructing double-counting ...")
dc = config.doublecounting_from_cfg(cfg, ineq_list, mylattice, atom_list,
                                    udd_full + udp_full + upp_full)

# Initialise solver
log("Initialising solver ...")
if cfg["General"]["DMFTsteps"]==0:
    sr = random.SystemRandom()
    cfg["QMC"]["NSeed"] = mpi_on_root(lambda: sr.randint(0, int(2e9)))
    log("Statistics gathering: setting seed to %d", cfg["QMC"]["NSeed"])

Nseed = cfg["QMC"]["NSeed"] + mpi_rank

#FIXME: hack
cfg["QMC"]["FTType"]=cfg["General"]["FTType"]

Uw = cfg["General"]["Uw"]
Uw_Mat = cfg["General"]["Uw_Mat"]

solver = impurity.CtHybSolver(cfg, Nseed, Uw, Uw_Mat, epsn, not use_mpi, mpi_comm)
#if use_mpi:
#    log("Using MPI-enabled solver")
#    solver = mpi.MPIStatisticalSolver(solver)

log("Initialising DMFT loop ...")
nftau = cfg["QMC"]["Nftau"]
GW = cfg["General"]["GW"]
GW_KAverage = cfg["General"]["GW_KAverage"]
dc_dp = cfg["General"]["dc_dp"]
dc_dp_orbitals = cfg["General"]["dc_dp_orbitals"]

restarted_run = cfg["General"]["readold"]

siw_mixer = mixing.FlatMixingDecorator(
    mixing.LinearMixer(cfg["General"]["mixing"]))
mu_mixer = mixing.LinearMixer(cfg["General"]["mu_mixing"])
if "diis" == cfg["General"]["mixing_strategy"]:
    siw_mixer = mixing.FlatMixingDecorator(mixing.RealMixingDecorator(
        mixing.DiisMixer(cfg["General"]["mixing"],
                         cfg["General"]["mixing_diis_history"],
                         cfg["General"]["mixing_diis_period"])))
if not restarted_run:
    siw_mixer = mixing.InitialMixingDecorator(1, mixing.NoMixingMixer(), siw_mixer)
dmft_step = selfcons.DMFTStep(
      beta, mylattice, ineq_list, niw, nftau, dc_dp, dc_dp_orbitals,
      GW, GW_KAverage, natoms, dc, udd_full, udp_full, upp_full,
      paramag, siw_mixer, mu_mixer, mpi_comm
      )

if restarted_run:
    fileold = cfg["General"]["fileold"]
    try:
        if not os.path.isfile(fileold):
            fileold = next(filter(lambda x: x != output.filename,
                                  sorted(glob.iglob(fileold), reverse=True)))
        iterold = restarted_run
        log("Reading old data from file %s, iteration %d ...", fileold, iterold)
        old_mu, siw_dd, smom_dd, dc_value, old_beta = output.load_old(fileold, iterold)
    except (StopIteration, OSError, TypeError):
        raise ValueError("When General.readold is set to non-zero to continue "
                         "from a previous calculation, General.fileold must "
                         "contain a file path or glob pattern referring to a "
                         "valid w2dynamics output file.")

    # ensure usable smom_dd for potential use in preprocessing
    if smom_dd is None:
        warnings.warn("Extracting old moments from old Sigma",
                      UserWarning, 2)
        smom_dd = [siw_block[:2].real for siw_block in siw_dd]

    # self-energy preprocessing (cropping/extending if necessary and
    # interpolating if desired)

    # step 1: if asked for in config file, continue to matsubara
    # frequencies not contained in the old grid by linear
    # interpolation and, toward omega = 0, extrapolation for positive
    # and negative half separately for at most as many frequencies as
    # needed
    if cfg["General"]["readold_mapsiw"] == 'iomega' and beta != old_beta:
        from scipy.interpolate import interp1d
        for isiw in range(len(siw_dd)):
            siw = siw_dd[isiw]

            niw_old = siw.shape[0]
            assert(niw_old % 2 == 0)
            if niw_old < 4:
                warnings.warn("No self-energy interpolation attempted for less"
                              " than two positive old Matsubara frequencies")
                break
            iwf_old = tf.matfreq(old_beta, 'fermi', niw_old)

            iwf_new = tf.matfreq(beta, 'fermi', niw)
            lastinc = np.searchsorted(iwf_new, iwf_old[0], side='left')
            iwf_new = iwf_new[lastinc:niw - lastinc]
            niw_new = len(iwf_new)

            siw_new = np.zeros([niw_new] + list(siw.shape[1:]), dtype=siw.dtype)
            for i in range(siw.shape[1]):
                for j in range(siw.shape[2]):
                    for k in range(siw.shape[3]):
                        for l in range(siw.shape[4]):
                            # Positive half (do halves separately to
                            # not force intrapolation Sigma(0) = 0)
                            siw_ipl = interp1d(iwf_old[niw_old//2:],
                                               siw[niw_old//2:, i, j, k, l],
                                               bounds_error=False,
                                               fill_value='extrapolate',
                                               assume_sorted=True)

                            siw_new[niw_new//2:,
                                    i, j, k, l] = siw_ipl(iwf_new[niw_new//2:])

                            # Negative half
                            siw_ipl = interp1d(iwf_old[:niw_old//2],
                                               siw[:niw_old//2, i, j, k, l],
                                               bounds_error=False,
                                               fill_value='extrapolate',
                                               assume_sorted=True)

                            siw_new[:niw_new//2,
                                    i, j, k, l] = siw_ipl(iwf_new[:niw_new//2])
            siw_dd[isiw] = siw_new

    # step 2: extension or cropping of old data to new frequency grid
    for isiw in range(len(siw_dd)):
        siw = siw_dd[isiw]
        niw_old = siw.shape[0]
        if niw_old != niw:
            assert(niw_old % 2 == 0)
            siw_new = np.empty([niw] + list(siw.shape[1:]), dtype=siw.dtype)
            siw_new[:, :, :, :, :] = smom_dd[isiw][np.newaxis, 0, :, :, :, :]
            if niw_old > niw:
                siw_new[:, :, :, :, :] = siw[niw_old//2 - niw//2:
                                             niw_old//2 + niw//2,
                                             :, :, :, :]
            else:
                siw_new[niw//2 - niw_old//2:niw//2 + niw_old//2,
                        :, :, :, :] = siw[:, :, :, :, :]
            siw_dd[isiw] = siw_new

    if fix_mu:
        if mu_method == "kappa":
            last_mudiff = cfg["General"]["initial_last_mudiff"]
            oldkappa_mu, oldkappa_n, oldkappa_kappa, oldkappa_mudiff = (
                output.load_old_kappa(fileold,
                                      iterold,
                                      kappa_sign=kappa_forced_sign)
            )
            if last_mudiff is not None:
                oldkappa_mudiff = (last_mudiff
                                   if last_mudiff != 0.0
                                   else None)
            kmuiter.initialize(oldkappa_mu, oldkappa_n,
                               oldkappa_kappa, oldkappa_mudiff)
            if kmuiter.has_mu():
                mu = kmuiter.next_mu()
            else:
                mu = old_mu
        else:
            mu = old_mu
    elif cfg["General"]["readold_mu"] == "always":
        mu = old_mu

    # Tell the user where the mu used for (at least) the first
    # iteration came from
    if fix_mu:
        if mu_method == "kappa" and kmuiter.has_mu():
            mu_source_string = "extrapolated"
        else:
            mu_source_string = "old"
    else:
        if cfg["General"]["readold_mu"] == "always":
            mu_source_string = "old"
        else:
            mu_source_string = "parameter-specified"

    log("Continuing old run at old self-energy for %s mu = %.6g ...",
        mu_source_string, mu)
    dmft_step.set_siws(siw_dd, smom_dd, init=True)
    dmft_step.set_mu(mu)
else:
    log("Starting new run without self-energy for %s mu = %.6g ...",
        ("parameter-specified", "LDA")[fix_mu], mu)
    dmft_step.set_siws(None, init=True, hartree_start=True)
    dmft_step.set_mu(mu)

dmft_iterations = cfg["General"]["DMFTsteps"]
stat_iterations = cfg["General"]["StatisticSteps"]
worm_iterations = cfg["General"]["WormSteps"]
total_iterations = dmft_iterations + stat_iterations + worm_iterations
#remember FourPnt for statisitc iters, disable for dmft iters
compute_fourpnt = cfg["QMC"]["FourPnt"]
siw_method = cfg["General"]["SelfEnergy"]
smom_method = cfg["General"]["siw_moments"]

if cfg["QMC"]["ReuseMCConfig"] != 0:
    mccfgs = []

# DMFT loop
for iter_no in range(total_iterations + 1):
    # figure out type of iteration
    iter_time1=time.time()
    if solver.abort:
        iter_type = "aborted"
        iter_no = None
        log("Cleaning up aborted run ...")
    elif iter_no < dmft_iterations:
        iter_type = "dmft"
        cfg["QMC"]["FourPnt"] = 0
        log("Starting DMFT iteration no. %d ...", iter_no+1)
    elif iter_no < dmft_iterations + stat_iterations:
        iter_type = "stat"
        iter_no -= dmft_iterations
        cfg["QMC"]["FourPnt"] = compute_fourpnt
        log("Starting statistics iteration no. %d ...", iter_no+1)
    elif iter_no < total_iterations:
        iter_type = "worm"
        iter_no -= dmft_iterations + stat_iterations
        cfg["QMC"]["FourPnt"] = compute_fourpnt
        log("Starting worm iteration no. %d ...", iter_no+1)    
    else:
        iter_type = "finish"
        iter_no = None
        log("Completing last self-consisting cycle ...")

    iter_start = time.time()

    if mpi_iamroot:
        output.open_iteration(iter_type, iter_no)

    if iter_type == "dmft" or ((iter_type == "stat" or iter_type == "worm") and iter_no == 0):
        if fix_mu and (mu_method == 'nloc'
                       or (mu_method == 'kappa' and not kmuiter.has_mu())):
            log("Computing lattice problem for old mu = %g ...", dmft_step.mu)
            dmft_step.siw2gloc()
            log("Total density = %g", dmft_step.densities.sum())

            log("Writing lattice quantities for old mu to output file ...")
            dmft_step.write_before_mu_search(output)

            log("Updating chemical potential ...")
            dmft_step.update_mu(totdens, epsn, 100.)
        elif fix_mu and mu_method == 'kappa':
            dmft_step.set_mu(kmuiter.next_mu())

    log("Computing lattice problem for mu = %g ...", dmft_step.mu)
    dmft_step.siw2gloc()
    log("Total density = %g", dmft_step.densities.sum())

    log("Writing lattice quantities to output file ...")
    dmft_step.write_lattice_problem(output)

    # The finish iteration is just there to get the final chemical potential
    # and lattice quantities. No more QMC run is performed, so we stop here.
    if iter_type == "finish" or iter_type == "aborted":
        break

    log("Generating impurity problems ...")
    dmft_step.gloc2fiw()
    dmft_step.write_imp_problems(output)

    giws = []
    siws = []
    smoms = []
    occs = []
    imp_electrons = 0.0
    imp_electrons_err = 0.0


    if iter_type == "dmft" or iter_type == "stat":
        if cfg["QMC"]["ReuseMCConfig"] != 0:
            if iter_no > 0 and cfg["QMC"]["Nwarmups2Plus"] >= 0:
                cfg["QMC"]["Nwarmups"] = cfg["QMC"]["Nwarmups2Plus"]

        for iimp, imp_problem in enumerate(dmft_step.imp_problems):
            log("Solving impurity problem no. %d ...", iimp+1)
            solver.set_problem(imp_problem, cfg["QMC"]["FourPnt"])
            if cfg["QMC"]["ReuseMCConfig"] != 0:
                if iter_no > 0:
                    mccfgcontainer = [mccfgs.pop(0)]
                else:
                    mccfgcontainer = []
                result = solver.solve(iter_no, mccfgcontainer)
                result_worm = None
                mccfgs.append(mccfgcontainer[0])
            elif cfg['QMC']['WormMeasGiw'] != 0 or cfg['QMC']['WormMeasGSigmaiw'] != 0 or cfg['QMC']['WormMeasQQ'] != 0:
                result, result_worm = solver.solve_worm(iter_no, log_function=log)
            else:
                result = solver.solve(iter_no)
                result_worm = None
            result.postprocessing(siw_method, smom_method)
            giws.append(result.giw.mean())
            siws.append(result.siw.mean())
            smoms.append(result.smom.mean())
            occs.append(result.other["occ"].mean())

            imp_electrons += (orbspin.trace(result.other["occ"].mean())
                              * len(ineq_list[iimp]))
            imp_electrons_err += (orbspin.trace(result.other["occ"].var())
                                  * len(ineq_list[iimp]))

            if cfg["QMC"]["WriteCovMeanGiw"]:
                output.write_impurity_result(iimp,
                                             {'giw-cov':
                                              diagonal_covariance(result.giw)})
            if cfg["QMC"]["WriteCovMeanSigmaiw"]:
                output.write_impurity_result(iimp,
                                             {'siw-cov':
                                              diagonal_covariance(result.siw)})

            output.write_impurity_result(iimp, result.other)
            # write quantities not contained in result.other after
            # pulling band and spin dimensions to the front
            output.write_impurity_result(
                iimp,
                {'siw-full': result.siw
                 .apply(lambda x: np.transpose(x, (1, 2, 3, 4, 0)), inplace=True),
                 'giw-full': result.giw
                 .apply(lambda x: np.transpose(x, (1, 2, 3, 4, 0)), inplace=True),
                 'smom-full': result.smom
                 .apply(lambda x: np.transpose(x, (1, 2, 3, 4, 0)), inplace=True),
                 })
            if result_worm is not None:
                output.write_impurity_result(iimp, result_worm.other)

        output.write_quantity("giw", giws)
        output.write_quantity("siw", siws)
        output.write_quantity("smom", smoms)
        output.write_quantity("siw-trial", dmft_step.siw_dd)

        if iter_type == "dmft": 
            log("Feeding back self-energies ...")
            dmft_step.set_siws(siws, smoms, giws=giws, occs=occs)

        if fix_mu and mu_method == "kappa":
            imp_electrons_err = np.sqrt(imp_electrons_err)
            kmuiter.step(dmft_step.mu, imp_electrons, imp_electrons_err)

    elif iter_type == "worm":
        for iimp, imp_problem in enumerate(dmft_step.imp_problems):
            log("Solving impurity problem no. %d ...", iimp+1)

            #empty mc configuration for first run
            mc_cfg_container = []
            worm_sector = worm.get_sector_index(cfg['QMC'])
            log("Sampling components of worm sector %d ...", worm_sector)

            #if WormComponents not specified -> sample all
            if not cfg["QMC"]["WormComponents"]:
                log("Sampling all components")
                if worm_sector in [2, 3, 10]:
                   component_list = range(1, imp_problem.nflavours**2 + 1)
                else:
                   component_list = range(1, imp_problem.nflavours**4 + 1)
            else:
                log("Sampling components from configuration file")
                #only integer components 1,... are considered
                component_list = [int(s) for s in cfg["QMC"]["WormComponents"] if int(s)>0]

            for icomponent in component_list:
                log("Sampling component %d", icomponent)
                solver.set_problem(imp_problem, cfg["QMC"]["FourPnt"])
                result_gen, result_comp = solver.solve_comp_stats(iter_no, worm_sector, icomponent, mc_cfg_container)

                # only write result if component-list is user-specified
                # or when eta>0, i.e. the component exists
                worm_eta = result_comp.other['worm-eta/{:05}'.format(icomponent)].mean()
                if (cfg["QMC"]["WormComponents"] or np.amax(worm_eta) > 0.):
                    output.write_impurity_result(iimp, result_comp.other)
                    output.write_impurity_result(iimp, result_gen.other)

                log("Done with component {}".format(icomponent))


    if mpi_iamroot:
        output.close_iteration(iter_type, iter_no)

    iter_time2=time.time()
    log("Time of iteration no. %d: %g sec",iter_no+1,iter_time2-iter_time1)

total_time2=time.time()
log("Total time of calculation: %g sec",total_time2-total_time1)
log("Finished calculation.")
