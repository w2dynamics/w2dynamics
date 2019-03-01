#!/usr/bin/env python
"""Program for the DMFT self-consistency loop"""
import os.path
import glob
import random
from subprocess import Popen, PIPE
import sys
import time  # necessary for tokyo cluster
import traceback
import optparse
import warnings

import numpy as np

import w2dyn.auxiliaries as aux
from w2dyn.auxiliaries import transform as tf
from w2dyn.auxiliaries import wien2k as wien2k
from w2dyn.auxiliaries import hdfout
from w2dyn.auxiliaries import config

from w2dyn.dmft import impurity
from w2dyn.dmft import lattice
from w2dyn.dmft import atoms
from w2dyn.dmft import interaction
from w2dyn.dmft import selfcons
from w2dyn.dmft import orbspin

def git_revision():
    try:
        return os.environ["DMFT_GIT_REVISION"]
    except KeyError:
        pass
    try:
        return Popen(["git", "rev-parse", "HEAD"], stdout=PIPE,
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
        mpi_comm.Abort()

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
    rerr = file(os.devnull, "w")

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

total_time1=time.time()

# Finished argument parsing, cfg now contains the configuration for the run

# creating output file
Output = hdfout.HdfOutput

# write important stuff to output file
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
fix_mu = epsn > 0
totdens = cfg["General"]["totdens"] * natoms
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
udd_full, udp_full, upp_full = atoms.construct_ufull(atom_list)

log("Constructing double-counting ...")
dc = config.doublecounting_from_cfg(cfg, ineq_list, mylattice, atom_list,
                                    udd_full + udp_full + upp_full)

# Initialise solver
log("Initialising solver ...")
if cfg["General"]["DMFTsteps"]==0:
    sr = random.SystemRandom()
    cfg["QMC"]["NSeed"] = sr.randint(0, int(2e9))
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

dmft_step = selfcons.DMFTStep(
      beta, mylattice, ineq_list, niw, nftau, dc_dp, dc_dp_orbitals, GW, GW_KAverage, natoms, dc, udd_full, udp_full, upp_full,
      paramag, selfcons.LinearMixingStrategy(cfg["General"]["mixing"]),
      selfcons.LinearMixingStrategy(cfg["General"]["mu_mixing"]), mpi_comm
      )

restarted_run = cfg["General"]["readold"]
if restarted_run:
    fileold = cfg["General"]["fileold"]
    if not os.path.isfile(fileold):
        fileold = filter(lambda x: x != output.filename,
                         sorted(glob.iglob(fileold), reverse=True))[0]
    iterold = restarted_run
    log("Reading old data from file %s, iteration %d ...", fileold, iterold)
    old_mu, siw_dd, smom_dd, dc_value = output.load_old(fileold, iterold)
    if fix_mu:
        mu = old_mu

    log("Continuing old run at old self-energy for %s mu = %.6g ...",
        ("parameter-specified", "old")[fix_mu], mu)
    dmft_step.set_siws(siw_dd, smom_dd, dc_value, init=True)
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
        if cfg["QMC"]["PercentageWormInsert"] != 0:
           log("Use WormSteps instead of DMFTsteps for worm sampling.")
           sys.exit()   
    elif iter_no < dmft_iterations + stat_iterations:
        iter_type = "stat"
        iter_no -= dmft_iterations + worm_iterations
        cfg["QMC"]["FourPnt"] = compute_fourpnt
        log("Starting statistics iteration no. %d ...", iter_no+1)
        if cfg["QMC"]["PercentageWormInsert"] != 0:
           log("Use WormSteps instead of StatisticSteps for worm sampling.")
           sys.exit()   
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
        output.next_iteration(iter_type, iter_no)

    if iter_type == "dmft" or ((iter_type == "stat" or iter_type == "worm") and iter_no == 0):
        if epsn > 0:
            log("Computing lattice problem for old mu = %g ...", dmft_step.mu)
            dmft_step.siw2gloc()
            log("Total density = %g", dmft_step.densities.sum())

            log("Writing lattice quantities for old mu to output file ...")
            dmft_step.write_before_mu_search(output)

            log("Updating chemical potential ...")
            dmft_step.update_mu(totdens, epsn, 100.)

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
                mccfgs.append(mccfgcontainer[0])
            else:
                result = solver.solve(iter_no)
            result.postprocessing(siw_method, smom_method)
            giws.append(result.giw)
            siws.append(result.siw)
            smoms.append(result.smom)
            occs.append(result.other["occ"])
            output.write_impurity_result(iimp, result.other)
       
        output.write_quantity("giw", giws)
        output.write_quantity("giw-full", giws)
        output.write_quantity("siw", siws)
        output.write_quantity("smom", smoms)
        output.write_quantity("siw-full", siws)
        output.write_quantity("smom-full", smoms)

        if iter_type == "dmft": 
            log("Feeding back self-energies ...")
            dmft_step.set_siws(siws, smoms, giws=giws, occs=occs)

    elif iter_type == "worm":
      for iimp, imp_problem in enumerate(dmft_step.imp_problems):
         log("Solving impurity problem no. %d ...", iimp+1)

         #empty mc configuration for first run
         mccfgcontainer = []

         #looping over all worm spaces
         maxsector = 9
         for isector in xrange(2, maxsector + 1):
            
            #skipping sectors if measurement is not enabled
            if isector == 2 and not (cfg["QMC"]["WormMeasGiw"] == 1 or cfg["QMC"]["WormMeasGtau"] == 1):
               log("Skipping worm sector %d ...", isector)
               continue
            if isector == 3 and not cfg["QMC"]["WormMeasGSigmaiw"] == 1:
               log("Skipping worm sector %d ...", isector)
               continue
            if isector == 4: 
               if not cfg["QMC"]["WormMeasG4iw"] == 1:
                  log("Skipping worm sector %d ...", isector)
                  continue
               if cfg["QMC"]["FourPnt"] != 8:
                  log("Set FourPnt to '8' to measure worm")
                  sys.exit()
            if isector == 5: 
               if not cfg["QMC"]["WormMeasH4iw"] == 1:
                  log("Skipping worm sector %d ...", isector)
                  continue
               if cfg["QMC"]["FourPnt"] != 8:
                  log("Set FourPnt to '8' to measure worm")
                  sys.exit()
            if isector == 6 and not (cfg["QMC"]["WormMeasP2iwPH"] == 1 or cfg["QMC"]["WormMeasP2tauPH"] == 1):
               log("Skipping worm sector %d ...", isector)
               continue
            if isector == 7 and not (cfg["QMC"]["WormMeasP2iwPP"] == 1 or cfg["QMC"]["WormMeasP2tauPP"] == 1): 
               log("Skipping worm sector %d ...", isector)
               continue
            if isector == 8 and not cfg["QMC"]["WormMeasP3iwPH"]:
               log("Skipping worm sector %d ...", isector)
               continue
            if isector == 9 and not cfg["QMC"]["WormMeasP3iwPP"]:
               log("Skipping worm sector %d ...", isector)
               continue

            log("Sampling components of worm sector %d ...", isector)

            #if WormComponents not specified -> sample all
            if not cfg["QMC"]["WormComponents"]:
               log("Sampling all components")
               if isector < 4:
                  component_list = xrange(1,imp_problem.nflavours**2+1)
               else:
                  component_list = xrange(1,imp_problem.nflavours**4+1)
            else:
               log("Sampling components from configuration file")
               #only integer components 1,... are considered
               component_list = [int(s) for s in cfg["QMC"]["WormComponents"] if int(s)>0]


            for icomponent in component_list:
               log("Sampling component %d", icomponent)
               solver.set_problem(imp_problem, cfg["QMC"]["FourPnt"])
               result, result_aux = solver.solve_component(iter_no,isector,icomponent,mccfgcontainer)

               #only write result if component returns ne 0
               if np.any([np.any(result.other[ky]["value"])  for ky in result.other.keys()]):
                  output.write_impurity_component(iimp, result.other)
                  try:
                     output.write_impurity_component(iimp, result_aux.other)
                  except ValueError:
                     sys.stderr.write(
                      "\nWARNING: Ignoring auxiliary entries for multiple worm estimators.\n\n")

            log("Done with component")

      #after worm-sampling we terminate all loops
      break 
    
    iter_time2=time.time()
    log("Time of iteration no. %d: %g sec",iter_no+1,iter_time2-iter_time1)

total_time2=time.time()
log("Total time of calculation: %g sec",total_time2-total_time1)
log("Finished calculation.")
