"""Module for generation of stuff from configurations"""

from warnings import warn
import os.path
import sys
import configobj
import validate
import numpy as np
import h5py as hdf5

import input as _input

from ..dmft import orbspin as orbspin
from ..dmft import lattice as lattice
from ..dmft import atoms as atoms
from ..dmft import interaction as interaction
from ..dmft import doublecounting as doublecounting
from ..dmft import dynamicalU as dynamicalU

def get_cfg(cfg_file_name="Parameters.in", kvargs={}, err=sys.stderr):
    """ Parse a config """
    def validate_convert(value, func, name=None):
        """Helper function for validate parameter conversion"""
        if value is None:
            return None
        try:
            return func(value)
        except (TypeError, ValueError):
            if name is None:
                raise validate.VdtValueError(value)
            else:
                raise validate.VdtParamError(name, value)

    def sci_int(val):
        """Converts type to int, honouring scientific notation if possible"""
        try:
            i = int(val)
        except ValueError:
            val = float(val)
            i = int(val)
        if type(val)(i) != val:
            raise ValueError("Not exactly representable as integer")
        return i

    def is_sci_integer(value, min=None, max=None):
        """Function for the validate module to support scientific notation"""
        value = validate_convert(value, sci_int)
        minval = validate_convert(min, sci_int, "min")
        maxval = validate_convert(max, sci_int, "max")
        if min is not None and value < minval:
            raise validate.VdtValueTooSmallError(value)
        if max is not None and value > maxval:
            raise validate.VdtValueTooBigError(value)
        return value

    def is_option_or_float_list(value, *options):
        try:
            return validate.is_option(value, *options)
        except validate.VdtTypeError:
            return validate.is_float_list(value)

    # config 
    configspec = file(os.path.dirname(__file__) + '/configspec','r')
    cfg = configobj.ConfigObj(infile=file(cfg_file_name, 'r'),
                              configspec=configspec.readlines(),
                              indent_type="\t")

    # update command line parameters
    for key, value in kvargs.iteritems():
        groups = key.split(".")
        parent = cfg
        for group in groups[:-1]:
            parent = parent[group]
        parent[groups[-1]] = value

    validator = validate.Validator()
    validator.functions["integer"] = is_sci_integer
    validator.functions["option_or_list"] = is_option_or_float_list
    valid = cfg.validate(validator, copy=True)

    try:
        pairs = configobj.get_extra_values(cfg)
    except AttributeError:
        print >> err, "WARNING: cannot check unknown entries in config"
        pairs = False

    if pairs:
        print >> err, "error: unknown entries in config: %s" % cfg_file_name
        print >> err, ">>>", ", ".join(".".join(e[0] + e[1:]) for e in pairs)
        raise CfgException()

    if valid is not True:
        print >> err, "error: invalid entries in config: %s" % cfg_file_name
        print >> err, ">>>", ", ".join(str(entry) for entry, ok 
                                       in flat_items(valid) if not ok)
        raise CfgException()

    return cfg

def lattice_from_cfg(cfg):
    latt_type = cfg["General"]["DOS"]
    beta = cfg["General"]["beta"]
    if latt_type == 'ReadIn' or latt_type == 'ReadInSO':
        has_spin_orbit = latt_type == 'ReadInSO'
        hkfile = file(cfg["General"]["HkFile"], "r")
        hk, _ = _input.read_hamiltonian(hkfile, has_spin_orbit)
        mylattice = lattice.KspaceHamiltonian(beta, hk)
    elif latt_type == 'Bethe' or latt_type == 'semicirc':
        # TODO: crystal field
        half_bw = np.array(cfg["General"]["half-bandwidth"], np.float)
        half_bw = half_bw[:, np.newaxis].repeat(2, 1)
        mylattice = lattice.Bethe(beta, half_bw)
    elif latt_type == 'nano':
    #GS:
        niw = 2*cfg["QMC"]["Niw"]
        readleads = cfg["General"]["readleads"]
        leadsfile = file(cfg["General"]["leadsfile"], "r")
        dos_deltino = cfg["General"]["dos_deltino"]
        beta = cfg["General"]["beta"]    # beta is needed for the Matsubara frequencies, in case leadsiw needs to be constructed
        hkfile = file(cfg["General"]["HkFile"], "r")
        has_spin_orbit = False
        hk, _ = _input.read_hamiltonian(hkfile, has_spin_orbit)
        norbitals = hk.shape[1]
        if readleads:
#AV:
#            leadsw, w_hyb, nleads = _input.read_ImHyb(leadsfile, norbitals,
#                                                      has_spin_orbit)
            leadsw, w_hyb, nleads = _input.read_Delta(leadsfile, norbitals,
                                                      has_spin_orbit)
        else:
            leadsw, w_hyb, nleads = None, None, 0
        mylattice = lattice.NanoLattice(hk, beta, niw, leadsw, w_hyb, nleads,
                                        deltino=dos_deltino)

    elif latt_type == 'EDcheck' or latt_type == 'readDelta':
        raise ValueError("Single-shot: use with `cthyb' instead")
    else:
        raise NotImplementedError("unknown lattice")
    return mylattice

def symmoves_from_cfg(atom_cfg):
    nspins = 2
    symmetries = []
    nds = atom_cfg["Nd"]
    for isymm in range(atom_cfg["NSymMove"]):
        perm = atom_cfg["SymMove%02d" % (isymm+1)]
        perm = np.array(perm.split(), int)
        perm -= 1
        if perm.size != nds * nspins:
            raise ValueError("invalid symmetry move: %s" % (perm+1))
        if (np.sort(perm) != np.arange(nds * nspins)).any():
            raise ValueError("invalid permutation: %s" % (perm+1))
        # remodel perm to be easier accessible
        symmetries.append(perm)
    return tuple(symmetries)

def interaction_from_cfg(atom_cfg, cfg):
    int_type = atom_cfg["Hamiltonian"]
    norbitals = int(atom_cfg["Nd"])
    zero_umatrix = np.zeros((norbitals*2,norbitals*2,norbitals*2,norbitals*2))

    Uw = cfg["General"]["Uw"]
    Uw_Mat = cfg["General"]["Uw_Mat"]

    if int_type == "Density":

        if Uw == 1:
          _Retarded2Shifts = dynamicalU.Replace_Retarded_Interaction_by_a_Shift_of_instantaneous_Potentials(0, zero_umatrix, Uw_Mat, 1)
          if _Retarded2Shifts.Uw == False:
            _Retarded2Shifts.U_shift = [0,0]
          if _Retarded2Shifts.Jw == False:
            _Retarded2Shifts.J_shift = [0,0]
          if _Retarded2Shifts.Vw == False:
            _Retarded2Shifts.V_shift = [0,0]

          result = interaction.Density(norbitals, float(atom_cfg["Udd"] + _Retarded2Shifts.U_shift[0]),
                        float(atom_cfg["Vdd"] + _Retarded2Shifts.V_shift[0]), float(atom_cfg["Jdd"] + _Retarded2Shifts.J_shift[0])
                        )
          print ' ==> The shifted interaction potentials in the U(w) implementation are:'
          print '       U=', float(atom_cfg["Udd"]), ' > U=' , float(atom_cfg["Udd"] + _Retarded2Shifts.U_shift[0])
          print '       J=', float(atom_cfg["Jdd"]), ' > J=' , float(atom_cfg["Jdd"] + _Retarded2Shifts.J_shift[0])
          print '       V=', float(atom_cfg["Vdd"]), ' > V=' , float(atom_cfg["Vdd"] + _Retarded2Shifts.V_shift[0])
          print "         ****************************"
          print "         **** U(w) Config Message ***"
          print "         ****************************"
          print " "

        else:
          result = interaction.Density( norbitals, float(atom_cfg["Udd"]), float(atom_cfg["Vdd"]), float(atom_cfg["Jdd"]) )

    elif int_type == "Kanamori":
        if Uw == 1:
          print 'Config.py: Invalid U(w) Configuration: Kanamori. Please use Density for now.'
          exit()
        result = interaction.Kanamori(norbitals, float(atom_cfg["Udd"]),
                        float(atom_cfg["Vdd"]), float(atom_cfg["Jdd"])
                        )

    elif int_type == "Coulomb":
        if Uw == 1:
          print 'Config.py: Invalid U(w) Configuration: Coulomb. Please use Density for now.'
          exit()
        if norbitals == 5:
           # this is a hack for the case, that F6 is not set in parameters file.
           # If it was set, then it is overwritten here. Nevermind for Nd=5
           # below the same is done for s and p shells
           atom_cfg["F6"] = 0.0
        elif norbitals == 3:
            atom_cfg["F4"] = 0.0
            atom_cfg["F6"] = 0.0
        elif norbitals == 1:
            atom_cfg["F2"] = 0.0
            atom_cfg["F4"] = 0.0
            atom_cfg["F6"] = 0.0

        result = interaction.Coulomb(norbitals, float(atom_cfg["F0"]),
                        float(atom_cfg["F2"]), float(atom_cfg["F4"]),
                        float(atom_cfg["F6"])
                        )

    elif int_type == "ReadUmatrix":
        if Uw == 1:
          print 'Config.py: Invalid U(w) configuration: ReadUmatrix. Please use Density for now.'
          exit()
        u_matrix = _input.read_u_matrix(atom_cfg["umatrix"], True)
        result = interaction.CustomFull(u_matrix)

    elif int_type == "ReadNormalUmatrix":
        orb_u_matrix = _input.read_u_matrix(atom_cfg["umatrix"], False)
        if Uw == 1:
          _Retarded2Shifts = dynamicalU.Replace_Retarded_Interaction_by_a_Shift_of_instantaneous_Potentials(0, zero_umatrix, Uw_Mat,1)
          orb_u_matrix = _Retarded2Shifts.shift_instantaneous_densitydensity_potentials_from_ReadNormalUmatrix_config(orb_u_matrix)
          print ' ==> The orbital interaction potentials (umatrix entries) were shifted'
          print "         ****************************"
          print "         **** U(w) Config Message ***"
          print "         ****************************"
          print " "
        result = interaction.CustomSU2Invariant(orb_u_matrix)

    else:
        raise ValueError('Unrecognized interaction: %s' % int_type)

    if atom_cfg["QuantumNumbers"]:
        manual_qn = tuple(atom_cfg["QuantumNumbers"].split())
        if manual_qn != result.quantum_numbers:
            warn("Overriding default quantum numbers %s with %s.\n"
                 "This may very well eat your run." %
                 (result.quantum_numbers, manual_qn))
            result.quantum_numbers = manual_qn

    return result

def atomlist_from_cfg(cfg, norbitals=None):
    start = 0
    atom_list = []
    if norbitals is None:
        norbitals = sum(acfg["Nd"] + acfg["Np"] for acfg in cfg["Atoms"])
    for atom_cfg in cfg["Atoms"].itervalues():
        nds = atom_cfg["Nd"]
        nps = atom_cfg["Np"]
        nlig =  atom_cfg["Nlig"]
        se_shift = atom_cfg["se-shift"]
        occ_dc_number = atom_cfg["occdcnumber"]
        typ = atom_cfg["type"]
        dd_int = interaction_from_cfg(atom_cfg, cfg)
        # read out symmetry moves
        symmetry_moves = symmoves_from_cfg(atom_cfg)

        # construct crystal field
        crystalfield = atom_cfg["crystalfield"]
        if crystalfield is not None:
            if len(crystalfield) == nds:
                raise ValueError("expecting number of d-orbitals")
            crystalfield = np.diag(crystalfield)
        
        # add d atom
        if np == 0:
            atom_list.append(atoms.Atom(norbitals, start, nds, dd_int,
                                        se_shift, symmetry_moves, typ,
                                        crystalfield=crystalfield,
                                        occ_dc_number=occ_dc_number))
        else:
            if nlig == 0: nlig = 1
            dp_int = atoms.DpInteraction(atom_cfg["Udp"], atom_cfg["Jdp"])
            pp_int = atoms.PpIntraLigand(atom_cfg["Upp"], atom_cfg["Vpp"],
                                         atom_cfg["Jpp"], atom_cfg["Uppod"],
                                         atom_cfg["Jppod"])
            atom_list.append(atoms.Atom(norbitals, start, nds, dd_int, se_shift,
                                        symmetry_moves, typ, nps, nlig, pp_int,
                                        dp_int, crystalfield,
                                        occ_dc_number=occ_dc_number))
        start += nds + nps
    if start != norbitals:
        raise RuntimeError("lattice dimensions and config are inconsistent")
    return atom_list

def doublecounting_from_cfg(cfg, ineq_list, mylattice, atom_list, u_full):
    dc_type = cfg["General"]["dc"]

    # No double counting needed (d-only and max. 1 ineq atom and not
    # user-enforced)
    if len(ineq_list) == 1 and not np.any([atom.np for atom in atom_list]) \
       and not isinstance(dc_type, list):
        return doublecounting.Zero(mylattice.norbitals, mylattice.nspins)

    lda_dens = mylattice.densities
    if isinstance(dc_type, list):
        # user-supplied DC
        dc_user_values = np.array(dc_type, float)
        if dc_user_values.size != mylattice.norbitals:
            raise ValueError("Expected %d elements" % mylattice.norbitals)
        dc_user_values = np.repeat(dc_user_values[:,None], mylattice.nspins, 1)
        dc_user_values = orbspin.promote_diagonal(dc_user_values)
        return doublecounting.UserSupplied(dc_user_values)
    elif dc_type == 'fll' or dc_type == 'anisimov':
        return doublecounting.FullyLocalisedLimit(lda_dens, atom_list, u_full)
    elif dc_type == 'amf':
        return doublecounting.AroundMeanField(lda_dens, atom_list, u_full)
    elif dc_type == 'siginfbar':
        return doublecounting.SigInfBar(ineq_list, lda_dens, atom_list, u_full)
    elif dc_type == 'sigzerobar':
        return doublecounting.SigZeroBar(ineq_list, lda_dens, atom_list, u_full)
    elif dc_type == 'trace':
        return doublecounting.Trace(ineq_list, lda_dens, atom_list, u_full)
    else:
        raise NotImplemented("Unknown double counting scheme")

# --------- helper routines and classes -----------

class CfgException(Exception):
    pass

def find_file(path):
    """ Locates the file named path relative to directories in sys.path """ 
    for directory in sys.path:
        candidate = os.path.join(directory, path)
        if os.path.isfile(candidate):
            return candidate
    raise IOError("Not found on any directory on path: " + path)

def flat_items(d, prefix=""):
    items = []
    for k, v in d.iteritems():
        fullk = prefix + k
        if isinstance(v, dict): 
            items.extend(flat_items(v, fullk + "."))
        else:
            items.append((fullk, v))
    return items

def parse_pairs(argv):
    """Extracts all key=value pairs from the argument vector"""
    kv_opts = {}
    clean_args = []
    for iarg, arg in enumerate(argv):
        # rules to be consistent with getopt/optparse's option parsing rules
        if arg[0] in ("-", "+") or arg.find("=") < 1:
            clean_args.append(arg)
        elif arg == "--":
            clean_args.extend(arg[iarg:])
        else: # a key/value-pair
            eqpos = arg.find("=")
            kv_opts[arg[:eqpos].strip()] = arg[eqpos+1:].strip()
    return kv_opts, clean_args

class Hdf5Type(object):         # OOP can be a joy :-)
    """Factory for creating h5py.File object types

    Instances of Hdf5Type are typically passed as type= arguments to
    the ArgumentParser.add_argument method.  Derived from
    argparser.FileType.

    Keyword arguments are passed on to h5py.File().
    """

    def __init__(self, mode=None, driver=None, libver=None,
                 notify=open(os.devnull, 'w'), **kwargs):
        self.notify = notify
        self.mode   = mode
        self.driver = driver
        self.libver = libver
        self.kwargs = kwargs

    def __call__(self, fname):
        if fname == 'latest':
            files = [f for f in os.listdir('.') if f.endswith('.hdf5')]

            if not files:
                raise Exception("no HDF5 file in current directory to use 'latest' with")

            fname = max(files, key=lambda f: os.stat(f).st_mtime)

            print >> self.notify, \
                "Using latest file from current directory: `%s'" % fname

        return hdf5.File(fname, self.mode, self.driver, self.libver, **self.kwargs)

    def __repr__(self):
        args = ('mode='   + repr(self.mode),
                'driver=' + repr(self.driver),
                'libver=' + repr(self.libver),
                'notify=' + repr(self.notify),
                repr(self.kwargs))

        return '%s(%s)' % (type(self).__name__, ', '.join(args))

