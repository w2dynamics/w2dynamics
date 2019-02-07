"""Package abstracting partition of the Hamiltonian into atoms as well as
   up- and downfolding between the impurity and the lattice picture.
"""
import numpy as np
import interaction

class Atom:
    def __init__(self, ntotal_orbitals, start_orbital, nd, dd_interaction,
                 se_shift=0., symmetry_moves=(), typ=None,
                 np=0, nlig=1, pp_interaction=None, dp_interaction=None,
                 crystalfield=None, occ_dc_number=None):
        self.ntotal = ntotal_orbitals
        self.start = start_orbital
        self.nd = nd
        self.np = np
        self.norbitals = nd + np
        self.nlig = nlig
        self.se_shift = se_shift
        self.occ_dc_number = occ_dc_number
        self.dslice = slice(start_orbital, start_orbital+nd)
        self.pslice = slice(start_orbital+nd, start_orbital+nd+np)
        self.dd_int = dd_interaction
        self.dp_int = dp_interaction
        self.pp_int = pp_interaction
        self.symmetry_moves = symmetry_moves
        self.typ = typ
        self.crystalfield = crystalfield

        if np > 0:
            if np % nlig:
                raise ValueError("np must be dividable by nlig")
            start_orbital += nd
            self.ligslices = []
            self.npplig = np // nlig
            for pstart in range(start_orbital, start_orbital+np, self.npplig):
                self.ligslices.append(slice(pstart, pstart + self.npplig))
            self.pp_int = pp_interaction

    def d_setpart(self, x, part):
        x[..., self.dslice, :, self.dslice, :] = part

    def d_downfold(self, x):
        return x[..., self.dslice, :, self.dslice, :]

    def d_downfolded_shape(self, shape):
        df_shape = list(shape)
        df_shape[-4] = self.nd
        df_shape[-2] = self.nd
        return tuple(df_shape)

    def d_upfold(self, x):
        res = np.zeros(x.shape[:-4] + (self.ntotal, x.shape[-3],
                                       self.ntotal, x.shape[-1]), x.dtype)
        res[..., self.dslice, :, self.dslice, :] = x
        return res

def check_equivalence(equiv_prev, atoms, criterion):
    """Check equivalence of atoms using some criterion."""
    # The i-th entry of the equivalence array contains i if the atom is
    # a 'template', and some other value j < i if it equivalent to the
    # template atom j.
    natoms = len(atoms)
    if equiv_prev is None:
        # initialise to all atoms being (potentially) equivalent
        equiv = np.zeros(natoms, np.uint)
    else:
        equiv = np.array(equiv_prev)
        if equiv.shape != (natoms,):
            raise ValueError("illegal equivalence array")

    # Check equivalence by this criterion
    for iat1, atom1 in enumerate(atoms):
        # I have already detected that atom to be equivalent to another
        if equiv[iat1] != iat1: continue

        # We only have to iterate over the "top triangle", since other
        # combinations are treated previously.
        first_inequiv = None
        for idiff, atom2 in enumerate(atoms[iat1+1:]):
            iat2 = idiff + iat1 + 1
            # first check that the two atoms belong to the same equivalence
            # class so that we can verify its equivalence.
            if iat1 != equiv[iat2]:
                continue

            # The first non-equivalent object marks the start of a new
            # equivalence class, while subsequent ones are also assigned to
            # this class for checking later.
            if not criterion(atom1, atom2):
                if first_inequiv is None:
                    first_inequiv = iat2
                equiv[iat2] = first_inequiv

    # logical AND of inequivalence with previous equivalence is taken care
    # by the condition inside the loop, so we can return equiv
    return equiv

def construct_ufull(atom_list):
    ntotal = atom_list[0].ntotal
    udd_full = np.zeros((ntotal, 2, ntotal, 2))
    udp_full = udd_full.copy()
    upp_full = udd_full.copy()

    # set density part of intra-atom U (d-d/p-p)
    for iatom1, atom1 in enumerate(atom_list):
        udd_full[atom1.dslice, :, atom1.dslice, :] = atom1.dd_int.get_udens()
        udp_full[atom1.dslice, :, atom1.pslice, :] = \
                                    atom1.dp_int.get_udens(atom1, atom1)
        upp_full[atom1.pslice, :, atom1.pslice, :] = \
                                    atom1.pp_int.get_udens(atom1, atom1)/2.
        for atom2 in atom_list[iatom1+1:]:
            # set upper triangle of inter-atom U (d-p), but use double the value
            # such that the symmetrisation below yields the right thing
            udp_full[atom1.dslice, :, atom2.pslice, :] = \
                                    atom1.dp_int.get_udens(atom1, atom2)
            upp_full[atom1.pslice, :, atom2.pslice, :] = \
                                    atom1.pp_int.get_udens(atom1, atom2)

    # set lower triangle: in d-p, the diagonal terms are obviously zero, in p-p
    # the diagonal is halved above.
    udp_full = udp_full + udp_full.transpose(2,3,0,1)
    upp_full = upp_full + upp_full.transpose(2,3,0,1)
    return udd_full, udp_full, upp_full

class InequivalentAtom(Atom):
    def _verify_clones(self):
        pass

    def __init__(self, atoms):
        self.template = atoms[0]
        self.clones = atoms[1:]

    def d_setpart(self, x, part):
        self.template.d_setpart(x, part)
        for clone in self.clones:
            clone.d_setpart(x, part)

    def d_downfold(self, x):
        return self.template.d_downfold(x)

    def d_upfold(self, x):
        res = self.template.d_upfold(x)
        for clone in self.clones:
            clone.d_setpart(res, x)
        return res

    def __getattr__(self, attr):
        return getattr(self.template, attr)

class DpInteraction:
    """Class for d-p density-density interaction"""
    def __init__(self, udp=0., jdp=0.):
        self.udp = udp
        self.jdp = jdp

    def get_udens(self, d_atom, p_atom):
        uinter = np.zeros((d_atom.nd, 2, p_atom.np, 2))
        uinter[...] = self.udp
        uinter[:, 0, :, 0] -= self.jdp
        uinter[:, 1, :, 1] -= self.jdp
        return uinter

class PpIntraLigand:
    """Class for p-p density-density interaction"""
    def __init__(self, upp=0., vpp=0., jpp=0., upp_offdiag=0., jpp_offdiag=0.):
        self.upp = upp
        self.vpp = vpp
        self.jpp = jpp
        # extra-ligand off-diagonal terms
        self.upp_offdiag = upp_offdiag
        self.jpp_offdiag = jpp_offdiag

    def get_udens(self, atom1, atom2):
        upart = np.zeros((atom1.np, 2, atom2.np, 2))
        if atom1.start != atom2.start:
            return upart
        if atom1.np == 0:
            return upart
        npplig = atom1.np // atom1.nlig
        for ligstart in range(0, atom1.np, npplig):
            ligslice = slice(ligstart, ligstart + npplig)
            upart[ligslice, :, ligslice, :] = interaction.udensity_values(
                          npplig, self.upp, self.vpp, self.jpp)
        # FIXME ligand off-diagonal stuff
        return upart
