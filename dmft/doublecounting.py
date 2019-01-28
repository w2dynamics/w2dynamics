"""Classes for double-counting correction"""

from __future__ import division
import numpy as np

#import dmft.lattice as latt
import orbspin as orbspin
import atoms as atoms
import lattice as lattice
from ..auxiliaries import transform as tf
from itertools import islice

class DoubleCounting:
    """Base class for Double counting correction"""
    def __init__(self):
        self.dc_value = None
    def get(self, impurity_results=None):
        """Gets the value of the double counting correction.

        This method is called *after* the solution to the impurity problems in
        iteration `iter_no` (zero-based) and is expected to yield a correction
        to the self-energy.
        """
        return self.dc_value
    def set(self,dc_value):
        """Sets the value of the double counting correction internally."""
        self.dc_value = dc_value

class Zero(DoubleCounting):
    """Trivial double counting that always yields zero"""
    def __init__(self, nbands, nspins):
        self.self_cons = False
        self.dc_value = np.zeros((nbands, nspins, nbands, nspins))

class UserSupplied(DoubleCounting):
    """Double counting as supplied by user"""
    def __init__(self, dc_user_supplied):
        DoubleCounting.__init__(self)
        self.self_cons = False
        self.dc_value = dc_user_supplied

class ShellAveraged(DoubleCounting):
    r"""Base class for l-shell-averaged double counting.

    The most common double-counting schemes work with fillings and interactions
    that are averaged over the orbitals of l-subshells of a single atom or
    ligand, e.g. averaged over 5 d-bands or 3 p-bands [1].

    Therefore, this class provides some subshell-averaged quantities:

      - `nshells`:   total number of shells over every atom
      - `slices`:    list of slices for selection of the different shells
      - `dens_sum`:  total filling of a shell in LDA
      - `dens_mean`: average filling of an orbital in a shell in LDA
      - `ubar`:      the mean value of the interaction between shells

    Averaged interaction
    --------------------
    The code computes the mean U and J for every sub-shell. According to Czyzyk
    and Sawatzky [1] and Parragh [2], these are given as:

        Ubar        = 1/N^2    \sum_{mn} U_{mnmn}
        Ubar - Jbar = 1/N(N-1) \sum_{mn} (U_{mnmn} - U_{mnnm})

    (the second line is divided by N(N-1) because the term for m=n is zero by
    construction). Comparing this with the formula for the density-density part
    of the U matrix, which is spin-dependent (s,t):

        U_{ms,nt} = U_{mnmn} - U_{mnnm} \delta_{st}

    we see that for one shell we can write this as (spin-dependent) mean U:

        ubar' = [ Ubar-Jbar Ubar      ] = sum U_{ms,nt} (/) [ N(N-1) N^2    ]
                [ Ubar      Ubar-Jbar ]    mn               [ N^2    N(N-1) ]

    where the bracketed slash `(/)` denotes element-wise division, which is
    omitted for N = 1 in order to avoid 0/0).

    [1]: Czyzyk and Sawatzky, Phys. Rev. B. 49, 14211 (1994)
    [2]: Nicolaus Parragh, Ph.D. thesis
    """
    def __init__(self, densities, atom_list, u_full):
        self.self_cons = False
        if atom_list is None:
            ndorb = u_full.shape[0]
            atom_list = [atoms.Atom(ndorb, 0, ndorb, None)]

        self.atom_list = atom_list

        # split problem into sub-shells.
        self.slices = []
        self.sizes = []
        for atom in self.atom_list:
            self.slices.append(atom.dslice)
            self.sizes.append(atom.nd)
            if atom.np > 0:
                self.slices.extend(atom.ligslices)
                self.sizes.extend([atom.npplig] * atom.nlig)

        self.nshells = len(self.slices)
        self.sizes = np.array(self.sizes)
        self.norbitals, self.nspins = u_full.shape[:2]

        # compute mean density for the different manifolds
        # spin-dependent LDA occupations averaged over the orbitals in the
        # subspace.
        self.dens_mean = np.array([densities[sl].mean() for sl in self.slices])
        self.dens_sum = np.array([densities[sl].sum() for sl in self.slices])
 
        # Compute the mean U and J for every sub-shell.
        ubar = np.zeros((self.nshells, self.nshells))
        jbar = np.zeros((self.nshells, self.nshells))
        for ishell1,shell1 in enumerate(self.slices):
            for ishell2,shell2 in enumerate(self.slices):
                ubar[ishell1, ishell2] = u_full[shell1,0,shell2,1].mean(1).mean(0)
                jbar[ishell1, ishell2] = u_full[shell1,0,shell2,0].sum(1).sum(0)

        # jbar normalization avoiding 0 division.
        for ipos, isize in enumerate(self.sizes):
            for jpos, jsize in enumerate(self.sizes):
                jbar[ipos, jpos] /= isize
                if ipos == jpos:
                    if isize == 1:
                        jbar[ipos, jpos] = 0
                    else:
                        jbar[ipos, jpos] /= isize - 1
                else:
                    jbar[ipos, jpos] /= jsize

        # We actually computed Ubar - Jbar as Jbar
        jbar = ubar - jbar
        self.ubar = ubar
        self.jbar = jbar

        self.dc_shell = None
        self.dc_value = None

    def from_shell(self, dc_shell):
        dc_diag = np.zeros((self.norbitals, self.nspins))
        for isl, sl in enumerate(self.slices):
            dc_diag[sl] = dc_shell[isl, np.newaxis]

        # note the minus sign!
        self.dc_shell = dc_shell
        self.dc_value = -orbspin.promote_diagonal(dc_diag)

class FullyLocalisedLimit(ShellAveraged):
    r"""Double counting in the fully localised limit.

    The fully localised limit double counting (also known as Anisimov's double
    counting) assumes the correlation included in the non-interacting problem
    to represent the atomic limit of the lattice problem:

        DC_d = Ubar_{dd} (N_d - 1/2) + Ubar_{dp} N_p
        DC_p = Ubar_{pp} (N_p - 1/2) + Ubar_{dp} N_d

    When the indices I, J run over "subshells" corresponding either to d-atom
    or a p-atom/ligand, one can use the definition of `ubar` (see documentation
    of `ShellAveraged`) to recast the formula into:

        DC_{Is} = \sum_{Jt} Ubar_{Is,Jt} N_{Jt} - 1/2 \sum_t Ubar_{Is,It}

    Note that the subtraction of 1/2 only occurs for sub-shell local terms.
    """
    def __init__(self, *args):
        ShellAveraged.__init__(self, *args)

        dc_shell = self.ubar.dot(self.dens_sum) - self.jbar.dot(self.dens_sum/2)
        for I in range(self.nshells):
            dc_shell[I] -= self.ubar[I, I] * 0.5
            dc_shell[I] += self.jbar[I, I] * 0.5

        self.from_shell(dc_shell)

class AroundMeanField(ShellAveraged):
    r"""Double counting in the fully localised limit.

    The around mean field double counting can be seen as a correction to the
    fully localised limit for systems away from half filling:

        DC_d = Ubar_{dd} (N_d - n_d) + Ubar_{dp} N_p
        DC_p = Ubar_{pp} (N_p - n_d) + Ubar_{dp} N_d

    When the indices I, J run over "subshells" corresponding either to d-atom
    or a p-atom/ligand, one can use the definition of `ubar` (see documentation
    of `ShellAveraged`) to recast the formula into:

        DC_{Is} = \sum_{Jt} Ubar_{Is,Jt} N_{Jt} - \sum_t Ubar_{Is,It} nbar_{It}

    Note that the subtraction of `nbar` only occurs for sub-shell local terms.
    """
    def __init__(self, *args):
        ShellAveraged.__init__(self, *args)
        
        dc_shell = self.ubar.dot(self.dens_sum) - self.jbar.dot(self.dens_sum/2)
        for I in range(self.nshells):
            dc_shell[I] -= self.ubar[I, I] * self.dens_mean[I]
            dc_shell[I] += self.jbar[I, I] * self.dens_mean[I]

        self.from_shell(dc_shell)

class SelfConsistent(DoubleCounting):
    def __init__(self, ineq_list, densities, atom_list, u_full):
        self.ineq_list = ineq_list
        self.shell_av = {}
        self.shell_av['densities'] = densities
        self.shell_av['atom_list'] = atom_list
        self.shell_av['u_full'] = u_full
        if atom_list is None:
            ndorb = u_full.shape[0]
            atom_list = [atoms.Atom(ndorb, 0, ndorb, None)]

        self.atom_list = atom_list
        # split problem into sub-shells.
        # we only want the slices for the d manifold of each atom
        self.d_slices = []
        for atom in self.atom_list:
            self.d_slices.append(atom.dslice)

        self.norbitals, self.nspins = u_full.shape[:2]

        # compute mean density for the different d manifolds       
        # spin-dependent LDA occupations averaged over the orbitals in the
        # subspace.
        self.d_dens_sum = np.array([densities[sl].sum() for sl in self.d_slices])
        self.dc_value = np.zeros((self.norbitals, self.nspins, self.norbitals, self.nspins))
        self.self_cons = True

class SigZeroBar(SelfConsistent):
    def get(self, siws=None, smoms=None, giws = None, occs = None):
        self.dc_value[...] = 0
        if siws is None:
            return self.dc_value
        for ineq, siw in zip(self.ineq_list, siws):
            niw = siw.shape[0]
            szero = siw[niw//2] # lowest Matsubara point of self-energy, 
                                # since ReS not very w-dependent. might not 
                                # be good enough for high T. Should do 
                                # (simple) extrapolation instead... 

            szerobar = szero.trace(0,1,3).trace().real
            szerobar /= ineq.nd * self.nspins
            szerobar = szerobar * np.eye(ineq.nd * self.nspins).reshape(
                                    ineq.nd, self.nspins, ineq.nd, self.nspins)
            ineq.d_setpart(self.dc_value, -szerobar)
        return self.dc_value

class SigInfBar(SelfConsistent):
    def get(self, siws=None, smoms=None, giws = None, occs = None):
        self.dc_value[...] = 0
        if smoms is None:
            return self.dc_value
            
        for ineq, smom in zip(self.ineq_list, smoms):
            sinf = smom[0]     # 0-th moment
            sinfbar = sinf.trace(0,1,3).trace()
            sinfbar /= ineq.nd * self.nspins
            sinfbar = sinfbar * np.eye(ineq.nd * self.nspins).reshape(
                                    ineq.nd, self.nspins, ineq.nd, self.nspins)
            ineq.d_setpart(self.dc_value, -sinfbar)
        return self.dc_value

class Trace(SelfConsistent):
    def get(self, siws=None, smoms=None, giws = None, occs = None):

        if giws is None:
            # initializing with FLL double counting... better than 0
            dc_fll = FullyLocalisedLimit(self.shell_av['densities'],
                                         self.shell_av['atom_list'],
                                         self.shell_av['u_full'])
            self.dc_value = dc_fll.get()
            return self.dc_value

        for ineq, giw, occ, occ0 in zip(self.ineq_list, giws, occs, self.d_dens_sum):
            # DOS as "approximated" by first Matsubara frequency
            iw0 = giw.shape[0]//2
            dos_fermi_level = -1/np.pi * orbspin.trace(giw[iw0].imag)

            # impurity density from double occupations.
            # TODO: this is tied to diagonal hybr
            dens_impurity = orbspin.promote_diagonal(
                        orbspin.extract_diagonal(occ))
            dens_impurity = np.array(dens_impurity).astype(np.float)
            dens_impurity = np.sum(dens_impurity)

            # we do not use the non interacting density from Weiss field
            # but that from the DFT Hamiltonian.
            # that is a lot more stable

            # impurity non-interacting density (from Weiss field)
            #iw = tf.matfreq(beta, 'fermi', imp_result.problem.niwf)
            #eye = np.eye(imp_result.problem.norbitals * imp_result.problem.nspins) \
            #        .reshape(imp_result.problem.norbitals, imp_result.problem.nspins,
            #                 imp_result.problem.norbitals, imp_result.problem.nspins)
            #g0iw = orbspin.invert(imp_result.problem.g0iwinv)
            #g0iw_model = orbspin.invert(1j * iw * eye + imp_result.problem.muimp)
            #dens_model = lattice.fermi(-imp_result.problem.muimp) # FIXME offdiag
            #dens0_impurity = (g0iw - g0iw_model).real.sum(0) / imp_result.problem.beta \
            #                 + dens_model
            if ineq.occ_dc_number is not None:
                dens0_impurity = ineq.occ_dc_number
            else:
                dens0_impurity = occ0
            
            # If the DOS at the Fermi level is large, we do not want to adjust
            # too much for stability reasons (fudging)
            if dos_fermi_level > 1.:
                delta_dc = (dens0_impurity - dens_impurity)/dos_fermi_level
            else:
                delta_dc = dens0_impurity - dens_impurity

            # Fudging factor
            delta_dc *= 0.2

            delta_dc = delta_dc * np.eye(ineq.nd * self.nspins).reshape(
                                    ineq.nd, self.nspins, ineq.nd, self.nspins)
            ineq.d_setpart(self.dc_value, ineq.d_downfold(self.dc_value) - delta_dc )
        return self.dc_value



# Experimental non-benchmarked DC Method!
class Fixed_dp_Distance():
  def get(self, dc, siw_dd, smom_dd, iwf, natoms, fixed_orbitals=None):
    print "    *** FIXING t2g-p ORBITAL DISTANCE ***"
    
    if fixed_orbitals == None:
      print '    ==> Fixing standard orbitals: 1, 3, 5'
      fixed_orbitals = np.array([0,2,4])
    
    # SINGLE ORBITAL SELECTION
    if np.isscalar(fixed_orbitals):
      print '    Fixing to orbital:', fixed_orbitals
      no_iwn = np.asarray(siw_dd).shape[1]
      siw_t2g_real_extrapolated = 0
      #print '    ==> Using EXTRAPOLATION!'
      #siw_t2g_real_extrapolated = np.real(siw_dd[0][int(no_iwn/2),fixed_orbitals,1,fixed_orbitals,1]) \
            #- ( np.real(siw_dd[0][int(no_iwn/2)+1,fixed_orbitals,1,fixed_orbitals,1]) - np.real(siw_dd[0][int(no_iwn/2),fixed_orbitals,1,fixed_orbitals,1]) ) \
            #/ ( iwf[int(no_iwn/2)+1] - iwf[int(no_iwn/2)]) * iwf[int(no_iwn/2)]
      #siw_t2g_real_avg = siw_t2g_real_extrapolated

      ## Alternative p+d DC Method
      #print '    ==> Using first Sigma value!'
      #siw_t2g_real_avg = np.real(siw_dd[0][int(no_iwn/2),fixed_orbitals,1,fixed_orbitals,1])

      ## Alternative p+d DC Method
      #print '    ==> Using Sigma oo value!'
      #siw_t2g_real_avg = np.real(siw_dd[0][int(0),fixed_orbitals,1,fixed_orbitals,1])
      
    # MULTIPLE ORBITAL SELECTION
    else:
      print '    Fixing to orbitals:', fixed_orbitals
      no_iwn = np.asarray(siw_dd).shape[1]

      # Variant 1
      #print '    ==> Using EXTRAPOLATION!'
      #siw_t2g_real_extrapolated = np.zeros(np.asarray(fixed_orbitals).shape[0])
      #for b in range(0,np.asarray(fixed_orbitals).shape[0]):
        #siw_t2g_real_extrapolated[b] = np.real(siw_dd[0][int(no_iwn/2),fixed_orbitals[b],1,fixed_orbitals[b],1]) \
              #- ( np.real(siw_dd[0][int(no_iwn/2)+1,fixed_orbitals[b],1,fixed_orbitals[b],1]) - np.real(siw_dd[0][int(no_iwn/2),fixed_orbitals[b],1,fixed_orbitals[b],1]) ) \
              #/ ( iwf[int(no_iwn/2)+1] - iwf[int(no_iwn/2)]) * iwf[int(no_iwn/2)]
      #siw_t2g_real_avg = 0
      #for b in range(0,np.asarray(fixed_orbitals).shape[0]):
        #siw_t2g_real_avg = siw_t2g_real_avg + siw_t2g_real_extrapolated[b]
      #siw_t2g_real_avg = siw_t2g_real_avg/np.asarray(fixed_orbitals).shape[0]    

      ## Alternative p+d DC Method
      #siw_t2g_real_avg = 0
      #for b in range(0,np.asarray(fixed_orbitals).shape[0]):
      #  print '    ==> Using first Sigma value!'
      #  siw_t2g_real_avg = siw_t2g_real_avg + np.real(siw_dd[0][int(no_iwn/2),fixed_orbitals[b],1,fixed_orbitals[b],1])
      #siw_t2g_real_avg = siw_t2g_real_avg/np.asarray(fixed_orbitals).shape[0]

      ## Alternative p+d DC Method
      #siw_t2g_real_avg = 0
      #print '    ==> Using Sigma oo value!'
      #for b in range(0,np.asarray(fixed_orbitals).shape[0]):
       #siw_t2g_real_avg = siw_t2g_real_avg + np.real(siw_dd[0][int(0),fixed_orbitals[b],1,fixed_orbitals[b],1])
      #siw_t2g_real_avg = siw_t2g_real_avg/np.asarray(fixed_orbitals).shape[0]

      # Z Factor Method
      siw_t2g_real_avg = 0
      for b in range(0,np.asarray(fixed_orbitals).shape[0]):
       print '    ==> Using first Sigma value times Z-factor!'
       siw_t2g_real_avg = siw_t2g_real_avg + np.real(siw_dd[0][int(no_iwn/2),fixed_orbitals[b],1,fixed_orbitals[b],1])
      siw_t2g_real_avg = siw_t2g_real_avg/np.asarray(fixed_orbitals).shape[0]

      for b in range(0,np.asarray(fixed_orbitals).shape[0]):
        alpha = -( np.real(siw_dd[0][int(no_iwn/2+1),fixed_orbitals[b],1,fixed_orbitals[b],1]) - np.real(siw_dd[0][int(no_iwn/2),fixed_orbitals[b],1,fixed_orbitals[b],1]) ) / ( iwf[int(no_iwn/2+1)] - iwf[int(no_iwn/2)] )
        Z = 1/(1+alpha)
        print '    ==> alpha:', alpha, ', Z:', Z
      siw_t2g_real_avg = Z * siw_t2g_real_avg


    a=0
    for a in range(0,natoms):
      for b in range(0,np.asarray(siw_dd).shape[2]):
        idx = np.int(b + a*dc.shape[0]/natoms)
        #print '              DC Bands:', idx, a
        dc[idx,0,idx,0] = siw_t2g_real_avg
        dc[idx,1,idx,1] = siw_t2g_real_avg
    
    
    return dc
