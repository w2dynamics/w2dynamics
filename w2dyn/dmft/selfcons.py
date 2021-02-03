from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from warnings import warn
import numpy as np
import scipy.optimize as opt
import sys

#import w2dyn.dmft.lattice as latt
import w2dyn.dmft.impurity as impurity
import w2dyn.dmft.doublecounting as doublecounting
import w2dyn.dmft.orbspin as orbspin
import w2dyn.dmft.gw as gw
import w2dyn.dmft.mixing as mixing

import w2dyn.auxiliaries.transform as tf

def iw_to_tau_fast(aiw, ntau, beta, axis=-1):
    """Fast Fourier transform from fermionic Matsubara to tau grid [0, beta]

    This functions takes a quantity on the fermionic Matsubara axis `aiw`,
    where all frequencies from n = {-n, ..., n-1} have to be given, and
    transforms it to the imaginary time axis, where both endpoints [0,beta]
    are included, by mapping it to a fast discrete Fourier transform (FFT).
    """
    if axis != -1:
        aiw = np.rollaxis(aiw, axis, aiw.ndim)
    niw = aiw.shape[-1]

    if ntau >= niw:
        fft_ntau = ntau
    elif ntau > 1:
        overprod_factor = int(np.ceil((niw - 1)/(ntau - 1)))
        fft_ntau = overprod_factor * (ntau - 1) + 1
    else:
        raise ValueError("tau axis must have more than 1 point")

    # We work over the doubled array. This corresponds to both fermionic
    # and bosonic frequencies in Fourier space and to periodic functions
    # on an interval of 2 * beta.  This allows us to use the standard
    # FFT algorithms, but we can then only use (1) half the result interval,
    # (2) every other result point, since we get the signal and its
    # mirror image.
    n = 2 * (fft_ntau - 1)
    a = np.zeros(aiw.shape[:-1] + (2*n,), aiw.dtype)
    a[..., n-niw+1:-n+niw:2] = aiw
    atau = np.fft.fft(a)
    atau = 1./beta * atau[...,:2*fft_ntau:2]
    if ntau < fft_ntau:
        atau = atau[..., 0::overprod_factor]

    assert atau.shape == aiw.shape[:-1] + (ntau,)
    if axis != -1:
        atau = np.rollaxis(atau, -1, axis)
    return atau


class FrequencyDistribution:
    def __init__(self, mpi_comm, niw):
        self.comm = mpi_comm
        self.rank = mpi_comm.Get_rank()
        self.size = mpi_comm.Get_size()
        self.niw = None
        self.set_niw(niw)

    def rdebug(self, fmt, *params):
        """Write debugging message to STDERR, but only on root"""
        if self.rank == 0:
            sys.stderr.write("MPI: %s\n" % (str(fmt) % params))

    def set_niw(self, niw):
        if self.niw == niw: return
        rank = self.comm.Get_rank()
        nranks = self.comm.Get_size()
        iw_per_rank = niw // nranks
        iw_excess = niw - iw_per_rank * nranks

        self.niw = niw
        self.sizes = iw_per_rank * np.ones(nranks, int)
        if iw_excess:
            self.sizes[-iw_excess:] += 1

        slice_ends = self.sizes.cumsum()
        self.slices = list(map(slice, slice_ends - self.sizes, slice_ends))

        self.myniw = self.sizes[rank]
        self.myslice = self.slices[rank]
        self.myshare = float(self.myniw) / self.niw
        self.rdebug("distributing %d frequencies (%d per rank)",
                    niw, iw_per_rank)

    def allgather(self, rank_result):
        if self.myniw != rank_result.shape[0]:
            raise ValueError("Illegal frequency gathering operation")
        tot_shape = (self.niw,) + rank_result.shape[1:]
        tot_result = np.empty(tot_shape, rank_result.dtype)
        tot_result[...] = np.nan
        other_dims = np.prod(rank_result.shape[1:])

        # The sizes argument needs the total number of elements rather than
        # just the first axis. The type argument is inferred.
        self.comm.Allgatherv(rank_result,
                             [tot_result, self.sizes*other_dims])
        return tot_result
    
    def shape(self, rank_result):
        return (self.niw,) + rank_result.shape[1:]
    
    def parts(self, rank_result):
        for rank, curr_slice in enumerate(list(self.slices)):
            curr_result = self.comm.bcast(rank_result, rank)
            yield curr_slice, curr_result

class DMFTStep:
    def __init__(self, beta, lattice, ineq_list, niwf, nftau, dc_dp, dc_dp_orbitals, GW, GW_KAverage, natoms, dc=None,
                 udd_full=None, udp_full=None, upp_full=None, paramag=False,
                 siw_mixer=None, mu_mixer=None, mpi_comm=None):

        if beta < 0: raise ValueError("beta must be positive")
        if niwf <= 0 or niwf % 2 != 0: raise ValueError("niwf must be even")
        if nftau <= 1: raise ValueError("nftau must be greater than 1")
        if dc is None: dc = doublecounting.Zero()
        if siw_mixer is None: siw_mixer = mixing.FlatMixingDecorator(mixing.LinearMixer())
        if mu_mixer is None: mu_mixer = mixing.LinearMixer()

        self.mpi_comm = mpi_comm

        self.beta = beta
        self.lattice = lattice
        self.ineq_list = ineq_list
        self.dc = dc
        self.udd_full = udd_full
        self.udp_full = udp_full
        self.upp_full = upp_full
        self.paramag = paramag

        if self.udp_full is None or upp_full is None:
            self.use_hartree = False
        else:
            self.use_hartree = udp_full.any() or upp_full.any()

        self.siw_mixer = siw_mixer
        self.mu_mixer = mu_mixer

        self._eye = lattice.eye

        self.niwf = niwf
        self.iwf = tf.matfreq(beta, 'fermi', niwf)
        
        if self.mpi_comm is not None:
            self.mpi_strategy = FrequencyDistribution(mpi_comm, self.niwf)
            self.my_slice = self.mpi_strategy.myslice
        else:
            self.my_slice = slice(None)
        
        self.my_iwf = self.iwf[self.my_slice]
        self.my_niwf = self.iwf.size
        self.nftau = nftau

        # build transformation yourself since it potentially only contains part
        # of the frequency axis
        self.tauf = np.arange(nftau)*self.beta/(nftau - 1.0)
        
        # Fixing d-p distance
        self.dc_dp = dc_dp
        self.dc_dp_orbitals = dc_dp_orbitals
        self.dc_dp_shift = None

        # GW Inclusion
        self.use_gw = GW
        self.use_gw_kaverage = GW_KAverage

        if self.use_gw == 1:
            try:
                dummy_try = self.lattice.nkpoints
            except:
                print('********* GW MODULE IS NOT AVAILABLE FOR BETHE LATTICE *********')
                print('...exiting')
                exit()
            self.gw = gw.GWInclusion(self.lattice.norbitals, natoms, self.lattice.nkpoints, self.my_niwf, self.paramag, self.use_gw_kaverage)

        self.siw_gw = None
        self.smom_gw = None
        self.natoms = natoms

    def set_siws(self, siw_dd=None, smom_dd=None, dc_full=None, init=False,
                 hartree_start=False,giws=None,occs=None):
        new_run = siw_dd is None and init

        nspins = self.lattice.nspins
        norbitals = self.lattice.norbitals

        if siw_dd is None:
            siw_dd = []
            smom_dd = []
            for ineq in self.ineq_list:
                siw_block = np.zeros((self.niwf, ineq.nd, nspins, ineq.nd, 
                                      nspins), np.complex)
                if not self.paramag:
                    siw_block[:,:,0,:,0] += .5 * ineq.se_shift
                    siw_block[:,:,1,:,1] -= .5 * ineq.se_shift

                siw_dd.append(siw_block)
                smom_dd.append(siw_block[:2].real)

        if smom_dd is None:
            warn("extracting moments from Sigma. Better use siw_mom",
                 UserWarning, 2)
            smom_dd = [siw_block[:2].real for siw_block in siw_dd]

        # Enforce paramagnetic option
        if self.paramag:
            siw_dd = [orbspin.symm_spins(siw_block) for siw_block in siw_dd]
            smom_dd = [orbspin.symm_spins(smom_block) for smom_block in smom_dd]

        # get double-counting if we do not force it to a value.
        # Note that the double-counting is essentially a correction to the
        # Hamiltonian and therefore present from the start.


        
        # fix for self consitent dc: If self.dc_full is present
        # (from last iteration or from init), then do not call the
        # standard dc call here but do it later
        # this implies that set_siws() can initialize dc with dc_full
        # exaclty once, which makes sense

        try:
            self.dc_full
        except AttributeError:
            if dc_full is None:
                self.dc_full = self.dc.get()
            else:
                self.dc_full = dc_full
                # We also have to set the double counting inside the
                # dc class so that trace double counting works
                self.dc.set(dc_full)

        if self.use_gw == 1:
            self.dc_full = self.dc_full * 0

        # Try to mimic the behaviour of the old code
        if new_run and hartree_start:
            self.sigma_hartree = -self.dc_full
            for ineq_no, ineq in enumerate(self.ineq_list):
                siw_dd[ineq_no][:] += ineq.d_downfold(self.sigma_hartree)
                ineq.d_setpart(self.sigma_hartree, 0)

        # Perform mixing of self-energy, its moments and
        # double-counting with the same mixer. This is important to
        # keep them "consistent", i.e., double-count correlations at
        # the same rate as switching them on, and with DIIS, the
        # mixing of all related quantities using the same mixer is
        # necessary in particular because the earlier values
        # themselves influence the ratios in which they are mixed in.
        self.dc_full, self.siw_dd, self.smom_dd = self.siw_mixer(self.dc_full, siw_dd, smom_dd)

        # Fix the distance between selected d-orbitals and the p-manifold to the original distance of the LDA/GW Hamiltonian
        if self.dc_dp == 1:
            self.dc_full = self.dc_full * 0
            dp_dc = doublecounting.Fixed_dp_Distance()
            self.dc_full = -dp_dc.get(self.dc_full, self.siw_dd, self.smom_dd, self.iwf, self.natoms, self.dc_dp_orbitals)

        # not doing mixing for self consistent dc
        # if we do want to mix... simply move this part above the mixing routine
        # if dc_full is not None, it comes from a restarted run. We do not want to
        # overwrite the restart value here!
        if self.dc.self_cons and ( dc_full is None ):
            # this currently should only work for trace dc
            # todo: fix siginfbar
            self.dc_full = self.dc.get(siws=siw_dd, smoms=smom_dd, giws=giws,occs=occs)
        # Upfold self-energy.
        self.siw_full = 0
        self.siw_moments = 0
        for ineq, siw_bl, siw_bl_mom in zip(self.ineq_list, self.siw_dd,
                                            self.smom_dd):
            # Here, we potentially encounter a memory problem (a set of
            # nflavours x nflavours x nfrequencies) object may very well
            # overflow, so we split the upfolded self-energy in frequency 
            # "chunks"
            self.siw_full += ineq.d_upfold(siw_bl[self.my_slice])
            self.siw_moments += ineq.d_upfold(siw_bl_mom)          

        # Add double-counting to the self-energy
        self.siw_full += self.dc_full
        self.siw_moments[0] += self.dc_full

        # This is an approximation: really, the partial densities and Hartree
        # self-energy an inter-dependent system. Here, we assume that the
        # partial densities won't change much adding  a static self-energy
        # shift...
        if not (hartree_start and new_run):
            if self.use_hartree:
                #self.siw2gloc() # computes densities as well
                hartree = np.einsum("isjt,jt->is", self.udp_full + self.upp_full,
                                    self.densities)
                self.sigma_hartree = orbspin.promote_diagonal(hartree)
            else:
                self.sigma_hartree = np.zeros(self.siw_full.shape[1:],
                                              dtype=self.siw_full.dtype)

        self.siw_full += self.sigma_hartree
        self.siw_moments[0] += self.sigma_hartree.real

        # Invalidate lattice quantities (call to siw2gloc necessary)
        self.glociw = None
        self.densities = None
        self.densmatrix = None

        # Only if the DMFT Self-Energy is non-zero, we add the GW Self-Energy
        if (self.use_gw == 1) and (np.sum(self.siw_full) != 0):        # siw_full is d+p siw_dd is d-only
            self.siw_gw, self.smom_gw = self.gw.get_GW_Sigma(self.beta, self.my_iwf)
            self.siw_gw = self.siw_gw.transpose(1,0,2,3,4,5)
            

    def update_mu(self, target_totdens, epsn, search_range):
        if epsn <= 0: raise ValueError("epsn must be greater than 0")
        if search_range <= 0: raise ValueError("search range must be > 0")

        # Again consider only frequencies of the own process
        trace_factory = self.lattice.trace(self.my_iwf, self.siw_full, self.siw_moments, 'local', self.use_gw, self.siw_gw, self.smom_gw)


        def root_func(mu):
            """f(mu) that has a root (zero) at the desired density"""
            trgloc = trace_factory.trgloc(mu)
            trgmodel, dmodel = trace_factory.trgmodel(mu)

            totdens = (trgloc - trgmodel).sum()/self.beta
            if self.mpi_comm is not None:
                totdens = self.mpi_comm.allreduce(totdens)
            totdens += dmodel
            
            if totdens.imag > epsn: warn("non-vanishing imaginary part")
            return totdens.real - target_totdens

        mu = opt.brentq(root_func, self.mu - search_range,
                        self.mu + search_range, xtol=epsn)
        self.set_mu(mu)

    def set_mu(self, mu, init=False):
        self.mu = self.mu_mixer(mu)
        self.trace = None

    def siw2gloc(self):
        self.glociw = self.lattice.gloc(self.my_iwf, self.mu, self.siw_full, self.siw_gw)
        self.gmodeliw, dmatrix = self.lattice.gmodel(self.my_iwf, self.mu, self.siw_moments,  self.smom_gw)

        self.densmatrix = (self.glociw - self.gmodeliw).sum(0)/self.beta
        if self.mpi_comm is not None:
            self.densmatrix = self.mpi_comm.allreduce(self.densmatrix)
        self.densmatrix += dmatrix
        
        self.densities = orbspin.extract_diagonal(self.densmatrix)
        # Destroy large GW siw for MC cycle 
        self.siw_GW_dyn = None
        #self.smom_GW = None


    def write_before_mu_search(self, output):
        # This quantity is potentially distributed over cores, so we need to
        # take that into account here.
        if self.mpi_comm is not None:
            output.write_distributed("glocold",
                                  self.mpi_strategy.parts(self.glociw),
                                  self.mpi_strategy.shape(self.glociw))
            glociw = self.mpi_strategy.allgather(orbspin.extract_diagonal(self.glociw))
            output.write_quantity("glocold-lattice", glociw)
        else:
            output.write_quantity("glocold", self.glociw)
            output.write_quantity("glocold-lattice", orbspin.extract_diagonal(self.glociw))

        output.write_quantity("gdensold", self.densmatrix)

    def write_lattice_problem(self, output):
        # This quantity is potentially distributed over cores, so we need to
        # take that into account here.
        if self.mpi_comm is not None:
            output.write_distributed("glocnew",
                                  self.mpi_strategy.parts(self.glociw),
                                  self.mpi_strategy.shape(self.glociw))
            glociw = self.mpi_strategy.allgather(orbspin.extract_diagonal(self.glociw))
            output.write_quantity("glocnew-lattice", glociw)
        else:
            output.write_quantity("glocnew", self.glociw)
            output.write_quantity("glocnew-lattice", orbspin.extract_diagonal(self.glociw))

        output.write_quantity("dc", self.dc_full)
        output.write_quantity("dc-latt", self.dc_full)
        if self.use_hartree:
            output.write_quantity("sigma-hartree", self.sigma_hartree)

        # FIXME: hack; should got into some write routine of ShellAveraged
        if isinstance(self.dc, doublecounting.ShellAveraged):
            output.write_quantity("ubar", self.dc.ubar)
            output.write_quantity("jbar", self.dc.jbar)

        # The rest is consistent on all nodes
        output.write_quantity("mu", self.mu)
        output.write_quantity("gdensnew", self.densmatrix)

    def gloc2fiw(self):
        self.imp_problems = []

        if (self.use_gw == 1) and (np.sum(self.siw_full) != 0):
            muimp = (self.mu * self._eye - self.dc_full - self.lattice.hloc.real - self.sigma_hartree - np.sum(self.smom_gw, axis=0)/self.lattice.nkpoints)
        else:
            muimp = (self.mu * self._eye - self.dc_full - self.lattice.hloc.real - self.sigma_hartree)
          
        atom = 0
        for ineq, siw_block, smom_block in zip(self.ineq_list, self.siw_dd, self.smom_dd):
            # extract d-d blocks and perform inversion.
            giw_block = ineq.d_downfold(self.glociw)
            g0inviw_block = orbspin.invert(giw_block)

            # Now we "undo" the frequency splitting in order to create
            # the same impurity problem for each core
            if self.mpi_comm is not None:
                g0inviw_block = self.mpi_strategy.allgather(g0inviw_block)

            # Remove self-energy using the Dyson equation.  Note that in 
            # d+p or general multi-site, the d-d block is extracted *before*
            # the inversion and therefore the impurity self-energy is used.
            g0inviw_block += siw_block

# giorgio: Symmetrisierung von G0inviw
            g0inviw_block.real[:,...] = 0.5*(g0inviw_block.real[:,...]+g0inviw_block.real[::-1,...])
            g0inviw_block.imag[:,...] = 0.5*(g0inviw_block.imag[:,...]-g0inviw_block.imag[::-1,...])

            nbands=g0inviw_block.shape[1]
            niw=g0inviw_block.shape[0]
            g0inviw_block=g0inviw_block.reshape(niw,nbands*2,nbands*2)

            g0_new=np.zeros_like(g0inviw_block,dtype=complex)
            for i in range(0,niw):
               tmp=g0inviw_block[i,:,:]
               g0_new[i,:,:]=0.5*(tmp.transpose(1,0)+tmp)

            g0_new=g0_new.reshape(niw,nbands,2,nbands,2)
            g0inviw_block=g0_new

# giorgio   ######

            # in principle, since the bare impurity propagator reads:
            # `1/G_0(iw) =  iw - muimp - F(iw)`, we have a freedom where to put
            # constant shifts.  However, constant terms in `F(iw)` translate to
            # an unresolvable delta peak at `F(tau=0)`. Therefore, we want the
            # zero-th moment of `F(iw)` to vanish, which in turn fixes `muimp`.
            eye_block = ineq.d_downfold(self._eye)

            # The first moment of the hybridisation function is given as
            # - P <H>^2 P + (P <H> P)^2, where P is the downfolding projector,
            # (see Eq. (4.63b) in Markus' thesis), so we cannot simply downfold
            # the second moment of the DOS.
            hloc_block = ineq.d_downfold(self.lattice.hloc)

            # GW Inclusion and Inclusion of GW 0th Moment  
            if (self.use_gw == 1) and (np.sum(self.siw_full) != 0):
                norbitals_per_atom = np.int(self.lattice.norbitals/self.natoms)
                orbits = slice( atom*norbitals_per_atom , (atom+1)*norbitals_per_atom )
                # => Model based
                if self.use_gw_kaverage == 0:    
                    smom_gw2 = np.sum(self.smom_gw[:,orbits,:,orbits,:] * self.smom_gw[:,orbits,:,orbits,:], axis=0)/self.lattice.nkpoints
                    fmom1_block = (np.einsum("isjt,jtku->isku", hloc_block, hloc_block) - ineq.d_downfold(self.lattice.hmom2)) #+ C1 - C2 - C3 - C4 - smom_gw2  # (-) CORRECT HERE!
                    # => K-average based
                else:                                 
                    smom_gw2 = np.sum(self.smom_gw[:,orbits,:,orbits,:] * self.smom_gw[:,orbits,:,orbits,:], axis=0)/self.lattice.nkpoints
                    fmom1_block = (np.einsum("isjt,jtku->isku", hloc_block, hloc_block) - ineq.d_downfold(self.lattice.hmom2)) - smom_gw2    # MINUS SIGN IS CORRECT HERE!
            # No GW in first iteration
            else:
                fmom1_block = (np.einsum("isjt,jtku->isku", hloc_block, hloc_block) - ineq.d_downfold(self.lattice.hmom2))
            
            muimp_block = ineq.d_downfold(muimp)

            # fiw is defined in the bath picture, which means that:
            #      G0(iw)^{-1} = iw + muimp - F(-iw);   F(-iw) = F*(iw)
            # TODO: check if that is still true with complex F(tau)
            fiw_block = -1j*self.iwf[:,None,None,None,None] * eye_block \
                            + muimp_block - g0inviw_block.conj()
            fiw_model = fmom1_block/(1j*self.iwf[:,None,None,None,None])

            # GW Inclusion 
            if (self.use_gw == 1) and (np.sum(self.siw_full) != 0) and (self.use_gw_kaverage == 0):
                extension_of = 1
                fiw_block_asymptotic, fiw_model_asymptotic, model_mom1 = self.gw.Reach_asymptotics_for_Hybridization(fiw_block, fiw_model, self.iwf, extension_of)
                ftau_block = iw_to_tau_fast(fiw_block_asymptotic - fiw_model_asymptotic, self.nftau, self.beta, axis=0) - model_mom1/2.
                fmom_block = np.asarray((np.zeros_like(model_mom1), model_mom1))
            else:
                ftau_block = iw_to_tau_fast(fiw_block - fiw_model, self.nftau, self.beta, axis=0) - fmom1_block/2.
                fmom_block = np.asarray((np.zeros_like(fmom1_block), fmom1_block))
            
            imp_problem = impurity.ImpurityProblem(
                    self.beta, g0inviw_block, fiw_block, fmom_block, ftau_block,
                    muimp_block, ineq.dd_int, None, None, ineq.symmetry_moves,
                    self.paramag)
            self.imp_problems.append(imp_problem)
            atom = atom + 1
        
    def write_imp_problems(self, output):
        g0iw = [orbspin.invert(p.g0inviw) for p in self.imp_problems]

        output.write_quantity("g0iw-full", g0iw)
        output.write_quantity("g0iw", g0iw)
        output.write_quantity("fiw", [p.fiw for p in self.imp_problems])
        output.write_quantity("fmom", [p.fmom for p in self.imp_problems])
        output.write_quantity("ftau", [p.ftau for p in self.imp_problems])
        output.write_quantity("ftau-full", [p.ftau for p in self.imp_problems])
        output.write_quantity("muimp", [p.muimp for p in self.imp_problems])
        output.write_quantity("muimp-full", [p.muimp for p in self.imp_problems])

class KappaMuExtrapolator:
    def __init__(self, totdens, last_mudiff=None, forced_kappa_sign=None):
        self.totdens = totdens
        if forced_kappa_sign is not None and np.sign(forced_kappa_sign) == 0:
            raise ValueError("Forced sign of kappa must not be 0")
        self.forced_kappa_sign = (np.sign(forced_kappa_sign)
                                  if forced_kappa_sign is not None
                                  else None)
        self.kappa = None
        self.mus = []
        self.ns = []
        self.last_n = None
        self.last_mudiff = last_mudiff
        self.last_kappa = None
        self.kappa_sign_volatile = False
        self.max_inc = 5.0
        self.max_dec = 10.0

    def initialize(self, last_mu, last_n, start_kappa=None, last_mudiff=None):
        self.last_mu = last_mu
        self.last_n = last_n
        self.mus.append(last_mu)
        self.ns.append(last_n)
        if (start_kappa is not None
           and (self.forced_kappa_sign is None
                or np.sign(start_kappa) == self.forced_kappa_sign)):
            self.kappa = start_kappa
        self.last_mudiff = last_mudiff

    def has_mu(self):
        return self.kappa is not None and self.last_n is not None

    def next_mu(self):
        diff = (self.totdens - self.last_n) / self.kappa

        if self.last_mudiff is None:
            self.last_mudiff = diff
            return self.last_mu + diff

        # heuristic: if the sign of kappa changed, force the first
        #            step with the new kappa to be smaller to prevent
        #            DMFT lag and error from erasing too much
        #            convergence progress
        max_inc = self.max_inc
        if self.forced_kappa_sign is None:
            if (self.last_kappa is not None
                    and np.sign(self.last_kappa) != np.sign(self.kappa)):
                if not self.kappa_sign_volatile:
                    max_inc = 1.0/3.0
                else:
                    max_inc = 3.0 * max_inc
                self.kappa_sign_volatile = not self.kappa_sign_volatile
            else:
                self.kappa_sign_volatile = False

        # heuristic: limit ratio of successive step size changes (to
        #            prevent wild swings due to outliers and freezing)
        sign = (np.sign(diff) if diff != 0.0 else 1.0)
        if abs(diff) > max_inc * abs(self.last_mudiff):
            diff = sign * max_inc * abs(self.last_mudiff)
        elif abs(diff) < abs(self.last_mudiff) / self.max_dec:
            diff = sign * abs(self.last_mudiff) / self.max_dec

        self.last_mudiff = diff

        return self.last_mu + diff

    def step(self, mu, n, nerr=None):
        # reduce max step size change factor if within error or after
        # first crossing the target value
        if nerr is not None and (abs(n - self.totdens) <= nerr):
            self.max_inc = 2.0
            self.max_dec = 5.0
        elif ((nerr is None or np.isnan(nerr)) and self.last_n is not None
              and (self.totdens - self.last_n) * (self.totdens - n) < 0.0):
            self.max_inc = 2.0
            self.max_dec = 5.0

        self.last_kappa = self.kappa
        if self.last_n is not None:
            kappa = (n - self.last_n)/(mu - self.last_mu)
            if (self.forced_kappa_sign is None
               or np.sign(kappa) == self.forced_kappa_sign):
                self.kappa = kappa

        self.last_mu, self.last_n = mu, n
        self.mus.append(self.last_mu)
        self.ns.append(self.last_n)

