""" @package postprocessing

Provides methods for annotation and post-processing of QMC output scripts.
"""
from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
import numpy as np
import inspect
from functools import wraps
import math

def divide(enum, denom, enum_err=None, denom_err=None):
    """divides and propagates the uncertainty, assuming no covariance"""
    has_enum_err = enum_err is not None
    has_denom_err = denom_err is not None
    if has_enum_err or has_denom_err:
        if not has_enum_err: enum_err = np.zeros_like(denom_err)
        if not has_denom_err: denom_err = np.zeros_like(enum_err)
        ratio_err = np.abs(enum_err/enum)**2 + np.abs(denom_err/denom)**2
        return enum/denom, ratio_err
    else:
        return enum/denom
    
def multiply(left, right, left_err=None, right_err=None):
    """multiplies and propagates the uncertainty, assuming no covariance"""
    has_left_err = left_err is not None
    has_right_err = right_err is not None
    if has_left_err or has_right_err:
        if not has_left_err: left_err = np.zeros_like(right_err)
        if not has_right_err: right_err = np.zeros_like(left_err)
        product_err = 0 
        return left*right, product_err
    else:
        return left*right
 
def fslice(x, minval, maxval, stride=None):
    """returns a slice of a sorted float interval"""
    if minval is not None: minval = np.searchsorted(x, minval)
    if maxval is not None: maxval = np.searchsorted(x, maxval)
    if stride is not None: stride = np.round(stride/(x[1]-x[0]))
    return slice(minval, maxval, stride)

def slice_center(subsize, size):
    """Returns a slice for the center elements of a bigger array"""
    if subsize < 0 or size < 0:
        raise ValueError("sizes must be bigger than zero")
    if subsize > size: 
        raise ValueError("subsize > size")
    if subsize % 2 != size % 2:
        raise ValueError("size and subsize are different MOD 2")
    if subsize == size:
        return slice(None)
    else:
        startidx = (size - subsize)//2
        cslice = slice(startidx, startidx + subsize)
        return cslice

def unwind_dim(dim):
    """Allows treatment of inequivalent atoms and such"""
    def unwind_decorator(func):
        @wraps(func)
        def unwind_wrapper(*args):
            adim = args[0].ndim
            if adim == dim - 1:
                return func(*args)
            if adim != dim:
                raise ValueError("expected %d or %d dimensions, got shape: %s"
                                 % (dim - 1, dim, str(args[0].shape)))
            dimsize = args[0].shape[0]
            unwound_args = []
            for arg in args:
                if not np.ndim(arg):
                    arg = (arg,) * dimsize
                unwound_args.append(arg)  

            res = [func(*fargs) for fargs in zip(*unwound_args)]
            # now res contains tuples or arrays
            if isinstance(res[0], tuple):
                nfields = len(res[0])
                res = tuple(np.asarray([elem[i] for elem in res])
                            for i in range(nfields))
                #print "shapes:", [elem.shape for elem in res]
                return res
            else:
                res = np.asarray(res)
                #print "shape:", np.asarray(res).shape
                return res
        return unwind_wrapper
    return unwind_decorator

def dia_slice(A, *dims):
    """ Returns a slicing object to be used for A that extracts diagonals """
    ndim = len(dims)//2
    dims = np.reshape(dims, (ndim, 2))
    slices = [slice(None)] * A.ndim
    for idim, (dim1, dim2) in enumerate(dims):
        shape = int(A.shape[dim1])
        if shape != A.shape[dim2]:
            raise NotImplementedError("not supported yet: %d/%d" % (dim1,dim2) 
                                      + str(A.shape))
        extshape = [1]*ndim
        extshape[idim] = -1
        indices = np.arange(shape).reshape(extshape)
        slices[dim1] = slices[dim2] = indices
    return tuple(slices)
 
# ----- derived quantities -----

def impurity_density(occ, occ_err=None):
    """density of the impurity model"""
    def get_diags(val):
        return np.diagonal(np.diagonal(val, 0, -4, -2), 0, -3, -2)
    if occ_err is not None:
        return get_diags(occ), get_diags(occ_err)
    else:
        return get_diags(occ)
     
def impurity_electrons(occ, occ_err=None):
    """total number of electrons in the impurity model"""
    def get_traces(val):
        return np.trace(np.trace(val, 0, -3, -1), 0, -2, -1)
    if occ_err is not None:
        return get_traces(occ), np.sqrt(get_traces(np.abs(occ_err)**2))
    else:
        return get_traces(occ)

def lattice_electrons(gdensnew, gdensnew_err=None):
    """number of electrons in the lattice model, total and per spin"""
    # (..., orb, sp, orb, sp) -> (..., sp)
    gdensnew = np.diagonal(np.trace(np.real(gdensnew), 0, -4, -2), 0, -2, -1)
    return (gdensnew[..., 0], gdensnew[..., 1], np.sum(gdensnew, axis=-1))

def sigmaiw_improved(gsigmaiw, giw, gsigmaiw_err=None, giw_err=None):
    """self energy in Matsubara expansion from improved estimators"""
    return divide(gsigmaiw, giw, gsigmaiw_err, giw_err)

@unwind_dim(5)
def moment(occ, occ_err=None):
    """spontaneuous local moment"""
    mu = occ[:,0,:,0] - occ[:,0,:,1] - occ[:,1,:,0] + occ[:,1,:,1]
    if occ_err is not None:
        mu_err = np.sqrt((occ_err**2).sum(-1).sum(1))
        return mu, mu_err
    else:
        return mu

def sigmaiw_short(iw, sigmaiw):
    """self energy in Matsubara expansion in the range 0 ... 20"""
    sl = fslice(iw, 0, 20)
    return iw[sl], sigmaiw[...,sl]

def get_tauint(data):
    """Bias-corrected logarithmic binning analysis"""
    # See arXiv:1810.05079v2, Sec. III.B
    ndat, nm = data.shape
    m = 2**np.arange(nm)
    lbap = m[None,:-1] * (4 * data[:,1:] - data[:,:-1]) / data[:,:1]
    tauint = (lbap - 1)/2
    return tauint

def get_spectrum(data, tau=None):
    """Perform logarithmic spectral analysis"""
    import scipy.optimize as sp_opt

    # See arXiv:1810.05079v2, Sec. IV.B
    def corr_wedge_conv(tau, m):
        alpha = np.exp(-1/tau)
        return alpha/(1 - alpha)**2 * (1 - alpha**m)**2/m

    # Cut away infinite variance points
    trusted = np.isfinite(data).all(0).nonzero()[0].max() + 1
    data = data[:,:trusted]
    # get grid: Eq. (33)
    ndat, nm = data.shape
    m = 2**np.arange(nm)
    if tau is None: tau = m[:-1].copy()
    # set up GLS: Eqs. (29) & (30)
    A = corr_wedge_conv(tau[None,:], m[:-1,None])
    b = m[None,:-1] * (2 * data[:,1:] - data[:,:-1])
    sigma_invsqrt = (m[:-1])**(-0.5)
    # GLS to OLS: Eq. (36)
    A = sigma_invsqrt[:,None] * A
    b = sigma_invsqrt[None,:] * b
    # Use NNLS: Eq. (36)
    x = np.array([sp_opt.nnls(A, bi)[0] for bi in b])
    # Divide by norm to pass from autocov to autocorr spectra
    x /= x.sum(1)[:,None]
    return x

def spec_tau_est(x, tau=None):
    """Compute tauint estimate from spectrum"""
    def corr_sum(tau):
        alpha = np.exp(-1/tau)
        return (1 + alpha)/(1 - alpha)

    ndat, ntau = x.shape
    if tau is None: tau = 2**np.arange(ntau)
    tauint_comp = corr_sum(tau)
    tauint = x.dot(tauint_comp)
    tauint = (tauint - 1) / 2
    return tauint

def gtau_lbap(gtau_blocks, gtau_blocks_err=None):
    """Bias-corrected logarithmic binning analysis for G(tau)"""
    # See arXiv:1810.05079v2, Sec. III.B
    nm = gtau_blocks.shape[-1]
    gtau_flat = gtau_blocks.reshape(-1, nm)
    lbap = get_tauint(gtau_flat)
    return lbap.reshape(gtau_blocks.shape[:-1] + (nm-1,))

def gtau_spec(gtau_blocks, gtau_blocks_err=None):
    """Logarithmic autocorrelation spectrum for G(tau)"""
    nm = gtau_blocks.shape[-1]
    gtau_flat = gtau_blocks.reshape(-1, nm)
    spectrum = get_spectrum(gtau_flat)
    return spectrum.reshape(gtau_blocks.shape[:-1] + spectrum.shape[-1:])

def gtau_tauint(gtau_blocks, gtau_blocks_err=None):
    """Integrated autocorrelation time estimate for G(tau)"""
    nm = gtau_blocks.shape[-1]
    gtau_flat = gtau_blocks.reshape(-1, nm)
    tauint = spec_tau_est(get_spectrum(gtau_flat))
    tauint = tauint.reshape(gtau_blocks.shape[:-1])
    return tauint

def meta_from_func(func, axes, fields=["value"], base=None):
    """Returns a meta dictionary based on a function"""
    meta = {}
    meta["axes"] = axes
    meta["fields"] = fields
    meta["func"] = func
    meta["desc"] = inspect.getdoc(func)
    if base is None:
        base = []
        error_mode = False
        for arg in inspect.getfullargspec(func).args:
            if arg.endswith("_err"): 
                error_mode = True
                continue
            if error_mode: raise ValueError("errors have to come last")
            base.append(arg.replace("_","-"))
    meta["base"] = tuple(base)
    return meta
  
@unwind_dim(4)
def get_ggstraight_ph(giw, niw4f):
    """Helper function for getting the straight part from GG
    
    The "straight" disconnected part of the two-particle Green's function is 
    given by:
    
       GG_AB(iv, iv', iw) = G_A(iv) G_B(iv') delta(iw,0)
       
    and is returned as six-dimensional array GG(A,B,iv,iv'), omitting the
    bosonic frequency for brevity.
    """
    nband, nspin, niw = giw.shape 
    iw4f_slice = slice_center(niw4f, niw)
    giw = giw.reshape(-1, niw)[:, iw4f_slice]
    gg_straight = np.tensordot(giw, giw, ((),()))  # i,iv,j,iv'
    gg_straight = gg_straight.transpose(0, 2, 1, 3)
    return gg_straight.reshape(nband, nspin, nband, nspin, niw4f, niw4f)

@unwind_dim(4)
def get_ggcross_ph(giw, niw4b, niw4f=None):
    """Helper function for getting the cross part from GG
    
    The "cross" disconnected part of the two-particle Green's function is 
    given by:
    
       GG_AB(iv, iv', iw) = G_A(iv) G_A(iv+iw) delta(iv,iv') delta(A,B)
    
    and is returned as four-dimensional array GG(A,iv,iw), omitting the index
    and fermionic frequency for the second particle for brevity.
    """ 
    def move_slice(sl, offset):  # Moves slice_center() by some offset
        return slice(sl.start + offset, sl.stop + offset)
    
    nband, nspin, niw = giw.shape 
    if niw4f is None: niw4f = niw - niw4b - 1
    iw4f_slice = slice_center(niw4f, niw)
    iw4b_range = range(-(niw4b//2), niw4b//2+1)
    
    giw = giw.reshape(-1, niw)
    gg_cross =  np.dstack(giw[:, iw4f_slice, None] * 
                          giw[:, move_slice(iw4f_slice, iwbos), None]
                          for iwbos in iw4b_range)  # i,iv,iw
    return gg_cross.reshape(nband, nspin, niw4f, niw4b)

@unwind_dim(4)
def get_g4iw_disconn_ph(giw, niw4b, niw4f):   
    nband, nspin, niw = giw.shape
    c = np.zeros(shape=(nband, nspin, nband, nspin, niw4f, niw4f,niw4b),dtype=complex)
    iw4b0 = niw4b//2
    c[...,iw4b0] = get_ggstraight_ph(giw, niw4f)
    diasl = dia_slice(c, 0, 2, 1, 3, 4, 5) # b,s,b,s,iv,iv',iw
    c[diasl] -= get_ggcross_ph(giw, niw4b, niw4f)
    return c

@unwind_dim(4)
def get_chi_ph(giw, g4iw, giw_err=None, g4iw_err=None):
    """generalised susceptibility (particle-hole channel)"""
    iw4b0 = g4iw.shape[-1]//2
    
    chi = g4iw.copy()
    chi[..., iw4b0] -= get_ggstraight_ph(giw, g4iw.shape[-2])
    return chi

@unwind_dim(4)
def get_g4iw_conn_ph(giw, g4iw, giw_err=None, g4iw_err=None):
    """connected part of susceptiblity (particle-hole channel)"""
    conn = get_chi_ph(giw, g4iw, giw_err, g4iw_err) 
    diasl = dia_slice(conn, 0, 2, 1, 3, 4, 5) # b,s,b,s,iv,iv',iw
    conn[diasl] += get_ggcross_ph(giw, g4iw.shape[-1], g4iw.shape[-2])
    return conn

@unwind_dim(4)
def get_ggstraight_pp(giw, g4iw_pp_shape):
    """Computes GG = G(iv)G(iv')delta(iw',-iv-iv')"""
    #print giw.shape, g4iw_pp_shape
    assert giw.shape[-3:] == g4iw_pp_shape[-5:-2], "Shape mismatch"
    dotnot = ((),())
    nneq = g4iw_pp_shape[0]
    N = g4iw_pp_shape[-3]
    K = g4iw_pp_shape[-1]
    KhpNm1 = + K//2 - N + 1
    # TODO slow
    chi0_pp = np.zeros(shape=g4iw_pp_shape, dtype=complex)
    #chi0_pp[...] = np.nan
    for m in range(N):
        for n in range(N):
            ktarg = KhpNm1 + m + n
            if 0 <= ktarg < K:
                chi0_pp[...,m,n,ktarg] = \
                           np.tensordot(giw[...,m], giw[...,n], dotnot)
    #
    return chi0_pp

@unwind_dim(4)
def get_chi_pp(giw, g4iw_pp, giw_err=None, g4iw_pp_err=None):
    """generalised susceptibility (particle-particle channel)"""
    iw4st = (giw.shape[-1] - g4iw_pp.shape[-2])//2
    iw4sl = slice(iw4st, -iw4st)
    chi0_pp = get_ggstraight_pp(giw[...,iw4sl], g4iw_pp.shape)
    return g4iw_pp - chi0_pp


def get_siw_mom(u_matrix, rho1, rho2):
   """Calculates the high frequency moments of the self-energy"""

   # reshape the arrays to contracted band and spin indices
   nbands=u_matrix.shape[0]//2

   rho1_contracted=rho1.reshape((nbands*2),(nbands*2))
   rho2_contracted=rho2.reshape((nbands*2),(nbands*2),(nbands*2),(nbands*2))

   # Use antisymmetrisized U-matrix
   u_antisym=-0.5*(u_matrix-np.swapaxes(u_matrix,0,1)-np.swapaxes(u_matrix,2,3)+np.swapaxes(np.swapaxes(u_matrix,0,1),2,3))

   sigma_inf=np.tensordot(u_antisym,rho1_contracted,axes=([1,2],[0,1]))

   # For formulas see Wang, Dang, Millis PRB 84, 073104
   L=-0.5*(np.tensordot(u_antisym,u_antisym,axes=([2,3],[0,1])))

   first_term=0.25*np.tensordot(u_antisym,u_antisym,axes=(1,2))
   second_term=-np.tensordot(u_antisym,u_antisym,axes=(3,0))

   k1=np.tensordot(first_term,rho2_contracted,axes=([1,2,3,4],[2,3,0,1]))
   k2=np.tensordot(second_term,rho2_contracted,axes=([1,2,3,4],[0,2,1,3]))

   K=np.add(k1,k2)

   sigma_1=-np.add(\
         np.add(K,\
         np.tensordot(L,rho1_contracted,axes=([1,2],[0,1]))),\
         -np.tensordot(sigma_inf,sigma_inf,axes=(1,0)))

   # collect computed moments
   smoms = np.asarray([sigma_inf, sigma_1])
   smoms = np.diagonal(smoms, axis1=1, axis2=2).copy()
   smoms = smoms.transpose().reshape(nbands, 2, -1)
   return smoms






def fix_the_legendre(legs, siw_mom_inf, siw_mom_1, nleg_order, beta, muimp):
   from scipy.misc import factorial

   np.set_printoptions(threshold='nan')
   nleg=legs.shape[2]
   nbands = siw_mom_inf.shape[0]//2

   #print "nleg_order"
   #print nleg_order

   #print "siw_mom_inf"
   #print "siw_mom_1"
   #print siw_mom_inf
   #print siw_mom_1

   for i in range(0,nbands):
      for j in range(0,2):
         for k in range(0,legs.shape[2]):
            if (k>nleg_order):
               #print "leg_geloescht"
               #print k
               legs[i][j][k]=0

   fakultaeten=np.arange(nleg+5,dtype=np.int64)
   fakultaeten[0]=1
   fakultaeten=factorial(fakultaeten)

   coefficient_t=np.zeros((4,nleg))

   for p in range(1,4):
      for l in range(0,nleg):
         coefficient_t[p][l]=(-1)**p*2*math.sqrt(2*l+1)*fakultaeten[l+p-1]/fakultaeten[p-1]/fakultaeten[l-p+1]
         if (l+p)%2==0:
            coefficient_t[p][l]=0
         if l>nleg_order:
            coefficient_t[p][l]=0

   #print "coefficient_t"
   #print coefficient_t

   mu_diagonal=np.zeros((muimp.shape[0]*2,muimp.shape[0]*2))

   for iB in range(muimp.shape[0]):
      for iS in range(0,2):
         mu_diagonal[2*iB+iS][2*iB+iS]=muimp[iB][iS]

   muimp_squared= np.tensordot(mu_diagonal,mu_diagonal,(1,0))
   siw_mom_inf_squared=np.tensordot(siw_mom_inf,siw_mom_inf,(1,0))
   siw_mom_inf_mu_contracted=np.tensordot(siw_mom_inf,mu_diagonal,(1,0))

   legs_copy=np.zeros_like(legs)
   legs_copy=legs

   ##############################
   # do the rescaling
   lam1=np.zeros((nbands,2))
   lam2=np.zeros((nbands,2))
   lam3=np.zeros((nbands,2))

   for i in range(0,nbands):
      for j in range(0,2):

         c1=-1
         c2=siw_mom_inf[2*i+j][2*i+j]+mu_diagonal[2*i+j][2*i+j]
         #c3=siw_mom_1[2*i+j][2*i+j]+siw_mom_inf_squared[2*i+j][2*i+j]-2*siw_mom_inf_mu_contracted[2*i+j][2*i+j]-muimp_squared[2*i+j][2*i+j]
         c3=siw_mom_1[2*i+j][2*i+j]-c2**2

         #print ""
         #print "check relations:"
         #print "c1"
         #print "c2"
         #print "c3"
         #print c1
         #print c2
         #print c3

         term1=beta**3*c3*np.dot(coefficient_t[1][:],coefficient_t[3][:])
         term2=np.dot(legs[i][j][:],coefficient_t[3][:])*np.dot(coefficient_t[1][:],coefficient_t[3][:])
         term3=beta*c1*np.dot(coefficient_t[3][:],coefficient_t[3][:])
         term4=np.dot(legs[i][j][:],coefficient_t[1][:])*np.dot(coefficient_t[3][:],coefficient_t[3][:])
         term5=np.dot(coefficient_t[1][:],coefficient_t[3][:])*np.dot(coefficient_t[1][:],coefficient_t[3][:])
         term6=np.dot(coefficient_t[1][:],coefficient_t[1][:])*np.dot(coefficient_t[3][:],coefficient_t[3][:])

         lam1[i][j] = -(-term1+term2+term3-term4)/(term5-term6)

         term1=beta**2*c2
         term2=np.dot(legs[i][j][:],coefficient_t[2][:])
         term3=np.dot(coefficient_t[2][:],coefficient_t[2][:])
            
         lam2[i][j]=-(-term1+term2)/term3

         term1=beta**3*c3*np.dot(coefficient_t[1][:],coefficient_t[1][:])
         term2=np.dot(legs[i][j][:],coefficient_t[3][:])*np.dot(coefficient_t[1][:],coefficient_t[1][:])
         term3=beta*c1*np.dot(coefficient_t[1][:],coefficient_t[3][:])
         term4=np.dot(legs[i][j][:],coefficient_t[1][:])*np.dot(coefficient_t[1][:],coefficient_t[3][:])
         term5=np.dot(coefficient_t[1][:],coefficient_t[3][:])*np.dot(coefficient_t[1][:],coefficient_t[3][:])
         term6=np.dot(coefficient_t[1][:],coefficient_t[1][:])*np.dot(coefficient_t[3][:],coefficient_t[3][:])

         lam3[i][j] = -(+term1-term2-term3+term4)/(term5-term6)

         legs[i][j][:]= legs_copy[i][j][:]+lam1[i][j]*coefficient_t[1][:]+lam2[i][j]*coefficient_t[2][:]+lam3[i][j]*coefficient_t[3][:]
               
         #print "np.dot(legs[i][j][:],coefficient_t[1][:])"
         #print "np.dot(legs[i][j][:],coefficient_t[2][:])"
         #print "np.dot(legs[i][j][:],coefficient_t[3][:])"
         #print np.dot(legs[i][j][:],coefficient_t[1][:])
         #print np.dot(legs[i][j][:],coefficient_t[2][:])
         #print np.dot(legs[i][j][:],coefficient_t[3][:])

         #print "check relation with formula from paper:"
         #cc1=0
         #for l in range(0,nleg):
            #ff=2*math.sqrt(2*l+1)/beta*legs[i][j][l]
            #if l%2==1:
               #ff=0
            #cc1=cc1-ff
         #print "cc1"
         #print cc1
         
         #cc2=0
         #for l in range(0,nleg):
            #ff=2*math.sqrt(2*l+1)/beta**2*legs[i][j][l]*l*(l+1)
            #if l%2==0:
               #ff=0
            #cc2=cc2+ff
         #print "cc2"
         #print cc2

         #cc3=0
         #for l in range(0,nleg):
            #ff=math.sqrt(2*l+1)/beta**3*legs[i][j][l]*(l+2)*(l+1)*l*(l-1)
            #if l%2==1:
               #ff=0
            #cc3=cc3-ff
         #print "cc3"
         #print cc3

   return legs




# test function to glue the high frequency part of the self energy; not in use
def klebe_siw(siw,niw,matsubara,siw_mom):

   def find_index_nearest(array,value):
      value_array=np.zeros_like(array)      
      value_array[:]=value
      return abs(array-value_array).argmin()


   print("siw_mom")
   print(siw_mom)

   print("siw[0][0][niw//2]") # this is the value of the self energy at the first matsubara frequency
   print(siw[0][0][niw//2])

   print("matsubara[niw//2]") # first positive matsubara frequency
   print(matsubara[niw//2])
   
   print("matsubara.shape")
   print(matsubara.shape)

   mat_50=find_index_nearest(matsubara,40)
   print("mat_50")
   print(mat_50)
   print(matsubara[mat_50])
   
   index_min=find_index_nearest(matsubara[niw//2:mat_50],5)
   index_max=find_index_nearest(matsubara[niw//2:mat_50],10)

   print("find_nearest(siw[0][0][:],5)")
   print(index_min)
   print(matsubara[niw//2+index_min])
   print(index_max)
   print(matsubara[niw//2+index_max])

   high_freq=np.zeros_like(siw)
   
   def weight_function(x,lower,upper):
      return (x-lower)/(upper-lower)

   #def weight_function(x,lower,upper):
      #x=(x-lower)/(lower-upper)/2-0.5
      #return 0.5*(1+math.cos(2*math.pi*x))

   for iB in range(0,siw_mom.shape[0]//2):
      for iS in range(0,siw_mom.shape[1]//2):

         lower=2*siw[iB][iS][niw//2:niw].argmin()
         upper=4*lower
         
         lower_mat=matsubara[niw//2+lower]
         upper_mat=matsubara[niw//2+upper]

         if siw[iB][iS][lower]<-10: 
            lower_mat=4
            lower_mat=15


         for i in matsubara:
            index=np.where(matsubara==i)
            if i<=-upper_mat:
               high_freq[iB][iS][index]=siw_mom[2*iB+iS][2*iB+iS]/i
            if i<-lower_mat and i>=-upper_mat:
               high_freq[iB][iS][index]=siw[iB][iS][index]*(1-weight_function(abs(i),lower_mat,upper_mat))+weight_function(abs(i),lower_mat,upper_mat)*siw_mom[2*iB+iS][2*iB+iS]/i
            if i>=-lower_mat and i<0:
               high_freq[iB][iS][index]=siw[iB][iS][index]
            if i<=lower_mat and i>0:
               high_freq[iB][iS][index]=siw[iB][iS][index]
            if i>lower_mat and i<=upper_mat:
               high_freq[iB][iS][index]=siw[iB][iS][index]*(1-weight_function(i,lower_mat,upper_mat))+weight_function(i,lower_mat,upper_mat)*siw_mom[2*iB+iS][2*iB+iS]/i
            if i>=upper_mat:
               high_freq[iB][iS][index]=siw_mom[2*iB+iS][2*iB+iS]/i

   return high_freq

def get_sztau_sz0_diag(ntau_n0, ntau_n0_err=None):
   "< S_z(tau) S_z(0) > diagonal in orbitals"
   nbands=ntau_n0.shape[0]
   ntau=ntau_n0.shape[-1]
   sztau_sz0=np.zeros((ntau))

   for b1 in range(0,nbands):

      # this is not the formula that Werner wrote in his spin-freezing PRL, but used to produce the plot
      sztau_sz0[:]=sztau_sz0[:]+(ntau_n0[b1,0,b1,0,:]+ntau_n0[b1,1,b1,1,:]-ntau_n0[b1,1,b1,0,:]-ntau_n0[b1,0,b1,1,:])

   if ntau_n0_err is not None:
       sztau_sz0_err = np.zeros_like(sztau_sz0)
       for b1 in range(0, nbands):
           sztau_sz0_err[:] = sztau_sz0_err[:] + (ntau_n0_err[b1,0,b1,0,:]**2
                                                  + ntau_n0_err[b1,1,b1,1,:]**2
                                                  + ntau_n0_err[b1,1,b1,0,:]**2
                                                  + ntau_n0_err[b1,0,b1,1,:]**2)
       sztau_sz0_err = np.sqrt(sztau_sz0_err)
       return sztau_sz0, sztau_sz0_err

   return sztau_sz0

def get_sztau_sz0(ntau_n0, ntau_n0_err=None):
   "< S_z(tau) S_z(0) > full"
   nbands=ntau_n0.shape[0]
   ntau=ntau_n0.shape[-1]
   sztau_sz0=np.zeros((ntau))

   for b1 in range(0,nbands):
      for b2 in range(0,nbands):

         sztau_sz0[:]=sztau_sz0[:]+(ntau_n0[b1,0,b2,0,:]+ntau_n0[b1,1,b2,1,:]-ntau_n0[b1,1,b2,0,:]-ntau_n0[b1,0,b2,1,:])

   if ntau_n0_err is not None:
       sztau_sz0_err = np.zeros_like(sztau_sz0)
       for b1 in range(0, nbands):
           for b2 in range(0, nbands):
               sztau_sz0_err[:] = sztau_sz0_err[:] + (ntau_n0_err[b1,0,b2,0,:]**2
                                                      + ntau_n0_err[b1,1,b2,1,:]**2
                                                      + ntau_n0_err[b1,1,b2,0,:]**2
                                                      + ntau_n0_err[b1,0,b2,1,:]**2)
       sztau_sz0_err = np.sqrt(sztau_sz0_err)
       return sztau_sz0, sztau_sz0_err

   return sztau_sz0

def get_Ntau_N0(ntau_n0, ntau_n0_err=None):
   "< N(tau) N(0) > full"
   if ntau_n0_err is not None:
       return np.sum(ntau_n0, axis=(0, 1, 2, 3)), np.sqrt(np.sum(ntau_n0_err**2, axis=(0, 1, 2, 3)))
   return np.sum(ntau_n0, axis=(0, 1, 2, 3))

def get_sztau_sz0_orb_resolved(ntau_n0, ntau_n0_err=None):
   "< S_z(tau) S_z(0) > orbital resolved"
   nbands=ntau_n0.shape[0]
   ntau=ntau_n0.shape[-1]
   sztau_sz0=np.zeros((nbands,ntau))

   for b1 in range(0,nbands):

      sztau_sz0[b1,:]=sztau_sz0[b1,:]+(ntau_n0[b1,0,b1,0,:]+ntau_n0[b1,1,b1,1,:]-ntau_n0[b1,1,b1,0,:]-ntau_n0[b1,0,b1,1,:])

   if ntau_n0_err is not None:
       sztau_sz0_err = np.zeros_like(sztau_sz0)
       for b1 in range(0, nbands):
           sztau_sz0_err[b1, :] = sztau_sz0_err[b1, :] + (ntau_n0_err[b1,0,b1,0,:]**2
                                                          + ntau_n0_err[b1,1,b1,1,:]**2
                                                          + ntau_n0_err[b1,1,b1,0,:]**2
                                                          + ntau_n0_err[b1,0,b1,1,:]**2)
       sztau_sz0_err = np.sqrt(sztau_sz0_err)
       return sztau_sz0, sztau_sz0_err

   return sztau_sz0

def get_sztau_sz0_offdiag(ntau_n0, ntau_n0_err=None):
   "< S_z(tau) S_z(0) > orbital offdiagonal components"
   nbands=ntau_n0.shape[0]
   ntau=ntau_n0.shape[-1]
   sztau_sz0=np.zeros((ntau))

   for b1 in range(0,nbands):
      for b2 in range(0,nbands):

         if b1!=b2:
            sztau_sz0[:]=sztau_sz0[:]+(ntau_n0[b1,0,b2,0,:]+ntau_n0[b1,1,b2,1,:]-ntau_n0[b1,1,b2,0,:]-ntau_n0[b1,0,b2,1,:])

   if ntau_n0_err is not None:
       sztau_sz0_err = np.zeros_like(sztau_sz0)
       for b1 in range(0, nbands):
           for b2 in range(0, nbands):
               if b1 != b2:
                   sztau_sz0_err[:] = sztau_sz0_err[:] + (ntau_n0_err[b1,0,b2,0,:]**2
                                                          + ntau_n0_err[b1,1,b2,1,:]**2
                                                          + ntau_n0_err[b1,1,b2,0,:]**2
                                                          + ntau_n0_err[b1,0,b2,1,:]**2)
       sztau_sz0_err = np.sqrt(sztau_sz0_err)
       return sztau_sz0, sztau_sz0_err

   return sztau_sz0

# Derived quantities for the iterations
derived_quantities = {
    "imp-density": meta_from_func(
                        impurity_density, 
                        ["ineq","band","spin"], ["value", "error"], ["occ"]
                        ),
    "imp-electrons": meta_from_func(
                        impurity_electrons,
                        ["ineq"], ["value", "error"]
                        ),
    "moment": meta_from_func(
                        moment,
                        ["ineq","band1","band2"], ["value","error"],
                        ["occ"]
                        ),
    "latt-electrons": meta_from_func(lattice_electrons,
                        [], ["spin-up", "spin-dn", "total"]
                        ),
    "chi-ph": meta_from_func(get_chi_ph,
                        ['ineq','band1','spin1','band2','spin2','iwf-g4','iwf-g4','iwb-g4'],
                        ['value'],
                        ['giw', 'g4iw']
                        ),
    "chi-pp": meta_from_func(get_chi_pp,
                        ['ineq','band1','spin1','band2','spin2','iwf-g4','iwf-g4','iwb-g4'],
                        ['value'],
                        ['giw', 'g4iw-pp'],
                        ),
    "g4iw-conn-ph": meta_from_func(get_g4iw_conn_ph,
                        ['ineq','band1','spin1','band2','spin2','iwf-g4','iwf-g4','iwb-g4'],
                        ['value'],
                        ['giw', 'g4iw']
                        ),
    "sztau-sz0": meta_from_func(get_sztau_sz0,
                        ['ineq', 'tausus'],
                        ['value', "error"],
                        ['ntau-n0']
                        ),
    "sztau-sz0-diag": meta_from_func(get_sztau_sz0_diag,
                        ['ineq', 'tausus'],
                        ['value', "error"],
                        ['ntau-n0']
                        ),
    "sztau-sz0-orb-resolved": meta_from_func(get_sztau_sz0_orb_resolved,
                        ['ineq', 'band', 'tausus'],
                        ['value', "error"],
                        ['ntau-n0']
                        ),
    "sztau-sz0-offdiag": meta_from_func(get_sztau_sz0_offdiag,
                        ['ineq', 'tausus',],
                        ['value', "error"],
                        ['ntau-n0']
                        ),
    "denstau-dens0": meta_from_func(get_Ntau_N0,
                        ['ineq', 'tausus'],
                        ['value', "error"],
                        ['ntau-n0']
                        ),
    "gtau-tauint": meta_from_func(
                        gtau_tauint,
                        ["ineq","band","spin","taubin"],
                        ["value"],
                        ["gtau-blocks"]
                        ),
    "gtau-spec": meta_from_func(
                        gtau_spec,
                        ["ineq","band","spin","taubin","m"],
                        ["value"],
                        ["gtau-blocks"]
                        ),
    "gtau-lbap": meta_from_func(
                        gtau_lbap,
                        ["ineq","band","spin","taubin","m"],
                        ["value"],
                        ["gtau-blocks"]
                        ),
    }
