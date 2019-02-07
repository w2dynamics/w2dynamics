r""" @package transform
Routines for transforming CT-QMC results between Legendre, Matsubara and
imaginary time representation.

Firstly, provides methods named like @c source2target returning transformations
between finite expansions.  @c source and @c target can be one of:
   - tau -- Imaginary time mesh @f$ \tau_i @f$
   - leg -- Legendre expansion @f$ P_l @f$
   - mat -- Matsubara expansion, either bosonic or fermionic @f$ i\omega_n @f$

The returned matrices can be tensor producted with the appropriate quantity to
yield the result.
@see leg2tau(), tau2leg(), leg2mat(), mat2tau(), tau2mat(), np.tensordot()

Secondly, provides methods working with the actual physical quantities as
retrieved by the CT-QMC algorithm.
@see matfreq(), transform(), g4leg2tau(), g4leg2mat(), ggleg()
"""
import warnings
import numpy as np
import scipy.special as sps
import itertools
from scipy import integrate, interpolate

_eps = 1e-8

def _trapez_weights(taus):
    """ Gets trapezoidal integration weights for some mesh """
    return (np.append(taus,2*[taus[-1]]) - np.append(2*[taus[0]],taus))[1:-1]/2

def _ensurebeta(beta):
    """ Helper for sanity check on beta """
    beta = float(beta)
    assert beta > _eps, "Beta must be positive"
    return beta

def _ensurearray(array, dim_pred):
    """ Helper for ensuring numpy array of the correct shape """
    array = np.asarray(array)
    assert dim_pred(array.ndim), "Invalid number of array dimensions"
    return array

def _tauarray(beta, taus):
    """ Return tau arrays also for number of tau points """
    try:
        ntaus = int(taus)
        assert taus == ntaus and ntaus > 0, "Illegal number of tau points"
        return np.arange(ntaus)*beta/(ntaus-1.0)
    except TypeError:
        taus = _ensurearray(taus, lambda dim: dim == 1)
        assert np.all(taus >= -_eps) and np.all(taus <= beta + _eps), \
               "Taus must be reduced to the interval [0, beta]"
        return taus

def _ensureweights(beta, taus, weights):
    """ Helper for integration functions """
    if weights is None:
        return _trapez_weights(taus)
    else:
        warr = np.asarray(weights)
        assert warr.shape == taus.shape, "Illegal weights (must match tau mesh)"
        assert abs(np.sum(warr) - beta) < _eps(), "Weights should add up to beta"
        return warr

def _ensureorder(order):
    """ Helper for ensuring that order is an integer """
    iorder = int(order)
    assert iorder == order and iorder > 0, "Order must be positive integer"
    return iorder

def tau_bins(beta, nbins, value='lower-bound'):
    r""" Get bin positions for a tau mesh

    @param beta    Thermodynamic beta
    @param nbins   Number of bins (equals the number of points returned)
    @param value   Position to return inside the bin, may be one of
                   'lower-bound', 'centre' or 'upper-bound'
    @return One-dimensional real numpy.array with \c nbins bin values
    """
    beta = _ensurebeta(beta)
    nbins = _ensureorder(nbins)
    loc = {'lower-bound': 0.,
           'centre':      0.5,
           'upper-bound': 1.}[value] * beta / nbins
    return np.linspace(loc, beta+loc, nbins, endpoint=False)

def matfreq(beta, mattype, matorder):
    r""" Get fermionic or bosonic Matsubara frequencies.

    @param beta      Thermodynamic beta
    @param mattype   Either 'bose' or 'fermi'
    @param matorder  (Even) number of fermionic or (odd) number of bosonic
                     Matsubara frequencies to be generated
    @return
        One-dimensional REAL numpy.array with matsubara frequencies
        @f$ w_i = \frac{2\pi}\beta (i + k) @f$, where @f$ k @f$ is zero
        for bosonic and 1/2 for fermionic frequencies. @f$ i @f$ runs
        symmetrically around 0.
    """
    beta = _ensurebeta(beta)
    matorder = _ensureorder(matorder)
    #
    numwhalf = int(matorder/2)   # number of matsubara functions
    if mattype == 'bose':
        assert matorder % 2 == 1, "Number of Bose frequencies must be odd"
        return (2*np.pi/beta)*np.arange(-numwhalf,numwhalf+1)
    elif mattype == 'fermi':
        assert matorder % 2 == 0, "Number of Fermi frequencies must be even"
        return (2*np.pi/beta)*(np.arange(-numwhalf,numwhalf) + 0.5)
    else:
        raise ValueError("mattype must be either 'fermi' or 'bose'")

def hybr_from_sites(epsk, vki, ffreq, tau, beta=None):
    """Constructs the Hybridisation function from bath sites

    Given a set of bath sites with energy levels `epsk` and Hybridisation
    strength to the impurity `vki`, computes the Hybridisation function
    on the Matsubara axis and on the imaginary time axis

    ..math ::  \Delta_{ij}(i\omega) =
                    \sum_k \frac {V^*_{ki} V_{kj}} {i\omega - \epsilon_k} \\
               \Delta_{ij}(\tau) = \sum_k V^*_{ki} V_{kj}
                    \frac {\exp \tau\epsilon_k} {1 + \exp \beta\epsilon_k}

    Parameters:
      - epsk     set of bath energies
      - vki      two-by-two Hybridisation matrix V_{k,i}, where the rows
                 correspond to bath levels and the columns to impurity sites
      - ffreq    fermionic Matsubara frequencies
      - tau      imagninary time grid
      - beta     inverse temperature
    """
    epsk = np.atleast_1d(epsk)
    vki = np.atleast_2d(vki)
    iwn = 1j*np.asarray(ffreq)
    tau = np.asarray(tau)
    if beta is None: beta = tau[-1]

    vijk = vki.T.conj()[:,None,:] * vki.T[None,:,:]
    hybriv = (vijk[:,:,None,:] / (iwn[:,None] - epsk)).sum(-1)
    hybrtau = ((vijk/(1 + np.exp(beta * epsk)))[:,:,None,:] *
               np.exp(tau[:,None] * epsk)).sum(-1)

    return hybriv, hybrtau

def leg2tau(beta, legorder, taus):
    r""" Get transformation matrix from a Legendre expansion to a tau mesh.

    Computes the transformation matrix @f$ L @f$ from a finite set of Legendre
    coefficients @f$ G_l @f$ to a imaginary time mesh @f$ (\tau_i)_i @f$:
    @f[
            G(\tau_i) = \sum_l L_{il} G(l)
    @f]
    given simply by the evaluation of the Legendre polynomials at the different
    tau points (see Boehnke et al., 2011, eq. 1):
    @f[
            L_{il} = \frac \sqrt{2l+1} \beta  P_l(x(\tau_i))
    @f]

    This matrix is isometric (L*L = 1) in the limit of a dense tau mesh,
    co-isometric (LL* = 1) in the limit of infinite Legendre order and unitary
    (LL* = L*L = 1) in the combined limit.

    @see tau2leg()
    @param beta      Thermodynamic beta
    @param legorder  Number of Legendre coefficients (highest order minus 1)
    @param taus      Either array (tau mesh) or number of tau points as integer
    """
    # Get array
    beta = _ensurebeta(beta)
    taus = _tauarray(beta, taus)
    red_taus = 2*taus/beta - 1
    legcoeffs = np.empty(shape=(taus.size,legorder))
    for l in range(0,legorder):
        legcoeffs[:,l] = np.sqrt(2*l+1)/beta * sps.eval_legendre(l,red_taus)

    return legcoeffs

def tau2leg(beta, taus, legorder, weights=None):
    r""" Get approximate transformation matrix from a tau mesh to a Legendre expansion.

    Computes an approximation to the transformation matrix @f$ L^{-1} @f$ from a
    (possibly equidistant) tau mesh @f$ (tau_i)_i @f$ to a finite set of Legendre
    coefficients @f$ G_l @f$ (see Boehnke et al. 2011, eq. 2):
    @f[
          G(l) =       \sqrt{2l+1} \int_0^\beta d\tau P_l(x(\tau)) G(\tau)
               \approx \sum_i L^{-1}_{li} G(\tau_i)
    @f]
    The integral is, therefore, approximated by a weighted sum over the sampling
    values, where the weights @f$ w_i @f$ depend on the integration method and
    the tau mesh. The explicit form of the matrix element then reads:
    @f[
              L^{-1}_{li} = \sqrt{2l+1} w_i P_l(x(\tau_i)) G(\tau_i)
    @f]

    This matrix is isometric (L*L = 1) in the limit of infinite Legendre order,
    co-isometric (LL* = 1) in the limit of a dense tau mesh and unitary
    (LL* = L*L = 1) in the combined limit.

    @see leg2tau()
    @param beta      Thermodynamic beta
    @param taus      Either array (tau mesh) or number of tau points as integer
    @param legorder  Number of Legendre coefficients (highest order minus 1)
    @param weights   Integration weights (by default given by Trapezoid rule)
    """
    legorder = _ensureorder(legorder)
    beta = _ensurebeta(beta)
    taus = _tauarray(beta, taus)
    weights = _ensureweights(beta, taus, weights)

    # Reduces taus
    red_taus = 2*taus/beta - 1
    coeffs = np.empty(shape=(legorder,taus.size))
    for l in range(0,legorder):
        coeffs[l,:] = np.sqrt(2*l+1) * ws * sps.eval_legendre(l,red_taus)

    return coeffs

def leg2mat(legorder, matorder, other=False):
    r""" Get transformation matrix from a Legendre expansion to a Matsubara expansion

    Compute the transformation matrix @f$ T @f$ from a finite set of Legendre
    coefficients @f$ G_l @f$ to a truncated expansion in terms of fermionic
    Matsubara frequencies @f$ \omega_n = (2n + 1)\frac \pi\beta @f$:
    @f[
            G(i\omega_n) = \sum_l T_{nl} G(l)
    @f]
    given in terms of the spherical Bessel functions @f$ j_l @f$ (see Boehnke et al.,
    2011, eq. E2):
    @f[
            T_{nl} = (-1)^n i^{l+1} \sqrt{2l+1} j_l(\omega_n \beta/2)
    @f]

    This matrix is independent of beta. It is isometric (T*T = 1) in the limit of
    infinite Matsubara order, co-isometric (TT* = 1) in the limit of infinite
    Legendre order and unitary (TT* = T*T = 1) in the combined limit.

    @see matfreq()
    @param legorder  Number of Legendre coefficients (highest order minus 1)
    @param matorder  Even number of fermionic Matsubara frequencies (symm around 0)
    """
    legorder = _ensureorder(legorder)

    # Sum over Legendre polynomials ("basis change" l -> iw_n), Boehnke 2011(E2)
    numwhalf = int(matorder/2)   # number of matsubara functions
    wshalf = matfreq(2, 'fermi', matorder) + 0j # returns 2pi/2*(n+1/2)
    legendres = np.arange(legorder)

    # Compute prefactor
    nfac = (-1.0)**np.arange(-numwhalf,numwhalf) * np.sqrt(np.pi/(2*wshalf))
    lfac = (1j)**(legendres+1) * np.sqrt(2*legendres+1)
    return nfac[np.newaxis].T * lfac * sps.jn(legendres+0.5, wshalf[np.newaxis].T)

def mat2tau(beta, mattype, matorder, taus):
    r""" Get transformation from truncated Matsubara expansion to tau mesh.

    Compute the truncated fourier series @f$ F^{-1} @f$ of a finite set
    of Matsubara frequencies @f$ \omega_n = (2n + k)\frac\pi\beta @f$, either
    fermionic (@f$ k=1 @f$) or bosonic (@f$ k=0 @f$) on a imaginary time mesh:
    @f[
            G(\tau_i) = \sum_n F^{-1}_{in} G(i\omega_n)
    @f]
    given in terms of phase factors:
    @f[
            F^{-1}_{in} = \frac 1\beta  \exp(-\imath\omega_n \tau)
    @f]

    @see tau2mat(), matfreq()
    @param beta      Thermodynamic beta
    @param mattype   Either 'bose', 'fermi' or 'all'
    @param matorder  Even number of fermionic Matsubara frequencies (symm around 0)
    @param taus      Either array (tau mesh) or number of tau points as integer
    """
    beta = _ensurebeta(beta)
    taus = _tauarray(beta, taus)
    miws = -1j*matfreq(beta, mattype, matorder)

    return np.exp(np.outer(taus, miws))/beta

def tau2mat(beta, taus, mattype, matorder, weights=None):
    r""" Get transformation from a tau mesh to a truncated Matsubara expansion.

    Compute an approximation to the fourier coefficients @f$ F @f$ of a
    finite set of Matsubara frequencies @f$ \omega_n = (2n + k)\frac\pi\beta @f$,
    either fermionic (@f$ k=1 @f$) or bosonic (@f$ k=0 @f$) based on a
    finite imaginary time mesh:
    @f[
          G(i\omega_n) =       \int_0^\beta d\tau \exp(-i\omega_n\tau) G(\tau)
                       \approx \sum_i F_{ni} G(\tau_i)
    @f]
    The integral is, therefore, approximated by a weighted sum over the sampling
    values, where the weights @f$ w_i @f$ depend on the integration method and
    the tau mesh. The explicit form of the matrix element then reads:
    @f[
          F_{ni} = \exp(+\imath\omega_n \tau) w_i
    @f]

    Please mind the poor stability of this function for a coarse tau mesh.

    @see mat2tau(), matfreq()
    @param beta      Thermodynamic beta
    @param taus      Either array (tau mesh) or number of tau points as integer
    @param mattype   Either 'bose', 'fermi' or 'all'
    @param matorder  Even number of fermionic Matsubara frequencies (symm around 0)
    @param weights   Integration weights (by default given by Trapezoid rule)
    """
    beta = _ensurebeta(beta)
    taus = _tauarray(beta, taus)
    weights = _ensureweights(beta, taus, weights)

    iws = 1j*matfreq(beta, mattype, matorder)
    # Since weights is a row vector, this works
    return np.exp(np.outer(iws, taus)) * weights


#GS:
def wre2mat(beta, wre, mattype, matorder, weights=None):
    r""" Get transformation from a real frequency mesh to Matsubara.

    @see mat2tau(), matfreq()
    @param beta      Thermodynamic beta
    @param taus      Either array (tau mesh) or number of tau points as integer
    @param mattype   Either 'bose', 'fermi' or 'all'
    @param matorder  Even number of fermionic Matsubara frequencies (symm around 0)
    @param weights   Integration weights (by default given by Trapezoid rule)
    """
    beta = _ensurebeta(beta)
#GS: do we need these weights?
#   weights = _ensureweights(beta, taus, weights)
#GS: attempt to create wre
    wre = _ensurearray(wre, lambda dim: dim == 1)
    dw = (wre[-1]-wre[0])/float(wre.shape[0])

    iws = 1j*matfreq(beta, mattype, matorder)
#GS: Broadcast works becase I am summing a (Niws,1) vector with a (Nwre,) one.
#GS: The result is a (Niws, Nwre) matrix. In function "transform" this will multiply
#GS: a vector of length Nwre (leadsw), this way performing the integral
    return dw/(iws.reshape(-1,1) - wre)


#GS interpolation + quadrature
def wre2mat_quad(beta, wre, leadsw, mattype, matorder):
    r""" Performs the following integral:
 
    int dwre leadsw(wre) / ( iws - wre ) = int dwre leadsw(wre)(-iws-wre)/(wre^2+ws^2)

    if leadsw(wre) is already -(1/pi)*Delta(w) the result is equal to Delta(iws)
    
    @see mat2tau(), matfreq()
    @param beta      Thermodynamic beta 
    @param taus      Either array (tau mesh) or number of tau points as integer  
    @param mattype   Either 'bose', 'fermi' or 'all'
    @param matorder  Even number of fermionic Matsubara frequencies (symm around 0)
    """

    print "leadsw.shape", leadsw.shape
    leadsiw=np.zeros(leadsw.shape[:-1]+(matorder,), dtype=complex)
    print "leadsiw.shape", leadsiw.shape
    print "leadsiw.dtype", leadsiw.dtype

    beta = _ensurebeta(beta)
#GS: create wre
    wre = _ensurearray(wre, lambda dim: dim == 1)

    iws = 1j*matfreq(beta, mattype, matorder)

#GS: mega loop over all components of leadsiw
    for ilead, iorb, isp, jorb, jsp in np.ndindex(leadsiw.shape[:-1]):
  
#GS: interpolate leadsw: is interp1d the best choice?
        finter = interpolate.interp1d(wre, leadsw[ilead, iorb, isp, jorb, jsp,:], kind='cubic')

        for iiw, iw in enumerate(iws):
            re, _ = integrate.quad(lambda x: finter(x)*x/(x**2+iw.imag**2), wre[0], wre[-1], epsabs=1.5e-06, epsrel=1.5e-06)
            im, _ = integrate.quad(lambda x: finter(x)*iw.imag/(x**2+iw.imag**2), wre[0], wre[-1], epsabs=1.5e-06, epsrel=1.5e-06)
            leadsiw[ilead, iorb, isp, jorb, jsp, iiw] = -re - 1j*im

    return leadsiw


#GS:
def wre_int(beta, wre, weights=None):
    r""" Get the integral on real frequencies

    @see mat2tau(), matfreq()
    @param beta      Thermodynamic beta
    @param taus      Either array (tau mesh) or number of tau points as integer
    @param mattype   Either 'bose', 'fermi' or 'all'
    @param matorder  Even number of fermionic Matsubara frequencies (symm around 0)
    @param weights   Integration weights (by default given by Trapezoid rule)
    """
    beta = _ensurebeta(beta)
#GS: do we need these weights?
#   weights = _ensureweights(beta, taus, weights)
#GS: attempt to create wre
    wre = _ensurearray(wre, lambda dim: dim == 1)
    dw = (wre[-1]-wre[0])/float(wre.shape[0])

#GS:iws = 1j*matfreq(beta, mattype, matorder)
#GS: I am creating a (1,Nwre) array. In function "transform" this will multiply
#GS: a vector of length Nwre (leadsw), this way performing the integral
    return dw*np.ones(wre.shape[0]).reshape(1,-1)


def errorprop(tfmat, full=False):
    r""" Get error propagation tensor for a linear transformation.

    Let @f$ A @f$ be a linear transformation (@f$ Y = A X @f$), compute the
    error propagation operator @f$ T(e) @f$ through @f$ A @f$ for a given set of
    uncorrelated statistical variables @f$ X @f$. This operator maps a set of
    variances @f$ e_i = Var[X_i] @f$ to the corresponding covariance matrix @f$
    \Sigma_Y @f$ in @f$ Y @f$ by propagation of uncertainty:
    @f[
         \Sigma_Y^{ij} = (A \Sigma_X A^+)^{ij} = \sum_k A^{ik} (A^*)^{jk} e_k
    @f]

    The covariance matrix is a positive-semidefinite operator, which generalises
    the variance to multiple variables @f (X_i) @f:
    @f[
                 \Sigma^{ij} =  E[ (X_i - E[X_i])(X_j - E[X_j])* ],
    @f]
    where @f$ E @f$ denotes the expectation value. The error (uncertainties) are
    given by the diagonal elements. Off-diagonal terms yield the covariance and
    vanish for uncorrelated variables.

    @see transform(), leg2mat(), leg2tau(), mat2tau()
    @return
        A numpy array of the shape @f$ \dim Y, \dim X @f$ (same as @c tfmat),
        relating the error in @f$ X @f$ to the variance in @f$ Y @f$. If @c full
        is @c True, instead a numpy array of the shape @f$ \dim Y,\dim Y,\dim X
        @f$ is returned, which gives the full covariance matrix @f$ \Sigma_Y @f$
        for any error.
    """
    tfmat = _ensurearray(tfmat, lambda dims: dims == 2)
    if not full:
        # just return A_ik A^+_ki
        return np.square(np.abs(tfmat))
    else:
        # construct the full projectors
        prop = np.empty(shape=tfmat.shape[0:1]+tfmat.shape, dtype=tfmat.dtype)
        for k in range(tfmat.shape[1]):
            prop[:,:,k] = np.outer(tfmat[:,k], tfmat[:,k].conj())
        return prop

def _moveslices(myslices, dim):
    """ Helper for slicing out a part of an array at a specific dimension """
    try:
        if dim >= 0:
            return (slice(None),)*dim + myslices + (Ellipsis,)
        else:
            return (Ellipsis,) + myslices + (slice(None),)*(-dim-1)
    except:
        raise TypeError("Dimension must be integer")

class Truncate:
    """ Defines a set of actions on axis truncation for use in transform() """
    @staticmethod
    def tofirst(dim,order,warn=False):
        """ Return slice for the first elements (Legendre expansion) """
        if warn:
            warnings.warn("Truncating %d to first %d elements" % (dim,order))
        return slice(order)

    @staticmethod
    def tocentre(dim,order,warn=False):
        """ Return slice for the centre elements (Matsubara expansion) """
        if warn:
            warnings.warn("Truncating %d to centre %d elements" % (dim,order))
        assert dim % 2 == order % 2, "Centre truncation must be symmetric"
        cut = (dim-order)/2
        return slice(cut,-cut)

    @staticmethod
    def error(dim,order,warn=False):
        """ Do not truncate, but raise a ValueError """
        raise ValueError("Error: dimension mismatch: %d != %d" % (dim,order))

def transform(tfmat, data, error=None, dim=-1, fullerror=False,
              onmismatch=Truncate.error, warn=True):
    """ Apply linear transformation to parts of data and error matrix.

    @see leg2mat(), leg2tau(), mat2tau(), errorprop(), value_error(), #Truncate
    @param tfmat       Transformation matrix (dimension min 1, exactly if error)
    @param data        Data tensor to apply the transformation on
    @param error       Optional error estimates for the data tensor
    @param dim         Dimension to transform
    @param fullerror   If True, calculate full covariance matrix
    @param onmismatch  Action on shape mismatch of transformation and data matrix
    @param warn        Issue warning on truncation
    @return
        A numpy array of the same shape as @c tfmat, except that the dimension
        @c dim has been removed. Instead the transformed dimension(s) are
        *appended* as new dimensions (tensordot behaviour). If @error is given,
        instead a tuple @c (tdata,terror) of numpy arrays is returned.
    """
    tfmat = _ensurearray(tfmat, lambda dims: dims >= 1)
    data = _ensurearray(data, lambda dims: dims > (dim,-dim-1)[dim>=0])
    #
    tforder = tfmat.shape[-1]
    dimextent = data.shape[dim]
    tfpart = Ellipsis
    dpart = Ellipsis
    if tforder > dimextent: # Transformation is "too big"
        tfpart = (Ellipsis, onmismatch(tforder,dimextent,warn))
    elif tforder < dimextent:  # Data matrix is "too big"
        dpart = _moveslices((onmismatch(dimextent,tforder,warn),), dim)
    #
    tdata = np.tensordot(data[dpart], tfmat[tfpart], (dim,-1))
    if dim != -1:
        tdata = np.rollaxis(tdata, -1, dim)
    if error is None:
        return tdata
    else:
        error = np.asarray(error)
        assert error.shape == data.shape, "Data and error shape must agree"
        terror = np.tensordot(error[dpart], errorprop(tfmat[tfpart], fullerror),
                              (dim,-1))
        if dim != -1:
            terror = np.rollaxis(terror, -1, dim)
        return tdata, terror

def _g4chk(g, legorder):
    """ Helper routine ensuring correct format of the four-point Green's function """
    assert len(g.shape) >= 3, "Expecting array at least in l,lp,iw"
    assert g.shape[-3] >= legorder, "Demanded Legendre order exceeds provided l"
    assert g.shape[-2] >= legorder, "Demanded Legendre order exceeds provided lp"

def occexpand(occ):
    occ = _ensurearray(occ, lambda dims: dims == 3)
    assert occ.shape[1] % 2 == 0, "Shape1 must be even"
    assert occ.shape[1] == occ.shape[2], "Indices must agree"
    nbands = occ.shape[1]/2
    #
    # split the combined index by reshaping the array when viewed in column-major
    # then permute the new columns to yield (band1,spin1,band2,spin2) sequence
    res = np.reshape(occ, (occ.shape[0], nbands, 2, nbands, 2), order='C')
    return res

def g4expand(g):
    """ Expands the compressed 4-point Green's function returned by the code due
        to the Fortran 90 dimension limit of seven.

    @see g4compress()
    @param g
        Four point Green's function as 6-dimensional complex numpy.array, where
        the final dimensions correspond to:
           - 0-th dim -- Non-equivalent atoms in the cell
           - 1-st dim -- Combined band index (band1 * Nbands + band2)
           - 2-nd dim -- Combined spin index (spin1 * 2 + spin2)
           - 3-rd and 4-th dim -- Legendre indices l and lp
           - 5-th dim -- Bosonic Matsubara frequency index n
    @return
        Same Four point Green's function as 8-dimensional complex numpy.array,
        where the 1-st and 2-nd dimension are expanded into two dimensions each:
           - 0-th dim -- Non-equivalent atoms in the cell
           - 1-st and 2-nd dim -- Band and spin index of first particle
           - 3-rd and 4-th dim -- Band and spin index of second particle
           - 5-rd and 6-th dim -- Legendre indices l and lp
           - 7-th dim -- Bosonic Matsubara frequency index n
    """
    g = _ensurearray(g, lambda dims: dims == 6)
    nbands = int(np.sqrt(g.shape[1]))
    #
    # split the combined index by reshaping the array when viewed in column-major
    # then permute the new columns to yield (band1,spin1,band2,spin2) sequence
    res = np.reshape(g, (g.shape[0], nbands, nbands, 2, 2) + g.shape[3:], order='C')
    res = np.transpose(res, (0,1,3,2,4,5,6,7))
    return res

def g4leg2tau(g, beta, legorder, taus):
    """ Transforms 4-point Green's function from Legendre basis to a tau mesh.

    Converts the four-point Green's function array in the mixed Legendre/Matsubara
    basis to the three corresponding tau differences (tau12, tau34 and tau14 as
    outlined in Boehnke 2011).

    @param g
        Four point Green's function as at least 3-dimensional complex numpy.array,
        where the last three dimensions correspond to:
           1. Legendre index l
           2. Legendre index lp
           3. Bosonic matsubara frequency index n
    @param taus
        One dimensional numpy array containing the desired tau mesh. The values
        must be normalized to the range 0 to beta.
    @param legorder
        Desired order+1 of the Legendre expansion, modelling a band-pass filter.
        The third and fourth dimension of g are truncated to this order.
    @return
        Green's function ascomplex numpy array of the same shape as g, with all but
        the last three dimensions unchanged. The latter correspond to the tau
        indices of the discretisation of tau12, tau34, and tau14, respectively.
    """
    _g4chk(g, legorder)
    nmats = g.shape[-1]
    beta = _ensurebeta(beta)
    taus = _tauarray(beta, taus)
    expiwt = mat2tau(beta, 'bose', nmats, taus)
    lcoeffs = leg2tau(beta, legorder, taus)
    # arange yields a row vector, which is then multiplied column by column with
    # every row in lcoeffs, so OK.
    lpcoeffs = lcoeffs * (-1)**np.arange(1,legorder+1)
    #
    # Perform the transformation on the parts of g by tensor multiplying the respective
    # indices. Note that tensordot orders its result first by all the remaining indices
    # of the first argument, therefore we operate on dimension actdim all the time.
    actdim = len(g.shape) - 3
    res = g[..., 0:legorder, 0:legorder, :]
    res = np.tensordot(res, lcoeffs, ([actdim],[1]))
    res = np.tensordot(res, lpcoeffs, ([actdim],[1]))
    res = np.tensordot(res, expiwt, ([actdim],[1]))
    assert res.shape == (g.shape[0:actdim] + 3*(taus.size,)), "Internal shape error"
    return res

def g4leg2mat(g, legorder, matorder):
    """ Transforms 4-point Green's function from Legendre basis to Matsubara basis.

    Converts the four-point Green's function array in the mixed Legendre/Matsubara
    basis to the Fourier transform in two Fermionic and one bosonic Matsubara
    frequency

    @param g
        Four point Green's function as at least 3-dimensional complex numpy.array,
        where the last three dimensions correspond to:
           1. Legendre index l
           2. Legendre index lp
           3. Bosonic matsubara frequency index n
    @param matorder
        Desired number of the Fermionic Matsubara frequencies (size of the 3-rd and
        4-th dimension of the output array). The Matsubara frequencies are then
        in the range -n/2 to n/2-1.
    @param legorder
        Desired order+1 of the Legendre expansion, modelling a band-pass filter.
        The third and fourth dimension of g are truncated to this order.
    @return
        Green's function ascomplex numpy array of the same shape as g, with all but
        the Legendre dimensions unchanged. The latter correspond to the indices of
        fermionic Matsubara frequencies of the FT Legendre expansion.
    """
    _g4chk(g, legorder)
    #
    tnl = leg2mat(legorder, matorder)
    # arange yields a row vector, which is then multiplied column by column with
    # every row in lcoeffs, so OK.
    tnlp = tnl * (-1)**np.arange(1,legorder+1)
    #
    # Perform the transformation on the parts of g by tensor multiplying the respective
    # indices. Note that tensordot orders its result first by all the remaining indices
    # of the first argument, therefore we operate on dimension actdim all the time and use
    # transpose at the end to shift the matsubara axis (actdim-th dim.) again to the back.
    actdim = len(g.shape) - 3
    res = g[..., 0:legorder, 0:legorder, :]
    res = np.tensordot(res, tnl, axes=([actdim],[1]))
    res = np.tensordot(res, tnlp, axes=([actdim],[1]))
    res = np.transpose(res, range(actdim) + [actdim+1, actdim+2, actdim])
    assert res.shape == g.shape[0:actdim] + 2*(matorder,) + (g.shape[-1],), "Shape error"
    return res

def ggleg(g2leg, beta, nmats, g2err=None):
    """ Build the disconnected (GG) part of the two-particle Green's function.

    Build the disconnected (GG) part of the two-particle Green's function by
    products of one-particle Green's functions.

    @param g2leg
        Single particle (two point) Green's function as at least 3-dimensional
        complex numpy array, where the last three dimensions correspond to
           1. Band index
           2. Spin index (length 2)
           3. Legendre index l
    @param g2err
        Error array of the same shape as g2leg
    @return
        Disconnected part of the four-point Green's function as complex numpy array,
        with all but the last three dimensions of g2leg unchanged. The other seven
        dimension correspond to:
           - 1-st and 2-nd dim -- Band and spin index of first particle
           - 3-rd and 4-th dim -- Band and spin index of second particle
           - 5-rd and 6-th dim -- Legendre indices l and lp
           - 7-th dim -- Bosonic Matsubara frequency index n

        If g2err was supplied, it instead returns a pair of numpy arrays of the
        above shape, where the first element corresponds to the mean and the second
        element corresponds to the error estimate.
    """
    beta = _ensurebeta(beta)
    assert len(g2leg.shape) == 4, "Expected at least array in (band, spin, l)"
    assert g2leg.shape[-2] == 2, "Dimension of spin index must be 2"
    assert nmats%2 == 0, "Number of matsubara frequencies must be even"
    assert g2err is None or g2err.shape == g2leg.shape, "Shapes disagree"
    #
    # Disconnected part is given by:
    #          N_{AABB}(l,l',n) = \beta (-1)^{l'+1} G_A(l) G_B(l') \delta_{n0},
    nbands = g2leg.shape[1]
    mzero = nmats/2                    # zeroth matsubara function
    lpfact = float(beta) * (-1)**np.arange(1,g2leg.shape[3]+1)    # beta * (-1)^(l'+1)
    res = np.zeros(shape=(g2leg.shape[0],)+g2leg.shape[1:3]*2+(g2leg.shape[3],)*2+(nmats,),
                   dtype=g2leg.dtype)
    if g2err is not None:
        err = np.zeros(shape=res.shape, dtype=g2err.dtype)
    #
    # Extract the individual band/spin indices out of the compound ones of the 4-point
    # Green's function and build the direct product of the two GF
    for neq in range(g2leg.shape[0]):
        for b1, b2 in itertools.product(range(nbands), range(nbands) ):
            for sp1, sp2 in itertools.product(range(2), range(2)):
                res[neq,b1,sp1,b2,sp2,:,:,mzero] = \
                        np.outer(g2leg[neq,b1,sp1,:], lpfact*g2leg[neq,b2,sp2,:])

                if g2err is not None:
                    # Calculate error in GG by using straight-forward error propagation
                    # in the above calculation:
                    #  \sigma_{GG} = \beta \sqrt{ G_A^2(l) \sigma_B^2(l') +
                    #                             \sigma_A^2(l) G_B^2(l') }
                    err[neq,b1,sp1,b2,sp2,:,:,mzero] = float(beta) * np.sqrt( \
                        np.outer(np.abs(g2leg[neq,b1,sp1,:]), g2err[neq,b2,sp2,:])**2 +\
                        np.outer(g2err[neq,b1,sp1,:], np.abs(g2leg[neq,b2,sp2,:]))**2 )
    #
    if g2err is not None:
        return res,err
    else:
        return res
