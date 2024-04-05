#!/usr/bin/env python3
"""Allows for manipulation of the Green's function in real frequency space"""
from __future__ import print_function
import numpy as np
import scipy.integrate
import scipy.optimize
import itertools as itt
from warnings import warn

_CONVERR_MESSAGE = """
sigmaw.py: error: failed to converge Sigma(w = %.9g), aborting.

You can try using a higher deltino with the -D option (currently -D %.2g).
Also, you can restart from a different frequency using --wmin and/or approach
from the other side using --wreversed.
"""

class ConvergenceError(RuntimeError):
    pass

def kkt(aws, ws, integrator=scipy.integrate.simps):
    """ Does a Kramers-Kronig transform for the spectral function:
                Re G(w) = p.v. integral  A(w')/(w - w') dw'
    """
    realgw = np.empty_like(aws)
    for iw, w in enumerate(ws):
        denom = w - ws
        denom[iw] = 1  # to avoid division by zero
        integrand = (aws - aws[:,iw:iw+1])/denom
        realgw[:,iw] = integrator(integrand, ws)
    return realgw

def aws_to_gws(aws, ws):
    """Takes A(w) and computes the (complex) G(w) using a KKT."""
    gws = np.empty(shape=aws.shape, dtype=np.cdouble)
    gws.real = kkt(aws, ws)
    gws.imag = -np.pi*aws
    return gws

class spopt_cache:
    """ Decorator that caches the result for the last argument.
    
    This class works around two BRAINDEAD flaws in scipy.optimize:
      - the triple invocation of func for the same x in some solvers
      - the split of f and its Jacobian for scipy < 0.11
    """
    def __init__(self, func):
        self.func = func
        self.last_arg = None
        
    def __call__(self, *args, **kwargs):
        if np.all(args[0] == self.last_arg):
            #print "Cache hit:", args[0]
            return self.last_result
        else:
            #print "Cache miss:", args[0]
            self.last_arg = np.copy(args[0])
            self.last_result = self.func(*args, **kwargs)
            return self.last_result

    def f(self, *args, **kwargs):
        return self.__call__(*args, **kwargs)[0]

    def fprime(self, *args, **kwargs):
        return self.__call__(*args, **kwargs)[1]

def spopt_complex(func):
    """ Decorator that makes a complex function accessible to scipy.optimize """
    def wrapper(*args, **kwargs):
        x = np.asarray(args[0])
        x = x[:x.size//2] + 1j*x[x.size//2:]
        res = func(x, *args[1:], **kwargs)
        if isinstance(res,tuple):
            print("   diff = %8.2e" % np.linalg.norm(res[0]))
            return (np.concatenate((res[0].real, res[0].imag)),
                    np.vstack((np.hstack((res[1].real, -res[1].imag)),
                               np.hstack((res[1].imag,  res[1].real)))))
        else:
            return np.concatenate((res.real, res.imag))
    
    return wrapper

def get_gw_jac(sigma, wplusmu, hk, gtarget=None, jac=False):
    """ Returns G for a diagonal Sigma plus the Jacobian dG/dSigma for one w """
    print(" sigma(w) =", cformat("%10g%+10gj ", sigma), end=' ')
    denom = np.diag(wplusmu - sigma)
    gw = np.zeros(shape=(sigma.size,sigma.size), dtype=np.cdouble)
    
    # do the k summation
    if jac:
        jacobian = gw.copy()
        for ham in hk:
            inv = np.linalg.inv(denom - ham)
            gw += inv
            jacobian += np.dot(inv,inv)
        jacobian /= hk.shape[0]
    else:
        for ham in hk: 
            gw += np.linalg.inv(denom - ham)    
    
    # check that off-diagonal terms vanish
    gw /= hk.shape[0]
    gwdiag = gw.diagonal().copy()
    if gtarget is not None:
        gwdiag -= gtarget
    #
    if jac:
        return gwdiag, jacobian
    else:
        return gwdiag

get_gw_jacr = spopt_complex(get_gw_jac)
get_gw_jacr.__doc__ = "Version of get_gw_jac that works with real arrays"

class wrapper_folding:
    """Consider only some bands and set rest of self-energies to small value"""
    def __init__(self, func, folding, const_sigmas):
        self.func = func
        self.mapping = np.asarray(folding, dtype=int)
        self.const_sigmas = np.asarray(const_sigmas, dtype=np.cdouble)
        unique, inverse = np.unique(self.mapping, return_index=True)
        self.inv_mapping = inverse[(unique[0] == -1):]
        assert self.mapping.size == self.const_sigmas.size

    def __call__(self, *args, **kwargs):
        firstarg = np.where(self.mapping == -1, self.const_sigmas,
                            args[0][self.mapping])
        fourtharg = args[3]
        if fourtharg is not None:
            fourtharg = np.where(self.mapping == -1, 0, fourtharg[self.mapping])
        res = self.func(firstarg, args[1], args[2], fourtharg, *args[4:], **kwargs)
        if isinstance(res, tuple):
            return (res[0][self.inv_mapping],
                    res[1][self.inv_mapping[:,None], self.inv_mapping[None,:]])
        else:
            return  res[0][self.inv_mapping]

def find_sigma(gw, w, hk, mu=0, get_gw=get_gw_jacr, sigma_guess=None):
    """ Find self energy Sigma(w) for a Green's function G(w) for one omega """
    if sigma_guess is None:
       sigma_guess = np.zeros_like(gw)
    nbands = gw.size
    get_gw = spopt_cache(get_gw)

    # to keep compatibility with scipy <= 0.10, we use fsolve instead of root
    x, infodict, success, message = scipy.optimize.fsolve(
                func=get_gw.f, fprime=get_gw.fprime,
                x0=np.concatenate((sigma_guess.real, sigma_guess.imag)),
                args=(w, hk, gw, True), full_output=True)
    #
    #jac_r = np.zeros((x.size,x.size))
    #jac_r[np.triu_indices_from(jac_r)] = infodict["r"]
    #jac_q = np.array(infodict["fjac"]).T
    if success != 1:
        raise ConvergenceError(message)
    return x[:nbands] + 1j * x[nbands:]

def hk_decompose(hks):
    """ Decomposes the Hamiltonian into its eigenbasis for any k-point.
    
    If E, V is the result, then the following holds for any ik:
    >>> E, V = hk_decompose(Hk)
    >>> np.dot(np.dot(V[k], np.diag(E[k])), V[k].T) == Hk[k] 
    """
    eigenvals = np.empty(shape=(hks.shape[:2]))
    eigenbasis = np.empty_like(hks)
    for ik, hk in enumerate(hks):
        eigenvals[ik], eigenbasis[ik] = np.linalg.eigh(hk)
    return eigenvals, eigenbasis

def intmask(array):
    """ For a sorted array, returns a mask for the nearest integers """
    nint = np.rint(array)
    diff = np.hstack((np.nan, np.abs(array - nint), np.nan))
    nint = np.hstack((np.iinfo(int).min, nint, np.iinfo(int).max))
    return np.logical_and(
                np.logical_or(nint[1:-1] != nint[:-2], diff[1:-1] <= diff[:-2]),
                np.logical_or(nint[1:-1] != nint[2:],  diff[1:-1] < diff[2:])
                )
    
def floatmask(arr, sl):
    """Takes a slice of floats and returns a mask for an array"""
    mask = np.empty(shape=arr.shape, dtype=np.bool8)
    mask[...] = True
    if sl.step is not None and sl.step < 0:
        sl = slice(sl.stop, sl.start, -sl.step)
    if sl.start is not None:
        mask = np.logical_and(mask, arr >= sl.start)
    if sl.stop is not None:
        mask = np.logical_and(mask, arr < sl.stop)
    if sl.step is not None:
        mask = np.logical_and(mask, intmask((arr - arr[0])/sl.step))
    return mask

def read_aw_files(*filenames):
    """ Reads A(w) files produces by maxent.py, returns w, A(w) """
    if not filenames:
        raise ValueError("expected at least one filename")
    ws = None
    for iband, awfile in enumerate(filenames):
        data = np.fromfile(awfile, sep=" ").reshape(-1,2)
        if ws is None:
            ws = data[:,0]
            nw = ws.size
            aws = np.empty(shape=(len(filenames),nw))
        elif not np.allclose(ws, data[:,0]):
            raise RuntimeError("Inconsistent frequencies")
        aws[iband,:] = data[:,1]
    return aws, ws

def cformat(fmt, iterable):
    """ Formats a list of complex numbers to a format specification """
    return "".join(map(lambda c: fmt % (c.real, c.imag), iterable))

if __name__ == "__main__":
    import optparse
    from w2dyn.auxiliaries.input import read_hamiltonian
    import sys, re
    
    parser = optparse.OptionParser(version="1.0",
                usage="%prog [-H <ham_file>] [-M <mu>] " +
                      "([file:]<aw_file> | symm:<ref> | sigma{<n>}:<const>) ...",
                description="""Computes the self energy from the spectral function.""")
    parser.add_option("-H", metavar="FILE", dest="hk_file",
                      help="Hamiltonian file H(k) to use")
    parser.add_option("-M", type="float", dest="mu", 
                      default=0., help="chemical potential")
    parser.add_option("-D", type="float", dest="deltino",
                      default=0., help="small imaginary part (>= 0, optional)")
    parser.add_option("--wmin", type="float", dest="wmin",
                      help="beginning of frequency window to use")
    parser.add_option("--wmax", type="float", dest="wmax",
                      help="end of frequency window to use (exclusive)")
    parser.add_option("--wstep", type="float", dest="wstep",
                      help="minimum frequency step size (for quicker runs)")
    parser.add_option("--reversed", action="store_true", dest="wreversed",
                      default=False, help="reverse the frequency iteration")
    parser.add_option("--gw", metavar="FILE", dest="gw_file",
                      default="gw.dat", help="file to store G(w) in")
    parser.add_option("--sigmaw", metavar="FILE", dest="sigmaw_file", 
                      default="sigmaw.dat", help="file to store Sigma(w) in")
    options, args = parser.parse_args()

    if options.deltino < 0: parser.error("deltino must be positive")

    aw_files = []
    sigmas = []
    folding = []
    arg_re = re.compile(r"""^(?: ([a-z]+)            # prefix name
                                 (?: { ([^}]*) } )?  # prefix arg
                              :)? (.+) $             # argument
                         """, re.X)
    try:
        for arg in args:
            prefix, prefix_arg, arg = arg_re.match(arg).groups()
            sigma = 0
            repeats = 1
            if not prefix or prefix == "file":
                ref = len(aw_files)
                aw_files.append(arg)
            elif prefix == "symm":
                ref = int(arg) - 1
            elif prefix == "sigma":
                if prefix_arg: repeats = int(prefix_arg)
                ref = -1
                sigma = complex(arg)
            else:
                raise ValueError("invalid prefix")
            folding += [ref] * repeats
            sigmas += [sigma - 1j*options.deltino] * repeats
    except Exception as e:
        parser.error("error parsing argument %s: %s" % (arg, e))

    print("Reading spectral functions A(w) ...")
    aws, ws = read_aw_files(*aw_files)

    if options.hk_file is not None:
        print("Reading Hamiltonian H(k) ...")
        hk = read_hamiltonian(file(options.hk_file,"r"))[0]
    else:
        print("Warning: setting Hamiltonian H(k) to zero (use -H)")
        if len(args) < 1: parser.error("expecting A(w) file as argument")
        hk = np.zeros((1,len(args),len(args)))

    nbands = hk.shape[1]
    if len(folding) != nbands:
        parser.error("expected %d arguments, got %d", nbands, len(folding))

    # Build the function that is used as argument for the root finding here.
    # In principle there are four levels added in reverse order
    #   3. spopt_cache - caches previously obtained G(w), Sigma(w) (find_sigma)
    #   2. spopt_complex - assembles complex Sigma(w) from real/imaginary part
    #   1. wrapper_folding - adds redundant/non-interacting bands to Sigma(w)
    #   0. get_gw_jac - computes G(w) for some full Sigma(w)
    get_gw = get_gw_jac
    if len(aw_files) != len(folding):
        get_gw = wrapper_folding(get_gw, folding, sigmas)
        print("Default values: %s" % (sigmas))
        print("Folding map %s (inverse %s)" % (get_gw.mapping, get_gw.inv_mapping))
    get_gw = spopt_complex(get_gw)

    print("Using Kramers-Kronig to get real part of Green's function...")
    gws = aws_to_gws(aws, ws)

    if options.deltino:
        print("Adding deltino=-%gj to Green's function ..." % options.deltino)
        gws = 1/(1/gws + 1j*options.deltino)

    # selecting only parts 
    wmask = floatmask(ws, slice(options.wmin, options.wmax, options.wstep))
    if options.wmin is not None:
        options.gw_file += ".from%g" % options.wmin
        options.sigmaw_file += ".from%g" % options.wmin

    # use ws
    ws = ws.compress(wmask)
    gws = gws.compress(wmask,-1)

    # revert the frequency order?
    if options.wreversed:
        ws = ws[...,::-1]
        gws = gws[...,::-1]
        options.gw_file += ".rev"
        options.sigmaw_file += ".rev"

    nw = ws.size
    print("Chosen %d out of %d frequencies ..." % (ws.size, nw))

    # write to file
    f = open(options.gw_file,"w")
    for iw, w in enumerate(ws):
        f.write("%12.6f" % w + cformat("%12.6f %12.6f   ", gws[:,iw]) + "\n")
    f.close()
                                             
    print("Finding Sigma(w) (this will take some time)...")
    f = open(options.sigmaw_file,"w")
    f.write("# w, Re Sigma(w), Im Sigma(w) (for each band)\n")
    sigma = np.zeros_like(gws[:,0])
    sigma[:] = -1j * options.deltino
    try:
        for iw, w in enumerate(ws):
            sys.stdout.flush()
            print("w = %9g (frequency %d of %d)" % (w, iw+1, nw))
            sigma = find_sigma(gws[:,iw], w, hk, options.mu, get_gw, sigma)
            if (sigma.imag > 0).any():
                print("WARNING: positive imaginary part of Sigma(w)!")
            f.write("%12.6f" % w + cformat("%12.6f %12.6f   ", sigma) + "\n")
            f.flush()
    except ConvergenceError:
        sys.stderr.write(_CONVERR_MESSAGE  % (w, options.deltino))
        sys.exit(1)
    finally:
        f.close()
    sys.exit(0)
