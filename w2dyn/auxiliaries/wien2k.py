""" @package wien2k
Provides methods for charge selfconsistency with Wien2k
"""
from __future__ import absolute_import, print_function, unicode_literals
import os
import numpy        as np
import numpy.linalg as la

def init(cfg, kpoints, notify=open(os.devnull, 'w')):
    global klist

    if cfg['General']['KListFile'] is None:
        cfg['General']['KListFile'] = os.path.basename(os.getcwd()) + '.klist'

    klist_file = open(cfg['General']['KListFile'], 'r')

    klist = read_klist(klist_file)

    klist_file.close()

    print("Read %g k-points from `%s' for charge selfconsistency"
          % (len(klist), cfg['General']['KListFile']), file=notify)

    klist = check_klist(klist, kpoints, cfg)


def delta_N(beta, w, mu_lda, mu_dmft, Hk, Siw, DC, nat, nneq):
    """DeltaN(k) for charge selfconsistency

    @return
         dnk[nbands, nbands, 2, nk_wien2k]: DeltaN(k) as per
         Eq. A6, Lechermann et al. PRB 74, 125120
    """

    nbands = Hk.shape[1]
    nw     = len(w)
    try:
        nk = len(klist)
    except NameError:
        raise Exception("`klist' not defined -- did you call wien2k.init()?")
 
    #TODO: Hk and gks need to be changed once supporting spin-orbit
    assert Hk .shape[1:4:2] == (nbands, nbands)
    assert Siw.shape      == (nbands, 2, len(w))
    assert DC .shape      == (nbands, 2)

    mu  = mu_dmft*np.eye(nbands, dtype=complex)
    mu0 = mu_lda *np.eye(nbands, dtype=complex)

    dnk = np.zeros((nbands, nbands, 2, nk), dtype=complex, order='F')
    
    for (iw, w) in enumerate(w):
        sigma = Siw[..., iw] + DC

        w = 1j*w*np.eye(nbands, dtype=complex)

        for ik in range(nk):
            #assuming Hk without spin-orbit
            h = Hk[ik, :, 0, :, 0]

            # gks: Kohn-Sham Green function
            # g  : interacting Green function
            #
            # We must resort to lattice.hk and an inversion because
            # we need this k-resolved

            gks = la.inv(w + mu0 - h)
            g   = np.dstack(( la.inv(w + mu - h - np.diag(sigma[:, 0])),
                              la.inv(w + mu - h - np.diag(sigma[:, 1])) ))

            dnk[..., ik] += \
                np.einsum(gks, [0,1], sigma, [1,2], g, [1,3,2], [0,3,2])/beta

    return dnk

def check_klist(klist, kpoints, cfg):
    r""" Check `klist' (from a Wien2k `klist' file) against `kpoints'
    (from H(k))

    @return
          a list of indices into `kpoints' corresponding to `klist'
    """

    # klist should be a subset of the k-points of H(k), and we expect
    # the same order for the two
    import collections
     
    kwien = collections.deque(klist)
    kham  = collections.deque(kpoints)
    rest  = collections.deque()
    found = collections.deque()
    
    w = kwien.popleft()
    ih = -1
    try:
        # First, sort out all k-vectors which appear in the same order
        # in both arrays.  Usually, that will be all of them.
        while kham:
            h = kham.popleft()
            ih += 1

            if np.allclose(w, h):
                found.append(ih)
                w = kwien.popleft()
            else:
                rest.append((ih, h))

        # If we got here, the orderings did not agree and we have to
        # go through the rest.
        while True:
            for (i, r) in enumerate(rest):
                if np.allclose(w, r[1]):
                    rest.rotate(i)
                    rest.popleft()
                    found.append(r[0])
                    break
                else:
                    raise Exception(
                        """K-vector %s present in `%s' not found in `%s'.
For ComputeDeltaN, the k-vectors in KListFile must be a subset of those in HkFile""" \
                            % (w, cfg['General']["KListFile"], cfg['General']["HkFile"])
                        )

                w = kwien.popleft()
    except IndexError:
        pass              # all found

    return np.array(found)


def read_klist(file):
    r""" Reads a Wien2k ``klist'' file
    
    @return 
        The k-points in a numpy.array((#k, 3), dtype=float)
    """

    from collections import deque

    klist = deque()

    ended = False
    for line in file:
        # The k-list must be terminated by "END"
        if line[0:3] == "END":
            ended = True
            break

        # Fortran line format is (A10, 4I10, ...) [Wien2k UG].  We
        # read this using split, which might conceivably lead to an
        # error (if someone is crazy enough to give 10-digit
        # k-coordinates).

        tok = line[10:].split(None, 4)

        klist.append(np.array(tok[:3], dtype=float) / int(tok[3]))

    assert ended, "klist file must be terminated by `END'"

    return np.array(klist)
