""" @package input
Provides methods for interacting with input files
"""
from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
import re
from warnings import warn

import numpy as np
import h5py
import sys
h5ustrs = h5py.special_dtype(vlen=(str
                                   if sys.version_info >= (3, )
                                   else unicode))


def read_u_matrix(u_file, spin=False):
    r""" Reads a full four-index U matrix f$ U_{ijkl} f$ from a text file.

    Expects a text file with a header of the following form:

        <no of orbitals> BANDS       # optional comment

    followed by a list of entries of the U matrix.  For a spin-independent U
    (SU2-invariant) matrix (spin is False), expects only the orbitals and a
    value, e.g.:

        1   1   1   1    2.50
        1   2   1   2    1.25

    For the full spin-dependent (spin is True), it also expects a spin index,
    which is either "u" (spin-up) or "d" (spin-down), e.g.:

        1 u 1 d 1 u 1 d  2.50
        1 d 1 u 1 d 1 u  2.50

    Returns the U-matrix as either four- or eight-dimensional numpy array.
    """
    if not hasattr(u_file, "readline"):
        u_file = open(u_file, "r")

    # reading header of the form "3 BANDS", possibly preceded by comments
    header_re = re.compile(r"\s*(?:([\d]+)\s+bands?\s*)?(?:\#.*)?$",
                           re.I | re.X)
    norbitals = None
    while norbitals is None:
        header = u_file.readline()
        match = header_re.match(header)
        if not match:
            raise ValueError("expecting U matrix header `n BANDS`:\n%s"
                             % header)
        norbitals = match.group(1)
    norbitals = int(norbitals)

    # now read the remaining lines
    if spin:
        line_re = r"\s*" + r"(\d+)\s*([ud])\s+" * 4 + r"([^\s]*)\s*$"
    else:
        line_re = r"\s*" + r"(\d+)\s+" * 4 + r"([^\s]*)\s*$"
    line_re = re.compile(line_re, re.I | re.X)

    entries = []
    u_file.flush()
    for line in u_file:
        def get_orbital(orb):
            orb = int(orb)
            if orb < 1 or orb > norbitals:
                raise ValueError("invalid orbital: %d\n%s" % (orb, line))
            return orb - 1
        def get_spin(sp):
            return "du".index(sp)

        match = line_re.match(line)
        if not match:
            raise ValueError("invalid U matrix entry:\n%s" % line)
        entry = list(match.groups())
        entry[-1] = float(entry[-1])

        if spin:
            entry[:-1:2] = map(get_orbital, entry[:-1:2])
            entry[1::2] = map(get_spin, entry[1::2])
        else:
            entry[:-1] = map(get_orbital, entry[:-1])

        entries.append(entry)

    # now re-model the entries as U matrix
    entries = np.array(entries)
    if spin:
        u_matrix = np.zeros((norbitals, 2) * 4, float)
    else:
        u_matrix = np.zeros((norbitals,) * 4, float)

    indices = tuple(np.array(entries[:, :-1].T, int))
    values = entries[:, -1]
    u_matrix[indices] = values
    return u_matrix


def write_u_matrix(out, u_matrix, comment=None):
    """Writes a full four-index U matrix f$ U_{ijkl} f$ to a text file.

    Expects a file object out that can be written to, u_matrix as
    either a four-dimensional array with (flavor, flavor, flavor,
    flavor) axes or an eight-dimensional array with a (band, spin)
    pair of axes per flavor and writes the data to the file in the
    text format that can be read using read_u_matrix, i.e.

        <no of orbitals> BANDS

    followed by a list of entries of the U matrix.  For a
    spin-independent U (SU2-invariant) matrix (automatically
    determined from u_matrix argument), these consist of orbital
    indices starting at one followed by the value of the matrix
    element:

        1   1   1   1    2.50
        1   2   1   2    1.25

    and for a spin-dependent matrix, a spin index written as "d" for
    the first and "u" for the second element along the axis follows
    each orbital index:

        1 u 1 d 1 u 1 d  2.50
        1 d 1 u 1 d 1 u  2.50

    A string optionally passed in comment is written first as a
    comment (will be prefixed with #); if it contains newlines,
    comment characters for subsequent lines are not taken care of.

    """
    if comment is not None:
        out.write("# {}\n".format(comment))

    # reshape flavor-axes array to (band, spin)-axes array
    if len(u_matrix.shape) == 4:
        u_matrix = np.reshape(u_matrix,
                              sum(tuple((dim//2, 2)
                                        for dim
                                        in u_matrix.shape),
                                  ()))
    elif len(u_matrix.shape) != 8:
        raise ValueError("u_matrix argument has wrong shape")

    out.write("{} BANDS\n".format(u_matrix.shape[0]))

    # check whether the interaction fits spin-independent format
    absoluteu = np.abs(u_matrix)
    atol = 1.0e-8 * np.amin(absoluteu[absoluteu > 0.0])
    if all((np.allclose(u_matrix[:, 0, :, 0, :, 0, :, 0],
                        u_matrix[:, s1, :, s2, :, s3, :, s4],
                        atol=atol)
            if (s1, s2, s3, s4) in ((0, 0, 0, 0),
                                    (0, 1, 0, 1),
                                    (1, 0, 1, 0),
                                    (1, 1, 1, 1))
            else np.allclose(0.0,
                             u_matrix[:, s1, :, s2, :, s3, :, s4],
                             atol=atol))
           for s1, s2, s3, s4 in np.ndindex(u_matrix.shape[1:8:2])):
        orbmatrix = u_matrix[:, 0, :, 0, :, 0, :, 0]
        nonzero_indices = np.nonzero(orbmatrix)
        np.savetxt(out,
                   np.concatenate((np.transpose(nonzero_indices) + 1,
                                   orbmatrix[nonzero_indices][:, np.newaxis]),
                                  axis=1),
                   fmt="%d %d %d %d %.18e")
    else:
        nonzero_indices = np.transpose(np.nonzero(u_matrix))

        for indices in nonzero_indices:
            value = u_matrix[tuple(indices)]
            # 1-based indices for bands and 'd' and 'u' for spins
            for icol, val in enumerate(indices):
                if icol % 2 == 0:
                    indices[icol] += 1
                else:
                    indices[icol] = (
                        ord('d') + (ord('u') - ord('d')) * val)
            out.write(("{} {:c} " * 4 + "{u_val:.18e}\n").format(*indices, u_val=value))


def read_epsk_vk_file(epsk_file, vk_file):
    r"""Reads an epsk and a vk file.

    Returns a triple:
      - diag:  Whether the hybridisation is spin/orbital diagonal
      - epsk:  the bath energies as 1-d nbath array
      - vki:   the hoppings from impurity to bath (as nbath x nimp array)
    """
    epsk = np.loadtxt(epsk_file)
    vki = np.loadtxt(vk_file)

    if epsk.ndim == 1: epsk = epsk[:,None]
    if vki.ndim == 1: vki = vki[:,None]

    if epsk.shape != vki.shape:
        raise ValueError("epsk and Vk files must agree in shape")

    # Add a spin dimension
    epsk = epsk[:,None].repeat(2, -1)
    vki = vki[:,:,None].repeat(2, -1).reshape(vki.shape[0], -1)

    epsk = epsk.reshape(-1)
    vki = np.vstack([np.diag(vki_row) for vki_row in vki])
    assert epsk.shape[0] == vki.shape[0]
    return True, epsk, vki

def read_muimp(muimp_file):
    r"""Reads a muimp (diagonal entries) file.

    Returns:
      - muimp:  the impurity chemical potentials as 4-d nbands x 2 x nbands x 2 array
    """
    muimp = np.loadtxt(muimp_file)

    assert muimp.size % 2 == 0
    Norb = muimp.size // 2
    muimp = np.reshape(muimp, (Norb, 2))
    fullmuimp = np.zeros((Norb, 2, Norb, 2), dtype=muimp.dtype)
    for i in range(Norb):
        for j in range(2):
            fullmuimp[i][j][i][j] = muimp[i][j]
    return fullmuimp

def read_Delta_iw_tau(deltaiw_file, deltatau_file, beta):
   r"""Reads the imaginary time and Matsubara hybridization function.
   
   Expects two text files:
   
   1. The matsubara hyb function (deltaiw_file) in the format:
   
      iv-val re_up_1 im_up_1 re_down_1 im_down_1 ...
   
   where the last index is for the orbital. 
   The matsubara hybridization function from file needs to be defined
   for positive matsubaras only.
   
   2. The imaginary time hyb function (deltatau_file) in the format:
   
      tau-val up_1 down_1 ...
   
   Returns:
      - hybriv: hybridization function in matsubara frequencies
      - hybrtau: hybridization function in imaginary time
      - iwf: fermionic frequency array
      - nftau: number of tau points
      - tauf: tau values
   
   FIXME: currently only implemented for diagonaly hybridization functions
   FIXME: function name conflicts with nano hybridization read-in
   """
   from .transform import matfreq

   #ignoring explicit matsubara and tau values in file
   hybriv = np.loadtxt(deltaiw_file)[:,1:]
   hybrtau = np.loadtxt(deltatau_file)[:,1:]
   
   #convert float array to complex
   hybriv = hybriv[:,0::2] + 1j*hybriv[:,1::2]
   
   #generate conjugate
   hybriv = np.resize(hybriv,(2*hybriv.shape[0],hybriv.shape[1]))
   hybriv[:int(hybriv.shape[0]/2),:] = np.conj(hybriv[int(hybriv.shape[0]/2):,:])[::-1,:]

   #generating tau and matsubara arrays
   #tau interval goes from [0,beta]
   niw = hybriv.shape[0]
   iwf = matfreq(beta, 'fermi', niw)
   nftau = hybrtau.shape[0]
   tauf = np.linspace(0., beta, num=nftau, endpoint=True)
   
   return hybriv, hybrtau, niw, iwf, nftau, tauf

def read_Delta_iw_tau_full(deltaiw_file, deltatau_file, beta, norbitals):
   r"""Reads the imaginary time and Matsubara hybridization function.
   
   Expects two text files:
   
   1. The matsubara hyb function (deltaiw_file) in the format:

      flavour1, flavour2, iv, Re(value), Im(value)

      e.g.

      # flav1 flav2     iv           Re             Im 
          0     0  -313.84511   3.869071e-07   -0.001049714
          0     0  -313.21679    3.88461e-07    -0.00105182
          0     0  -312.58847   3.900242e-07   -0.001053934
          0     0  -311.96015   3.915969e-07   -0.001056057
   
      The matsubara hybridization function from file needs to be defined
      for negative and positive matsubaras.
   
   2. The imaginary time hyb function (deltatau_file) in the format:

      flavour1, flavour2, tau, value

      e.g.

      # flav1 flav2    tau          value
          0     0    0.00000      0.1819913
          0     0    0.00500      0.1817078
          0     0    0.01000       0.181425
          0     0    0.01500      0.1811428

   Returns:
      - hybriv: hybridization function in matsubara frequencies
      - hybrtau: hybridization function in imaginary time
      - iwf: fermionic frequency array
      - nftau: number of tau points
      - tauf: tau values
   
   FIXME: function name conflicts with nano hybridization read-in
   """
   from .transform import matfreq

   # load the files
   hybriv = np.loadtxt("fiw.dat")
   hybrtau = np.loadtxt("ftau.dat")

   # extract how many points
   niw=int(hybriv.shape[0]/(2*norbitals)**2)
   nftau=int(hybrtau.shape[0]/(2*norbitals)**2)

   # check integrity of the files
   if abs(niw-hybriv.shape[0]/(2*norbitals)**2) > 1e-10:
      #print "niw", niw
      #print "hybriv.shape[0]/(2*norbitals)**2 ",  hybriv.shape[0]/(2*norbitals)**2
      exit(-1)

   if abs(nftau-hybrtau.shape[0]/(2*norbitals)**2) > 1e-10:
      #print "nftau", nftau
      #print "hybrtau.shape[0]/(2*norbitals)**2 ",  hybrtau.shape[0]/(2*norbitals)**2
      exit(-1)

   # create arrays to store the functions in
   ftau=np.zeros(shape=(nftau,norbitals*2,norbitals*2))
   fiw=np.zeros(shape=(niw,norbitals*2,norbitals*2),dtype=complex)

   # loop over flavour**2
   for i in range(0,(2*norbitals)**2):

      # line in full array, where a new (flavour,flavour) block begins
      firstpoint=i*nftau+1
      # the flavours that are currently read
      flav1=int(hybrtau[firstpoint,0])
      flav2=int(hybrtau[firstpoint,1])

      # here comes the tau 
      ftau[:,flav1,flav2]=hybrtau[i*nftau:i*nftau+nftau,3]
      tau=hybrtau[i*nftau:i*nftau+nftau,2]

      # and here comes the Matsubara 
      firstpoint=i*niw+1
      #print "firstpoint", firstpoint
      flav1=int(hybriv[firstpoint,0])
      flav2=int(hybriv[firstpoint,1])
      iv=hybriv[firstpoint,2]
      value=hybriv[firstpoint,3]

      fiw[:,flav1,flav2]=hybriv[i*niw:i*niw+niw,3]+1.0j*hybriv[i*niw:i*niw+niw,4]
      iv=hybriv[i*niw:i*niw+niw,2]
      #print "iv", iv
      #print "tau", tau

   return fiw, ftau, niw, iv, nftau, tau

def read_hamiltonian(hk_file, spin_orbit=False):
    r""" Reads a Hamiltonian f$ H_{bb'}(k) f$ from a text file.

    Expects a text file with white-space separated values in the syntax as
    generated by wannier90:  the first line is a header with three integers,
    optionally followed by '#'-prefixed comment:

        <no of k-points> <no of wannier functions> <no of bands (ignored)>

    For each k-point, there is a header line with the x, y, z coordinates of the
    k-point, followed by <nbands> rows as lines of 2*<nbands> values each, which
    are the real and imaginary part of each column.

    Returns: A pair (hk, kpoints):
     - hk       three-dimensional complex array of the Hamiltonian H(k),
                where the first index corresponds to the k-point and the other
                dimensions mark band indices.
     - kpoints  two-dimensional real array, which contains the
                components (x,y,z) of the different kpoints.
    """
 
    def nextline():
        line = hk_file.readline()
        return line[:line.find('#')].split()

    # parse header
    header = nextline()
    if header[0] == 'VERSION':
        warn("Version 2 headers are obsolete (specify in input file!)")
        nkpoints, natoms = map(int, nextline())
        lines = np.array([nextline() for _ in range(natoms)], int)
        nbands = np.sum(lines[:,:2])
        del lines, natoms
    elif len(header) != 3:
        warn("Version 1 headers are obsolete (specify in input file!)")
        header = list(map(int, header))
        nkpoints = header[0]
        nbands = header[1] * (header[2] + header[3])
    else:
        nkpoints, nbands, _ = map(int, header)
    del header

    # nspins is the spin dimension for H(k); G(iw), Sigma(iw) etc. will always
    # be spin-dependent
    if spin_orbit:
        if nbands % 2: raise RuntimeError("Spin-structure of Hamiltonian!")
        nbands //= 2
        nspins = 2    
    else:
        nspins = 1
#GS: inside read_hamiltonian nspins is therefore equal to 1 if spin_orbit=0
#GS: outside nspins is however always set to 2. Right?
#GS: this also means that we have to set nspins internally in read_ImHyb too

    hk_file.flush()

    # parse data
    hk = np.fromfile(hk_file, sep=" ")
    hk = hk.reshape(-1, 3 + 2 * nbands**2 * nspins**2)
    kpoints_file = hk.shape[0]
    if kpoints_file > nkpoints:
        warn("truncating Hk points")
    elif kpoints_file < nkpoints:
        raise ValueError("problem! %d < %d" % (kpoints_file, nkpoints))
    kpoints = hk[:nkpoints, :3]

    # TODO: check with Martin if this the actual spin structure ... 
    hk = hk[:nkpoints, 3:].reshape(nkpoints, nspins, nbands, nspins, nbands, 2)
    hk = hk[...,0] + 1j * hk[...,1]
    if not np.allclose(hk, hk.transpose(0,3,4,1,2).conj()):
        warn("Hermiticity violation detected in Hk file")

    # go from wannier90/convert_Hamiltonian structure to our Green's function
    # convention
    hk = hk.transpose(0, 2, 1, 4, 3)
    return hk, kpoints

#GS:
def read_ImHyb(hyb_file, norbitals, spin_orbit):
    r""" Reads the imaginary part of ImDelta(omega) for each lead from hyb_file

    Expects a file with a similar structure as hk_file, 
    with real-frequency values instead of k-points headers. 
  
    Each frequency header is followed by nleads square 
    
         <nbands>*<nspins> times <nbands>*<nspins> 

    matrices containing the imaginary part of the Hybridization functions for
    each lead. We assume that we do not need explicit lead-lead hybridizations.
    """

#GS: nspins is the spin dimension for H(k); G(iw), Sigma(iw) etc. will always be spin-dependent
#GS: outside nspins is therefore always 2 (assuming not to have any Nambu....)    
#GS: inside here, nbands is the dimension of the Hyb-matrix in hyb_file. 
#GS: Set it to the right size using the (internal) nspins, similarly to read_hk:
    if spin_orbit:
        nspins = 2
    else:
        nspins = 1

    nbands = norbitals*nspins

#GS(stupid check): 
#    line=hyb_file.readline()
#    print "line = ", line
#    line=line[:line.find('#')].split()
#    print "line[:line...] = ", line[0], line[1]
#    numw_hyb = int(line[0])
#    nleads = int(line[1])
#   
#    print "numw_hyb = ", numw_hyb
#    print "nleads = ", nleads
#GS: 
    def nextline():
        line = hyb_file.readline()
        return line[:line.find('#')].split()
# 
    header = nextline()
#   print "header = ", header
#GS:
#GS: and then we need a check that this is compatible with the number of columns contained in hyb_file
    numw_hyb = int(header[0])
    nleads = int(header[1])
    del header


    hyb_file.flush()

#   print "numw_hyb = ", numw_hyb
#   print "nleads = ", nleads
#   print "nbands = ", nbands
#   print "nspins (read_ImHyb internal) = ", nspins

    # parse data
    leadsw = np.fromfile(hyb_file, sep=" ")
    leadsw = leadsw.reshape(-1, 1 + nleads * nbands**2 * nspins**2)
    numw_in_file = leadsw.shape[0]
    if numw_in_file > numw_hyb:
        warn("truncating frequencies in Hybridization to leads")
    elif numw_in_file < numw_hyb:
        raise ValueError("problem! %d < %d" % (numw_in_file, numw_hyb))
    w_hyb = leadsw[:numw_hyb, 0]
#GS: Assuming that the leads are listed as matrices sequentially for each frequency, the following:
#GS: leadsw = leadsw[:numw_hyb, 1:].reshape(numw_hyb, nspins, nbands, nspins, nbands, nleads)
#GS: was wrong, right? 
    leadsw = leadsw[:numw_hyb, 1:].reshape(numw_hyb, nleads, nspins, nbands, nspins, nbands)

#GS: do we have to check any symmetry property of leadsw before going on?

    # go from wannier90/convert_Hamiltonian structure to our standard Green's function structure
#GS: In order to use it as argument for the function transform, the frequency must be the last argument:
#   leadsw = leadsw.transpose(0, 5, 2, 1, 4, 3)   # w, lead, band, spin, band, spin
    leadsw = leadsw.transpose(1, 3, 2, 5, 4, 0)   # lead, band, spin, band, spin, w

#GS: if the file contains ImDelta(w), I have to perform the integral for -(1/pi)*ImDelta(w):
    leadsw = -leadsw/np.pi
  
    return leadsw, w_hyb, nleads



#AV: slightly modified version of read_ImHyb
def read_Delta(hyb_file, norbitals, spin_orbit):
    r""" Reads the real and imaginary part of Delta(omega) for each lead from hyb_file

    Expects a file with a similar structure as hk_file, 
    with real-frequency values instead of k-points headers. 
  
    Each frequency header is followed by nleads square 
    
         <nbands>*<nspins> times <nbands>*<nspins> 

    matrices containing the real and imaginary part of the Hybridization functions for
    each lead. We assume that we do not need explicit lead-lead hybridizations.
    """

#GS: nspins is the spin dimension for H(k); G(iw), Sigma(iw) etc. will always be spin-dependent
#GS: outside nspins is therefore always 2 (assuming not to have any Nambu....)    
#GS: inside here, nbands is the dimension of the Hyb-matrix in hyb_file. 
#GS: Set it to the right size using the (internal) nspins, similarly to read_hk:
    if spin_orbit:
        nspins = 2
    else:
        nspins = 1

    nbands = norbitals*nspins

    def nextline():
        line = hyb_file.readline()
        return line[:line.find('#')].split()
# 
    header = nextline()
#   print "header = ", header
#GS:
#GS: and then we need a check that this is compatible with the number of columns contained in hyb_file
    numw_hyb = int(header[0])
    nleads = int(header[1])
    del header

    hyb_file.flush()

    # parse data
    leadsw = np.fromfile(hyb_file, dtype=float, sep=" ")
#AV: reshape to take complex nature (i.e., 2 floats per entry) of the data into account
    leadsw = leadsw.reshape(-1, 1 + 2 * nleads * nbands**2 * nspins**2)
    numw_in_file = leadsw.shape[0]
    if numw_in_file > numw_hyb:
        warn("truncating frequencies in Hybridization to leads")
    elif numw_in_file < numw_hyb:
        raise ValueError("problem! %d < %d" % (numw_in_file, numw_hyb))
    w_hyb = leadsw[:numw_hyb, 0]
    leadsw = leadsw[:numw_hyb, 1:].reshape(numw_hyb, nleads, nspins, nbands, nspins, nbands, 2)
#AV: extract real and imaginary part separately and recombine into a complex vector
    leadsw = leadsw[...,0] + 1j * leadsw[...,1]

#GS: do we have to check any symmetry property of leadsw before going on?

    # go from wannier90/convert_Hamiltonian structure to our standard Green's function structure
#GS: In order to use it as argument for the function transform, the frequency must be the last argument:
#   leadsw = leadsw.transpose(0, 5, 2, 1, 4, 3)   # w, lead, band, spin, band, spin
    leadsw = leadsw.transpose(1, 3, 2, 5, 4, 0)   # lead, band, spin, band, spin, w

#AV: we are supposed to read ReDelta and ImDelta, so we should not do this...
#GS: if the file contains ImHyb(w), I have to perform the integral for -(1/pi)*ImDelta(w):
#    leadsw = -leadsw/np.pi
  
    return leadsw, w_hyb, nleads
   

  
def read_uomega(uw_file):
    def nextline():
        line = uw_file.readline()
        return line[:line.find('#')].split()

    nw = int(nextline())  # parse header
    uw_file.flush()
    uw = np.loadtxt(uw_file)  # load data - loadtxt okay because not too large

    nw_in_file = uw.shape[0]
    if nw < nw_in_file:
        warn("truncating U(w) frequencies", RuntimeWarning, 2)
        uw = uw[:nw_in_file]
    elif nw > nw_in_file:
        raise ValueError("U(w) file seems to be truncated")

    w = uw[:, 0]
    uw = uw[:, 1:]
    return w, uw

def read_SGWnl_iw_internal(inputfilename,nd,nw,nwout,nk,iw):
    """
    reads unformatted Fortran output of the non-local complex GW self-energy
    all k-points and orbitals for a GIVEN FREQUENCY
    
    nd number of orbitals. in case of presence of GW correction to np orbitals
      this has to be generalized.
    nw number of frequencies per record
    nwout number of frequencies in ctqmc, needs to be smaller/equal than nw
    nk is number of k-points
    iw is frequency index (iw=0,... nwout)
    """
    if nwout > nw:
        raise RuntimeError("not enough frequencies in the GW... increase nw"
                           "in inp.gw to at least %d or decrease Niw in "
                           "Parameters.in" % nwout)

    f = open(inputfilename,'rb')
    f.seek(16*nk*nd*nd*iw)
#        stuff = np.zeros(shape=(nd,nd,nk), dtype=complex,order='F')
    stuff = np.fromfile(f,dtype=complex,count=nd*nd*nk)
    stuff = np.reshape(stuff,(nd,nd,nk))
    f.close()
    
    #have to change order of columns to match with hk[nk,nband,nband,ns]
    stuff = stuff.transpose(2, 0, 1)
    return stuff

#JMT addition START
def read_SGWnl_iw(sgwfile,nd,nk,gwnw,niw,iw):
    """ Reads a frequency and momentum dependent self-energy.
    This is supposed to stem from a GW calculation and be already restricted to
    fully non-local components for the nd block of dmft interacting orbitals
    for all k and a GIVEN FREQUENCY iw.
    Needs to be generalized for correction to np orbitals.
    The GW has to be run with the same Hamiltonian.

    access to direct unformatted fortran output of GW

    sgwfile is the filename to read from
    nd number of bands
    nk number of k-points
    gwnw number of frequencies in the GW
    niw number of frequencies in the ctqmc (needs to be at least as large as in GW)
    iw is the specific frequency that is being read.

    """
    if niw > gwnw:
        raise RuntimeError("not enough frequencies in the GW... increase nw"
                           "in inp.gw to at least %d or or decrease Niw in "
                           "Parameters.in" % niw)

    siwnl = np.empty(shape=(nk,nd,nd), dtype=complex,order='F')
    # this is segmented so that one can parallize over FREQUENCIES...
    siwnl[...] = read_SGWnl_iw_internal(sgwfile,nd,gwnw,niw,nk,iw)
    return siwnl
#JMT addition END




def hdf_dataset(parent, name, value, error):
    """ Creates a dataset and stores some metadata along with it """
    value = np.asarray(value)
    grp = parent.create_group(name)
    if error is not None:
        error = np.asarray(error)
        assert value.shape == error.shape, "Shapes must agree"
        grp.create_dataset("error", data=error)
    grp.create_dataset("value", data=value)

def hdf_qmeta(hdf_qmeta_node,name,axes=None,desc=None):
    if desc is None: desc = ""
    qnode = hdf_qmeta_node.create_group(name)
    qnode.attrs["desc"] = desc
    if axes is not None:
        try:
            qnode.attrs["axes"] = axes
        except TypeError:
            qnode.attrs.create("axes", axes, dtype=h5ustrs)
    return qnode
