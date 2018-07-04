#!/usr/bin/env python
"""
Calculates local spin suceptibility from an HDF5 file

As for now, output is written to sdtout and needs to be piped to file.  The 
output format contains the data fit for every bosonic frequency in a block
seperated by an empty line, followed by a block of the form udd, omega, 
chi_loc in order to plot omega chi_loc over omega behavior.

As for now chi_loc needs to be devided by beta, and the freuqencies need to
be shifted around zero in the plot.

Author: Patrik Gunacker
Author: Markus Wallerberger
"""
import numpy as np
import h5py as hdf5
import sys
import re
import scipy.optimize
import itertools as itt
import optparse

import auxiliaries.postprocessing as pp
import auxiliaries.quantities as qttys

__version__ = 1, 0

#hdf5 file name to enter by user, otherwise exit
parser = optparse.OptionParser(usage = "%prog [OPTIONS] <hdf-file>",
                               version = ".".join(map(str,__version__)))
parser.add_option("-a", type="int", dest="order", default=1,
                  help="order of analytic approximation (0 or [1])")
parser.add_option("-f", type="int", dest="fit", default=1,
		  help="type of fit")
parser.add_option("-n", type="string", dest="name",
                  help="filename infix")

options, args = parser.parse_args()
if len(args) != 1:
    parser.error("expected HDF5 file as single argument")

filename = args[0]
hf = hdf5.File(filename, "r")
file_version = tuple(hf.attrs["outfile-version"])

if file_version >= (2,0):
    beta = hf[".config"].attrs["general.beta"]
else:
    beta = float(next(e[:e.find("#")].split("=")[1]
                      for e in hf["configfile"].value 
                      if e.startswith("beta")))

print >> sys.stderr, "beta =", beta

#extracting one and two particle green function and errors from hdf5
axes_node = hf[(None, "axes",".axes")[file_version[0]]]
iw4f = axes_node["iwf-g4"].value
iw4b = axes_node["iwb-g4"].value

iw4f_pos = iw4f.size//2
iw4b_zero = iw4b.size//2

iter = hf["stat-last"]
qtty_nodes = itt.izip(qttys.ineq_quantity(iter, "giw"), 
                      qttys.ineq_quantity(iter, "g4iw"),
                      qttys.ineq_quantity(iter, "occ"))

if options.order == 1:
    # exactly sum the GG crossing term, fit 1/x^2
    gg_exact = True
elif options.order == 0:
    # just do the sum plus the fitting
    gg_exact = False
else:
    parser.error("order must currently be either 0 or 1")

if options.fit == 0:
    fit_func = None
elif options.fit == 1:
    fit_func = lambda x, a, b, c: a + b/x + c/(x**2)
elif options.fit == 2:
    fit_func = lambda x, a, b: a + b/(x**2)
else:
    parser.error("fit order must currently be either 1 or 2")

suffix = "dat"
if options.name: suffix = options.name + "." + suffix
fitf = file("chi_asymp.%s" % suffix, "w")
chi0f = file("chi0.%s" % suffix, "w")
fit0f = file("chi0_asymp.%s" % suffix, "w")
tot0f = file("chi0_tot.%s" % suffix, "w")
chif = file("chi.%s" % suffix, "w")
totf = file("chi_tot.%s" % suffix, "w")
diaf = file("chi_diag.%s" % suffix, "w")
offf = file("chi_off.%s" % suffix, "w")
momf = file("moments.%s" % suffix, "w")

print >> chi0f, ("#%2s %3s %9s %12s %12s" % 
                 ("at", "bi", "iwbos", "ReX0_loc","ImX0_loc")),
print >> fitf, ("#%2s %3s %3s %9s %9s %12s %12s" % 
                ("at", "bi","bj","iwbos","ivmax","ReX_loc","ImX_loc")),
print >> chif, (("#%2s %3s %3s %9s"+4*" %12s"+"   further_fit_parameters") %
                ("at","bi","bj","iwbos","ReX_loc","err(ReX)","ImX_loc","err(ImX)")),
print >> totf, ("#%2s %9s %12s %12s" % 
                ("at", "iwbos", "ReX_loc","ImX_loc"))
print >> tot0f, ("#%2s %9s %12s %12s" % 
                ("at", "iwbos", "ReX_loc","ImX_loc"))
print >> diaf, ("#%2s %9s %12s %12s" % 
                ("at", "iwbos", "ReX_loc","ImX_loc"))
print >> offf, ("#%2s %9s %12s %12s" % 
                ("at", "iwbos", "ReX_loc","ImX_loc"))
print >> momf, ("#%2s %3s %3s %12s %12s %12s" % 
                ("at", "bi", "bj", "Re(mu)","Im(mu)","mu(QMC)"))

for iineq, ((giw, giw_err), (g4iw, g4iw_err), 
            (occ, occ_err)) in enumerate(qtty_nodes):

    # chi0 = -GG, therefore the minus (that is why in get_g4iw_conn_ph,
    # there is a plus when subtracting chi0). The bubble is a nice scale
    # to compare stuff to.
    chi0 = -pp.get_ggcross_ph(giw, iw4b.size).sum(2) # b,s,(iv),iw
    chi0 = np.sum(chi0, axis=1)   # paramagnetic
    for (bi, iiw), val in np.ndenumerate(chi0/beta):
        if not iiw:
            print >> chi0f
            print >> chi0f
        print >> chi0f, ("%3i %3i %9.3f %12.6g %12.6g" % 
                         (iineq, bi, iw4b[iiw], val.real, val.imag))
        
    chi0x = -pp.get_ggcross_ph(giw, iw4b.size)
    iwf_pos = chi0x.shape[-2]//2
    sl_neg = slice(iwf_pos-1, None, -1)
    sl_pos = slice(iwf_pos, None, 1)
    chi0x = (chi0x[...,sl_neg,:] + chi0x[...,sl_pos,:])
    chi0x = chi0x.cumsum(2).transpose(0,1,3,2).sum(1)
     
    for (bi, iiw, iv), val in np.ndenumerate(chi0x/beta):
        if not iv:
            print >> fit0f
            print >> fit0f
        print >> fit0f, ("%3i %3i %9.3f %5i %12.6g %12.6g" % 
                         (iineq, bi, iw4b[iiw], iv, val.real, val.imag))

    # subtract first order
    if gg_exact:
        chi = pp.get_g4iw_conn_ph(giw, g4iw, giw_err, g4iw_err)
    else:
        chi = pp.get_chi_ph(giw, g4iw, giw_err, g4iw_err)
    
    # do superposition in the spin channel
    chi = chi[:,1,:,1] - chi[:,1,:,0] - chi[:,0,:,1] + chi[:,0,:,0]
    
    # add negative frequencies to positive ones in order to find a good
    # frequency box (0, ... v,)
    sl_neg = slice(iw4f_pos-1, None, -1)
    sl_pos = slice(iw4f_pos, None, 1)
    chi = (chi[...,sl_neg,sl_neg,:] + chi[...,sl_neg,sl_pos,:] +
           chi[...,sl_pos,sl_neg,:] + chi[...,sl_pos,sl_pos,:])
    
    # This does the cumulative sum in fermionic dimensions, allowing for a
    # fit in v, v'. Since we only do a fit in the combined frequency right
    # now, we take the diagonal, which corresponds to an expanding quadratic
    # frequency box. chi_loc is a (band, band, iW, ivmax) array
    chi = chi.cumsum(-3).cumsum(-2)
    chi_loc = chi.diagonal(0, -3, -2)
    #chi_loc = chi[...,-1,:].transpose(0,1,3,2)
    ivmax = iw4f[iw4f_pos:]

    for (bi, bj, iiw, iiv), val in np.ndenumerate(chi_loc):
        if not iiv: 
            print >> fitf
            print >> fitf
        print >> fitf, ("%3i %3i %3i %9.3f %9.3f %12.6g %12.6g" % 
                        (iineq,bi,bj,iw4b[iiw],ivmax[iiv], val.real,val.imag))
        
    # begin the fit here
    chi_fit = np.zeros_like(chi_loc[...,0])
    ivbegin = np.round(iw4f_pos/2)
    
    for fit_idx in np.ndindex(*chi_fit.shape):
        x = ivmax[ivbegin:] 
        y = chi_loc[fit_idx][ivbegin:]
	if fit_func is not None:
	    rfit, rcov = scipy.optimize.curve_fit(fit_func, x, y.real)
            ifit, icov = scipy.optimize.curve_fit(fit_func, x, y.imag)
        else:
            rfit = np.asanyarray(y[-1].real).reshape(1)
            ifit = np.asanyarray(y[-1].imag).reshape(1)
            rcov, icov = np.zeros((2,1,1))
            
        rerr = np.sqrt(rcov.diagonal())
        ierr = np.sqrt(icov.diagonal())
        
        # the band indices have to be equal in order for the crossing term to
        # contribute.
        if gg_exact and fit_idx[0] == fit_idx[1]:
            rfit += chi0[fit_idx[1:]].real
            ifit += chi0[fit_idx[1:]].imag
            
        rfit /= beta
        ifit /= beta
        rerr /= beta
        ierr /= beta

        chi_fit[fit_idx] = rfit[0] + 1j*ifit[0]
        
        if not fit_idx[-1]: 
            print >> chif
            print >> chif
        print >> chif, "%3i %3i %3i %9.3f" % (iineq,fit_idx[0],fit_idx[1],
                                              iw4b[fit_idx[2]]),
        print >> chif, " ".join(map(lambda x: "%12.6g" % x,
                                    np.transpose((rfit, rerr, ifit, ierr)).flat))

    # PADE
    np.savez("chi.%s.npz" % options.name, chi=chi_fit, chi0=chi0, iw=iw4b)

    # computing diagonal, off-diagonal and total part
    chi_tot = chi_fit.sum(1).sum(0)
    chi_dia = chi_fit.diagonal(0,0,1).sum(-1)
    chi_off = chi_tot - chi_dia
    chi0_tot = chi0.sum(0)

    for iiw, iw in enumerate(iw4b):
        fmt = "%3i %9.3f %12.6g %12.6g"
        print >> totf, fmt % (iineq, iw, chi_tot[iiw].real, chi_tot[iiw].imag)
        print >> diaf, fmt % (iineq, iw, chi_dia[iiw].real, chi_dia[iiw].imag)
        print >> offf, fmt % (iineq, iw, chi_off[iiw].real, chi_off[iiw].imag)
        print >> tot0f, fmt % (iineq, iw, chi0_tot[iiw].real, chi0_tot[iiw].imag)

    moments = (1./beta) * chi_fit.sum(-1)
    moments_qmc = pp.moment(occ)

    for (bi, bj), val in np.ndenumerate(moments):
        print >> momf, ("%3i %3i %3i %12.6g %12.6g %12.6g" % 
                        (iineq, bi, bj, val.real, val.imag, moments_qmc[bi,bj]))
    
print >> sys.stderr, "Success."

