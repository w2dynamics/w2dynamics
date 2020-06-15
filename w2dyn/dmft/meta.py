from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from textwrap import dedent
QUANTITIES = {
    "g0iw-full": dict(
        axes=["ineq", "band1", "spin1", "band2", "spin2", "iw"],
        desc="Weiss field including spin/orbital off-diagonal terms"
        ),
    "siw-full": dict(
        axes=["ineq", "band1", "spin1", "band2", "spin2", "iw"],
        desc="Full self-energy in matsubara expansion (with jackknife error)"
        ),
    "siw-trial": dict(
        axes=["ineq", "band1", "spin1", "band2", "spin2", "iw"],
        desc="Full trial self-energy in matsubara expansion"
        ),
    "smom-full": dict(
        axes=["ineq", "band1", "spin1", "band2", "spin2", "moment"],
        desc="Moments of the self-energy including all terms"
        ),
    "g0iw": dict(
        axes=["ineq", "band", "spin", "iw"],
        desc="Weiss field"
        ),
    "leadsiw-full": dict(
        axes=["band1", "spin1", "band2", "spin2", "lead", "iw"],
        desc="Weiss field including spin/orbital off-diagonal terms"
        ),
    "fmom": dict(
        axes=["ineq", "band", "spin", "order"],
        desc="moments of the Hybridisation function in iw"
        ),
    "fiw": dict(
        axes=["ineq", "band", "spin", "iw"],
        desc="hybridisation function in Matsubara frequencies"
        ),
    "ftau": dict(
        axes=["ineq", "band", "spin", "tauf"],
        desc="hybridisation function in imaginary time"
        ),
    "ftau-full": dict(
        axes=["ineq", "band1", "spin1", "band2", "spin2", "tauf"],
        desc="hybridisation function in imaginary time"
        ),
    "muimp-full": dict(
        axes=["ineq", "band1", "spin1", "band2", "spin2"],
        desc="impurity chemical potential"
        ),
    "muimp": dict(
        axes=["ineq", "band", "spin"],
        desc="impurity chemical potential"
        ),
    "gtau-full": dict(
        axes=["ineq", "band1", "spin1", "band2", "spin2", "taubin"],
        desc="impurity Green's function on the imaginary time axis, full offdiagonal"
        ),
    "gtau": dict(
        axes=["ineq", "band", "spin", "taubin"],
        desc="impurity Green's function on the imaginary time axis"
        ),
    "gtau-blocks": dict(
        axes=["ineq", "band", "spin", "taubin", "m"],
        desc="blocking analysis for the impurity Green's function"
        ),
    "gtau-mean-step": dict(
        axes=["ineq", "band", "spin", "iNmeas"],
        desc="impurity Green's function averaged over tau, resolved in time (= Monte carlo steps)"
        ),
    "gtau-mid-step": dict(
        axes=["ineq", "band", "spin", "iNmeas"],
        desc="impurity Green's function averaged over tau between 0.4 * beta and 0.6 * beta, resolved in time (= Monte carlo steps)"
        ),
    "sign-step": dict(
        axes=["ineq", "iNmeas"],
        desc="Sign of a configuration's total weight, resolved in time (= Monte carlo steps)"
        ),
    "gleg": dict(
        axes=["ineq", "band", "spin", "l"],
        desc="impurity Green's function in Legendre expansion"
        ),
    "gleg-full": dict(
        axes=["ineq", "band1", "spin1", "band2", "spin2", "l"],
        desc="full impurity Green's function in Legendre expansion"
        ),
    "giw-full": dict(
        axes=["ineq", "band1", "spin1", "band2", "spin2", "iw"],
        desc="impurity Green's function used as input for Dyson equation"
        ),
    "giw": dict(
        axes=["ineq", "band", "spin", "iw"],
        desc="impurity Green's function used as input for Dyson equation"
        ),
    "giw-cov": dict(
        axes=["ineq", "band", "spin", "pos-iw", "part", "pos-iw", "part"],
        desc="covariance of diagonal impurity Green's function"
        ),
    "giw-meas": dict(
        axes=["ineq", "band", "spin", "iw"],
        desc="Impurity Green's function in Matsubara"
        ),
    "giw-worm": dict(
        axes=["ineq", "iw"],
        desc="Impurity Green's function in Matsubara from worm sampling"
        ),
    "gtau-worm": dict(
        axes=["ineq", "taubin"],
        desc="Impurity Green's function in imaginary time from worm sampling"
        ),
    "gsigmaiw-worm": dict(
        axes=["ineq", "iw"],
        desc="Worm Improved estimators in Matsubara basis"
        ),
    'quddag-worm': dict(
        axes=['ineq', 'iw'],
        desc='Worm Improved estimators in Matsubara basis'
        ),
    "gsigmaiw": dict(
        axes=["ineq", "band", "spin", "iw"],
        desc="Improved estimators in Matsubara basis"
        ),
    "g2iw": dict(
        axes=["ineq", "band", "spin", "iwf-g2", "iwf-g2"],
        desc="Two-particle Green's function in Matsubara basis"
        ),
    "p2iw-worm": dict(
        axes=["ineq", "iwb-p2"],
        desc="Two legged two-particle Green's function ph-convention (1 bosonic)"
        ),
    "p2iwpp-worm": dict(
        axes=["ineq", "iwb-p2"],
        desc="Two legged two-particle Green's function pp-convention (1 bosonic)"
        ),
    "p2tau-worm": dict(
        axes=["ineq", "taubin"],
        desc="Two legged two-particle Green's function in imaginary time"
        ),
    "p2taupp-worm": dict(
        axes=["ineq", "taubin"],
        desc="Two legged two-particle Green's function pp in imaginary time"
        ),
    "p3iw-worm": dict(
        axes=["ineq", "iwf-p3", "iwb-p3"],
        desc="Three legged two-particle Green's function ph-convention (1 fermionic, 1 bosonic)"
        ),
    "p3iwpp-worm": dict(
        axes=["ineq", "iwf-p3", "iwb-p3"],
        desc="Three legged two-particle Green's function pp-convention (1 fermionic, 1 bosonic)"
        ),
    "g4iw-worm": dict(
        axes=["ineq", "iwf-g4", "iwf-g4", "iwb-g4"],
        desc=dedent("""\
        Two-particle Green's function in particle-hole Matsubara frequencies from worm sampling
        There are two conventions:
        0: (v+w) tau_1 - v tau_2 + v' tau_3 - (v'+w) tau_4
        1: v tau_1 - (v-w) tau_2 + (v'-w) tau_3 - v' tau_4""")
        ),
    "g4iwpp-worm": dict(
        axes=["ineq", "iwf-g4", "iwf-g4", "iwb-g4"],
        desc=dedent("""\
        Two-particle Green's function in particle-particle Matsubara frequencies from worm sampling
        Convention: v tau_1 - (w-v') tau_2 + (w-v) tau_3 - v' tau_4""")
        ),
    "h4iw-worm": dict(
        axes=["ineq","wf-g4","iwf-g4","iwb-g4"],
        desc="Two-particle improved estimator in Matsubara basis from worm sampling"
        ),
    "qqiw-worm": dict(
        axes=["ineq", "iw"],
        desc="Worm Symmetric Improved 1P estimator in Matsubara basis"
        ),
    "qqtau-worm": dict(
        axes=["ineq", "taubin"],
        desc="Worm Symmetric Improved 1P estimator in imaginary time"
        ),
    "nqqdag-worm": dict(
        axes=["ineq", "iwf-p3", "iwb-p3"],
        desc="Part of Symmetric Improved 2P estimator in Matsubara frequencies"
        ),
    "qqdd-worm": dict(
        axes=["ineq", "iwf-p3", "iwb-p3"],
        desc="Part of Symmetric Improved 2P estimator in Matsubara frequencies"
        ),
    "qqqq-worm": dict(
        axes=["ineq", "iwf-g4","iwf-g4","iwb-g4"],
        desc="Worm Symmetric Improved 2P estimators in Matsubara basis"
        ),
    "ucaca-worm": dict(
        axes=["ineq", "iwb-p2"],
        desc="Worm Symmetric Improved 2P estimators in Matsubara basis"
        ),
    "ucacatau-worm": dict(
        axes=["ineq", "taubin"],
        desc="Worm Symmetric Improved 2P estimators in Matsubara basis"
        ),
    "uccaa-worm": dict(
        axes=["ineq", "iwb-p2"],
        desc="Worm Symmetric Improved 2P estimators in Matsubara basis"
        ),
    "uccaatau-worm": dict(
        axes=["ineq", "taubin"],
        desc="Worm Symmetric Improved 2P estimators in Matsubara basis"
        ),
    "g4tau": dict(
        axes=["ineq", "band1", "spin1", "band2", "spin2", "tau-g4", "tau-g4", "tau-g4"],
        desc="Two-particle Green's function in tau12, tau34, tau14"
        ),
    "g4leg": dict(
        axes=["ineq", "band1", "spin1", "band2", "spin2", "l", "l", "iwb-g4"],
        desc="Two-particle Green's function in Legendre/Matsubara basis"
        ),
    "g4iw": dict(
        axes=["ineq", "band1", "spin1", "band2", "spin2", "iwf-g4", "iwf-g4", "iwb-g4"],
        desc="Two-particle Green's function in Matsubara basis (particle-hole channel)"
        ),
    "g4iw-pp": dict(
        axes=["ineq", "band1", "spin1", "band2", "spin2", "iwf-g4", "iwf-g4", "iwb-g4"],
        desc="Two-particle Green's function in Matsubara basis (particle-particle channel)"
        ),
    "occ": dict(
        axes=["ineq", "band1", "spin1", "band2", "spin2"],
        desc="Occupancies"
        ),
    "rho2": dict(
        axes=["ineq", "band1", "spin1", "band2", "spin2", "band3", "spin3", "band4", "spin4"],
        desc="< c^+ c^+ c c >"
        ),
    "rho1": dict(
        axes=["ineq", "band1", "spin1", "band2", "spin2"],
        desc="< c^+ c >"
        ),
    "smom": dict(
        axes=["ineq", "band", "spin", "order"],
        desc="Moments of the self-energy (in Matsubara)"
        ),
    "ntau-n0": dict(
        axes=["ineq", "band1", "spin1", "band2", "spin2", "tausus"],
        desc="density-density correlation function in tau bins"
        ),
    "densitymatrix": dict(
        axes=["ineq", "state", "state"],
        desc="Density matrix in occupation number basis measured at beta/2"
        ),
    "expresdensitymatrix": dict(
        axes=["ineq", "band" , "spin", "traceorder", "state", "state"],
        desc="density matrix resolved by expansion order"
        ),
    "time-qmc": dict(
        axes=["ineq"],
        desc="Mean CPU seconds spent in CT-QMC per process"
        ),
    "time": dict(
        desc="Total elapsed walltime of the run"
        ),
    "time-warmup": dict(
        axes=["ineq"],
        desc="CPU time used for the QMC warmup"
        ),
    "time-simulation": dict(
        axes=["ineq"],
        desc="CPU time used for the QMC simulation"
        ),
    "time-sampling": dict(
        axes=["ineq"],
        desc="CPU time used for the QMC simulation excluding measurements"
        ),
    "sign": dict(
        axes=["ineq"],
        desc="Mean sign"
        ),
    "single-occ": dict(
        axes=["ineq", "band1", "spin1"],
        desc="Single occupancies <n_i>"
        ),
    "double-occ": dict(
        axes=["ineq", "band1", "spin1", "band2", "spin2"],
        desc="Double occupancies <n_i n_j>"
        ),
    "accept-ins": dict(
        axes=["ineq", "space"],
        desc="Acceptance rate for insertions into local trace"
        ),
    "accept-ins4": dict(
        axes=["ineq"],
        desc="Acceptance rate for insertions of 4 operators into local trace"
        ),
    "accept-rem": dict(
        axes=["ineq", "space"],
        desc="Acceptance rate for removals of local trace operators"
        ),
    "accept-rem4": dict(
        axes=["ineq"],
        desc="Acceptance rate for removal of 4 operators in the local trace"
        ),
    "accept-pair-tau": dict(
        axes=["ineq", "band", "spin", "taubin"],
        desc="Total number of accepted pair insertions and removals in tau bins"
        ),
    "accept-glob": dict(
        axes=["ineq"],
        desc="Acceptance rate for global symmetry moves"
        ),
    "accept-shift": dict(
        axes=["ineq"],
        desc="Acceptance rate for tau-shift moves"
        ),
    "accept-flavourchange": dict(
        axes=["ineq"],
        desc="Acceptance rate for flavourchange moves"
        ),
    "accept-worm-ins": dict(
        axes=["ineq", "space"],
        desc="Acceptance rate for worm insertions into local trace"
        ),
    "accept-worm-rem": dict(
        axes=["ineq", "space"],
        desc="Acceptance rate for worm removals from local trace"
        ),
    "accept-worm-rep": dict(
        axes=["ineq", "space"],
        desc="Acceptance rate for worm replacement moves"
        ),
    "steps-worm-partition": dict(
        axes=["ineq", "space"],
        desc="Time spent in worm space and partition space"
        ),
    "worm-eta": dict(
        axes=["ineq", "space"],
        desc="Worm balancing factor"
        ),
    "hist-seg": dict(
        axes=["ineq", "order"],
        desc="number of empty flavours in trace"
        ),
    "hist": dict(
        axes=["ineq", "band", "spin", "traceorder"],
        desc="Histogram for expansion orders"
        ),
    "ssts-states": dict(
        axes=["ineq", "state"],
        desc="Assignment of states to superstates by state index"
        ),
    "energies-eigenstates": dict(
        axes=["ineq", "state"],
        desc="Energies of eigenstates by state index"
        ),
    "occbasis-mapping": dict(
        axes=["ineq", "state", "band", "spin"],
        desc="Occupation number eigenvalues by state index as used in densitymatrix"
        ),
    "hist-sst": dict(
        axes=["ineq", "superstate"],
        desc="Histogram for outer superstates"
        ),
    "hist-state": dict(
        axes=["ineq", "state"],
        desc="Histogram for outer states"
        ),
    "sign-sst": dict(
        axes=["ineq", "superstate"],
        desc="Summed total sign per outer superstate"
        ),
    "sign-state": dict(
        axes=["ineq", "state"],
        desc="Summed total sign per outer state"
        ),
    "contrib-sst": dict(
        axes=["ineq", "superstate"],
        desc="Trace contribution per outer superstate"
        ),
    "contrib-state": dict(
        axes=["ineq", "state"],
        desc="Trace contribution per outer state"
        ),
    "lhist": dict(
        axes=["ineq", "lsteps"],
        desc="Histogram for Lanczos steps in the time evolution"
        ),
    "rhist": dict(
        axes=["ineq", "rsteps"],
        desc="Histogram for additional Lanczos steps needed for reversible trace"
        ),
    "time-giw": dict(
        axes=["ineq", "traceorder"],
        desc="Time spent on measuring G(iw)"
        ),
    "time-g4iw-ft": dict(
        axes=["ineq", "traceorder"],
        desc="Time spent on fourier-transforming G4(iw,iw,iW)"
        ),
    "time-g4iw-add": dict(
        axes=["ineq", "traceorder"],
        desc="Time spent on constructing and adding G4(iw,iw,iW)"
        ),
    # === LATTICES ===
    "nneq": dict(
        desc="number of non-equivalent atoms in unit cell"
        ),
    "neq2at": dict(
        axes=["inequiv-index"],
        desc="atom index first occurrence of inequivalent atom"
        ),
    "at2neq": dict(
        axes=["atom"],
        desc="corresponding inequivalent index to atom index"
        ),
    "mu": dict(
        desc="chemical potential (lattice model)"
        ),
    "lda-mu": dict(
        desc="chemical potential (non-interacting)"
        ),
    "siw": dict(
        axes=["ineq", "band", "spin", "iw"],
        desc="Band-spin-diagonal self-energy in matsubara expansion"
        ),
    "siw-cov": dict(
        axes=["ineq", "band", "spin", "pos-iw", "part", "pos-iw", "part"],
        desc="covariance of diagonal self-energy in matsubara expansion"
        ),
    "sigma-hartree": dict(
        axes=["band1", "spin1", "band2", "spin2"],
        desc="Hartree self-energy"
        ),
    "glocold": dict(
        axes=["ineq", "band", "spin", "iw"],
        desc="local Green's function in Matsubara (old self-energy)"
        ),
    "glocold-lattice": dict(
        axes=["lda-band", "spin", "iw"],
        desc="local Green's function in Matsubara (old self-energy), diagonal part for all lda-bands"
        ),
    "glocnew": dict(
        axes=["ineq", "band", "spin", "iw"],
        desc="local Green's function in Matsubara (new self-energy)"
        ),
    "glocnew-lattice": dict(
        axes=["lda-band", "spin", "iw"],
        desc="local Green's function in Matsubara (old self-energy), diagonal part for all lda-bands"
        ),
    "dc": dict(
        axes=["ineq", "band", "spin"],
        desc="double counting correction"
        ),
    "dc-latt": dict(
        axes=["band1", "spin1", "band2", "spin2"],
        desc="double counting correction on the lattice"
        ),
    "lda-dens": dict(
        axes=["lda-band", "spin"],
        desc="LDA orbital-resolved densities"
        ),
    "h-mean": dict(
        axes=["lda-band1", "spin1", "lda-band2", "spin2"],
        desc="Mean over k-points of Hamiltonian H(k)"
        ),
    "h-mom2": dict(
        axes=["lda-band1", "spin1", "lda-band2", "spin2"],
        desc="Second moment of the density of states"
        ),
    "hk": dict(
        axes=["kpoint", "lda-band1", "spin1", "lda-band2", "spin2"],
        desc="Full Hamiltonian H(k)"
        ),
    "hk-mean": dict(
        axes=["ineq", "band", "spin"],
        desc="Mean of diagonals of Hamiltonian H(k)"
        ),
    "lda-dos": dict(
        axes=["lda-band", "spin", "w-dos"],
        desc="LDA density of states"
        ),
    "gdensold": dict(
        axes=["lda-band1", "spin1", "lda-band2", "spin2"],
        desc="densities with old self-energy before adjustment of mu"
        ),
    "gdensnew": dict(
        axes=["lda-band1", "spin1", "lda-band2", "spin2"],
        desc="densities with new self-energy after adjustment of mu"
        ),
    "ubar": dict(
        axes=["shell1","shell2"],
        desc="shell-averaged U values"
        ),
    "jbar": dict(
        axes=["shell1","shell2"],
        desc="shell-averaged J values"
        ),
}
