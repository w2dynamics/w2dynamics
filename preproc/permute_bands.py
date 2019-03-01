#!/usr/bin/env python
# example:
#  ./permute_bands.py Hk.dat $(seq 1 5) $(seq 11 28) $(seq 6 10) $(seq 29 46)

import numpy as np
import w2dyn.auxiliaries.input as inp
import sys

args = sys.argv[1:]
try:
    hk_file_name = args.pop(0)
except IndexError:
    print >> sys.stderr, "Usage:"
    print >> sys.stderr, "   ./permute_bands.py HK_FILE BAND1 BAND2 ... BANDn"
    sys.exit(1)
    

Hk, kpoints = inp.read_hamiltonian(file(hk_file_name, "r"))
nbands = Hk.shape[-1]

tp = np.array(map(int, args)) - 1
if tp.size != nbands:
    raise RuntimeError("not enough bands specified")
if np.any(np.sort(tp) != np.arange(nbands)):
    raise RuntimeError("did not specify permutation")

Hk_transposed = Hk[:, tp[:, None], tp[None, :]]
k_format = "%12.5e " * kpoints.shape[-1]
ham_format = "%12.5e %12.5e  " * Hk.shape[-1]

print "%d %d %d  # Generated-by permute_bands.py" % Hk.shape
for ik, Honek in enumerate(Hk_transposed):
    print k_format % tuple(kpoints[ik])
    np.savetxt(sys.stdout, Honek, fmt=ham_format)
