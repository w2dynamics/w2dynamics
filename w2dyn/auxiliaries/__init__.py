from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from warnings import warn

try:
    from . import CTQMC
except:
    warn("CTQMC binary module failed to load. Please make sure to correctly build the compiled components of w2dynamics using the cmake build system.")


BANNER = r"""
           ______
         _/XXXXXX\ ___ __  ____   /|   W2DYNAMICS - Wuerzburg/Wien strong
  |\    | |X/  \X| __ \\ \/ /  \ | |                coupling impurity solver
  \ \ |\| |____/X| | \ \\  /| \ \| |
   \ \\ \ /XXXXXX| |_/ // / | |\ \ |   Authors: M Wallerberger, A Hausoel,
    \__\__|X/____|____//_/  |_| \__|   P Gunacker, A Kowalski, N Parragh, F Goth,
         |XXXXXXXX|                              K Held, G Sangiovanni
                                       Version %s, %s
"""

CODE_VERSION = 1, 1, "1"
CODE_VERSION_STRING = ".".join(map(str,CODE_VERSION))
CODE_DATE = "February 2021"
OUTPUT_VERSION = 2, 2
