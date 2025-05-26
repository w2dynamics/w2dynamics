from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from warnings import warn

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

CODE_VERSION = 1, 1, "6"
CODE_VERSION_STRING = ".".join(map(str,CODE_VERSION))
CODE_DATE = "May 2025"
OUTPUT_VERSION = 2, 2
