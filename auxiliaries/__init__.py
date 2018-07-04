from warnings import warn

try:
    import CTQMC
except:
    warn("CTQMC.so binary module not found. Run make in ctqmcf")


BANNER = r"""
           ______
         _/XXXXXX\ ___ __  ____   /|   W2DYNAMICS - Wuerzburg/Wien strong
  |\    | |X/  \X| __ \\ \/ /  \ | |                coupling impurity solver
  \ \ |\| |____/X| | \ \\  /| \ \| |
   \ \\ \ /XXXXXX| |_/ // / | |\ \ |   Authors: M Wallerberger, A Hausoel,
    \__\__|X/____|____//_/  |_| \__|       P Gunacker, N Parragh, G Sangiovanni
         |XXXXXXXX|                    Version %s, %s
"""

CODE_VERSION = 1, 0, "alpha1"
CODE_VERSION_STRING = ".".join(map(str,CODE_VERSION))
CODE_DATE = "May 2016"
OUTPUT_VERSION = 2, 2
