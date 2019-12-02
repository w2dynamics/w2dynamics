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
    \__\__|X/____|____//_/  |_| \__|   P Gunacker, A Kowalski, N Parragh, F Goth,
         |XXXXXXXX|                              K Held, G Sangiovanni
                                       Version %s, %s
"""

CODE_VERSION = 1, 0, "0"
CODE_VERSION_STRING = ".".join(map(str,CODE_VERSION))
CODE_DATE = "July 2018"
OUTPUT_VERSION = 2, 2
