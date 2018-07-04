""" 
EXPLANATION:
The presence of retarded density-density interactions, can be transformed into 
(1) a shift of the density-density potentials and 
(2) a modified weight of all fermion operators.
In the present file dynamicalU.py the shift of all potentials 
is implemented, i.e. the umatrix is updated.

The full derivation behind this implementation (currently limited access) can be found here:
http://www.groupmatter.org/Progress/CondMat/index.php?title=DMFT_U%28w%29:_Multiband_Hubbard_Stratonovich

___________________________________________________________________________________________________________
 FIXME (optional): The multi-band version for retarded interactions should only be used 
 when the Hamiltonian option is set to ReadUmatrix. This is the only way to allow 
 multiple U_i and V_ij for a single atom. Check whether this is true and inform the user.
 
 FIXME (must): If only part of the bands have retarded interactions, 
 the present code will crash as it expects retardation for all bands.
____________________________________________________________________________________________________________
"""


import numpy as np
import os

class Replace_Retarded_Interaction_by_a_Shift_of_instantaneous_Potentials:

# ================================================================================================================================================================================================
# <== CLASS PREPARATION - prepares all required potential shifts (_*_) and makes them accessible (self) accross the entire class
# Replace_Retarded_Interaction_by_a_Shift_of_instantaneous_Potentials.
# ================================================================================================================================================================================================

  def __init__(self, beta, umatrix, matsubara_screening_switch):
    print '********* U(w) MODULE IS NOT AVAILABLE IN PRESENT W2DYNAMICS VERSION *********'
    print '...exiting'''
    exit()
