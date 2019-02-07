""" 
Package providing MaxEnt stuff and some fortran specific files

 MakeFile
Compiles MaximumEntropy.f90 using f2py
Fortran compiler is chosen by setting the FC flag (i.e. FC=gfortran)

MaxEnt.py: 
Interface to parse hdf5 files from w2dynamics and
simple dat files for the maximum entropy solver.

The file format for the dat files historically includes two
comment lines, which will be ignored (can be changed by setting the
offset parameter). The first column includes the tau/iw values. The
second column includes the data points. The third colum includes
the errorbars. In case of giw files the imaginary part is
included in the third column and the error in the fourth colum.

More help is displayed by typing "./MaxEnt.py --help"

Solver.py:
Base class for wrapping the MaxEnt solver with some sane default
parameters. There are two derived classes (Green2Aw,Chi2Aw) to
specify some further default parameters.

MaximumEntropy.F90:
Maximum Entropy code. In any case this part is "over-documented",
so details can be found in the code directly.

MersenneTwister.F90:
Just a note of caution: this random number generator works on a closed
interval [0,1] not on an open interval [0,1) as one would assume.

"""
