w2dynamics - Wien/Wuerzburg strong coupling solver
==================================================

w2dynamics is a hybridization-expansion continuous-time quantum Monte Carlo
package, developed jointly in Wien and Würzburg. 

**In any published papers arising from the use of w2dynamics, please cite:**

   M. Wallerberger, A. Hausoel, P. Gunacker, A. Kowalski, N. Paragh, F. Goth, K. Held, and G. Sangiovanni,  
   Comput. Phys. Commun. 235, 2 (2019)  
   <https://www.sciencedirect.com/science/article/pii/S0010465518303217>  
   arXiv:1801.10209 <https://arxiv.org/abs/1801.10209>  
   When using additional codes in  conjunction with w2dynamics, do not forget  
   to give credit to them as well.  

w2dynamics contains:

 - a multi-orbital quantum impurity solver for the Anderson impurity model
 - dynamical mean field theory self-consistency loop,
 - a maximum-entropy analytic continuation, as well as
 - coupling to density functional theory.

The w2dynamics package allows for calculating one- and two-particle quantities;
it includes worm and further novel sampling schemes. Details about its download,
installation, functioning and the relevant parameters are provided.

Maintainers and principal authors:

  - Markus Wallerberger
  - Andreas Hausoel
  - Patrik Gunacker
  - Alexander Kowalski
  - Nicolaus Parragh
  - Florian Goth
  - Karsten Held
  - Giorgio Sangiovanni


Installation
------------

Requirements:

  - Python (>= 2.4)
  - Fortran 90 compiler
  - C++11 compiler
  - cmake (>= 2.8.5)
  - HDF5 (preferably installed via your systems package manager)

Further dependencies (automatically installed if not found):

  - Python packages: numpy >= 1.4, scipy >= 0.10, h5py, mpi4py, configobj
  - NFFT3

To get the code use `git`:

    $ git clone git@github.com:w2dynamics/w2dynamics

To build the code follow the `cmake` build model:

    $ mkdir build
    $ cd build
    $ cmake .. [FURTHER_FLAGS_GO_HERE]
    $ make

To run the unit tests (optional), run the following in the build directory:

    $ make test

To install the code (optional), run the following in the build directory:

    $ cmake .. -DCMAKE_INSTALL_PREFIX=/path/to/prefix
    $ make install


Running the code
----------------

First, create a new run directory, prepare an input file, which is usually
named `Parameters.in`.  You can use the input files for the benchmarkus in
the benchmark repository as templates.

Next, run `DMFT.py` to use the self-consistency loop or `cthyb` if you do a
single-shot calculation. (If you have installed the code (see above), both
executable should be placed in your path.)

    $ man DMFT.py
    $ DMFT.py [Parameters.in]

The code will produce a file with a name like `RunIdentifier-Timestamp.hdf5`.
It is an archive of all quantities written by w2dynamics.  You can navigate this
file using any hdf5-compatible analysis tool, such as jupyter, matlab, etc.
For your convenience, we have also included the tool `hgrep`, which allows
quick analysis of the data:

    $ man hgrep
    $ hgrep -p latest siw 1 1 1 1


Files and directories
---------------------

  - `w2dyn/auxiliaries/`: auxiliary python routines (in/output, config files, etc.)
  - `w2dyn/dmft/`: python package for DMFT self-consistency loop
  - `w2dyn/maxent/`: Python wrapper for maximum entropy analytic continuation
  - `clusters/`: template submission scripts for different clusters
  - `cmake/`: cmake custom modules
  - `doc/`: documentation files
  - `documentation`: doxygen configuration
  - `Postproc/`: postprocessing scripts
  - `preproc/`: preprocessing scripts
  - `src/`: compiled modules loaded from python
    - `ctqmc_fortran`: Fortran 90 continuous-time quantum Monte Carlo solver
    - `maxent`: maximum entropy analytic continuation solver
    - `mtrng`: Mersenne twister pseudorandom number generator
  - `tests/`: minimal run-through tests
  - `testsuite/`: unit tests for the code

  - `cfg_converter.py`: small script converting old-style config files
  - `completions.sh`: file for bash completions
  - `cprun`: convenience script copying input files to different directory
  - `DMFT.py`: main entry point for DMFT self-consistency loop
  - `hgrep`: utility for extracting data from HDF5 file
  - `Maxent.py`: main entry point for maximum entropy code
  - `run_tests.sh`: run the run-through tests
  - `setup.py`: Python installation script

Citation
--------
If you use the w2dynamics package, please mention the following in the acknowledments:

The QMC simulations were carried out with the w2dynamics package available at https://github.com/w2dynamics/w2dynamics .


[Disclaimer](https://www.uni-wuerzburg.de/sonstiges/impressum/)

[Github Privacy Statement](https://help.github.com/articles/github-privacy-statement/)
