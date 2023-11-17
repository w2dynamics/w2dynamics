w2dynamics - Wien/Wuerzburg strong coupling solver
==================================================
w2dynamics is a hybridization-expansion continuous-time quantum Monte Carlo
package, developed jointly in Wien and Würzburg. 

**In any published papers arising from the use of w2dynamics, please cite:**

   M. Wallerberger, A. Hausoel, P. Gunacker, A. Kowalski, N. Parragh, F. Goth, K. Held, and G. Sangiovanni,  
   Comput. Phys. Commun. 235, 2 (2019)  
   <https://www.sciencedirect.com/science/article/pii/S0010465518303217>  
   arXiv:1801.10209 <https://arxiv.org/abs/1801.10209>  
   When using additional codes in  conjunction with w2dynamics, do not forget to give credit to them as well.  

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

  - [Python](https://www.python.org/) (>= 2.6)
  - Fortran 90 compiler, e.g. [GCC](https://gcc.gnu.org/)
  - C++11 compiler, e.g. [GCC](https://gcc.gnu.org/)
  - [CMake](https://cmake.org/) (>= 3.18)
  - [BLAS](https://www.netlib.org/blas/), ideally an optimized version as provided e.g. by [OpenBLAS](https://www.openblas.net/) or [Intel MKL](https://software.intel.com/mkl)
  - [LAPACK](https://www.netlib.org/lapack/), ideally an optimized version as provided e.g. by [OpenBLAS](https://www.openblas.net/) or [Intel MKL](https://software.intel.com/mkl)

Further dependencies (automatically installed if not found):

  - Python packages: [numpy](https://pypi.org/project/numpy/) >= 1.10, [scipy](https://pypi.org/project/scipy/) >= 0.10, [h5py](https://pypi.org/project/h5py/), [mpi4py](https://pypi.org/project/mpi4py/), [configobj](https://pypi.org/project/configobj/)
  - [NFFT3](https://www-user.tu-chemnitz.de/~potts/nfft/) (for automatic building, its dependency [FFTW3](http://www.fftw.org/) is required)
  - [HDF5](https://www.hdfgroup.org/solutions/hdf5) (as a dependency of h5py)

To get the code use `git`:

    $ git clone https://github.com/w2dynamics/w2dynamics.git

To build the code follow the `cmake` build model:

    $ mkdir build
    $ cd build
    $ cmake .. [FURTHER_FLAGS_GO_HERE]
    $ make

If CMake fails to find some native dependency automatically, command
line arguments of the form `-D<PACKAGE>_ROOT=/path/to/package/` can
usually be used to explicitly specify installation prefixes that
should be searched, for example `-DNFFT_ROOT=/path/to/nfft/` for
NFFT. For the Python packages, you should try to ensure that they can
be imported in the interpreter used by CMake. CMake can be instructed
to use a specific interpreter executable using
`-DPYTHON_EXECUTABLE=/path/to/python`. Consider especially that if a
Python 3 installation is found automatically, it will be preferred to
a Python 2 installation.

To run the unit tests (optional), run the following in the build directory:

    $ make test

To install the code (optional), run the following in the build directory:

    $ cmake .. -DCMAKE_INSTALL_PREFIX=/path/to/prefix
    $ make install

For example installation instructions for some operating systems and
computing clusters, you can also visit our wiki page on
[installation](https://github.com/w2dynamics/W2Dynamics/wiki/Installation).


Running the code
----------------

First, prepare a parameter file, which is usually named
`Parameters.in`. You can use the input files for the
[tutorials](https://github.com/w2dynamics/w2dynamics/wiki/Tutorials)
on our wiki as templates.

Next, run `DMFT.py` to use the self-consistency loop or `cthyb` if you
do a single-shot calculation. If you have installed the code (see
above), both executables should be placed in your path. If your
parameter file is not named `Parameters.in` or does not lie in your
current working directory, you have to specify it explicitly as a
command line argument. If you want to run your calculation using MPI
parallelization, use `mpiexec` or similar to execute the script.

    $ DMFT.py [Parameters.in]
    $ cthyb [Parameters.in]
    $ mpiexec -n 10 DMFT.py [Parameters.in]
    $ mpiexec -n 10 cthyb [Parameters.in]

The code will produce a file with a name like `RunIdentifier-Timestamp.hdf5`.
It is an archive of all quantities written by w2dynamics.  You can navigate this
file using any hdf5-compatible analysis tool, such as jupyter, matlab, etc.
For your convenience, we have also included the tool `hgrep`, which allows
quick analysis of the data:

    $ hgrep [options] (file|latest) quantity [[index] ...]
    $ hgrep latest siw 1 1 1 1

will print the self-energy on the Matsubara axis from the first
iteration for the first inequality, orbital and spin as tabular
data. Add option `-p` for automatic plotting or have a look at the man
page `hgrep.man` or the examples in our
[tutorials](https://github.com/w2dynamics/w2dynamics/wiki/Tutorials)
for details.


Files and directories
---------------------

  - `w2dyn/auxiliaries/`: auxiliary python routines (in/output, config files, etc.)
  - `w2dyn/dmft/`: python package for DMFT self-consistency loop
  - `w2dyn/maxent/`: Python wrapper for maximum entropy analytic continuation
  - `clusters/`: template submission scripts for different clusters
  - `cmake/`: cmake custom modules
  - `docs/`: documentation (github wiki)
  - `Postproc/`: postprocessing scripts
  - `preproc/`: preprocessing scripts
  - `src/`: compiled modules loaded from python
    - `ctqmc_fortran`: Fortran 90 continuous-time quantum Monte Carlo solver
    - `maxent`: maximum entropy analytic continuation solver
    - `mtrng`: Mersenne twister pseudorandom number generator
  - `testsuite/`: unit tests for the code

  - `cfg_converter.py`: small script converting old-style config files
  - `completions.sh`: file for bash completions
  - `cprun`: convenience script copying input files to different directory
  - `DMFT.py`: main entry point for DMFT self-consistency loop
  - `hgrep`: utility for extracting data from HDF5 file
  - `Maxent.py`: main entry point for maximum entropy code
  - `setup.py`: Python installation script

Citation
--------
If you use the w2dynamics package, please mention the following in the acknowledments:

The QMC simulations were carried out with the w2dynamics package available at https://github.com/w2dynamics/w2dynamics .


[Disclaimer](https://www.uni-wuerzburg.de/sonstiges/impressum/)

[Github Privacy Statement](https://help.github.com/articles/github-privacy-statement/)
