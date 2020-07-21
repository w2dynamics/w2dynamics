### Singularity definition files

Using a singularity container provides an alternative to setting up
the needed environment and building w2dynamics on the machine where
you want to use it while still allowing you to access files on the
host and to use MPI parallelization almost as easily as if you did
that (under certain conditions) with usually negligible overhead.

In this directory, you can find example definition files
`w2dynamics_distribution.def` that you can use to create a singularity
image file `w2dynamics.sif` by running `singularity build
w2dynamics.sif w2dynamics_distribution.def` with administrative
privileges on some machine. Singularity will then automatically
download a base image from the Singularity container library (Ubuntu,
Debian) or use an installation of pacman on the host system (Arch) to
set up a base environment in the container and then execute the
commands in the `%post` section of the definition file inside the
container to install all dependencies, clone the w2dynamics GitHub
repository into `/w2dynamics` in the container and build the native
code extension modules using GCC. Consult the singularity
documentation for more details on how to build an image.

You can then run the main executable scripts, e.g. `DMFT.py`, using a
command such as `singularity run --app DMFT w2dynamics.sif [OPTIONS]`
(predefined are `DMFT`, `cthyb`, `hgrep`, and `Maxent`), which would
run a DMFT calculation reading from `Parameters.in` in your working
directory and writing to the file specified using `FileNamePrefix` in
the configuration just as usual provided that your working directory
and the output path are both correctly mounted in the container, which
should for example be the case if both are in your home directory. For
MPI parallelization, you need an installation of MPI on the host
system as well and can then just run the same command using `mpiexec`
or similar. For more details on running commands in the container,
mounting host directories in the container, and MPI parallelization
you can, again, consult the singularity documentation.

We do currently not provide prebuilt images. If you want to use the
example definition files unmodified, consider that the Ubuntu one
currently produces the smallest image and the Arch one requires an
installation of pacman on the host system and takes longer to build
because OpenBLAS and NFFT are compiled during the build as
well. Further, there could be compatibility issues if the MPI
installations in the container (OpenMPI) and on the host are different
which you might need to fix yourself. In that case, it may also be
necessary to build the `mpi4py` Python package in the container from
source.
