### Apptainer / Singularity definition files

Using a singularity container provides an alternative to setting up
the needed environment and building w2dynamics on the machine where
you want to use it while still allowing you, under certain conditions,
to access files on the host and to use MPI parallelization almost as
easily as if you did that with usually negligible overhead.

In this directory, you can find definition files
`w2dynamics_distribution.def` that you can use to build a singularity
image file `w2dynamics.sif` by running `singularity build
w2dynamics.sif w2dynamics_distribution.def` with administrative
privileges on some machine. When you do this, Singularity will
download a base image from Docker Hub and perform all necessary
operations to compile and install w2dynamics (into paths contained in
the resulting image) automatically. Consult the singularity
documentation for more details on how to build an image.

You can then run the main executable scripts, e.g. `DMFT.py`, using a
command such as `singularity run --app DMFT.py w2dynamics.sif
[OPTIONS]` (you can use `DMFT.py`, `cthyb`, `hgrep`, or `Maxent.py` as
app arguments). Without further options, this example command would
run a DMFT calculation reading from `Parameters.in` in your working
directory and writing to the file specified using `FileNamePrefix` in
the configuration just as usual provided that your working directory
and the output path are both correctly mounted in the container, which
should for example be the case if both are in your home directory. For
MPI parallelization, you need an installation of MPI on the host
system as well and can then just run the same command using `mpiexec`
or similar. For more details on running commands in the container,
mounting host directories in the container, and MPI parallelization
you can consult the singularity documentation.

We do currently not provide prebuilt images. You should be able to use
any of the provided definition files unmodified, but there could be
compatibility issues if the MPI installations in the container (some
version of Open MPI in all cases) and on the host are different which
you might need to fix yourself. In that case, it may also be necessary
to build the `mpi4py` Python package in the container from source.
