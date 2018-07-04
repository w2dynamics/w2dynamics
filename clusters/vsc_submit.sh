#!/bin/bash
#$ -N CTQMC
#$ -pe mpich 48
#$ -V
#$ -l h_rt=24:00:00

trap '' USR1 USR2
mpirun -v -x LD_LIBRARY_PATH -x PYTHONUSERBASE -x DMFT_GIT_REVISION \
       -machinefile $TMPDIR/machines -np 48 DMFT.py
