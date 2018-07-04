#!/bin/bash
#PBS -V
#PBS -l nodes=2:ppn=12
#PBS -l walltime=01:00:00

echo "Starting script..."
BIN=$PBS_O_HOME/opt/bin

cd $PBS_O_WORKDIR

NPROCS=`wc -l < $PBS_NODEFILE`

trap '' USR1 USR2
mpirun -v -x LD_LIBRARY_PATH -x PYTHONUSERBASE -x DMFT_GIT_REVISION \
       --mca yield_when_idle 1 -machinefile $PBS_NODEFILE -np $NPROCS \
       python $BIN/DMFT.py >ctqmc-$PBS_JOBID.log 2>&1

echo "Script complete."
