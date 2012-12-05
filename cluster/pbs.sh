#!/bin/sh
#PBS -l walltime=10:00:00
#PBS -l nodes=1:ppn=16
#PBS -l mem=2048

cd $PBS_O_WORKDIR

# Make a temporary directory on the local node
tmp_dir=`mktemp -d`
cp -r ./* $tmp_dir
cd $tmp_dir

# Start the execution
export OMP_NUM_THREADS=16
./habitat.exe

# Copy the results back
mkdir -p $PBS_O_WORKDIR/$PBS_JOBID
cp -r * $PBS_O_WORKDIR/$PBS_JOBID

rm -rf $tmp_dir

