#!/bin/sh
#SBATCH --time=10:00:00
#SBATCH --nodes=1   #asdf
#SBATCH --ntasks-per-node=16
#SBATCH --mem-per-cpu=2048
#SBATCH --job-name=habitat
#SBATCH --error=job.%J.err
#SBATCH --output=job.%J.out

# Make a temporary directory
orig_dir=`cwd`
tmp_dir=`mktemp -d`
cp -r ./* $tmp_dir
cd $tmp_dir

# Start execution
export OMP_NUM_THREADS=16
./habitat.exe

# Copy the results back
cp -r * $orig_dir
rm -rf $tmp_dir

