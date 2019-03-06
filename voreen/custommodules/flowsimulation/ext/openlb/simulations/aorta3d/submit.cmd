#!/bin/bash

# set the number of nodes
#SBATCH --nodes=1

# MPI/OMP Hybrid config
#SBATCH --ntasks=18
#SBATCH --ntasks-per-node=18
#SBATCH --cpus-per-task=4

# set the number of CPU cores per node
#SBATCH --exclusive

# How much memory is needed (per node)
#SBATCH --mem=16G

# set a partition
#SBATCH --partition normal

# set max wallclock time
#SBATCH --time=00:10:00

# set name of job
#SBATCH --job-name=olb-simulation-aorta3d

# mail alert at start, end and abortion of execution
#SBATCH --mail-type=ALL

# set an output file
#SBATCH --output output.dat

# send mail to this address
#SBATCH --mail-user=s_leis06@uni-muenster.de

# run the application
OMP_NUM_THREADS=4 mpirun /home/s/s_leis06/OpenLB-intel/simulations/aorta3d/aorta3d
