# Example build config for OpenLB using CUDA and OpenMPI
#
# Tested using CUDA 11.4 and OpenMPI 4.1 (CUDA aware)
#
# Usage:
#  - Copy this file to OpenLB root as `config.mk`
#  - Run `make clean; make`
#  - Switch to example directory, e.g. `examples/laminar/cavity3dBenchmark`
#  - Run `make`
#  - Start the simulation using `mpirun -np 2 ./cavity3d` (All processes share default GPU, not optimal)
#
# Usage on a multi GPU system: (recommended when using MPI, use non-MPI version on single GPU systems)
#  - Run `mpirun -np 4 bash -c 'export CUDA_VISIBLE_DEVICES=${OMPI_COMM_WORLD_LOCAL_RANK}; ./cavity3d'
#    (for a 4 GPU system, further process mapping advisable, consult cluster documentation)
#
# CXXFLAGS and LDFLAGS may need to be adjusted depending on the specific MPI installation.
# Compare to `mpicxx --showme:compile` and `mpicxx --showme:link` when in doubt.

CXX             := nvcc
CC              := nvcc

CXXFLAGS        := -O3
CXXFLAGS        += -std=c++17

PARALLEL_MODE   := MPI

MPIFLAGS        := -lmpi_cxx -lmpi

PLATFORMS       := CPU_SISD GPU_CUDA

# for e.g. RTX 30* (Ampere), see table in `rules.mk` for other options
CUDA_ARCH       := 86

USE_EMBEDDED_DEPENDENCIES := ON
