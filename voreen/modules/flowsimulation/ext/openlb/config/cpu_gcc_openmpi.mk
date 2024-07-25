# Example build config for OpenLB using GNU C++ and OpenMPI
#
# Usage:
#  - Copy this file to OpenLB root as `config.mk`
#  - Run `make clean; make`
#  - Switch to example directory, e.g. `examples/laminar/cavity3dBenchmark`
#  - Run `make`
#  - Start the simulation using `mpirun ./cavity3d`
#
# Usage of the Intel C++ compiler is recommended for Intel CPU clusters.
# See `config/cpu_simd_intel_mpi.mk` for guidance.

CXX             := mpic++
CC              := gcc

# The `march=native` flag enables AVX2 / AVX-512 instructions if available.
# However, actually using them requires adding the `CPU_SIMD` platform.
#
# Note that on some clusters the head node for compilation may differ from
# the compute nodes, necessitating manual selection of the correct
# architecture / SIMD flags. Alternatively, compilation at the start of
# the HPC jobs is a common option.
CXXFLAGS        := -O3 -Wall -march=native -mtune=native
CXXFLAGS        += -std=c++17

# HYBRID mode is also possible but more complex to run correctly
PARALLEL_MODE   := MPI

# optional MPI and OpenMP flags
OMPFLAGS        := -fopenmp

# SIMD support may optionally be enabled by adding the `CPU_SIMD` platform
PLATFORMS       := CPU_SISD # CPU_SIMD

USE_EMBEDDED_DEPENDENCIES := ON
