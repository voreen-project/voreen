# Example build config for OpenLB with SIMD using Intel C++ and MPI library
#
# Recommended for CPU-only clusters with AVX2 / AVX-512 capability
#
# Usage:
#  - Copy this file to OpenLB root as `config.mk`
#  - Run `make clean; make`
#  - Switch to example directory, e.g. `examples/laminar/cavity3dBenchmark`
#  - Run `make`
#  - Start the simulation using `mpirun ./cavity3d`

CXX             := mpicxx
CC              := icc

CXXFLAGS        := -O3 -Wall -xHost -ipo -axMIC-AVX512,CORE-AVX2
CXXFLAGS        += -std=c++17

# HYBRID mode is also possible but more complex to run correctly
PARALLEL_MODE   := MPI

OMPFLAGS        := -fopenmp

PLATFORMS       := CPU_SISD CPU_SIMD

USE_EMBEDDED_DEPENDENCIES := ON
