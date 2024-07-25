# Example build config for OpenLB using CUDA on single GPU systems
#
# Tested using CUDA 11.4
#
# Usage:
#  - Copy this file to OpenLB root as `config.mk`
#  - Run `make clean; make`
#  - Switch to example directory, e.g. `examples/laminar/cavity3dBenchmark`
#  - Run `make`
#  - Start the simulation using `./cavity3d`

CXX             := nvcc
CC              := nvcc

CXXFLAGS        := -O3
CXXFLAGS        += -std=c++17

PARALLEL_MODE   := NONE

PLATFORMS       := CPU_SISD GPU_CUDA

# for e.g. RTX 30* (Ampere), see table in `rules.mk` for other options
CUDA_ARCH       := 86

USE_EMBEDDED_DEPENDENCIES := ON
