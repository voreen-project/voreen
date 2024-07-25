{ pkgs, ... }:

{
  gcc = {
    sequential = {
      cxx   = "g++";
      cc    = "gcc";
      flags = "-O3 -Wall -march=native -mtune=native";
      parallel_mode = "OFF";
      buildInputs = with pkgs; [
        gcc
      ];
    };

    debug = {
      cxx   = "g++";
      cc    = "gcc";
      flags = "-g -Wall -DOLB_DEBUG -march=native";
      parallel_mode = "OFF";
      buildInputs = with pkgs; [
        gcc
      ];
    };

    openmpi = {
      cxx   = "mpic++";
      cc    = "gcc";
      flags = "-O3 -Wall -march=native -mtune=native";
      parallel_mode = "MPI";
      buildInputs = with pkgs; [
        gcc
        openmpi
      ];
    };

    openmp = {
      cxx   = "g++";
      cc    = "gcc";
      flags = "-O3 -Wall -march=native -mtune=native";
      parallel_mode = "OMP";
      buildInputs = with pkgs; [
        gcc
      ];
    };
  };

  clang = {
    sequential = {
      cxx   = "clang++";
      cc    = "clang";
      flags = "-O3 -Wall -march=native -mtune=native";
      ldflags = "-fuse-ld=lld";
      parallel_mode = "OFF";
      buildInputs = with pkgs; [
        clang
        llvmPackages.bintools-unwrapped
      ];
    };

    openmp = {
      cxx   = "clang++";
      cc    = "clang";
      flags = "-O3 -Wall -march=native -mtune=native";
      ldflags = "-fuse-ld=lld";
      parallel_mode = "OMP";
      buildInputs = with pkgs; [
        clang
        llvmPackages.bintools-unwrapped
        llvmPackages.openmp
      ];
    };
  };

  nvcc = {
    sequential = {
      cxx   = "nvcc";
      cc    = "nvcc";
      flags = "-O3";
      parallel_mode = "OFF";
      platform = "GPU_CUDA";
      buildInputs = with pkgs; [
        gcc
        cudatoolkit_11
      ];
    };

    openmpi = {
      cxx   = "nvcc";
      cc    = "nvcc";
      flags = "-O3";
      parallel_mode = "MPI";
      ldflags = "-lmpi_cxx -lmpi";
      platform = "GPU_CUDA";
      buildInputs = with pkgs; [
        gcc
        cudatoolkit_11
        (openmpi.override { cudaSupport = true; cudatoolkit = cudatoolkit_11; })
      ];
    };
  };
}
