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
        llvmPackages.bintools
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
        llvmPackages.bintools
        llvmPackages.openmp
      ];
    };
  };
  intel = {
    sequential = {
      cxx   = "icpc -D__aligned__=ignored";
      cc    = "icc";
      flags = "-O3 -Wall -xHost";
      arprg = "xiar";
      parallel_mode = "OFF";
      buildInputs = [ ];
    };
  };
}
