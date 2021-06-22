{ cxx, cc, flags, ldflags ? "", arprg ? "ar", parallel_mode, ... }:

''
CXX             := ${cxx}
CC              := ${cc}

CXXFLAGS        := ${flags}
CXXFLAGS        += -std=c++14

LDFLAGS         := ${ldflags}

ARPRG           := ${arprg}

PARALLEL_MODE   := ${parallel_mode}

MPIFLAGS        :=
OMPFLAGS        := -fopenmp

BUILDTYPE       := generic

FEATURES        :=
''
