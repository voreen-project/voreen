{ cxx, cc, flags, ldflags ? "", parallel_mode, platform ? "", ... }:

''
CXX             := ${cxx}
CC              := ${cc}

CXXFLAGS        := ${flags}
CXXFLAGS        += -std=c++17

LDFLAGS         := ${ldflags}

PARALLEL_MODE   := ${parallel_mode}

MPIFLAGS        :=
OMPFLAGS        := -fopenmp

PLATFORMS       := CPU_SISD ${platform}

FEATURES        :=

USE_EMBEDDED_DEPENDENCIES := ON
''
