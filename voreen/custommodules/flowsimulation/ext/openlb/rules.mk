###########################################################################
## conditional settings

ifeq ($(PARALLEL_MODE), MPI)
   CXXFLAGS := -DPARALLEL_MODE_MPI $(MPIFLAGS) $(CXXFLAGS)
   LDFLAGS  := $(MPIFLAGS) $(LDFLAGS)
endif

ifeq ($(PARALLEL_MODE), OMP)
   CXXFLAGS := -DPARALLEL_MODE_OMP $(OMPFLAGS) $(CXXFLAGS)
   LDFLAGS  := $(OMPFLAGS) $(LDFLAGS)
endif

ifeq ($(PARALLEL_MODE), HYBRID)
   CXXFLAGS := -DPARALLEL_MODE_OMP -DPARALLEL_MODE_MPI $(OMPFLAGS) $(MPIFLAGS) $(CXXFLAGS)
   LDFLAGS  := $(MPIFLAGS) $(OMPFLAGS) $(LDFLAGS)
endif

LDFLAGS += $(if $(filter $(FEATURES), OPENBLAS),-lopenblas)
LDFLAGS += -lz -ltinyxml

ifneq ($(filter CPU_SIMD,$(PLATFORMS)),)
	LDFLAGS += -lrt
endif

ifneq ($(filter GPU_CUDA,$(PLATFORMS)),)
## | CUDA Architecture | Version    |
## |-------------------+------------|
## | Fermi             | 20         |
## | Kepler            | 30, 35, 37 |
## | Maxwell           | 50, 52, 53 |
## | Pascal            | 60, 61, 62 |
## | Volta             | 70, 72     |
## | Turing            | 75         |
## | Ampere            | 80, 86, 87 |
	CUDA_ARCH ?= 60

	LDFLAGS += -lcuda -lcudadevrt -lcudart
	CXXFLAGS += --generate-code=arch=compute_$(CUDA_ARCH),code=[compute_$(CUDA_ARCH),sm_$(CUDA_ARCH)]
	CXXFLAGS += --extended-lambda --expt-relaxed-constexpr -x cu
	CXXFLAGS += -Xcudafe "--diag_suppress=implicit_return_from_non_void_function --display_error_number --diag_suppress=20014 --diag_suppress=20011"
endif

CXXFLAGS += $(foreach platform,$(PLATFORMS),-DPLATFORM_$(platform))

CXXFLAGS += $(foreach feature,$(FEATURES),-DFEATURE_$(feature))
