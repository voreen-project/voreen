#  This file is part of the OpenLB library
#
#  Copyright (C) 2007, 2017 Markus Mohrhard, Mathias Krause
#  E-mail contact: info@openlb.net
#  The most recent release of OpenLB can be downloaded at
#  <http://www.openlb.net/>
#
#  This program is free software; you can redistribute it and/or
#  modify it under the terms of the GNU General Public License
#  as published by the Free Software Foundation; either version 2
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public 
#  License along with this program; if not, write to the Free 
#  Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
#  Boston, MA  02110-1301, USA.

###########################################################################
###########################################################################

mkfile_path := $(abspath $(lastword $(MAKEFILE_LIST)))

include $(dir $(mkfile_path))/config.mk

###########################################################################
## conditional settings

ifeq ($(BUILDTYPE), precompiled)
   CXXFLAGS := -DOLB_PRECOMPILED $(CXXFLAGS)
endif

ifeq ($(PARALLEL_MODE), MPI)
   CXXFLAGS := -DPARALLEL_MODE_MPI $(MPIFLAGS) $(CXXFLAGS)
endif

ifeq ($(PARALLEL_MODE), OMP)
   CXXFLAGS := -DPARALLEL_MODE_OMP $(OMPFLAGS) $(CXXFLAGS)
   LDFLAGS  := $(OMPFLAGS) $(LDFLAGS)
endif

ifeq ($(PARALLEL_MODE), HYBRID)
   CXXFLAGS := -DPARALLEL_MODE_OMP -DPARALLEL_MODE_MPI $(OMPFLAGS) $(MPIFLAGS) $(CXXFLAGS)
   LDFLAGS  := $(OMPFLAGS) $(LDFLAGS)
endif

###########################################################################
## defines shell

SHELL           := /bin/sh

###########################################################################
## dependencies, object, library directory and library name

DEPENDDIR       := build/$(BUILDTYPE)/dep
OBJDIR          := build/$(BUILDTYPE)/obj
LIBDIR          := build/$(BUILDTYPE)/lib
LIB             := olb
LIBS            := -l$(LIB) -lz

###########################################################################
## search directories

SUBDIRS         := src/boundary \
                   src/communication \
                   src/dynamics \
                   src/core \
                   src/geometry \
                   src/external/tinyxml \
                   src/external/zlib \
                   src/functors \
                   src/functors/analytical \
                   src/functors/analytical/indicator \
                   src/functors/lattice \
                   src/functors/lattice/indicator \
                   src/functors/lattice/integral \
                   src/io \
                   src/particles \
                   src/particles/forces \
                   src/particles/boundaries \
                   src/utilities

EXAMPLEDIRS     := examples/aorta3d \
                   examples/bifurcation3d/eulerEuler \
                   examples/bifurcation3d/eulerLagrange \
                   examples/bstep2d \
                   examples/bstep3d \
                   examples/cavity2d/sequential \
                   examples/cavity2d/parallel \
                   examples/cavity3d/sequential \
                   examples/cavity3d/parallel \
                   examples/cylinder2d \
                   examples/cylinder3d \
                   examples/multiComponent2d \
                   examples/multiComponent3d \
                   examples/nozzle3d \
                   examples/phaseSeparation2d \
                   examples/phaseSeparation3d \
                   examples/poiseuille2d \
                   examples/poiseuille3d \
                   examples/porousPoiseuille2d \
                   examples/powerLaw2d \
                   examples/tgv3d \
                   examples/thermalFlows/porousPlate2d \
                   examples/thermalFlows/porousPlate3d \
                   examples/thermalFlows/rayleighBenard2d \
                   examples/thermalFlows/rayleighBenard3d \
                   examples/thermalFlows/squareCavity2d \
                   examples/thermalFlows/squareCavity3d \
                   examples/venturi3d

INCLUDEDIRS     := src \
                   src/ \
                   src/external \
                   src/external/zlib

BUILDTYPEDIRS   := build/precompiled \
                   build/generic

IDIR            := $(foreach d,$(INCLUDEDIRS),-I$(ROOT)/$(d))
