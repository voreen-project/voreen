#  This file is part of the OpenLB library
#
#  Copyright (C) 2017 Markus Mohrhard, Mathias Krause
#  E-mail contact: info@openlb.net
#  The most recent release of OpenLB can be downloaded at
#  <http://www.openlb.net/>
#
#  This program is free software; you can redistribute it and/or
#  modify it under the terms of the GNU General Public License
#  as published by the Free Software Foundation; either version 2
#  of the License, or (at your option) any later version.
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

CXX             := g++
#CXX             := icpc -D__aligned__=ignored
#CXX             := mpiCC
#CXX             := mpic++

CC              := gcc                                          # necessary for zlib, for Intel use icc

OPTIM           := -O3 -Wall -march=native -mtune=native        # for gcc
#OPTIM           := -O3 -Wall -xHost                            # for Intel compiler
#OPTIM           := -O3 -Wall -xHost -ipo                       # optional for Intel compiler
DEBUG           := -g -Wall -DOLB_DEBUG

CXXFLAGS        := $(OPTIM)
#CXXFLAGS        := $(DEBUG)

# to enable std::shared_ptr in functor arithmetik
# works in gcc 4.3 and later, source https://gcc.gnu.org/projects/cxx0x.html
#CXXFLAGS        += -std=c++0x
# works in gcc 4.7 and later (recommended)
#CXXFLAGS        += -std=c++11
# Voreen modi
CXXFLAGS        += -std=c++11 -fPIC #fPIC is necessary!!!

# optional for gcc 4.9 and newer
#CXXFLAGS        += -fdiagnostics-color=auto
#CXXFLAGS        += -std=gnu++14

ARPRG           := ar
#ARPRG           := xiar                  # mandatory for intel compiler

LDFLAGS         :=

#PARALLEL_MODE   := OFF
#PARALLEL_MODE   := MPI
PARALLEL_MODE   := OMP
#PARALLEL_MODE   := HYBRID

MPIFLAGS        :=
OMPFLAGS        := -fopenmp

BUILDTYPE       := precompiled
#BUILDTYPE       := generic
