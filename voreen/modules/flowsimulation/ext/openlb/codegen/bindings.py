# This file is part of the OpenLB library
#
# Copyright (C) 2021 Adrian Kummerlaender
# E-mail contact: info@openlb.net
# The most recent release of OpenLB can be downloaded at
# <http://www.openlb.net/>
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public
# License along with this program; if not, write to the Free
# Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
# Boston, MA  02110-1301, USA.

import cppyy

cppyy.add_include_path('src')
cppyy.add_include_path('external/tinyxml')

cppyy.cppdef('#define DISABLE_CSE')
cppyy.cppdef('#define PLATFORM_CPU_SISD')
cppyy.cppdef('#define any_platform')

cppyy.include('utilities/expr.h')
cppyy.include('olb2D.h')
cppyy.include('olb2D.hh')
cppyy.include('olb3D.h')
cppyy.include('olb3D.hh')

from cppyy.gbl import olb

# Using FORCE field as dummy due to instantiation failure on empty variadic pack in cppyy
descriptors = {
    'D2Q5': olb.descriptors.D2Q5[olb.descriptors.FORCE],
    'D2Q9': olb.descriptors.D2Q9[olb.descriptors.FORCE],
    'D3Q7': olb.descriptors.D3Q7[olb.descriptors.FORCE],
    'D3Q19': olb.descriptors.D3Q19[olb.descriptors.FORCE],
    'D3Q27': olb.descriptors.D3Q27[olb.descriptors.FORCE],
}

mrtDescriptors = {
    'D2Q5': olb.descriptors.D2Q5[olb.descriptors.tag.MRT,olb.descriptors.FORCE],
    'D2Q9': olb.descriptors.D2Q9[olb.descriptors.tag.MRT,olb.descriptors.FORCE],
    'D3Q7': olb.descriptors.D3Q7[olb.descriptors.tag.MRT,olb.descriptors.FORCE],
    'D3Q19': olb.descriptors.D3Q19[olb.descriptors.tag.MRT,olb.descriptors.FORCE],
}

rtlbmDescriptors = {
    'D3Q7': olb.descriptors.D3Q7[olb.descriptors.tag.RTLBM],
    'D3Q15': olb.descriptors.D3Q15[olb.descriptors.tag.RTLBM],
    'D3Q27': olb.descriptors.D3Q27[olb.descriptors.tag.RTLBM],
}

freeEnergyDescriptors = {
    'D2Q9': olb.descriptors.D2Q9[olb.descriptors.CHEM_POTENTIAL,olb.descriptors.FORCE],
    'D3Q19': olb.descriptors.D3Q19[olb.descriptors.CHEM_POTENTIAL,olb.descriptors.FORCE],
}
