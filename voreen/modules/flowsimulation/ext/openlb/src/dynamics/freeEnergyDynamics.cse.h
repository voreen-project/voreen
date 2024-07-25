/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2021 Adrian Kummerlaender
 *  E-mail contact: info@openlb.net
 *  The most recent release of OpenLB can be downloaded at
 *  <http://www.openlb.net/>
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2
 *  of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public
 *  License along with this program; if not, write to the Free
 *  Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
 *  Boston, MA  02110-1301, USA.
*/

#ifndef DYNAMICS_FREE_ENERGY_DYNAMICS_CSE_H
#define DYNAMICS_FREE_ENERGY_DYNAMICS_CSE_H


#ifndef DISABLE_CSE

#include "latticeDescriptors.h"

namespace olb {

namespace collision {

template <typename... FIELDS>
struct FreeEnergy::type<descriptors::D2Q9<FIELDS...>,momenta::FreeEnergyBulkTuple,equilibria::FreeEnergy> {

template <typename CELL, typename PARAMETERS, typename V=typename CELL::value_t>
CellStatistic<V> apply(CELL& cell, PARAMETERS& parameters) any_platform
{
auto x12 = parameters.template get<olb::descriptors::OMEGA>();
auto x10 = cell.template getFieldComponent<olb::descriptors::FORCE>(0);
auto x13 = parameters.template get<olb::collision::FreeEnergy::GAMMA>();
auto x9 = cell.template getFieldComponent<olb::descriptors::CHEM_POTENTIAL>(0);
auto x11 = cell.template getFieldComponent<olb::descriptors::FORCE>(1);
auto x14 = x12 + V{-1};
auto x15 = cell[0] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + cell[7] + cell[8];
auto x16 = x15 + V{1};
auto x17 = -V{3}*x13*x9 + x16;
auto x18 = x15 + V{1};
auto x19 = x11*x11;
auto x20 = V{1.5}*x19;
auto x21 = x10*x10;
auto x22 = V{1.5}*x21;
auto x23 = x22 + V{-1};
auto x24 = x20 + x23;
auto x25 = V{0.0277777777777778}*x12;
auto x26 = V{3}*x11;
auto x27 = -x10 + x11;
auto x28 = V{3}*x10;
auto x29 = x24 + x28;
auto x30 = x17*x25;
auto x31 = V{0.111111111111111}*x12;
auto x32 = -x28;
auto x33 = V{1} - x20;
auto x34 = V{3}*x21 + x33;
auto x35 = x17*x31;
auto x36 = -x35;
auto x37 = x10 + x11;
auto x38 = V{4.5}*(x37*x37);
auto x39 = V{3}*x19;
auto x40 = x10 - x11;
auto x41 = -x22 + x26;
cell[0] = V{0.555555555555556}*x12*x17 - V{0.444444444444444}*x12*(x18*x24 + V{1}) - x14*cell[0];
cell[1] = -x14*cell[1] - x25*(x18*(-x26 + x29 - V{4.5}*x27*x27) + V{1}) - x30;
cell[2] = -x14*cell[2] + x31*(x18*(x32 + x34) + V{-1}) + x36;
cell[3] = -x14*cell[3] - x25*(x18*(x26 + x29 - x38) + V{1}) - x30;
cell[4] = -x14*cell[4] - x31*(x18*(x23 + x26 - x39) + V{1}) - x35;
cell[5] = -x14*cell[5] - x25*(x18*(x24 + x26 + x32 - V{4.5}*x40*x40) + V{1}) - x30;
cell[6] = -x14*cell[6] + x31*(x18*(x28 + x34) + V{-1}) + x36;
cell[7] = -x14*cell[7] + x25*(x18*(x28 + x33 + x38 + x41) + V{-1}) - x30;
cell[8] = -x14*cell[8] + x31*(x18*(x39 + x41 + V{1}) + V{-1}) + x36;
return { x16, x19 + x21 };
}

};

template <typename... FIELDS>
struct FreeEnergy::type<descriptors::D3Q19<FIELDS...>,momenta::FreeEnergyBulkTuple,equilibria::FreeEnergy> {

template <typename CELL, typename PARAMETERS, typename V=typename CELL::value_t>
CellStatistic<V> apply(CELL& cell, PARAMETERS& parameters) any_platform
{
auto x19 = cell.template getFieldComponent<olb::descriptors::CHEM_POTENTIAL>(0);
auto x20 = cell.template getFieldComponent<olb::descriptors::FORCE>(0);
auto x21 = cell.template getFieldComponent<olb::descriptors::FORCE>(1);
auto x22 = cell.template getFieldComponent<olb::descriptors::FORCE>(2);
auto x24 = parameters.template get<olb::collision::FreeEnergy::GAMMA>();
auto x23 = parameters.template get<olb::descriptors::OMEGA>();
auto x25 = x23 + V{-1};
auto x26 = cell[0] + cell[10] + cell[11] + cell[12] + cell[13] + cell[14] + cell[15] + cell[16] + cell[17] + cell[18] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + cell[7] + cell[8] + cell[9];
auto x27 = x26 + V{1};
auto x28 = -V{3}*x19*x24 + x27;
auto x29 = x26 + V{1};
auto x30 = x20*x20;
auto x31 = V{1.5}*x30;
auto x32 = x21*x21;
auto x33 = V{1.5}*x32;
auto x34 = x22*x22;
auto x35 = V{1.5}*x34;
auto x36 = x33 + x35 + V{-1};
auto x37 = x31 + x36;
auto x38 = V{0.0555555555555556}*x23;
auto x39 = V{3}*x20;
auto x40 = V{3}*x30;
auto x41 = x28*x38;
auto x42 = V{3}*x21;
auto x43 = V{3}*x32;
auto x44 = x31 + V{-1};
auto x45 = V{3}*x22;
auto x46 = V{3}*x34;
auto x47 = V{0.0277777777777778}*x23;
auto x48 = x20 + x21;
auto x49 = V{4.5}*(x48*x48);
auto x50 = x37 + x39;
auto x51 = x28*x47;
auto x52 = -x42;
auto x53 = -x20;
auto x54 = x21 + x53;
auto x55 = x20 + x22;
auto x56 = V{4.5}*(x55*x55);
auto x57 = -x45;
auto x58 = x22 + x53;
auto x59 = x21 + x22;
auto x60 = V{4.5}*(x59*x59);
auto x61 = x37 + x42;
auto x62 = -x21;
auto x63 = x22 + x62;
auto x64 = -x33;
auto x65 = V{1} - x35;
auto x66 = x64 + x65;
auto x67 = x39 + x66;
auto x68 = -x41;
auto x69 = -x31;
auto x70 = x42 + x69;
auto x71 = x45 + x69;
auto x72 = -x51;
auto x73 = -x39;
auto x74 = x20 + x62;
auto x75 = -x22;
auto x76 = x20 + x75;
auto x77 = x37 + x45;
auto x78 = x21 + x75;
cell[0] = V{0.666666666666667}*x23*x28 - V{0.333333333333333}*x23*(x29*x37 + V{1}) - x25*cell[0];
cell[1] = -x25*cell[1] - x38*(x29*(x36 + x39 - x40) + V{1}) - x41;
cell[2] = -x25*cell[2] - x38*(x29*(x35 + x42 - x43 + x44) + V{1}) - x41;
cell[3] = -x25*cell[3] - x38*(x29*(x33 + x44 + x45 - x46) + V{1}) - x41;
cell[4] = -x25*cell[4] - x47*(x29*(x42 - x49 + x50) + V{1}) - x51;
cell[5] = -x25*cell[5] - x47*(x29*(x50 + x52 - V{4.5}*x54*x54) + V{1}) - x51;
cell[6] = -x25*cell[6] - x47*(x29*(x45 + x50 - x56) + V{1}) - x51;
cell[7] = -x25*cell[7] - x47*(x29*(x50 + x57 - V{4.5}*x58*x58) + V{1}) - x51;
cell[8] = -x25*cell[8] - x47*(x29*(x45 - x60 + x61) + V{1}) - x51;
cell[9] = -x25*cell[9] - x47*(x29*(x57 + x61 - V{4.5}*x63*x63) + V{1}) - x51;
cell[10] = -x25*cell[10] + x38*(x29*(x40 + x67) + V{-1}) + x68;
cell[11] = -x25*cell[11] + x38*(x29*(x43 + x65 + x70) + V{-1}) + x68;
cell[12] = -x25*cell[12] + x38*(x29*(x46 + x64 + x71 + V{1}) + V{-1}) + x68;
cell[13] = -x25*cell[13] + x47*(x29*(x49 + x67 + x70) + V{-1}) + x72;
cell[14] = -x25*cell[14] - x47*(x29*(x61 + x73 - V{4.5}*x74*x74) + V{1}) - x51;
cell[15] = -x25*cell[15] + x47*(x29*(x56 + x67 + x71) + V{-1}) + x72;
cell[16] = -x25*cell[16] - x47*(x29*(x73 + x77 - V{4.5}*x76*x76) + V{1}) - x51;
cell[17] = -x25*cell[17] + x47*(x29*(x45 + x60 + x66 + x70) + V{-1}) + x72;
cell[18] = -x25*cell[18] - x47*(x29*(x52 + x77 - V{4.5}*x78*x78) + V{1}) - x51;
return { x27, x30 + x32 + x34 };
}

};



}

}

#endif

#endif