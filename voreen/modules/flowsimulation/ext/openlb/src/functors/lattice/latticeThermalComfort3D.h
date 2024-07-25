/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012 Lukas Baron, Tim Dornieden, Mathias J. Krause,
 *  Albert Mink
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

#ifndef LATTICE_THERMAL_COMFORT_3D_H
#define LATTICE_THERMAL_COMFORT_3D_H

#include<vector>

#include "superBaseF3D.h"
#include "superCalcF3D.h"
#include "functors/analytical/indicator/indicatorBaseF3D.h"

#include "blockBaseF3D.h"
#include "geometry/blockGeometry.h"
#include "functors/analytical/indicator/indicatorBaseF3D.h"
#include "indicator/blockIndicatorBaseF3D.h"
#include "dynamics/smagorinskyBGKdynamics.h"
#include "dynamics/porousBGKdynamics.h"


/* Note: Throughout the whole source code directory genericFunctions, the
 *  template parameters for i/o dimensions are:
 *           F: S^m -> T^n  (S=source, T=target)
 */

namespace olb {

/// functor to get pointwise PMV and PPD on local lattices to evaluate thermal comfort
template <typename T, typename DESCRIPTOR, typename TDESCRIPTOR>
class SuperLatticeThermalComfort3D final : public SuperLatticeThermalPhysF3D<T,DESCRIPTOR,TDESCRIPTOR> { // templatename before <
public:
  SuperLatticeThermalComfort3D(SuperLattice<T,TDESCRIPTOR>& sLattice,
                               ThermalUnitConverter<T,DESCRIPTOR,TDESCRIPTOR> const& converter);
};

/// BlockLatticeThermalComfort3D returns pointwise PMV and PPD on local lattice
template <typename T, typename DESCRIPTOR, typename TDESCRIPTOR>
class BlockLatticeThermalComfort3D final : public BlockLatticeThermalPhysF3D<T,DESCRIPTOR,TDESCRIPTOR> {
public:
  BlockLatticeThermalComfort3D(BlockLattice<T,TDESCRIPTOR>& blockLattice,
                               ThermalUnitConverter<T,DESCRIPTOR,TDESCRIPTOR> const& converter);
  bool operator() (T output[], const int input[]);
private:
  /// defining the constants for one specific setup
  const T al = 1.0;               // activity level [met]
  const T rel_Hum = 50.0;         // relative humidity [%]
  const T clo  = 1.0;             // clothing value [-]
  const T mech_Work = 0.0;        // mechanical Work [W/m²]
  const T temp_Mrt = 20.4;        // mean radiant temperature [°C]
  // ATTENTION: the mean radiant temperature is hardcoded here and not calculated
  const T temp_Bod = 34.0;        // body temperature [°C]

  const T met_Rate = al * 58.15;  // activity level or metabolic rate per M/A_Du [W/m²]
  const T i_clo  = clo*0.155;     // thermal resistance of the clothing [m² °C/W]

  const int iterMax = 150;
  const T eps = 0.0001;
  const T mw = met_Rate - mech_Work;
};

}
#endif
