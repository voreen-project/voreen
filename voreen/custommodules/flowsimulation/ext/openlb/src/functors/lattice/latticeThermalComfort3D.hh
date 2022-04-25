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

#ifndef LATTICE_THERMAL_COMFORT_3D_HH
#define LATTICE_THERMAL_COMFORT_3D_HH

#include<vector>    // for generic i/o
#include<cmath>     // for lpnorm
#include<math.h>

#include "latticeThermalComfort3D.h"
#include "superBaseF3D.h"
#include "functors/analytical/indicator/indicatorBaseF3D.h"
#include "indicator/superIndicatorF3D.h"
#include "dynamics/lbm.h"  // for computation of lattice rho and velocity
#include "geometry/superGeometry.h"
#include "blockBaseF3D.h"
#include "communication/mpiManager.h"
#include "utilities/vectorHelpers.h"

namespace olb {

template <typename T, typename DESCRIPTOR, typename TDESCRIPTOR>
SuperLatticeThermalComfort3D<T,DESCRIPTOR,TDESCRIPTOR>::SuperLatticeThermalComfort3D(
  SuperLattice<T,TDESCRIPTOR>& sLattice, ThermalUnitConverter<T,DESCRIPTOR,TDESCRIPTOR> const& converter)
  : SuperLatticeThermalPhysF3D<T,DESCRIPTOR,TDESCRIPTOR>(sLattice, converter, 2)
{
  this->getName() = "thermalComfort";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back(new BlockLatticeThermalComfort3D<T,DESCRIPTOR,TDESCRIPTOR>(this->_sLattice.getBlock(iC), this->_converter));
  }
}

template <typename T, typename DESCRIPTOR, typename TDESCRIPTOR>
BlockLatticeThermalComfort3D<T,DESCRIPTOR,TDESCRIPTOR>::BlockLatticeThermalComfort3D
(BlockLattice<T,TDESCRIPTOR>& blockLattice, ThermalUnitConverter<T,DESCRIPTOR,TDESCRIPTOR> const& converter)
  : BlockLatticeThermalPhysF3D<T,DESCRIPTOR,TDESCRIPTOR>(blockLattice,converter,2)
{
  this->getName() = "thermalComfort";
}


template <typename T, typename DESCRIPTOR, typename TDESCRIPTOR>
bool BlockLatticeThermalComfort3D<T,DESCRIPTOR,TDESCRIPTOR>::operator() (T output[], const int input[])
{
  OstreamManager clout(std::cout,"thermalComfort");

  // importing the local velocity and temperature for one specific lattice
  T latticeTemp = this->_blockLattice.get( input[0], input[1], input[2]).computeRho();
  T physTemp = this->_converter.getPhysTemperature(latticeTemp);

  auto latticeVel = this->_blockLattice.get(input[0], input[1], input[2]).template getFieldPointer<descriptors::VELOCITY>();
  T physVel[3];
  physVel[0]= this->_converter.getPhysVelocity(latticeVel[0]);
  physVel[1]= this->_converter.getPhysVelocity(latticeVel[1]);
  physVel[2]= this->_converter.getPhysVelocity(latticeVel[2]);

  T pmv;   // Predicted Mean Vote
  T ppd;   // Predicted Percantage of Dissatisfied

  T temp_Air;                // air temperature [°C]
  T vel_Air;                 // local relative air velocity [m/s]
  T press_Vapor;             // pressure of water vapour in ambient air [Pa]
  T h_c;                     // convective heat transfer coefficient [W/(m² K)]
  T f_clo;                   // ratio of surface area clothed/nude [-]
  T temp_Clo = temp_Bod;     // clothing temperature [°C]
  T temp_Clo_old = temp_Clo;

  // calculation of the air temperature in °C and the magnitude of the velocity
  temp_Air = physTemp - 273.15;
  vel_Air = util::sqrt(physVel[0]*physVel[0] + physVel[1]*physVel[1] + physVel[2]*physVel[2]);

  /// calculation of f_clo
  // ratio of surface area clothed/nude [-]
  // 1.15 = typical business suit
  if (i_clo <= 0.078) {
    f_clo = 1.0 + 1.290 * i_clo;
  }
  else {
    f_clo = 1.05 + 0.645 * i_clo;
  }

  /// calculation of press_Vapor with multiple equations to choose from
  // however, there is no significant difference between these equations
  // calculation with simplified equation from a comparable algorithm
  T ppk = 673.4 - 1.8 * temp_Air;
  T ppa = (3.2437814 + 0.00326014*ppk) + (2.00658 * 0.000000001 * ppk * ppk * ppk);
  T ppb = (1165.09 - ppk) * (1.0 + 0.00121547 * ppk);
  press_Vapor = rel_Hum/100.0 * 22105.8416/util::exp(2.302585 * ppk * ppa / ppb) * 1000.0;

  // calculation like EN ISO 7730
  //press_Vapor = rel_Hum * 10 * util::exp(16.6536 - 4030.183/(temp_Air+235.0));
  // calculation with Magnus Equation
  //press_Vapor = rel_Hum/100.0 * 6.112 * 100 * util::exp((17.62*temp_Air)/(243.12+temp_Air));

  /// iterative calculation of h_c and t_cl
  for (int iter = 0; iter < iterMax; iter++) {
    temp_Clo = 0.8 * temp_Clo_old + 0.2 * temp_Clo;

    h_c = 12.1 * util::sqrt(vel_Air);                         // pure forced convection

    if (2.38 * util::pow(util::fabs(temp_Clo - temp_Air), 0.25) > h_c) {
      h_c = 2.38 * util::pow(util::fabs(temp_Clo - temp_Air), 0.25);      // pure free convection
    }

    temp_Clo_old = temp_Clo;
    temp_Clo = 35.7 - 0.028 * mw -
               i_clo * (3.96e-8 * f_clo * (util::pow(temp_Clo + 273.0, 4) - util::pow(temp_Mrt + 273.0, 4)) +
                        f_clo * h_c * (temp_Clo - temp_Air));

    if ((util::fabs((temp_Clo_old - temp_Clo) / temp_Clo)) < eps) {
      break;
    }
    if (iter == (iterMax-1)) {
      clout << "WARNING: Maximal iteration step limit: PMV calculation not possible!" << std::endl;
    }
  }

  // calculation of PMV
  pmv = (0.303 * util::exp(-0.036 * mw) + 0.028) * (mw                                    // heat gain by internal metabolic process
        -3.05e-3*(5733.0 - 6.99*mw - press_Vapor)                                  // heat loss by skin diffusion
        -0.42*(mw - 58.15)                                                         // heat loss by evaporation of sweat secretion
        -1.7e-5 * mw * (5867.0 - press_Vapor) - 0.0014 * mw * (34.0 - temp_Air)    // heat loss by latent respiration and dry respiration
        -3.96e-8 * f_clo *(util::pow(temp_Clo + 273.0, 4) - util::pow(temp_Mrt + 273.0, 4))    // heat loss by radiation
        -f_clo * h_c * (temp_Clo - temp_Air));                                     // heat loss by convection

  // calculation of PPD
  ppd = 100.0 - 95.0*util::exp(-0.03353* util::pow(pmv,4.0) - 0.2179*util::pow(pmv,2.0));

  // return of PMV and PPD
  output[0] = pmv;
  output[1] = ppd;

  return true;
}

}
#endif
