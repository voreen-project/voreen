/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2017 Max Gaedtke, Albert Mink
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

/** \file
 * Unit conversion handling -- header file.
 */

#ifndef RADIATIVEUNITCONVERTER_H
#define RADIATIVEUNITCONVERTER_H


#include "io/ostreamManager.h"
#include "core/unitConverter.h"


/// All OpenLB code is contained in this namespace.
namespace olb {

double getThetaRefracted(double const& thetaIncident, double const& refractiveRelative);
double getFresnelFunction(double const& theta, double const& refractiveRelative);
double R_phi_diff (double const& theta, double const& refractiveRelative);
double R_j_diff (double const& theta, double const& refractiveRelative);
double getRefractionFunction(const double& refractiveRelative);
double getRefractionFunction(const double& refractiveRelative);
double getPartialBBCoefficient(double const& latticeDiffusionCoefficient, double const& relativeRefractiveIndex );

// wrapper for above function
template <typename T, template<typename U> class Lattice>
double getPartialBBCoefficient(RadiativeUnitConverter<T,Lattice> const& converter)
{
  return getPartialBBCoefficient( converter.getLatticeDiffusionCoefficient(), converter.getRefractiveRelative() );
};

// wrapper for above function
template <typename T, template<typename U> class Lattice>
double getRefractionFunction(RadiativeUnitConverter<T,Lattice> const& converter)
{
  return getRefractionFunction(converter.getRefractiveRelative());
};

/** Conversion between physical and lattice units, as well as discretization.
* Be aware of the nomenclature:
* We distingish between physical (dimensioned) and lattice (dimensionless) values.
* A specific conversion factor maps the two different scopes,
* e.g. __physLength = conversionLength * latticeLength__
*
*/
template <typename T, template<typename U> class Lattice>
class RadiativeUnitConverter : public UnitConverterFromResolutionAndRelaxationTime<T,Lattice> {
public:
  /** Documentation of constructor:
    *   \param resolution   is number of voxel per 1 meter
    *   \param latticeRelaxationTime    see class UnitConverterFromResolutionAndRelaxationTime
    *   \param physAbsorption   physical absorption in 1/meter
    *   \param physScattering   physical scattering in 1/meter
    */
  constexpr RadiativeUnitConverter( int resolution, T latticeRelaxationTime, T physAbsorption, T physScattering, T refractiveMedia=1, T refractiveAmbient=1 )
    : UnitConverterFromResolutionAndRelaxationTime<T, Lattice>( resolution, latticeRelaxationTime, T(1), T(1), T(1), T(1) ),
      clout(std::cout, "RadiativeUnitConverter"),
      _physAbsorption(physAbsorption),
      _physScattering(physScattering),
      _extinctionCoefficient( physAbsorption+physScattering ),
      _scatteringAlbedo( physScattering/(physAbsorption+physScattering) ),
      _physDiffusionCoefficient( 1.0 / (3.0*(physAbsorption+physScattering)) ),
      _refractiveRelative(refractiveMedia/refractiveAmbient),
      _latticeAbsorption( physAbsorption/double(resolution) ),
      _latticeScattering( physScattering/double(resolution) ),
      _latticeDiffusionCoefficient(_physDiffusionCoefficient*resolution)
    { };

  constexpr T getPhysAbsorption() const {
    return _physAbsorption;
  };

  constexpr T getPhysScattering() const {
    return _physScattering;
  };

  constexpr T getExtinctionCoefficient() const {
    return _extinctionCoefficient;
  };

  constexpr T getScatteringAlbedo() const {
    return _scatteringAlbedo;
  };

  constexpr T getPhysDiffusionCoefficient() const {
    return _physDiffusionCoefficient;
  };

  constexpr T getLatticeAbsorption() const {
    return _latticeAbsorption;
  };

  constexpr T getLatticeScattering() const {
    return _latticeScattering;
  };

  constexpr T getLatticeDiffusionCoefficient() const {
    return _latticeDiffusionCoefficient;
  };

  constexpr T getRefractiveRelative() const {
    return _refractiveRelative;
  };

  void print() const override;

private:
  mutable OstreamManager clout;

  double _physAbsorption;
  double _physScattering;
  double _extinctionCoefficient;
  double _scatteringAlbedo;
  double _physDiffusionCoefficient;

  double _refractiveRelative;

  double _latticeAbsorption;
  double _latticeScattering;
  double _latticeDiffusionCoefficient;
};

template <typename T, template<typename U> class Lattice>
void RadiativeUnitConverter<T, Lattice>::print() const
{
  clout << "----------------- UnitConverter information -----------------" << std::endl;
  clout << "-- Parameters:" << std::endl;
  clout << "Resolution:                       N=              " << this->getResolution() << std::endl;
  clout << "Lattice relaxation frequency:     omega=          " << this->getLatticeRelaxationFrequency(  ) << std::endl;
  clout << "Lattice relaxation time:          tau=            " << this->getLatticeRelaxationTime() << std::endl;
  clout << "Characteristical length(m):       charL=          " << this->getCharPhysLength() << std::endl;
  clout << "Phys. Density(kg/m^d):            charRho=        " << this->getPhysDensity() << std::endl;
  clout << "Physical absorption:              absorption=     " << getPhysAbsorption() << std::endl;
  clout << "Physical scattering:              scattering=     " << getPhysScattering() << std::endl;
  clout << "Extinction coefficient:           extinction=     " << getExtinctionCoefficient() << std::endl;
  clout << "Singl scattering albedo:          albedo=         " << getScatteringAlbedo() << std::endl;
  clout << "Physical diffusion coefficient:   D=              " << getPhysDiffusionCoefficient() << std::endl;
  clout << "Singls scattering albedo:         albedo=         " << getScatteringAlbedo() << std::endl;

  clout << std::endl;
  clout << "Lattice diffusion coefficient:    D^*=            " << getLatticeDiffusionCoefficient() << std::endl;
  clout << "C_R: " << getRefractionFunction(getRefractiveRelative()) << std::endl;
  clout << "r_F: " << getPartialBBCoefficient(getLatticeDiffusionCoefficient(),getRefractiveRelative()) << std::endl;

  clout << std::endl;
  clout << "-- Conversion factors:" << std::endl;
  clout << "Voxel length(m):                  physDeltaX=     " << this->getConversionFactorLength() << std::endl;
  clout << "Time step(s):                     physDeltaT=     " << this->getConversionFactorTime() << std::endl;
  clout << "Density factor(kg/m^3):           physDensity=    " << this->getConversionFactorDensity() <<  std::endl;
  clout << "-------------------------------------------------------------" << std::endl;

}


}  // namespace olb

#endif
