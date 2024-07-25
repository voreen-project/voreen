/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2019 Davide Dapelo
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

/* Models for Lagrangian back-coupling methods -- generic implementation.
 */

#ifndef LB_BACK_COUPLING_MODELS_HH
#define LB_BACK_COUPLING_MODELS_HH

namespace olb {


////////////////////// Class BaseBackCouplingModel ////////////////////////

template<typename T, typename Lattice, template<typename V> class Particle>
BaseBackCouplingModel<T,Lattice,Particle>::BaseBackCouplingModel (
  UnitConverter<T, Lattice>& converter,
  SuperLattice<T, Lattice>& sLattice,
  SuperGeometry<T,3>& sGeometry,
  int overlap )
  : _converter(converter),
    _sGeometry(sGeometry),
    _sLattice(sLattice),
    _commPopulation(_sLattice)
{
  _commPopulation.template requestField<olb::descriptors::POPULATION>();
  _commPopulation.requestOverlap(overlap);
  _commPopulation.exchangeRequests();

  _zeroAnalytical = std::make_shared<AnalyticalConst3D<T, T> > (T());
  _zeroField = std::make_shared<AnalyticalComposed3D<T, T> > (*_zeroAnalytical, *_zeroAnalytical, *_zeroAnalytical);
}

template<typename T, typename Lattice, template<typename V> class Particle>
void BaseBackCouplingModel<T,Lattice,Particle>::resetExternalField(int material)
{
  // resets external field
  this->_sLattice.template defineField<descriptors::FORCE>(this->_sGeometry, material, *_zeroField);

  // NECESSARY to communicate values before using them in operator()
  _commPopulation.communicate();
}


////////////////////// Class BaseLocalBackCouplingModel ////////////////////////

template<typename T, typename Lattice, template<typename V> class Particle>
BaseLocalBackCouplingModel<T,Lattice,Particle>::BaseLocalBackCouplingModel (
  UnitConverter<T, Lattice>& converter,
  SuperLattice<T, Lattice>& sLattice,
  SuperGeometry<T,3>& sGeometry,
  int overlap )
  : BaseBackCouplingModel<T,Lattice,Particle>(converter, sLattice, sGeometry, overlap)
{}

template<typename T, typename Lattice, template<typename V> class Particle>
void BaseLocalBackCouplingModel<T,Lattice,Particle>::communicate()
{
  this->_commPopulation.communicate();
}


////////////////////// Class BaseNonLocalBackCouplingModel ////////////////////////

template<typename T, typename Lattice, template<typename V> class Particle>
BaseNonLocalBackCouplingModel<T,Lattice,Particle>::BaseNonLocalBackCouplingModel (
  UnitConverter<T, Lattice>& converter,
  SuperLattice<T, Lattice>& sLattice,
  SuperGeometry<T,3>& sGeometry,
  int overlap )
  : BaseBackCouplingModel<T,Lattice,Particle>(converter, sLattice, sGeometry, overlap),
    _commForce(sLattice)
{
  _commForce.template requestField<olb::descriptors::FORCE>();
  _commForce.requestOverlap(overlap);
  _commForce.exchangeRequests();
}

template<typename T, typename Lattice, template<typename V> class Particle>
void BaseNonLocalBackCouplingModel<T,Lattice,Particle>::communicate()
{
  this->_commPopulation.communicate();
  _commForce.communicate();
}


////////////////////// Class LocalBackCouplingModel ////////////////////////

template<typename T, typename Lattice, template<typename V> class Particle>
LocalBackCouplingModel<T,Lattice,Particle>::LocalBackCouplingModel (
  UnitConverter<T, Lattice>& converter,
  SuperLattice<T, Lattice>& sLattice,
  SuperGeometry<T,3>& sGeometry,
  int overlap )
  : BaseLocalBackCouplingModel<T,Lattice,Particle>(converter, sLattice, sGeometry, overlap)
{}

template<typename T, typename Lattice, template<typename V> class Particle>
bool LocalBackCouplingModel<T,Lattice,Particle>::operator() (Particle<T>* p, int globic, int material, int subCycles)
{
  int locIC = this->_sLattice.getLoadBalancer().loc(globic);

  // reading the force from the value stored inside the particle
  std::vector<T> physForceP = p->getStoreForce(); // physical force acting on the particle

  T latticeForceP[3] = { physForceP[0] / this->_converter.getConversionFactorForce(),
                         physForceP[1] / this->_converter.getConversionFactorForce(),
                         physForceP[2] / this->_converter.getConversionFactorForce()
                       }; // dimensionless force acting on the particle

  T physPosP[3] = { (p->getPos()[0]),
                    (p->getPos()[1]),
                    (p->getPos()[2])
                  }; // particle's physical position

  // particle's dimensionless position, rounded at neighbouring voxel
  int latticeRoundedPosP[3] = {0, 0, 0};
  this->_sLattice.getCuboidGeometry().get(globic).getLatticeR( latticeRoundedPosP, physPosP );

  if (this->_sGeometry.getBlockGeometry(locIC).getMaterial( latticeRoundedPosP[0], latticeRoundedPosP[1], latticeRoundedPosP[2] ) == material) {
    FieldD<T,Lattice,descriptors::FORCE> F{
      -latticeForceP[0] / (T)(subCycles),
      -latticeForceP[1] / (T)(subCycles),
      -latticeForceP[2] / (T)(subCycles)
    };  // dimensionless smoothed force

    this->_sLattice.getBlock(locIC).get(
      latticeRoundedPosP[0],
      latticeRoundedPosP[1],
      latticeRoundedPosP[2]
    ).template getFieldPointer<descriptors::FORCE>() += F;
  }

  return true;
}


////////////////////// Class CubicDeltaBackCouplingModel ////////////////////////

template<typename T, typename Lattice, template<typename V> class Particle>
CubicDeltaBackCouplingModel<T,Lattice,Particle>::CubicDeltaBackCouplingModel (
  UnitConverter<T, Lattice>& converter,
  SuperLattice<T, Lattice>& sLattice,
  SuperGeometry<T,3>& sGeometry,
  int overlap )
  : BaseNonLocalBackCouplingModel<T,Lattice,Particle>(converter, sLattice, sGeometry, overlap)
{
  _cubicDeltaFunctional = std::make_shared<SuperLatticeSmoothDiracDelta3D<T, Lattice> > (
                            this->_sLattice, this->_converter, this->_sGeometry );
}

template<typename T, typename Lattice, template<typename V> class Particle>
bool CubicDeltaBackCouplingModel<T,Lattice,Particle>::operator() (Particle<T>* p, int globic, int material, int subCycles)
{
  int locIC = this->_sLattice.getLoadBalancer().loc(globic);

  // reading the force from the value stored inside the particle
  std::vector<T> physForceP = p->getStoreForce(); // physical force acting on the particle

  T latticeForceP[3] = { physForceP[0] / this->_converter.getConversionFactorForce(),
                         physForceP[1] / this->_converter.getConversionFactorForce(),
                         physForceP[2] / this->_converter.getConversionFactorForce()
                       }; // dimensionless force acting on the particle

  T physPosP[3] = { p->getPos()[0],
                    p->getPos()[1],
                    p->getPos()[2]
                  }; // particle's physical position

  // particle's dimensionless position, rounded at neighbouring voxel
  int latticeRoundedPosP[3] = {0, 0, 0};
  this->_sLattice.getCuboidGeometry().get(globic).getLatticeR( latticeRoundedPosP, physPosP );

  // smooth Dirac delta
  this->_cubicDeltaFunctional->operator() (_delta, physPosP, globic);

  T tempDelta = T();
  FieldD<T,Lattice,descriptors::FORCE> F; // dimensionless smoothed force

  for (int i = -_range; i <= _range; ++i) {
    for (int j = -_range; j <= _range; ++j) {
      for (int k = -_range; k <= _range; ++k) {
        if (this->_sGeometry.getBlockGeometry(locIC).getMaterial(
              latticeRoundedPosP[0] + i, latticeRoundedPosP[1] + j,
              latticeRoundedPosP[2] + k) == material) {

          tempDelta = _delta[i + _range][j + _range][k + _range];

          F[0] = -latticeForceP[0] * tempDelta / (T)(subCycles);
          F[1] = -latticeForceP[1] * tempDelta / (T)(subCycles);
          F[2] = -latticeForceP[2] * tempDelta / (T)(subCycles);

          this->_sLattice.getBlock(locIC).get(
            latticeRoundedPosP[0] + i,
            latticeRoundedPosP[1] + j,
            latticeRoundedPosP[2] + k
          ).template getFieldPointer<descriptors::FORCE>() += F;
        }
      }
    }
  }
  return true;
}


////////////////////// Class NonLocalBaseBackCouplingModel ////////////////////////

template<typename T, typename Lattice, template<typename V> class Particle>
NonLocalBaseBackCouplingModel<T,Lattice,Particle>::NonLocalBaseBackCouplingModel (
  UnitConverter<T, Lattice>& converter,
  SuperLattice<T, Lattice>& sLattice,
  SuperGeometry<T,3>& sGeometry,
  std::shared_ptr<SmoothingFunctional<T, Lattice>> smoothingFunctional,
  int overlap )
  : BaseNonLocalBackCouplingModel<T,Lattice,Particle>(converter, sLattice, sGeometry, overlap),
    _smoothingFunctional(smoothingFunctional)
{}

template<typename T, typename Lattice, template<typename V> class Particle>
bool NonLocalBaseBackCouplingModel<T,Lattice,Particle>::operator() (Particle<T>* p, int globic, int material, int subCycles)
{
  int locIC = this->_sLattice.getLoadBalancer().loc(globic);

  // reading the force from the value stored inside the particle
  std::vector<T> physForceP = p->getStoreForce(); // physical force acting on the particle
  T latticeForceP[3] = { physForceP[0] / this->_converter.getConversionFactorForce(),
                         physForceP[1] / this->_converter.getConversionFactorForce(),
                         physForceP[2] / this->_converter.getConversionFactorForce()
                       }; // dimensionless force acting on the particle

  // Updating force through voxels within kernel smoothing length from the bubble's position
  for (auto&& i : this->_smoothingFunctional->getData()) {

    // Updating iterated voxel
    if ( this->_sGeometry.getBlockGeometry(locIC).getMaterial(i.latticePos[0], i.latticePos[1], i.latticePos[2]) == material ) {

      // Weighted force acting on the iterated voxel
      FieldD<T,Lattice,descriptors::FORCE> F{
        -latticeForceP[0] * i.weight / (T)(subCycles),
        -latticeForceP[1] * i.weight / (T)(subCycles),
        -latticeForceP[2] * i.weight / (T)(subCycles)
      }; // dimensionless smoothed force

      this->_sLattice.getBlock(locIC).get(
        i.latticePos[0], i.latticePos[1], i.latticePos[2]
      ).template getFieldPointer<descriptors::FORCE>() += F;
    }
  }
  return true;
}

}

#endif
