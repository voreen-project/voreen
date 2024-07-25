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

/* Models for Lagrangian forward-coupling methods -- generic implementation.
 */

#ifndef LB_FORWARD_COUPLING_MODELS_HH
#define LB_FORWARD_COUPLING_MODELS_HH

namespace olb {


////////////////////// Class ForwardCouplingModel ////////////////////////

template<typename T, template<typename V> class Particle>
ForwardCouplingModel<T,Particle>::ForwardCouplingModel (
  SuperGeometry<T,3>& sGeometry,
  std::shared_ptr<DragModel<T, Particle>> dragModel )
  : _sGeometry(sGeometry),
    _dragModel(dragModel)
{}


////////////////////// Class BaseForwardCouplingModel ////////////////////////

template<typename T, typename Lattice, template<typename V> class Particle>
BaseForwardCouplingModel<T,Lattice,Particle>::BaseForwardCouplingModel (
  UnitConverter<T, Lattice>& converter,
  SuperLattice<T, Lattice>& sLattice,
  SuperGeometry<T,3>& sGeometry,
  std::shared_ptr<DragModel<T, Particle>> dragModel )
  : ForwardCouplingModel<T,Particle>(sGeometry, dragModel),
    _converter(converter),
    _sLattice(sLattice)
{
  this->_interpLatticeDensity = std::make_shared<SuperLatticeInterpDensity3Degree3D<T, Lattice> > (
                                  this->_sLattice, this->_sGeometry, this->_converter );
  this->_interpLatticeVelocity = std::make_shared<SuperLatticeInterpPhysVelocity3D<T, Lattice> > (
                                   this->_sLattice, this->_converter );
}


////////////////////// Class LocalBaseForwardCouplingModel ////////////////////////

template<typename T, typename Lattice, template<typename V> class Particle>
LocalBaseForwardCouplingModel<T,Lattice,Particle>::LocalBaseForwardCouplingModel (
  UnitConverter<T, Lattice>& converter,
  SuperLattice<T, Lattice>& sLattice,
  SuperGeometry<T,3>& sGeometry,
  std::shared_ptr<DragModel<T, Particle>> dragModel )
  : BaseForwardCouplingModel<T,Lattice,Particle>(converter, sLattice, sGeometry, dragModel)
{}

template<typename T, typename Lattice, template<typename V> class Particle>
bool LocalBaseForwardCouplingModel<T,Lattice,Particle>::operator() (Particle<T>* p, int globic)
{
  /// Getting the particle and its containing cell's position
  T physPosP[3] = { p->getPos()[0],
                    p->getPos()[1],
                    p->getPos()[2]
                  }; // particle's physical position
  // particle's dimensionless position, rounded at neighbouring voxel
  int latticeRoundedPosP[3] = {0, 0, 0};
  this->_sLattice.getCuboidGeometry().get(globic).getLatticeR (
    latticeRoundedPosP, physPosP );
  // { globic, latticeRoundedP[0, ..., 2] }
  int globicFull[4] = { globic,
                        latticeRoundedPosP[0],
                        latticeRoundedPosP[1],
                        latticeRoundedPosP[2]
                      };

  // Particle's velocity
  T physVelP[3] = { p->getVel()[0],
                    p->getVel()[1],
                    p->getVel()[2]
                  }; // Physical
  // Lattice
  T latticeVelP[3] = { this->_converter.getLatticeVelocity(physVelP[0]),
                       this->_converter.getLatticeVelocity(physVelP[1]),
                       this->_converter.getLatticeVelocity(physVelP[2])
                     }; // particle's dimensionless velocity

  // Lattice's velocity at particle's location
  T physVelF[3] = {T(), T(), T()}; // Physical
  this->_interpLatticeVelocity->operator() (physVelF, physPosP, globic);
  // Lattice
  T latticeVelF[3] = { this->_converter.getLatticeVelocity(physVelF[0]),
                       this->_converter.getLatticeVelocity(physVelF[1]),
                       this->_converter.getLatticeVelocity(physVelF[2])
                     }; // Lattice's dimensionless velocity at particle's location

  // Computing fluid-particle momentum transfer
  T gF[3] = {T(), T(), T()}; // force density gF
  this->_momentumExchange->operator() (gF, latticeVelF, latticeVelP, physPosP, latticeRoundedPosP, globic);

  // Computing drag coefficient
  T Cd = this->_dragModel->operator()(p, latticeVelF, latticeVelP, globicFull);
#ifdef VERBOSE
  std:: cout << Cd
             << " physPosP=(" << physPosP[0] << ", " << physPosP[1] << ", " << physPosP[2] << ") "
             << " physVelP=(" << physVelP[0] << ", " << physVelP[1] << ", " << physVelP[2] << ") "
             << " physVelF=(" << physVelF[0] << ", " << physVelF[1] << ", " << physVelF[2] << ") "
             << std::endl;
#endif
  /*
  if (Cd > 100.)
    throw std::range_error ( "LocalBaseForwardCouplingModel::operator(). Cd="+std::to_string(Cd)
                      + "\nphysVelP=(" + std::to_string(physVelP[0]) + ", " + std::to_string(physVelP[1]) + ", " + std::to_string(physVelP[2]) + ")"
                      + "\nphysVelF=(" + std::to_string(physVelF[0]) + ", " + std::to_string(physVelF[1]) + ", " + std::to_string(physVelF[2]) + ")"
                      + "\n");
  */

  /// Computing drag force in dimensionless units
  T latticePRad = p->getRad() / this->_converter.getConversionFactorLength();
  T latticeForceP[3] = { .5 * Cd * M_PI*util::pow(latticePRad,2) * gF[0] * (latticeVelF[0] - latticeVelP[0]),
                         .5 * Cd * M_PI*util::pow(latticePRad,2) * gF[1] * (latticeVelF[1] - latticeVelP[1]),
                         .5 * Cd * M_PI*util::pow(latticePRad,2) * gF[2] * (latticeVelF[2] - latticeVelP[2])
                       }; // dimensionless force acting on the particle

  /// Computing physical drag force
  std::vector<T> physForceP(3, T()); // physical force acting on the particle
  physForceP[0] = latticeForceP[0] * this->_converter.getConversionFactorForce();
  physForceP[1] = latticeForceP[1] * this->_converter.getConversionFactorForce();
  physForceP[2] = latticeForceP[2] * this->_converter.getConversionFactorForce();

  /// Updating the particle
  p->setStoreForce(physForceP);
  return true;
}

/*
template<typename T, typename Lattice, template<typename V> class Particle>
bool LocalBaseForwardCouplingModel<T,Lattice,Particle>::operator() (SuperParticleSystem3D<T,Particle>& spSys)
{
  return true;
}
*/


////////////////////// Class NaiveForwardCouplingModel ////////////////////////

template<typename T, typename Lattice, template<typename V> class Particle>
NaiveForwardCouplingModel<T,Lattice,Particle>::NaiveForwardCouplingModel (
  UnitConverter<T, Lattice>& converter,
  SuperLattice<T, Lattice>& sLattice,
  SuperGeometry<T,3>& sGeometry,
  std::shared_ptr<DragModel<T, Particle>> dragModel )
  : LocalBaseForwardCouplingModel<T,Lattice,Particle>(converter, sLattice, sGeometry, dragModel)
{
  this->_momentumExchange = std::make_shared<NaiveMomentumExchange<T, Lattice> > (
                              this->_converter, this->_sLattice, this->_interpLatticeDensity );
}


////////////////////// Class LaddForwardCouplingModel ////////////////////////

template<typename T, typename Lattice, template<typename V> class Particle>
LaddForwardCouplingModel<T,Lattice,Particle>::LaddForwardCouplingModel (
  UnitConverter<T, Lattice>& converter,
  SuperLattice<T, Lattice>& sLattice,
  SuperGeometry<T,3>& sGeometry,
  std::shared_ptr<DragModel<T, Particle>> dragModel )
  : LocalBaseForwardCouplingModel<T,Lattice,Particle>(converter, sLattice, sGeometry, dragModel)
{
  this->_momentumExchange = std::make_shared<LaddMomentumExchange<T, Lattice> > (
                              this->_converter, this->_sLattice,
                              this->_interpLatticeDensity, this->_interpLatticeVelocity );
}


////////////////////// Class NonLocalBaseForwardCouplingModel ////////////////////////

template<typename T, typename Lattice, template<typename V> class Particle>
NonLocalBaseForwardCouplingModel<T,Lattice,Particle>::NonLocalBaseForwardCouplingModel (
  UnitConverter<T, Lattice>& converter,
  SuperLattice<T, Lattice>& sLattice,
  SuperGeometry<T,3>& sGeometry,
  std::shared_ptr<DragModel<T, Particle>> dragModel,
  std::shared_ptr<SmoothingFunctional<T, Lattice>> smoothingFunctional )
  : BaseForwardCouplingModel<T,Lattice,Particle>(converter, sLattice, sGeometry, dragModel),
    _smoothingFunctional(smoothingFunctional)
{}

template<typename T, typename Lattice, template<typename V> class Particle>
bool NonLocalBaseForwardCouplingModel<T,Lattice,Particle>::operator() (Particle<T>* p, int globic)
{
  /// Getting the particle and its containing cell's position
  T physPosP[3] = { p->getPos()[0],
                    p->getPos()[1],
                    p->getPos()[2]
                  }; // particle's physical position
  // particle's dimensionless position, rounded at neighbouring voxel
  int latticeRoundedPosP[3] = {0, 0, 0};
  this->_sLattice.getCuboidGeometry().get(globic).getLatticeR (
    latticeRoundedPosP, physPosP );
  // { globic, latticeRoundedP[0, ..., 2] }
  int globicFull[4] = { globic,
                        latticeRoundedPosP[0],
                        latticeRoundedPosP[1],
                        latticeRoundedPosP[2]
                      };

  // Particle's velocity
  T physVelP[3] = { p->getVel()[0],
                    p->getVel()[1],
                    p->getVel()[2]
                  }; // Physical
  // Lattice
  T latticeVelP[3] = { this->_converter.getLatticeVelocity(physVelP[0]),
                       this->_converter.getLatticeVelocity(physVelP[1]),
                       this->_converter.getLatticeVelocity(physVelP[2])
                     }; // particle's dimensionless velocity

  // Update the smoothing functional
  if ( ! this->_smoothingFunctional->update(physPosP, globic)) {
    std::cout << "ERROR: no lattice point enclosed in particle's kernel length!" << std::endl;
    return false;
  }

  /// Computing dimensionless drag force through smoothed average within the kernel smoothing length
  T latticeForceP[3] = {T(), T(), T()}; // dimensionless force acting on the particle
  for (auto&& i : this->_smoothingFunctional->getData()) {

    // Position of the iterated voxel
    int iLatticePosF[3] = {i.latticePos[0], i.latticePos[1], i.latticePos[2]};
    // Physical
    T iPhysPosF[3] = {T(), T(), T()};
    this->_sLattice.getCuboidGeometry().get(globic).getPhysR (
      iPhysPosF, iLatticePosF );

    // Fluid velocity at the iterated voxel
    T iPhysVelF[3] = {T(), T(), T()}; // Physical
    this->_interpLatticeVelocity->operator() (iPhysVelF, iPhysPosF, globic);
    // Lattice
    T iLatticeVelF[3] = { this->_converter.getLatticeVelocity(iPhysVelF[0]),
                          this->_converter.getLatticeVelocity(iPhysVelF[1]),
                          this->_converter.getLatticeVelocity(iPhysVelF[2])
                        }; // Lattice's dimensionless velocity at particle's location

    // Computing fluid-particle momentum transfer
    T gF[3] = {T(), T(), T()}; // force density gF
    this->_momentumExchange->operator() ( gF, iLatticeVelF, latticeVelP, physPosP, iLatticePosF, globic);

    // Computing drag coefficient
    T Cd = this->_dragModel->operator()(p, iLatticeVelF, latticeVelP, globicFull);

    /// Computing drag force in dimensionless units
    T latticePRad = p->getRad() / this->_converter.getConversionFactorLength();
    latticeForceP[0] += .5 * Cd * M_PI*util::pow(latticePRad,2) * gF[0] * (iLatticeVelF[0] - latticeVelP[0]) * i.weight;
    latticeForceP[1] += .5 * Cd * M_PI*util::pow(latticePRad,2) * gF[1] * (iLatticeVelF[1] - latticeVelP[1]) * i.weight;
    latticeForceP[2] += .5 * Cd * M_PI*util::pow(latticePRad,2) * gF[2] * (iLatticeVelF[2] - latticeVelP[2]) * i.weight;
  }

  /// Computing physical drag force
  std::vector<T> physForceP(3, T()); // physical force acting on the particle
  physForceP[0] = latticeForceP[0] * this->_converter.getConversionFactorForce();
  physForceP[1] = latticeForceP[1] * this->_converter.getConversionFactorForce();
  physForceP[2] = latticeForceP[2] * this->_converter.getConversionFactorForce();

  /// Updating the particle
  p->setStoreForce(physForceP);
  return true;
}


////////////////////// Class NaiveNonLocalForwardCouplingModel ////////////////////////

template<typename T, typename Lattice, template<typename V> class Particle>
NaiveNonLocalForwardCouplingModel<T,Lattice,Particle>::NaiveNonLocalForwardCouplingModel (
  UnitConverter<T, Lattice>& converter,
  SuperLattice<T, Lattice>& sLattice,
  SuperGeometry<T,3>& sGeometry,
  std::shared_ptr<DragModel<T, Particle>> dragModel,
  std::shared_ptr<SmoothingFunctional<T, Lattice>> smoothingFunctional )
  : NonLocalBaseForwardCouplingModel<T,Lattice,Particle>(converter, sLattice, sGeometry, dragModel, smoothingFunctional)
{
  this->_momentumExchange = std::make_shared<NaiveMomentumExchange<T, Lattice> > (
                              this->_converter, this->_sLattice, this->_interpLatticeDensity );
}


////////////////////// Class vanWachemForwardCouplingModel ////////////////////////

template<typename T, typename Lattice, template<typename V> class Particle>
vanWachemForwardCouplingModel<T,Lattice,Particle>::vanWachemForwardCouplingModel (
  UnitConverter<T, Lattice>& converter,
  SuperLattice<T, Lattice>& sLattice,
  SuperGeometry<T,3>& sGeometry,
  std::shared_ptr<DragModel<T, Particle>> dragModel,
  std::shared_ptr<SmoothingFunctional<T, Lattice>> smoothingFunctional, int nVoxelInterpPoints )
  : BaseForwardCouplingModel<T,Lattice,Particle>(converter, sLattice, sGeometry, dragModel),
    _smoothingFunctional(smoothingFunctional)
{
  this->_momentumExchange = std::make_shared<NaiveMomentumExchange<T, Lattice> > (
                              this->_converter, this->_sLattice, this->_interpLatticeDensity );
}

template<typename T, typename Lattice, template<typename V> class Particle>
bool vanWachemForwardCouplingModel<T,Lattice,Particle>::operator() (Particle<T>* p, int globic)
{
  /// Getting the particle and its containing cell's position
  T physPosP[3] = { p->getPos()[0],
                    p->getPos()[1],
                    p->getPos()[2]
                  }; // particle's physical position
  // particle's dimensionless position, rounded at neighbouring voxel
  int latticeRoundedPosP[3] = {0, 0, 0};
  this->_sLattice.getCuboidGeometry().get(globic).getLatticeR (
    latticeRoundedPosP, physPosP );
  // { globic, latticeRoundedP[0, ..., 2] }
  int globicFull[4] = { globic,
                        latticeRoundedPosP[0],
                        latticeRoundedPosP[1],
                        latticeRoundedPosP[2]
                      };

  // Particle's velocity
  T physVelP[3] = { p->getVel()[0],
                    p->getVel()[1],
                    p->getVel()[2]
                  }; // Physical
  // Lattice
  T latticeVelP[3] = { this->_converter.getLatticeVelocity(physVelP[0]),
                       this->_converter.getLatticeVelocity(physVelP[1]),
                       this->_converter.getLatticeVelocity(physVelP[2])
                     }; // particle's dimensionless velocity

  // Lattice's velocity at particle's location
  T physVelF[3] = {T(), T(), T()}; // Physical
  this->_interpLatticeVelocity->operator() (physVelF, physPosP, globic);
  // Lattice
  T latticeVelF[3] = { this->_converter.getLatticeVelocity(physVelF[0]),
                       this->_converter.getLatticeVelocity(physVelF[1]),
                       this->_converter.getLatticeVelocity(physVelF[2])
                     }; // Lattice's dimensionless velocity at particle's location

  // Computing fluid-particle momentum transfer
  T gF[3] = {T(), T(), T()}; // force density gF
  this->_momentumExchange->operator() (gF, latticeVelF, latticeVelP, physPosP, latticeRoundedPosP, globic);

  // Update the smoothing functional
  if ( ! this->_smoothingFunctional->update(physPosP, globic)) {
    std::cout << "ERROR: no lattice point enclosed in particle's kernel length!" << std::endl;
    return false;
  }

  // Computing continuous phase fraction across kernel
  T continuousPhaseFraction = T();
  for (auto&& i : this->_smoothingFunctional->getData()) {
    continuousPhaseFraction += i.continuousPhaseFraction * i.weight;
  }

  // Computing drag coefficient
  T Cd = this->_dragModel->operator()(p, latticeVelF, latticeVelP, globicFull, continuousPhaseFraction);
#ifdef VERBOSE
  std:: cout << Cd
             << " physPosP=(" << physPosP[0] << ", " << physPosP[1] << ", " << physPosP[2] << ") "
             << " physVelP=(" << physVelP[0] << ", " << physVelP[1] << ", " << physVelP[2] << ") "
             << " physVelF=(" << physVelF[0] << ", " << physVelF[1] << ", " << physVelF[2] << ") "
             << std::endl;
#endif
  /*
  if (Cd > 100.)
    throw std::range_error ( "LocalBaseForwardCouplingModel::operator(). Cd="+std::to_string(Cd)
                      + "\nphysVelP=(" + std::to_string(physVelP[0]) + ", " + std::to_string(physVelP[1]) + ", " + std::to_string(physVelP[2]) + ")"
                      + "\nphysVelF=(" + std::to_string(physVelF[0]) + ", " + std::to_string(physVelF[1]) + ", " + std::to_string(physVelF[2]) + ")"
                      + "\n");
  */

  /// Computing drag force in dimensionless units
  T latticePRad = p->getRad() / this->_converter.getConversionFactorLength();
  T latticeForceP[3] = { .5 * Cd * M_PI*util::pow(latticePRad,2) * gF[0] * (latticeVelF[0] - latticeVelP[0]) / util::pow(continuousPhaseFraction, 1.65),
                         .5 * Cd * M_PI*util::pow(latticePRad,2) * gF[1] * (latticeVelF[1] - latticeVelP[1]) / util::pow(continuousPhaseFraction, 1.65),
                         .5 * Cd * M_PI*util::pow(latticePRad,2) * gF[2] * (latticeVelF[2] - latticeVelP[2]) / util::pow(continuousPhaseFraction, 1.65)
                       };

  /// Computing physical drag force
  std::vector<T> physForceP(3, T()); // physical force acting on the particle
  physForceP[0] = latticeForceP[0] * this->_converter.getConversionFactorForce();
  physForceP[1] = latticeForceP[1] * this->_converter.getConversionFactorForce();
  physForceP[2] = latticeForceP[2] * this->_converter.getConversionFactorForce();

  /// Updating the particle
  p->setStoreForce(physForceP);
  return true;
}

}

#endif
