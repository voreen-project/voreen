/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012 Lukas Baron, Mathias J. Krause, Albert Mink
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

#ifndef BLOCK_LATTICE_LOCAL_F_2D_HH
#define BLOCK_LATTICE_LOCAL_F_2D_HH

#include<vector>
#include<cmath>

#include "blockLatticeLocalF2D.h"
#include "blockBaseF2D.h"
#include "functors/genericF.h"
#include "functors/analytical/analyticalF.h"
#include "indicator/indicatorF2D.h"
#include "core/blockLattice2D.h"
#include "communication/mpiManager.h"
#include "core/blockLatticeStructure2D.h"
#include "dynamics/lbHelpers.h"  // for computation of lattice rho and velocity

namespace olb {


template <typename T, template <typename U> class DESCRIPTOR>
BlockLatticeDissipation2D<T,DESCRIPTOR>::BlockLatticeDissipation2D
(BlockLatticeStructure2D<T,DESCRIPTOR>& blockLattice, const UnitConverter<T,DESCRIPTOR>& converter)
  : BlockLatticeF2D<T,DESCRIPTOR>(blockLattice,1), _converter(converter)
{
  this->getName() = "dissipation";
}

// todo: get functor working
template <typename T, template <typename U> class DESCRIPTOR>
bool BlockLatticeDissipation2D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  //  int globIC = input[0];
  //  int locix= input[1];
  //  int lociy= input[2];
  //  int lociz= input[3];
  //  if ( this->_blockLattice.get_load().rank(globIC) == singleton::mpi().getRank() ) {
  //    // local coords are given, fetch local cell and compute value(s)
  //    T rho, uTemp[DESCRIPTOR<T>::d], pi[util::TensorVal<DESCRIPTOR<T> >::n];
  //    int overlap = this->_blockLattice.getOverlap();
  //    this->_blockLattice.getExtendedBlockLattice(this->_blockLattice.get_load().loc(globIC)).get(locix+overlap, lociy+overlap, lociz+overlap).computeAllMomenta(rho, uTemp, pi);

  //    T PiNeqNormSqr = pi[0]*pi[0] + 2.*pi[1]*pi[1] + pi[2]*pi[2];
  //    if (util::TensorVal<DESCRIPTOR<T> >::n == 6)
  //      PiNeqNormSqr += pi[2]*pi[2] + pi[3]*pi[3] + 2.*pi[4]*pi[4] +pi[5]*pi[5];

  //    T nuLattice = converter.getLatticeNu();
  //    T omega = converter.getOmega();
  //    T finalResult = PiNeqNormSqr*nuLattice*pow(omega*DESCRIPTOR<T>::invCs2,2)/rho/2.;

  //    return std::vector<T>(1,finalResult);
  //  } else {
  //    return std::vector<T>(); // empty vector
  //  }

  return false;
}

template <typename T, template <typename U> class DESCRIPTOR>
BlockLatticePhysDissipation2D<T,DESCRIPTOR>::BlockLatticePhysDissipation2D
(BlockLatticeStructure2D<T,DESCRIPTOR>& blockLattice, const UnitConverter<T,DESCRIPTOR>& converter)
  : BlockLatticePhysF2D<T,DESCRIPTOR>(blockLattice,converter,1)
{
  this->getName() = "physDissipation";
}

template <typename T, template <typename U> class DESCRIPTOR>
bool BlockLatticePhysDissipation2D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  T rho, uTemp[DESCRIPTOR<T>::d], pi[util::TensorVal<DESCRIPTOR<T> >::n];
  this->_blockLattice.get( input[0], input[1] ).computeAllMomenta(rho, uTemp, pi);

  T PiNeqNormSqr = pi[0]*pi[0] + 2.*pi[1]*pi[1] + pi[2]*pi[2];
  if (util::TensorVal<DESCRIPTOR<T> >::n == 6) {
    PiNeqNormSqr += pi[2]*pi[2] + pi[3]*pi[3] + 2.*pi[4]*pi[4] + pi[5]*pi[5];
  }

  T nuLattice = this->_converter.getLatticeViscosity();
  T omega = 1. / this->_converter.getLatticeRelaxationTime();
  T dt = this->_converter.getConversionFactorTime();
  output[0] = PiNeqNormSqr*nuLattice*pow(omega*DESCRIPTOR<T>::invCs2/rho,2)/2.*this->_converter.getPhysViscosity()/this->_converter.getLatticeViscosity()/dt/dt;

  return true;
}

template <typename T, template <typename U> class DESCRIPTOR>
BlockLatticeDensity2D<T,DESCRIPTOR>::BlockLatticeDensity2D
(BlockLatticeStructure2D<T,DESCRIPTOR>& blockLattice)
  : BlockLatticeF2D<T,DESCRIPTOR>(blockLattice,1)
{
  this->getName() = "density";
}


template <typename T, template <typename U> class DESCRIPTOR>
bool BlockLatticeDensity2D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  output[0] = this->_blockLattice.get( input[0] , input[1] ).computeRho();
  return true;
}


template <typename T, template <typename U> class DESCRIPTOR>
BlockLatticeVelocity2D<T,DESCRIPTOR>::BlockLatticeVelocity2D
(BlockLatticeStructure2D<T,DESCRIPTOR>& blockLattice)
  : BlockLatticeF2D<T,DESCRIPTOR>(blockLattice,2)
{
  this->getName() = "velocity";
}

template <typename T, template <typename U> class DESCRIPTOR>
bool BlockLatticeVelocity2D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  T rho;
  this->_blockLattice.get(input[0], input[1]).computeRhoU(rho,output);
  return true;
}

/*template <typename T, template <typename U> class DESCRIPTOR>
BlockLatticeStrainRate2D<T,DESCRIPTOR>::BlockLatticeStrainRate2D
  (BlockLatticeStructure2D<T,DESCRIPTOR>& blockLattice, const UnitConverter<T>& converter)
  : BlockLatticePhysF2D<T,DESCRIPTOR>(blockLattice,converter,4)
{ this->getName() = "strainRate"; }

template <typename T, template <typename U> class DESCRIPTOR>
std::vector<T> BlockLatticeStrainRate2D<T,DESCRIPTOR>::operator() (std::vector<int> input) {

  T strainRate[4];
  T rho, uTemp[DESCRIPTOR<T>::d], pi[util::TensorVal<DESCRIPTOR<T> >::n];
  this->_blockLattice.get( input[0], input[1] ).computeAllMomenta(rho, uTemp, pi);

  T omega = this->_converter.getOmega();

  strainRate[0] = -pi[0]*omega*DESCRIPTOR<T>::invCs2/rho/2.;
  strainRate[1] = -pi[1]*omega*DESCRIPTOR<T>::invCs2/rho/2.;
  strainRate[2] = -pi[1]*omega*DESCRIPTOR<T>::invCs2/rho/2.;
  strainRate[3] = -pi[2]*omega*DESCRIPTOR<T>::invCs2/rho/2.;

  //cout << pi[0] << " " << pi[1] << " " << pi[2] << " " << DESCRIPTOR<T>::invCs2 << endl;

  std::vector<T> output(strainRate, strainRate+4); // first adress, last adress
  return output;
}*/

template <typename T, template <typename U> class DESCRIPTOR>
BlockLatticePhysStrainRate2D<T,DESCRIPTOR>::BlockLatticePhysStrainRate2D
(BlockLatticeStructure2D<T,DESCRIPTOR>& blockLattice, const UnitConverter<T,DESCRIPTOR>& converter)
  : BlockLatticePhysF2D<T,DESCRIPTOR>(blockLattice,converter,4)
{
  this->getName() = "strainRate";
}

template <typename T, template <typename U> class DESCRIPTOR>
bool BlockLatticePhysStrainRate2D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  T rho, uTemp[DESCRIPTOR<T>::d], pi[util::TensorVal<DESCRIPTOR<T> >::n];
  this->_blockLattice.get( input[0], input[1] ).computeAllMomenta(rho, uTemp, pi);

  T omega = 1. / this->_converter.getLatticeRelaxationTime();
  T dt = this->_converter.getConversionFactorTime();

  output[0] = -pi[0]*omega*DESCRIPTOR<T>::invCs2/rho/2./dt;
  output[1] = -pi[1]*omega*DESCRIPTOR<T>::invCs2/rho/2./dt;
  output[2] = -pi[1]*omega*DESCRIPTOR<T>::invCs2/rho/2./dt;
  output[3] = -pi[2]*omega*DESCRIPTOR<T>::invCs2/rho/2./dt;

  return true;
}

template <typename T, template <typename U> class DESCRIPTOR>
BlockLatticeGeometry2D<T,DESCRIPTOR>::BlockLatticeGeometry2D
(BlockLatticeStructure2D<T,DESCRIPTOR>& blockLattice, BlockGeometryStructure2D<T>& blockGeometry, int material)
  : BlockLatticeF2D<T,DESCRIPTOR>(blockLattice,1), _blockGeometry(blockGeometry), _material(material)
{
  this->getName() = "geometry";
}

template <typename T, template <typename U> class DESCRIPTOR>
bool BlockLatticeGeometry2D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  int materialTmp = _blockGeometry.getMaterial( input[0], input[1] );

  if (_material != -1) {
    if (_material == materialTmp) {
      output[0] = T(1);
      return true;
    } else {
      output[0] = T();
      return true;
    }
  }
  output[0]=T(materialTmp);
  return false;
}


template <typename T, template <typename U> class DESCRIPTOR>
BlockLatticeRank2D<T,DESCRIPTOR>::BlockLatticeRank2D
(BlockLatticeStructure2D<T,DESCRIPTOR>& blockLattice)
  : BlockLatticeF2D<T,DESCRIPTOR>(blockLattice,1)
{
  this->getName() = "rank";
}

template <typename T, template <typename U> class DESCRIPTOR>
bool BlockLatticeRank2D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  output[0] = singleton::mpi().getRank() + 1;
  return false;
}


template <typename T, template <typename U> class DESCRIPTOR>
BlockLatticeCuboid2D<T,DESCRIPTOR>::BlockLatticeCuboid2D
(BlockLatticeStructure2D<T,DESCRIPTOR>& blockLattice, const int iC)
  : BlockLatticeF2D<T,DESCRIPTOR>(blockLattice,1), _iC(iC)
{
  this->getName() = "cuboid";
}

template <typename T, template <typename U> class DESCRIPTOR>
bool BlockLatticeCuboid2D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  output[0] = _iC + 1;
  return false;
}


template <typename T, template <typename U> class DESCRIPTOR>
BlockLatticePhysPressure2D<T,DESCRIPTOR>::BlockLatticePhysPressure2D
(BlockLatticeStructure2D<T,DESCRIPTOR>& blockLattice, const UnitConverter<T,DESCRIPTOR>& converter)
  : BlockLatticePhysF2D<T,DESCRIPTOR>(blockLattice,converter,1)
{
  this->getName() = "physPressure";
}


template <typename T, template <typename U> class DESCRIPTOR>
bool BlockLatticePhysPressure2D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  // lattice pressure = c_s^2 ( rho -1 )
  T latticePressure = ( this->_blockLattice.get( input[0] , input[1] ).computeRho() - 1.0 ) / DESCRIPTOR<T>::invCs2;
  output[0] = this->_converter.getPhysPressure(latticePressure);

  return true;
}


template <typename T, template <typename U> class DESCRIPTOR>
BlockLatticePhysVelocity2D<T,DESCRIPTOR>::BlockLatticePhysVelocity2D
(BlockLatticeStructure2D<T,DESCRIPTOR>& blockLattice, const UnitConverter<T,DESCRIPTOR>& converter)
  : BlockLatticePhysF2D<T,DESCRIPTOR>(blockLattice,converter,2)
{
  this->getName() = "physVelocity";
}

template <typename T, template <typename U> class DESCRIPTOR>
bool BlockLatticePhysVelocity2D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  T rho;
  this->_blockLattice.get( input[0], input[1] ).computeRhoU(rho,output);
  output[0]=this->_converter.getPhysVelocity(output[0]);
  output[1]=this->_converter.getPhysVelocity(output[1]);
  return true;
}

template <typename T, template <typename U> class DESCRIPTOR>
BlockLatticePhysExternalPorosity2D<T,DESCRIPTOR>::BlockLatticePhysExternalPorosity2D
(BlockLatticeStructure2D<T,DESCRIPTOR>& blockLattice, const UnitConverter<T,DESCRIPTOR>& converter)
  : BlockLatticePhysF2D<T,DESCRIPTOR>(blockLattice,converter,2)
{
  this->getName() = "ExtPorosityField";
}

template <typename T, template <typename U> class DESCRIPTOR>
bool BlockLatticePhysExternalPorosity2D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  this->_blockLattice.get(input[0],input[1]).computeExternalField(DESCRIPTOR<T>::ExternalField::porosityIsAt, 1, output);
  return true;
}

template <typename T, template <typename U> class DESCRIPTOR>
BlockLatticePhysExternalVelocity2D<T,DESCRIPTOR>::BlockLatticePhysExternalVelocity2D
(BlockLatticeStructure2D<T,DESCRIPTOR>& blockLattice, const UnitConverter<T,DESCRIPTOR>& converter)
  : BlockLatticePhysF2D<T,DESCRIPTOR>(blockLattice,converter,2)
{
  this->getName() = "ExtVelocityField";
}

template <typename T, template <typename U> class DESCRIPTOR>
bool BlockLatticePhysExternalVelocity2D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  this->_blockLattice.get(input[0],input[1]).computeExternalField(1, 2, output);
  //  this->_blockLattice.get(input[0],input[1]).computeExternalField(DESCRIPTOR<T>::ExternalField::localDragBeginsAt,
  //                                                                  DESCRIPTOR<T>::ExternalField::sizeOfLocalDrag, output);
  output[0]=this->_converter.getPhysVelocity(output[0]);
  output[1]=this->_converter.getPhysVelocity(output[1]);
  return true;
}

template <typename T, template <typename U> class DESCRIPTOR>
BlockLatticePhysExternalParticleVelocity2D<T,DESCRIPTOR>::BlockLatticePhysExternalParticleVelocity2D
(BlockLatticeStructure2D<T,DESCRIPTOR>& blockLattice, const UnitConverter<T,DESCRIPTOR>& converter)
  : BlockLatticePhysF2D<T,DESCRIPTOR>(blockLattice,converter,2)
{
  this->getName() = "ExtParticleVelocityField";
}

template <typename T, template <typename U> class DESCRIPTOR>
bool BlockLatticePhysExternalParticleVelocity2D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  T* foo = this->_blockLattice.get(input[0],input[1]).getExternal(0);

  if (foo[3] > std::numeric_limits<T>::epsilon()) {
    output[0]=this->_converter.getPhysVelocity(foo[1]/foo[3]);
    output[1]=this->_converter.getPhysVelocity(foo[2]/foo[3]);
    return true;
  }
  output[0]=this->_converter.getPhysVelocity(foo[1]);
  output[1]=this->_converter.getPhysVelocity(foo[2]);
  return true;
}

template <typename T, template <typename U> class DESCRIPTOR>
BlockLatticePhysWallShearStress2D<T,DESCRIPTOR>::BlockLatticePhysWallShearStress2D
(BlockLatticeStructure2D<T,DESCRIPTOR>& blockLattice, BlockGeometryStructure2D<T>& blockGeometry,
 int material, const UnitConverter<T,DESCRIPTOR>& converter)
  : BlockLatticePhysF2D<T,DESCRIPTOR>(blockLattice,converter,1),
    _blockGeometry(blockGeometry),_blockLattice(blockLattice), _material(material), _converter(converter)
{
  this->getName() = "physWallShearStress";
  T omega = 1. / _converter.getLatticeRelaxationTime();
  T dt = _converter.getConversionFactorTime();
  _physFactor = -omega * DESCRIPTOR<T>::invCs2 / dt * _converter.getPhysDensity() * _converter.getPhysViscosity();
  std::vector<int> discreteNormalOutwards(3, 0);
  for (int iX = 0; iX < _blockGeometry.getNx(); iX++) {
    _discreteNormal.resize(_blockGeometry.getNx());
    for (int iY = 0; iY < _blockGeometry.getNy(); iY++) {
      _discreteNormal[iX].resize(_blockGeometry.getNy());
      if (_blockGeometry.get(iX, iY) == _material) {
        discreteNormalOutwards = _blockGeometry.getStatistics().getType(iX, iY);
        _discreteNormal[iX][iY].resize(2);
        _discreteNormal[iX][iY][0] = -discreteNormalOutwards[1];
        _discreteNormal[iX][iY][1] = -discreteNormalOutwards[2];
      }
    }
  }
}

template <typename T, template <typename U> class DESCRIPTOR>
bool BlockLatticePhysWallShearStress2D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  output[0] = T();

  if (this->_blockGeometry.get(input[0],input[1])==_material) {
    T traction[2];
    T pi[3];
    T rho = _blockLattice.get(input[0] + _discreteNormal[input[0]][input[1]][0],
                              input[1] + _discreteNormal[input[0]][input[1]][1]).computeRho();
    _blockLattice.get(input[0] +   _discreteNormal[input[0]][input[1]][0],
                      input[1] +   _discreteNormal[input[0]][input[1]][1]).computeStress(pi);

    traction[0] = pi[0]*_physFactor/rho*_discreteNormal[input[0]][input[1]][0] +
                  pi[1]*_physFactor/rho*_discreteNormal[input[0]][input[1]][1];
    traction[1] = pi[1]*_physFactor/rho*_discreteNormal[input[0]][input[1]][0] +
                  pi[2]*_physFactor/rho*_discreteNormal[input[0]][input[1]][1];

    T tractionNormalComponent[3];
    // scalar product of traction and normal vector
    T traction_normal_SP = traction[0] * _discreteNormal[input[0]][input[1]][0] +
                           traction[1] * _discreteNormal[input[0]][input[1]][1];
    tractionNormalComponent[0] = traction_normal_SP * _discreteNormal[input[0]][input[1]][0];
    tractionNormalComponent[1] = traction_normal_SP * _discreteNormal[input[0]][input[1]][1];

    T WSS[3];
    WSS[0] = traction[0] - tractionNormalComponent[0];
    WSS[1] = traction[1] - tractionNormalComponent[1];
    // magnitude of the wall shear stress vector
    output[0] = sqrt(WSS[0]*WSS[0] + WSS[1]*WSS[1]);
    return true;
  } else {
    return true;
  }
}

template <typename T, template <typename U> class DESCRIPTOR>
BlockLatticePhysBoundaryForce2D<T,DESCRIPTOR>::BlockLatticePhysBoundaryForce2D
(BlockLatticeStructure2D<T,DESCRIPTOR>& blockLattice, BlockGeometryStructure2D<T>& blockGeometry,
 int material, const UnitConverter<T,DESCRIPTOR>& converter)
  : BlockLatticePhysF2D<T,DESCRIPTOR>(blockLattice,converter,2),
    _blockGeometry(blockGeometry), _material(material)
{
  this->getName() = "physBoundaryForce";
}

template <typename T, template <typename U> class DESCRIPTOR>
bool BlockLatticePhysBoundaryForce2D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  output[0] = T();
  output[1] = T();

  if (this->_blockGeometry.get(input[0],input[1])==_material) {
    for (int iPop = 1; iPop < DESCRIPTOR<T>::q; ++iPop) {
      // Get direction
      const int* c = DESCRIPTOR<T>::c[iPop];
      // Get next cell located in the current direction
      // Check if the next cell is a fluid node
      if (_blockGeometry.get(input[0] + c[0], input[1] + c[1]) == 1) {
        // Get f_q of next fluid cell where l = opposite(q)
        T f = this->_blockLattice.get(input[0] + c[0], input[1] + c[1])[iPop];
        // Get f_l of the boundary cell
        // Add f_q and f_opp
        f += this->_blockLattice.get(input[0], input[1])[util::opposite<DESCRIPTOR<T> >(iPop)];
        // Update force
        output[0] -= c[0] * f;
        output[1] -= c[1] * f;
      }
    }
    output[0] = this->_converter.getPhysForce(output[0]);
    output[1] = this->_converter.getPhysForce(output[1]);
    return true;
  } else {
    return true;
  }
}

template <typename T, template <typename U> class DESCRIPTOR>
BlockLatticePhysCorrBoundaryForce2D<T,DESCRIPTOR>::BlockLatticePhysCorrBoundaryForce2D
(BlockLatticeStructure2D<T,DESCRIPTOR>& blockLattice,BlockGeometry2D<T>& blockGeometry,
 int material, const UnitConverter<T,DESCRIPTOR>& converter)
  : BlockLatticePhysF2D<T,DESCRIPTOR>(blockLattice,converter,2),
    _blockGeometry(blockGeometry), _material(material)
{
  this->getName() = "physCorrBoundaryForce";
}

template <typename T, template <typename U> class DESCRIPTOR>
bool BlockLatticePhysCorrBoundaryForce2D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  //  int globIC = input[0];
  //  int locix= input[1];
  //  int lociy= input[2];
  //  int lociz= input[3];

  //  std::vector<T> force(3, T());
  //  if ( this->_blockLattice.get_load().rank(globIC) == singleton::mpi().getRank() ) {
  //    int globX = (int)this->_blockLattice.get_cGeometry().get(globIC).get_globPosX() + locix;
  //    int globY = (int)this->_blockLattice.get_cGeometry().get(globIC).get_globPosY() + lociy;
  //    int globZ = (int)this->_blockLattice.get_cGeometry().get(globIC).get_globPosZ() + lociz;
  //    if(BlockGeometry.getMaterial(globX, globY, globZ) == material) {
  //      for (int iPop = 1; iPop < DESCRIPTOR<T>::q ; ++iPop) {
  //        // Get direction
  //        const int* c = DESCRIPTOR<T>::c[iPop];
  //        // Get next cell located in the current direction
  //        // Check if the next cell is a fluid node
  //        if (BlockGeometry.getMaterial(globX + c[0], globY + c[1], globZ + c[2]) == 1) {
  //          int overlap = this->_blockLattice.getOverlap();
  //          // Get f_q of next fluid cell where l = opposite(q)
  //          T f = this->_blockLattice.getExtendedBlockLattice(this->_blockLattice.get_load().loc(globIC)).get(locix+overlap + c[0], lociy+overlap + c[1], lociz+overlap + c[2])[iPop];
  //          // Get f_l of the boundary cell
  //          // Add f_q and f_opp
  //          f += this->_blockLattice.getExtendedBlockLattice(this->_blockLattice.get_load().loc(globIC)).get(locix+overlap, lociy+overlap, lociz+overlap)[util::opposite<DESCRIPTOR<T> >(iPop)];
  //          // Update force
  //          force[0] -= c[0]*(f-2.*DESCRIPTOR<T>::t[iPop]);
  //          force[1] -= c[1]*(f-2.*DESCRIPTOR<T>::t[iPop]);
  //          force[2] -= c[2]*(f-2.*DESCRIPTOR<T>::t[iPop]);
  //        }
  //      }
  //      force[0]=this->_converter.physForce(force[0]);
  //      force[1]=this->_converter.physForce(force[1]);
  //      force[2]=this->_converter.physForce(force[2]);
  //      return force;
  //    }
  //    else {
  //      return force;
  //    }
  //  }
  //  else {
  //    return force;
  //  }
  return false;
}


template <typename T, template <typename U> class DESCRIPTOR>
BlockLatticePorosity2D<T,DESCRIPTOR>::BlockLatticePorosity2D
(BlockLatticeStructure2D<T,DESCRIPTOR>& blockLattice, BlockGeometryStructure2D<T>& blockGeometry,
 int material, const UnitConverter<T,DESCRIPTOR>& converter)
  : BlockLatticePhysF2D<T,DESCRIPTOR>(blockLattice,converter,1), _blockGeometry(blockGeometry), _material(material)
{
  this->getName() = "porosity";
}

// under construction
template <typename T, template <typename U> class DESCRIPTOR>
bool BlockLatticePorosity2D<T,DESCRIPTOR>::operator()(T output[], const int input[])
{
  //  int globIC = input[0];
  //  int locix= input[1];
  //  int lociy= input[2];
  //  int lociz= input[3];

  //  T* value = new T[1];
  //  int overlap = this->_blockLattice.getOverlap();
  //  this->_blockLattice.getExtendedBlockLattice(this->_blockLattice.get_load().loc(globIC)).get(locix+overlap, lociy+overlap, lociz+overlap).computeExternalField(0,1,value);
  //  std::vector<T> result(1,value[0]);
  //  delete value;
  //  return result;
  return false;
}


template <typename T, template <typename U> class DESCRIPTOR>
BlockLatticePhysPermeability2D<T,DESCRIPTOR>::BlockLatticePhysPermeability2D
(BlockLatticeStructure2D<T,DESCRIPTOR>& blockLattice, BlockGeometryStructure2D<T>& blockGeometry,
 int material, const UnitConverter<T,DESCRIPTOR>& converter)
  : BlockLatticePhysF2D<T,DESCRIPTOR>(blockLattice,converter,1), _blockGeometry(blockGeometry), _material(material)
{
  this->getName() = "permeability";
}

template <typename T, template <typename U> class DESCRIPTOR>
bool BlockLatticePhysPermeability2D<T,DESCRIPTOR>::operator()(T output[], const int input[])
{
  //  int globIC = input[0];
  //  int locix= input[1];
  //  int lociy= input[2];
  //  int lociz= input[3];

  //  T* value = new T[1];
  //  int overlap = this->_blockLattice.getOverlap();
  //  this->_blockLattice.getExtendedBlockLattice(this->_blockLattice.get_load().loc(globIC)).get(locix+overlap, lociy+overlap, lociz+overlap).computeExternalField(0,1,value);
  //  std::vector<T> result(1,this->_converter.physPermeability(value[0]));
  //  delete value;
  //  if (!(result[0]<42)&&!(result[0]>42)&&!(result[0]==42)) result[0] = 999999;
  //  if (isinf(result[0])) result[0] = 1e6;
  //  return result;
  return false;
}


template <typename T, template <typename U> class DESCRIPTOR>
BlockLatticePhysDarcyForce2D<T,DESCRIPTOR>::BlockLatticePhysDarcyForce2D
(BlockLatticeStructure2D<T,DESCRIPTOR>& blockLattice, BlockGeometry2D<T>& blockGeometry,
 int material, const UnitConverter<T,DESCRIPTOR>& converter)
  : BlockLatticePhysF2D<T,DESCRIPTOR>(blockLattice,converter,2),
    _blockGeometry(blockGeometry), _material(material)
{
  this->getName() = "alphaU";
}

template <typename T, template <typename U> class DESCRIPTOR>
bool BlockLatticePhysDarcyForce2D<T,DESCRIPTOR>::operator()(T output[], const int input[])
{
  BlockLatticePhysPermeability2D<T,DESCRIPTOR> permeability(this->_blockLattice,this->_blockGeometry,this->_material,this->_converter);
  BlockLatticeVelocity2D<T,DESCRIPTOR> velocity(this->_blockLattice);

  T nu = this->_converter.getPhysViscosity();
  permeability(output,input);
  T K = output[0];
  velocity(output,input);

  output[0] *= -nu/K;
  output[1] *= -nu/K;

  return true;
}


template <typename T, template <typename U> class DESCRIPTOR>
BlockLatticeAverage2D<T,DESCRIPTOR>::BlockLatticeAverage2D
(BlockLatticeF2D<T,DESCRIPTOR>& f, BlockGeometry2D<T>& blockGeometry,
 int material, T radius)
  : BlockLatticeF2D<T,DESCRIPTOR>(f.getBlockLattice(), f.getTargetDim()),
    _f(f), _blockGeometry(blockGeometry), _material(material), _radius(radius)
{
  this->getName() = "Average("+f.getName()+")";
}


template <typename T, template <typename U> class DESCRIPTOR>
bool BlockLatticeAverage2D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  //  CuboidGeometry2D<T>& cGeometry = f.getBlockLattice2D().get_cGeometry();
  //  loadBalancer& load = f.getBlockLattice2D().get_load();

  //  //create boolean indicator functor isInSphere
  //  std::vector<T> center(3,0);
  //  center[0] = (int)cGeometry.get(load.glob(input[0])).get_globPosX() + input[1];
  //  center[1] = (int)cGeometry.get(load.glob(input[0])).get_globPosY() + input[2];
  //  center[2] = (int)cGeometry.get(load.glob(input[0])).get_globPosZ() + input[3];
  //  SphereAnalyticalF2D<bool,T> isInSphere(center,radius);

  // iterate over all cuboids & points and test for material && isInSphere
  //  std::vector<T> tmp( this->_n, T() );
  //  int numVoxels(0);
  //  if (this->_blockGeometry.getMaterial(center[0],center[1],center[2]) == material) {
  //    for (int iC=0; iC<load.size(); iC++) {
  //      int nX = cGeometry.get(load.glob(iC)).getNx();
  //      int nY = cGeometry.get(load.glob(iC)).getNy();
  //      int nZ = cGeometry.get(load.glob(iC)).getNz();
  //      for (int iX=0; iX<nX; ++iX) {
  //        for (int iY=0; iY<nY; ++iY) {
  //          for (int iZ=0; iZ<nZ; iZ++) {
  //            std::vector<T> glob(3,0);
  //            glob[0] = (int)cGeometry.get(load.glob(iC)).get_globPosX() + iX;
  //            glob[1] = (int)cGeometry.get(load.glob(iC)).get_globPosY() + iY;
  //            glob[2] = (int)cGeometry.get(load.glob(iC)).get_globPosZ() + iZ;
  //            if (this->_blockGeometry.getMaterial(glob[0],glob[1],glob[2]) == material
  //                && isInSphere(glob)[0]==true) {
  //              for (unsigned iD=0; iD<f(load.glob(0),0,0,0).size(); iD++) {
  //                tmp[iD]+=f(load.glob(iC),iX,iY,iZ)[iD];
  //              }
  //              numVoxels++;
  //            }
  //          }
  //        }
  //      }
  //    }

  //#ifdef PARALLEL_MODE_MPI
  //    singleton::mpi().reduceAndBcast(numVoxels, MPI_SUM);
  //#endif
  //    for (int iD=0; iD<f.getTargetDim(); iD++) {
  //#ifdef PARALLEL_MODE_MPI
  //      singleton::mpi().reduceAndBcast(tmp[iD], MPI_SUM);
  //#endif
  //      if (numVoxels>0) {
  //        tmp[iD] /= numVoxels;
  //      }
  //    }
  //  }
  //  return tmp;

  return false;
}


template <typename T, template <typename U> class DESCRIPTOR>
BlockEuklidNorm2D<T,DESCRIPTOR>::BlockEuklidNorm2D(BlockF2D<T>& f)
  : BlockF2D<T>(f.getBlockStructure(),1), _f(f)
{
  this->getName() = "l2("+f.getName()+")";
}

template <typename T, template <typename U> class DESCRIPTOR>
bool BlockEuklidNorm2D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  output[0] = T();  // flash output, since this methods adds values.
  T data[_f.getTargetDim()];
  _f(data,input);
  for ( int i = 0; i < _f.getTargetDim(); ++i) {
    output[0] += data[i]*data[i];
  }
  output[0] = sqrt(output[0]);
  return true;
}

template<typename T, template<typename U> class DESCRIPTOR>
BlockLatticeMomentumExchangeForce2D<T, DESCRIPTOR>::BlockLatticeMomentumExchangeForce2D(
  BlockLatticeStructure2D<T, DESCRIPTOR>& blockLattice, BlockGeometryStructure2D<T>& blockGeometry,
  std::vector<ParticleIndicatorF2D<T,T>* >& indicator, const UnitConverter<T,DESCRIPTOR>& converter, int overlap)
  : BlockLatticePhysF2D<T, DESCRIPTOR>(blockLattice, converter, 4*indicator.size()), _blockGeometry(blockGeometry), _vectorOfIndicator(indicator), _overlap(overlap)
{
  this->getName() = "physMomentumExchangeForce";
}

template<typename T, template<typename U> class DESCRIPTOR>
bool BlockLatticeMomentumExchangeForce2D<T, DESCRIPTOR>::operator()(T output[], const int input[])
{
  // iterate over all particles in _indicator
  for (size_t iInd=0; iInd!=_vectorOfIndicator.size(); iInd++) {

    int numVoxels = 0;
    // check for intersection of cuboid and indicator
    if (_blockGeometry.getOrigin()[0] <= _vectorOfIndicator[iInd]->getCircumRadius()+_vectorOfIndicator[iInd]->getPos()[0]
        && _blockGeometry.getOrigin()[1] <= _vectorOfIndicator[iInd]->getCircumRadius()+_vectorOfIndicator[iInd]->getPos()[1]
        && _vectorOfIndicator[iInd]->getPos()[0]-_vectorOfIndicator[iInd]->getCircumRadius() <= _blockGeometry.getOrigin()[0] + _blockGeometry.getExtend()[0] * _blockGeometry.getDeltaR()
        && _vectorOfIndicator[iInd]->getPos()[1]-_vectorOfIndicator[iInd]->getCircumRadius() <= _blockGeometry.getOrigin()[1] + _blockGeometry.getExtend()[1] * _blockGeometry.getDeltaR()) {

      // compute size of intersection for iteration
      int start[2] = {0}, span[2] = {0};
      T invDeltaR = 1./_blockGeometry.getDeltaR();
      for (int k=0; k<2; k++) {
        start[k] = int( (_vectorOfIndicator[iInd]->getPos()[k]-_vectorOfIndicator[iInd]->getCircumRadius() - _blockGeometry.getOrigin()[k]) * invDeltaR );
        if (start[k] < 0) {
          start[k] = 0;
        }
        span[k] = int( (2.*_vectorOfIndicator[iInd]->getCircumRadius())*invDeltaR + 3 );
        if (span[k] + start[k] > _blockGeometry.getExtend()[k]) {
          span[k] = _blockGeometry.getExtend()[k] - start[k];
        }
      }

      // iterate over cells in the constructed intersection box
      for (int iX = start[0]; iX < start[0]+span[0]; iX++) {
        for (int iY = start[1]; iY < start[1]+span[1]; iY++) {
          // check if cell belongs to particle
          T inside[1] = {0.};
          T posIn[2] = {0.};
          _blockGeometry.getPhysR(posIn, iX, iY);
          (*(_vectorOfIndicator[iInd]))( inside, posIn);
          if ( !util::nearZero(inside[0]) && this->_blockGeometry.get(iX,iY)==1) {
            // compute momentum exchange force on particle
            T tmpForce[2] = {0.,0.};
            for (int iPop = 1; iPop < DESCRIPTOR<T>::q ; ++iPop) {
              // Get direction
              const int* c = DESCRIPTOR<T>::c[iPop];

              T posOut[2] = {0.};
              T inside2[1] = {0.};
              _blockGeometry.getPhysR(posOut, iX+c[0], iY+c[1]);
              (*(_vectorOfIndicator[iInd]))( inside2, posOut);
              if (util::nearZero(inside2[0])) {
                // Get f_q of next fluid cell where l = opposite(q)
                T f = this->_blockLattice.get(iX+_overlap + c[0], iY+_overlap + c[1])[iPop];
                // Get f_l of the boundary cell
                // Add f_q and f_opp
                f += this->_blockLattice.get(iX+_overlap, iY+_overlap)[util::opposite<DESCRIPTOR<T> >(iPop)];
                // Update force
                tmpForce[0] -= c[0]*f;
                tmpForce[1] -= c[1]*f;
              }
            }
            // convert force to SI units and compute torque
            numVoxels++;
            tmpForce[0] = this->_converter.getPhysForce(tmpForce[0]);
            tmpForce[1] = this->_converter.getPhysForce(tmpForce[1]);
            output[0+iInd*4] += tmpForce[0];
            output[1+iInd*4] += tmpForce[1];
            output[2+iInd*4] += (posIn[0]-_vectorOfIndicator[iInd]->getPos()[0])*tmpForce[1] - (posIn[1]-_vectorOfIndicator[iInd]->getPos()[1])*tmpForce[0];
          }

        }
      }
    }
    output[3+iInd*4] = numVoxels;
  }
  return true;
}

template <typename T, template <typename U> class DESCRIPTOR, template <typename V> class ThermalDESCRIPTOR>
BlockLatticePhysTemperature2D<T,DESCRIPTOR,ThermalDESCRIPTOR>::BlockLatticePhysTemperature2D
(BlockLatticeStructure2D<T,ThermalDESCRIPTOR>& blockLattice, ThermalUnitConverter<T,DESCRIPTOR,ThermalDESCRIPTOR> const& converter)
  : BlockLatticeThermalPhysF2D<T,DESCRIPTOR,ThermalDESCRIPTOR>(blockLattice,converter,1)
{
  this->getName() = "physTemperature";
}


template <typename T, template <typename U> class DESCRIPTOR, template <typename V> class ThermalDESCRIPTOR>
bool BlockLatticePhysTemperature2D<T,DESCRIPTOR,ThermalDESCRIPTOR>::operator() (T output[], const int input[])
{
  T latticeTemperature = this->_blockLattice.get( input[0] , input[1] ).computeRho();
  output[0] = this->_converter.getPhysTemperature(latticeTemperature);

  return true;
}
} // end namespace olb

#endif
