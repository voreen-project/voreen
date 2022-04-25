/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2013 Patrick Nathen, Mathias J. Krause
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

#ifndef TURBULENT_F_3D_HH
#define TURBULENT_F_3D_HH

#include<vector>
#include<cmath>

#include "turbulentF3D.h"
#include "dynamics/smagorinskyBGKdynamics.h"
#include "core/superLattice.h"
#include "core/finiteDifference.h"
#include "geometry/superGeometry.h"
#include "utilities/vectorHelpers.h"
#include "dynamics/lbm.h"  // for computation of lattice rho and velocity


namespace olb {


///////////////////////////// SuperLatticeYplus3D //////////////////////////////
template <typename T, typename DESCRIPTOR>
SuperLatticeYplus3D<T,DESCRIPTOR>::SuperLatticeYplus3D(SuperLattice<T,DESCRIPTOR>& sLattice,
    const UnitConverter<T,DESCRIPTOR>& converter, SuperGeometry<T,3>& superGeometry,
    IndicatorF3D<T>& indicator, const int material )
  : SuperLatticePhysF3D<T,DESCRIPTOR>(sLattice,converter,1),
    _superGeometry(superGeometry), _indicator(indicator), _material(material)
{
  this->getName() = "yPlus";
}

template <typename T, typename DESCRIPTOR>
bool SuperLatticeYplus3D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  int globIC = input[0];
  int locix = input[1];
  int lociy = input[2];
  int lociz = input[3];

  output[0]=T();
  if ( this->_sLattice.getLoadBalancer().rank(globIC) == singleton::mpi().getRank() ) {
    std::vector<T> normalTemp(3,T());
    std::vector<T> normal(3,T());       // normalized
    T counter = T();
    T distance = T();
    if (_superGeometry.get(input) == 1) {
      for (int iPop = 1; iPop < DESCRIPTOR::q; ++iPop) {
        if (_superGeometry.get(input[0],
                               input[1] + descriptors::c<DESCRIPTOR>(iPop,0),
                               input[2] + descriptors::c<DESCRIPTOR>(iPop,1),
                               input[3] + descriptors::c<DESCRIPTOR>(iPop,2)) == _material) {
          counter++;
          normalTemp[0] += descriptors::c<DESCRIPTOR>(iPop,0);
          normalTemp[1] += descriptors::c<DESCRIPTOR>(iPop,1);
          normalTemp[2] += descriptors::c<DESCRIPTOR>(iPop,2);
        }
      }
      if ( !util::nearZero(counter) ) {
        // get physical Coordinates at intersection

        std::vector<T> physR(3, T());
        _superGeometry.getCuboidGeometry().getPhysR(&(physR[0]), &(input[0]));

        T voxelSize = _superGeometry.getCuboidGeometry().get(globIC).getDeltaR();

        normal = util::normalize(normalTemp);

        std::vector<T> direction(normal);
        direction[0] = voxelSize*normal[0]*2.;
        direction[1] = voxelSize*normal[1]*2.;
        direction[2] = voxelSize*normal[2]*2.;

        // calculate distance to STL file
        if ( _indicator.distance(distance, physR, direction) ) {
          // call stress at this point
          T rho;
          T u[3];
          T pi[6];
          auto cell = this->_sLattice.getBlock(this->_sLattice.getLoadBalancer().loc(globIC)).get(locix, lociy, lociz);
          cell.computeAllMomenta(rho, u, pi);

          // Compute phys stress tau = mu*du/dx
          T omega = 1. / this->_converter.getLatticeRelaxationTime();
          T dt = this->_converter.getConversionFactorTime();
          T physFactor = -omega*descriptors::invCs2<T,DESCRIPTOR>()/rho/2./dt*this->_converter.getPhysDensity(rho)*this->_converter.getPhysViscosity();

          //  Totel Stress projected from cell in normal direction on obstacle
          T Rx = pi[0]*physFactor*normal[0] + pi[1]*physFactor*normal[1] + pi[2]*physFactor*normal[2];
          T Ry = pi[1]*physFactor*normal[0] + pi[3]*physFactor*normal[1] + pi[4]*physFactor*normal[2];
          T Rz = pi[2]*physFactor*normal[0] + pi[4]*physFactor*normal[1] + pi[5]*physFactor*normal[2];

          // Stress appearing as pressure in corresponding direction is calculated and substracted
          T R_res_pressure = normal[0]*pi[0]*physFactor*normal[0] + normal[0]*pi[1]*physFactor*normal[1] + normal[0]*pi[2]*physFactor*normal[2]
                             +normal[1]*pi[1]*physFactor*normal[0] + normal[1]*pi[3]*physFactor*normal[1] + normal[1]*pi[4]*physFactor*normal[2]
                             +normal[2]*pi[2]*physFactor*normal[0] + normal[2]*pi[4]*physFactor*normal[1] + normal[2]*pi[5]*physFactor*normal[2];

          Rx -= R_res_pressure *normal[0];
          Ry -= R_res_pressure *normal[1];
          Rz -= R_res_pressure *normal[2];

          T tau_wall = util::sqrt(Rx*Rx+Ry*Ry+Rz*Rz);
          T u_tau = util::sqrt(tau_wall/this->_converter.getPhysDensity(rho));
          //y_plus
          output[0] = distance*u_tau / this->_converter.getPhysViscosity();
        } // if 4
      }
    }
  }
  return true;
}

////////////////////////BlockFiniteDifference3D//////////////////////////////////
template <typename T>
BlockFiniteDifference3D<T>::BlockFiniteDifference3D
(BlockGeometry<T,3>& blockGeometry, BlockF3D<T>& blockFunctor, std::list<int>& matNumber)
  : BlockF3D<T>(blockFunctor.getBlockStructure(), 3*blockFunctor.getTargetDim()), _blockGeometry(blockGeometry), _blockFunctor(blockFunctor), _matNumber(matNumber)
{
  this->getName() = "FiniteDifference";
  _targetDim = _blockFunctor.getTargetDim();
  _n[0] = this-> _blockGeometry.getNx()-1;
  _n[1] = this-> _blockGeometry.getNy()-1;
  _n[2] = this-> _blockGeometry.getNz()-1;

}

template <typename T>
bool BlockFiniteDifference3D<T>::operator() (T output[], const int input[])
{
//  // derivation tensor
  std::vector<std::vector<T>> fdGrad;

  fdGrad.resize(_targetDim);
  for (int i = 0; i < _targetDim; i++) {
    fdGrad[i].resize(3);
  }

  for (int i = 0; i < 3; i++) {
    int fInput_p[3];
    fInput_p[0] = input[0];
    fInput_p[1] = input[1];
    fInput_p[2] = input[2];
    fInput_p[i]+=1;

    int fInput_2p[3];
    fInput_2p[0] = input[0];
    fInput_2p[1] = input[1];
    fInput_2p[2] = input[2];
    fInput_2p[i]+=2;

    int fInput_3p[3];
    fInput_3p[0] = input[0];
    fInput_3p[1] = input[1];
    fInput_3p[2] = input[2];
    fInput_3p[i]+=3;

    int fInput_4p[3];
    fInput_4p[0] = input[0];
    fInput_4p[1] = input[1];
    fInput_4p[2] = input[2];
    fInput_4p[i]+=4;

    int fInput_n[3];
    fInput_n[0] = input[0];
    fInput_n[1] = input[1];
    fInput_n[2] = input[2];
    fInput_n[i]-=1;

    int fInput_2n[3];
    fInput_2n[0] = input[0];
    fInput_2n[1] = input[1];
    fInput_2n[2] = input[2];
    fInput_2n[i]-=2;

    int fInput_3n[3];
    fInput_3n[0] = input[0];
    fInput_3n[1] = input[1];
    fInput_3n[2] = input[2];
    fInput_3n[i]-=3;

    int fInput_4n[3];
    fInput_4n[0] = input[0];
    fInput_4n[1] = input[1];
    fInput_4n[2] = input[2];
    fInput_4n[i]-=4;

    T fOutput[_targetDim];
    _blockFunctor(fOutput,input);

    if (input[i] < 3) {
      if (std::find(_matNumber.begin(), _matNumber.end(), _blockGeometry.get({fInput_2p[0], fInput_2p[1], fInput_2p[2]})) == _matNumber.end()) {
        T fOutput_p[_targetDim];
        _blockFunctor(fOutput_p,fInput_p);
        for (int j=0; j < _targetDim; j++) {
          fdGrad[j][i]= -fOutput[j] + fOutput_p[j];
        }
      }
      else {
        T fOutput_p[_targetDim];
        _blockFunctor(fOutput_p,fInput_p);
        T fOutput_2p[_targetDim];
        _blockFunctor(fOutput_2p,fInput_2p);
        for (int j=0; j < _targetDim; j++) {
          fdGrad[j][i]=fd::boundaryGradient(fOutput[j], fOutput_p[j], fOutput_2p[j]);
        }
      }
    }
    else if (input[i] > _n[i]-3) {
      if (std::find(_matNumber.begin(), _matNumber.end(), _blockGeometry.get({fInput_2n[0], fInput_2n[1], fInput_2n[2]})) == _matNumber.end()) {
        T fOutput_n[_targetDim];
        _blockFunctor(fOutput_n,fInput_n);
        for (int j=0; j < _targetDim; j++) {
          fdGrad[j][i]= -fOutput_n[j] + fOutput[j];
        }
      }
      else {
        T fOutput_n[_targetDim];
        _blockFunctor(fOutput_n,fInput_n);
        T fOutput_2n[_targetDim];
        _blockFunctor(fOutput_2n,fInput_2n);
        for (int j=0; j < _targetDim; j++) {
          fdGrad[j][i]=fd::boundaryGradient(-fOutput[j], -fOutput_n[j], -fOutput_2n[j]);
        }
      }
    }
    else {
      if ( std::find(_matNumber.begin(), _matNumber.end(), _blockGeometry.get({fInput_n[0], fInput_n[1], fInput_n[2]})) == _matNumber.end()  &&
          std::find(_matNumber.begin(), _matNumber.end(), _blockGeometry.get({fInput_p[0], fInput_p[1], fInput_p[2]})) == _matNumber.end() ) {
        for (int j=0; j < _targetDim; j++) {
          fdGrad[j][i]=0.;
        }
        // boundary treatment with Second-order asymmetric gradient
      }
      else if (std::find(_matNumber.begin(), _matNumber.end(), _blockGeometry.get({fInput_n[0], fInput_n[1], fInput_n[2]})) == _matNumber.end()) {
        if (std::find(_matNumber.begin(), _matNumber.end(), _blockGeometry.get({fInput_2p[0], fInput_2p[1], fInput_2p[2]})) == _matNumber.end()) {
          T fOutput_p[_targetDim];
          _blockFunctor(fOutput_p,fInput_p);
          for (int j=0; j < _targetDim; j++) {
            fdGrad[j][i]= -fOutput[j] + fOutput_p[j];
          }
        }
        else {
          T fOutput_p[_targetDim];
          _blockFunctor(fOutput_p,fInput_p);
          T fOutput_2p[_targetDim];
          _blockFunctor(fOutput_2p,fInput_2p);
          for (int j=0; j < _targetDim; j++) {
            fdGrad[j][i]=fd::boundaryGradient(fOutput[j], fOutput_p[j], fOutput_2p[j]);
          }
        }
      }
      else if (std::find(_matNumber.begin(), _matNumber.end(), _blockGeometry.get({fInput_p[0], fInput_p[1], fInput_p[2]})) == _matNumber.end() ) {
        if (std::find(_matNumber.begin(), _matNumber.end(), _blockGeometry.get({fInput_2n[0], fInput_2n[1], fInput_2n[2]})) == _matNumber.end()) {
          T fOutput_n[_targetDim];
          _blockFunctor(fOutput_n,fInput_n);
          for (int j=0; j < _targetDim; j++) {
            fdGrad[j][i]= -fOutput_n[j] + fOutput[j];
          }
        }
        else {
          T fOutput_n[_targetDim];
          _blockFunctor(fOutput_n,fInput_n);
          T fOutput_2n[_targetDim];
          _blockFunctor(fOutput_2n,fInput_2n);
          for (int j=0; j < _targetDim; j++) {
            fdGrad[j][i]=fd::boundaryGradient(-fOutput[j], -fOutput_n[j], -fOutput_2n[j]);
          }
        }
      }
      else {
        //inner domain 8th order central difference
        T fOutput_n[_targetDim];
        _blockFunctor(fOutput_n,fInput_n);

        T fOutput_2n[_targetDim];
        _blockFunctor(fOutput_2n,fInput_2n);

        T fOutput_3n[_targetDim];
        _blockFunctor(fOutput_3n,fInput_3n);

        T fOutput_4n[_targetDim];
        _blockFunctor(fOutput_4n,fInput_4n);

        T fOutput_p[_targetDim];
        _blockFunctor(fOutput_p,fInput_p);

        T fOutput_2p[_targetDim];
        _blockFunctor(fOutput_2p,fInput_2p);

        T fOutput_3p[_targetDim];
        _blockFunctor(fOutput_3p,fInput_3p);

        T fOutput_4p[_targetDim];
        _blockFunctor(fOutput_4p,fInput_4p);
        for (int j=0; j < _targetDim; j++) {
          //fdGrad[j][i]=fd::centralGradient(fOutput_p[j], fOutput_n[j]);
          fdGrad[j][i]=((T)672*(fOutput_p[j]-fOutput_n[j])+(T)168*(fOutput_2n[j]-fOutput_2p[j])
                        +(T)32*(fOutput_3p[j]-fOutput_3n[j])+(T)3*(fOutput_4n[j]-fOutput_4p[j])) / 840.;
        }
      }
    }
    for (int i=0; i < 3; i++) {
      for (int j=0; j < _targetDim; j++) {
        output[i*3+j] = fdGrad[i][j];
      }
    }
  }
  return true;
}

////////////////////////SuperFiniteDifference3D//////////////////////////////////
template <typename T>
SuperFiniteDifference3D<T>::SuperFiniteDifference3D
(SuperGeometry<T,3>& sGeometry, SuperF3D<T>& sFunctor, std::list<int>& matNumber) : SuperF3D<T>(sFunctor.getSuperStructure(),3*sFunctor.getTargetDim()),
  _sGeometry(sGeometry),_sFunctor(sFunctor), _matNumber(matNumber)
{
  this->getName() = "FiniteDifference";
  int maxC = this->_superStructure.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back(new BlockFiniteDifference3D<T> ( _sGeometry.getBlockGeometry(iC), _sFunctor.getBlockF(iC), _matNumber ));
  }
}

////////////////////////BlockPhysFiniteDifference3D//////////////////////////////////
template <typename T, typename DESCRIPTOR>
BlockPhysFiniteDifference3D<T,DESCRIPTOR>::BlockPhysFiniteDifference3D
(BlockF3D<T>& blockFinDiff, const UnitConverter<T,DESCRIPTOR>& converter)
  : BlockF3D<T>(blockFinDiff.getBlockStructure(), 3*blockFinDiff.getTargetDim()), _blockFinDiff(blockFinDiff), _converter(converter)
{
  this->getName() = "PhysFiniteDifference";
  _targetDim = _blockFinDiff.getTargetDim();

}

template <typename T, typename DESCRIPTOR>
bool BlockPhysFiniteDifference3D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  _blockFinDiff(output,input);
  for (int i = 0; i < _targetDim; i++) {
    output[i] /= _converter.getConversionFactorLength();
  }
  return true;
}

////////////////////////SuperPhysFiniteDifference3D//////////////////////////////////
template <typename T, typename DESCRIPTOR>
SuperPhysFiniteDifference3D<T,DESCRIPTOR>::SuperPhysFiniteDifference3D
(SuperGeometry<T,3>& sGeometry, SuperF3D<T>& sFunctor, std::list<int>& matNumber, const UnitConverter<T,DESCRIPTOR>& converter) : SuperF3D<T>(sFunctor.getSuperStructure(),3*sFunctor.getTargetDim()),
  _sFinDiff(sGeometry,sFunctor,matNumber),_converter(converter)
{
  this->getName() = "PhysFiniteDifference";
  int maxC = this->_superStructure.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back(new BlockPhysFiniteDifference3D<T,DESCRIPTOR> (_sFinDiff.getBlockF(iC), _converter ));
  }
}

////////////////////////BlockLatticeVelocityGradientFD3D//////////////////////////////////
template <typename T, typename DESCRIPTOR>
BlockLatticeVelocityGradientFD3D<T,DESCRIPTOR>::BlockLatticeVelocityGradientFD3D
(BlockLattice<T,DESCRIPTOR>& blockLattice, BlockF3D<T>& blockFinDiff)
  : BlockLatticeF3D<T,DESCRIPTOR>(blockLattice, 9), _blockFinDiff(blockFinDiff)
{
  this->getName() = "VelocityGradientFD";
}

template <typename T, typename DESCRIPTOR>
bool BlockLatticeVelocityGradientFD3D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  //1 dudx 2 dudy 3 dudz
  //4 dydx 5 dydy 6 dydz
  //7 dwdx 8 dwdy 9 dwdz
  _blockFinDiff(output,input);
  return true;
}

////////////////////////BlockLatticeExternalVelocityGradientFD3D//////////////////////////////////
template <typename T, typename DESCRIPTOR>
BlockLatticeExternalVelocityGradientFD3D<T,DESCRIPTOR>::BlockLatticeExternalVelocityGradientFD3D
(BlockLattice<T,DESCRIPTOR>& blockLattice, BlockF3D<T>& blockFinDiff)
  : BlockLatticeF3D<T,DESCRIPTOR>(blockLattice, 9), _blockFinDiff(blockFinDiff)
{
  this->getName() = "externalVelocityGradientFD";
}

template <typename T, typename DESCRIPTOR>
bool BlockLatticeExternalVelocityGradientFD3D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  //1 dudx 2 dudy 3 dudz
  //4 dydx 5 dydy 6 dydz
  //7 dwdx 8 dwdy 9 dwdz
  _blockFinDiff(output,input);
  return true;
}

////////////////////////SuperLatticeVelocityGradientFD3D//////////////////////////////////
template <typename T, typename DESCRIPTOR>
SuperLatticeVelocityGradientFD3D<T,DESCRIPTOR>::SuperLatticeVelocityGradientFD3D
(SuperGeometry<T,3>& sGeometry, SuperLattice<T,DESCRIPTOR>& sLattice, std::list<int>& matNumber) : SuperLatticeF3D<T,DESCRIPTOR>(sLattice,9),
  _sVelocity(sLattice), _sFinDiff(sGeometry, _sVelocity, matNumber)
{
  this->getName() = "VelocityGradientFD";
  int maxC = this->_superStructure.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back(new BlockLatticeVelocityGradientFD3D<T,DESCRIPTOR> (this->_sLattice.getBlock(iC), this->_sFinDiff.getBlockF(iC)));
  }
}

////////////////////////SuperLatticeExternalVelocityGradientFD3D//////////////////////////////////
template <typename T, typename DESCRIPTOR>
SuperLatticeExternalVelocityGradientFD3D<T,DESCRIPTOR>::SuperLatticeExternalVelocityGradientFD3D
(SuperGeometry<T,3>& sGeometry, SuperLattice<T,DESCRIPTOR>& sLattice, std::list<int>& matNumber) : SuperLatticeF3D<T,DESCRIPTOR>(sLattice,9),
  _sVelocity(sLattice), _sFinDiff(sGeometry, _sVelocity, matNumber)
{
  this->getName() = "externalVelocityGradientFD";
  int maxC = this->_superStructure.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back(new BlockLatticeExternalVelocityGradientFD3D<T,DESCRIPTOR> (this->_sLattice.getBlock(iC), this->_sFinDiff.getBlockF(iC)));
  }
}

////////////////////////BlockLatticePhysVelocityGradientFD3D//////////////////////////////////
template <typename T, typename DESCRIPTOR>
BlockLatticePhysVelocityGradientFD3D<T,DESCRIPTOR>::BlockLatticePhysVelocityGradientFD3D
(BlockLattice<T,DESCRIPTOR>& blockLattice, BlockF3D<T>& blockFinDiff, const UnitConverter<T,DESCRIPTOR>& converter)
  : BlockLatticeF3D<T,DESCRIPTOR>(blockLattice, 9), _blockFinDiff(blockFinDiff), _converter(converter)
{
  this->getName() = "PhysVelocityGradientFD";
}

template <typename T, typename DESCRIPTOR>
bool BlockLatticePhysVelocityGradientFD3D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  _blockFinDiff(output,input);
  return true;
}

////////////////////////SuperLatticePhysVelocityGradientFD3D//////////////////////////////////
template <typename T, typename DESCRIPTOR>
SuperLatticePhysVelocityGradientFD3D<T,DESCRIPTOR>::SuperLatticePhysVelocityGradientFD3D
(SuperGeometry<T,3>& sGeometry, SuperLattice<T,DESCRIPTOR>& sLattice, std::list<int>& matNumber, const UnitConverter<T,DESCRIPTOR>& converter) : SuperLatticeF3D<T,DESCRIPTOR>(sLattice,9),
  _sVelocity(sLattice, converter), _sFinDiff(sGeometry, _sVelocity, matNumber, converter), _converter(converter)
{
  this->getName() = "PhysVelocityGradientFD";
  int maxC = this->_superStructure.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back(new BlockLatticePhysVelocityGradientFD3D<T,DESCRIPTOR> (this->_sLattice.getBlock(iC), this->_sFinDiff.getBlockF(iC), this->_converter));
  }
}

////////////////////////BlockLatticeStrainRateFD3D//////////////////////////////////
template <typename T, typename DESCRIPTOR>
BlockLatticeStrainRateFD3D<T,DESCRIPTOR>::BlockLatticeStrainRateFD3D
(BlockLattice<T,DESCRIPTOR>& blockLattice, BlockF3D<T>& blockVeloGrad)
  : BlockLatticeF3D<T,DESCRIPTOR>(blockLattice, 9), _blockVeloGrad(blockVeloGrad)
{
  this->getName() = "StrainRateFD";
}

template <typename T, typename DESCRIPTOR>
bool BlockLatticeStrainRateFD3D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  T velograd[9];
  _blockVeloGrad(velograd,input);
  output[0] = velograd[0];
  output[1] = 0.5 * (velograd[1] + velograd[3]);
  output[2] = 0.5 * (velograd[2] + velograd[6]);
  output[3] = output[1];
  output[4] = velograd[4];
  output[5] = 0.5 * (velograd[5] + velograd[7]);
  output[6] = output[2];
  output[7] = output[5];
  output[8] = velograd[8];
  return true;
}

////////////////////////SuperLatticeStrainRateFD3D//////////////////////////////////
template <typename T, typename DESCRIPTOR>
SuperLatticeStrainRateFD3D<T,DESCRIPTOR>::SuperLatticeStrainRateFD3D
(SuperGeometry<T,3>& sGeometry, SuperLattice<T,DESCRIPTOR>& sLattice, std::list<int>& matNumber) : SuperLatticeF3D<T,DESCRIPTOR>(sLattice,9),
  _sVeloGrad(sGeometry, sLattice, matNumber)
{
  this->getName() = "StrainRateFD";
  int maxC = this->_superStructure.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back(new BlockLatticeStrainRateFD3D<T,DESCRIPTOR> (this->_sLattice.getBlock(iC), this->_sVeloGrad.getBlockF(iC), this->_converter));
  }
}

////////////////////////BlockLatticePhysStrainRateFD3D//////////////////////////////////
template <typename T, typename DESCRIPTOR>
BlockLatticePhysStrainRateFD3D<T,DESCRIPTOR>::BlockLatticePhysStrainRateFD3D
(BlockLattice<T,DESCRIPTOR>& blockLattice, BlockF3D<T>& blockVeloGrad, const UnitConverter<T,DESCRIPTOR>& converter)
  : BlockLatticeF3D<T,DESCRIPTOR>(blockLattice, 9), _blockVeloGrad(blockVeloGrad), _converter(converter)
{
  this->getName() = "PhysStrainRateFD";
}

template <typename T, typename DESCRIPTOR>
bool BlockLatticePhysStrainRateFD3D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  T velograd[9];
  _blockVeloGrad(velograd,input);
  output[0] = velograd[0];
  output[1] = 0.5 * (velograd[1] + velograd[3]);
  output[2] = 0.5 * (velograd[2] + velograd[6]);
  output[3] = output[1];
  output[4] = velograd[4];
  output[5] = 0.5 * (velograd[5] + velograd[7]);
  output[6] = output[2];
  output[7] = output[5];
  output[8] = velograd[8];

  return true;
}

////////////////////////SuperLatticePhysStrainRateFD3D//////////////////////////////////
template <typename T, typename DESCRIPTOR>
SuperLatticePhysStrainRateFD3D<T,DESCRIPTOR>::SuperLatticePhysStrainRateFD3D
(SuperGeometry<T,3>& sGeometry, SuperLattice<T,DESCRIPTOR>& sLattice, std::list<int>& matNumber, const UnitConverter<T,DESCRIPTOR>& converter) : SuperLatticeF3D<T,DESCRIPTOR>(sLattice,9),
  _sVeloGrad(sGeometry, sLattice, matNumber, converter), _converter(converter)
{
  this->getName() = "PhysStrainRateFD";
  int maxC = this->_superStructure.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back(new BlockLatticePhysStrainRateFD3D<T,DESCRIPTOR> (this->_sLattice.getBlock(iC), this->_sVeloGrad.getBlockF(iC), this->_converter));
  }
}

////////////////////////BlockLatticeDissipationFD3D//////////////////////////////////
template <typename T, typename DESCRIPTOR>
BlockLatticeDissipationFD3D<T,DESCRIPTOR>::BlockLatticeDissipationFD3D
(BlockLattice<T,DESCRIPTOR>& blockLattice, BlockF3D<T>& blockVeloGrad, const UnitConverter<T,DESCRIPTOR>& converter)
  : BlockLatticeF3D<T,DESCRIPTOR>(blockLattice, 1), _blockVeloGrad(blockVeloGrad), _converter(converter)
{
  this->getName() = "DissipationFD";
}

template <typename T, typename DESCRIPTOR>
bool BlockLatticeDissipationFD3D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  T velograd[9];
  _blockVeloGrad(velograd,input);
  output[0] = velograd[0] * velograd[0] + velograd[1] * velograd[1] + velograd[2] * velograd[2] +
              velograd[3] * velograd[3] + velograd[4] * velograd[4] + velograd[5] * velograd[5] +
              velograd[6] * velograd[6] + velograd[7] * velograd[7] + velograd[8] * velograd[8];
  output[0] *= _converter.getLatticeViscosity();

  return true;
}

////////////////////////SuperLatticeDissipationFD3D//////////////////////////////////
template <typename T, typename DESCRIPTOR>
SuperLatticeDissipationFD3D<T,DESCRIPTOR>::SuperLatticeDissipationFD3D
(SuperGeometry<T,3>& sGeometry, SuperLattice<T,DESCRIPTOR>& sLattice, std::list<int>& matNumber, const UnitConverter<T,DESCRIPTOR>& converter) : SuperLatticeF3D<T,DESCRIPTOR>(sLattice,1),
  _sVeloGrad(sGeometry, sLattice, matNumber, converter), _converter(converter)
{
  this->getName() = "DissipationFD";
  int maxC = this->_superStructure.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back(new BlockLatticeDissipationFD3D<T,DESCRIPTOR> (this->_sLattice.getBlock(iC), this->_sVeloGrad.getBlockF(iC), this->_converter));
  }
}

////////////////////////BlockLatticePhysDissipationFD3D//////////////////////////////////
template <typename T, typename DESCRIPTOR>
BlockLatticePhysDissipationFD3D<T,DESCRIPTOR>::BlockLatticePhysDissipationFD3D
(BlockLattice<T,DESCRIPTOR>& blockLattice, BlockF3D<T>& blockVeloGrad, const UnitConverter<T,DESCRIPTOR>& converter)
  : BlockLatticeF3D<T,DESCRIPTOR>(blockLattice, 1), _blockVeloGrad(blockVeloGrad), _converter(converter)
{
  this->getName() = "PhysDissipationFD";
}

template <typename T, typename DESCRIPTOR>
bool BlockLatticePhysDissipationFD3D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  T velograd[9];
  _blockVeloGrad(velograd,input);
  output[0] = velograd[0] * velograd[0] + velograd[1] * velograd[1] + velograd[2] * velograd[2] +
              velograd[3] * velograd[3] + velograd[4] * velograd[4] + velograd[5] * velograd[5] +
              velograd[6] * velograd[6] + velograd[7] * velograd[7] + velograd[8] * velograd[8];
  output[0] *= _converter.getPhysViscosity();

  return true;
}

////////////////////////SuperLatticePhysDissipationFD3D//////////////////////////////////
template <typename T, typename DESCRIPTOR>
SuperLatticePhysDissipationFD3D<T,DESCRIPTOR>::SuperLatticePhysDissipationFD3D
(SuperGeometry<T,3>& sGeometry, SuperLattice<T,DESCRIPTOR>& sLattice, std::list<int>& matNumber, const UnitConverter<T,DESCRIPTOR>& converter) : SuperLatticeF3D<T,DESCRIPTOR>(sLattice,1),
  _sVeloGrad(sGeometry, sLattice, matNumber, converter), _converter(converter)
{
  this->getName() = "PhysDissipationFD";
  int maxC = this->_superStructure.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back(new BlockLatticePhysDissipationFD3D<T,DESCRIPTOR> (this->_sLattice.getBlock(iC), this->_sVeloGrad.getBlockF(iC), this->_converter));
  }
}

////////////////////////BlockLatticeEffectiveDissipationFD3D//////////////////////////////////
template <typename T, typename DESCRIPTOR>
BlockLatticeEffectiveDissipationFD3D<T,DESCRIPTOR>::BlockLatticeEffectiveDissipationFD3D
(BlockLattice<T,DESCRIPTOR>& blockLattice, BlockF3D<T>& blockVeloGrad, const UnitConverter<T,DESCRIPTOR>& converter, LESDynamics<T, DESCRIPTOR>& LESdynamics)
  : BlockLatticeF3D<T,DESCRIPTOR>(blockLattice, 1), _blockVeloGrad(blockVeloGrad), _converter(converter), _LESdynamics(LESdynamics)
{
  this->getName() = "EffectiveDissipationFD";
}

template <typename T, typename DESCRIPTOR>
bool BlockLatticeEffectiveDissipationFD3D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  T velograd[9];
  _blockVeloGrad(velograd,input);
  output[0] = velograd[0] * velograd[0] + velograd[1] * velograd[1] + velograd[2] * velograd[2] +
              velograd[3] * velograd[3] + velograd[4] * velograd[4] + velograd[5] * velograd[5] +
              velograd[6] * velograd[6] + velograd[7] * velograd[7] + velograd[8] * velograd[8];

  auto cell = this->_blockLattice.get(input[0], input[1], input[2]);
  T omegaEff = _LESdynamics.getEffectiveOmega(cell);
  T nuEff = ((1./omegaEff)-0.5)/descriptors::invCs2<T,DESCRIPTOR>();
  output[0] *= nuEff;

  return true;
}

////////////////////////SuperLatticeEffectiveDissipationFD3D//////////////////////////////////
template <typename T, typename DESCRIPTOR>
SuperLatticeEffectiveDissipationFD3D<T,DESCRIPTOR>::SuperLatticeEffectiveDissipationFD3D
(SuperGeometry<T,3>& sGeometry, SuperLattice<T,DESCRIPTOR>& sLattice, std::list<int>& matNumber, const UnitConverter<T,DESCRIPTOR>& converter, LESDynamics<T, DESCRIPTOR>& LESdynamics)
  : SuperLatticeF3D<T,DESCRIPTOR>(sLattice,1), _sVeloGrad(sGeometry, sLattice, matNumber, converter), _converter(converter)
{
  this->getName() = "EffectiveDissipationFD";
  int maxC = this->_superStructure.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back(new BlockLatticeEffectiveDissipationFD3D<T,DESCRIPTOR> (this->_sLattice.getBlock(iC), this->_sVeloGrad.getBlockF(iC), this->_converter));
  }
}

////////////////////////BlockLatticePhysEffectiveDissipationFD3D//////////////////////////////////
template <typename T, typename DESCRIPTOR>
BlockLatticePhysEffectiveDissipationFD3D<T,DESCRIPTOR>::BlockLatticePhysEffectiveDissipationFD3D(
  BlockLattice<T,DESCRIPTOR>& blockLattice,
  BlockF3D<T>& blockVeloGrad,
  const UnitConverter<T,DESCRIPTOR>& converter,
  std::function<T(Cell<T,DESCRIPTOR>&)> effectiveOmegaF)
  : BlockLatticeF3D<T,DESCRIPTOR>(blockLattice, 1),
    _blockVeloGrad(blockVeloGrad),
    _converter(converter),
    _effectiveOmegaF(effectiveOmegaF)
{
  this->getName() = "PhysEffectiveDissipationFD";
}

template <typename T, typename DESCRIPTOR>
bool BlockLatticePhysEffectiveDissipationFD3D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  T velograd[9];
  _blockVeloGrad(velograd,input);
  output[0] = velograd[0] * velograd[0] + velograd[1] * velograd[1] + velograd[2] * velograd[2] +
              velograd[3] * velograd[3] + velograd[4] * velograd[4] + velograd[5] * velograd[5] +
              velograd[6] * velograd[6] + velograd[7] * velograd[7] + velograd[8] * velograd[8];

  auto cell = this->_blockLattice.get(input[0], input[1], input[2]);
  T omegaEff = _effectiveOmegaF(cell);
  T nuEff = ((1./omegaEff)-0.5)/descriptors::invCs2<T,DESCRIPTOR>();
  output[0] *= _converter.getPhysViscosity( nuEff );

  return true;
}

////////////////////////SuperLatticePhysEffectiveDissipationFD3D//////////////////////////////////
template <typename T, typename DESCRIPTOR>
SuperLatticePhysEffectiveDissipationFD3D<T,DESCRIPTOR>::SuperLatticePhysEffectiveDissipationFD3D
(SuperGeometry<T,3>& sGeometry, SuperLattice<T,DESCRIPTOR>& sLattice, std::list<int>& matNumber, const UnitConverter<T,DESCRIPTOR>& converter,
  std::function<T(Cell<T,DESCRIPTOR>&)> effectiveOmegaF)
  : SuperLatticeF3D<T,DESCRIPTOR>(sLattice,1), _sVeloGrad(sGeometry, sLattice, matNumber, converter), _converter(converter)
{
  this->getName() = "PhysEffectiveDissipationFD";
  int maxC = this->_superStructure.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back(new BlockLatticePhysEffectiveDissipationFD3D<T,DESCRIPTOR> (this->_sLattice.getBlock(iC), this->_sVeloGrad.getBlockF(iC), this->_converter, effectiveOmegaF));
  }
}

////////////////////////BlockLatticeVorticityFD3D//////////////////////////////////
template <typename T, typename DESCRIPTOR>
BlockLatticeVorticityFD3D<T,DESCRIPTOR>::BlockLatticeVorticityFD3D
(BlockLattice<T,DESCRIPTOR>& blockLattice, BlockF3D<T>& blockVeloGrad)
  : BlockLatticeF3D<T,DESCRIPTOR>(blockLattice, 3), _blockVeloGrad(blockVeloGrad)
{
  this->getName() = "VorticityFD";
}

template <typename T, typename DESCRIPTOR>
bool BlockLatticeVorticityFD3D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  T velograd[9];
  _blockVeloGrad(velograd,input);
  output[0] = velograd[7] - velograd[5];
  output[1] = velograd[2] - velograd[6];
  output[2] = velograd[3] - velograd[1];
  return true;
}

////////////////////////SuperLatticeVorticityFD3D//////////////////////////////////
template <typename T, typename DESCRIPTOR>
SuperLatticeVorticityFD3D<T,DESCRIPTOR>::SuperLatticeVorticityFD3D
(SuperGeometry<T,3>& sGeometry, SuperLattice<T,DESCRIPTOR>& sLattice, std::list<int>& matNumber) : SuperLatticeF3D<T,DESCRIPTOR>(sLattice,3),
  _sVeloGrad(sGeometry, sLattice, matNumber)
{
  this->getName() = "VorticityFD";
  int maxC = this->_superStructure.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back(new BlockLatticeVorticityFD3D<T,DESCRIPTOR> (this->_sLattice.getBlock(iC), this->_sVeloGrad.getBlockF(iC), this->_converter));
  }
}

////////////////////////BlockLatticePhysVorticityFD3D//////////////////////////////////
template <typename T, typename DESCRIPTOR>
BlockLatticePhysVorticityFD3D<T,DESCRIPTOR>::BlockLatticePhysVorticityFD3D
(BlockLattice<T,DESCRIPTOR>& blockLattice, BlockF3D<T>& blockVeloGrad, const UnitConverter<T,DESCRIPTOR>& converter)
  : BlockLatticeF3D<T,DESCRIPTOR>(blockLattice, 3), _blockVeloGrad(blockVeloGrad), _converter(converter)
{
  this->getName() = "PhysVorticityFD";
}

template <typename T, typename DESCRIPTOR>
bool BlockLatticePhysVorticityFD3D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  T velograd[9];
  _blockVeloGrad(velograd,input);
  output[0] = velograd[7] - velograd[5];
  output[1] = velograd[2] - velograd[6];
  output[2] = velograd[3] - velograd[1];
  return true;
}

////////////////////////SuperLatticePhysVorticityFD3D//////////////////////////////////
template <typename T, typename DESCRIPTOR>
SuperLatticePhysVorticityFD3D<T,DESCRIPTOR>::SuperLatticePhysVorticityFD3D
(SuperGeometry<T,3>& sGeometry, SuperLattice<T,DESCRIPTOR>& sLattice, std::list<int>& matNumber, const UnitConverter<T,DESCRIPTOR>& converter) : SuperLatticeF3D<T,DESCRIPTOR>(sLattice,3),
  _sVeloGrad(sGeometry, sLattice, matNumber, converter), _converter(converter)
{
  this->getName() = "PhysVorticityFD";
  int maxC = this->_superStructure.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back(new BlockLatticePhysVorticityFD3D<T,DESCRIPTOR> (this->_sLattice.getBlock(iC), this->_sVeloGrad.getBlockF(iC), this->_converter));
  }
}

////////////////////////BlockLatticePhysStressFD3D//////////////////////////////////
template <typename T, typename DESCRIPTOR>
BlockLatticePhysStressFD3D<T,DESCRIPTOR>::BlockLatticePhysStressFD3D
(BlockLattice<T,DESCRIPTOR>& blockLattice, BlockF3D<T>& blockStrainRate, const UnitConverter<T,DESCRIPTOR>& converter)
  : BlockLatticeF3D<T,DESCRIPTOR>(blockLattice, 9), _blockStrainRate(blockStrainRate), _converter(converter)
{
  this->getName() = "PhysStressFD";
}

template <typename T, typename DESCRIPTOR>
bool BlockLatticePhysStressFD3D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  _blockStrainRate(output,input);
  for (int i = 0; i < 9; i++) {
    output[i] /= _converter.getPhysViscosity() * _converter.getPhysDensity();
  }
  return true;
}

////////////////////////SuperLatticePhysStressFD3D//////////////////////////////////
template <typename T, typename DESCRIPTOR>
SuperLatticePhysStressFD3D<T,DESCRIPTOR>::SuperLatticePhysStressFD3D
(SuperGeometry<T,3>& sGeometry, SuperLattice<T,DESCRIPTOR>& sLattice, std::list<int>& matNumber, const UnitConverter<T,DESCRIPTOR>& converter) : SuperLatticeF3D<T,DESCRIPTOR>(sLattice,9),
  _sStrainRate(sGeometry, sLattice, matNumber, converter), _converter(converter)
{
  this->getName() = "PhysStressFD";
  int maxC = this->_superStructure.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back(new BlockLatticePhysStressFD3D<T,DESCRIPTOR> (this->_sLattice.getBlock(iC), this->_sStrainRate.getBlockF(iC), this->_converter));
  }
}

////////////////////////BlockIsotropicHomogeneousTKE//////////////////////////////////
template<typename T, typename DESCRIPTOR>
BlockIsotropicHomogeneousTKE3D<T, DESCRIPTOR>::BlockIsotropicHomogeneousTKE3D(BlockLattice<T,DESCRIPTOR>& blockLattice, BlockF3D<T>& blockVelocity)
  : BlockLatticeF3D<T,DESCRIPTOR>(blockLattice, 1), _blockVelocity(blockVelocity)
{
  this->getName() = "IsotropicHomogeneousTKE";
}

template<typename T, typename DESCRIPTOR>
bool BlockIsotropicHomogeneousTKE3D<T, DESCRIPTOR>::operator()(T output[], const int input[])
{
  output[0] = T();
  T data[_blockVelocity.getTargetDim()];
  _blockVelocity(data,input);
  for (int i = 0; i < _blockVelocity.getTargetDim(); ++i) {
    output[0] +=  data[i] * data[i];
  }
  output[0] = 0.5 * output[0];
  return true;
}

////////////////////////SuperIsotropicHomogeneousTKE//////////////////////////////////
template<typename T, typename DESCRIPTOR>
SuperIsotropicHomogeneousTKE3D<T, DESCRIPTOR>::SuperIsotropicHomogeneousTKE3D(
  SuperLattice<T,DESCRIPTOR>& sLattice, const  UnitConverter<T,DESCRIPTOR>& converter)
  : SuperLatticeF3D<T,DESCRIPTOR>(sLattice, 1), _sVelocity(sLattice, converter), _converter(converter)

{
  this->getName() = "IsotropicHomogeneousTKE";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back(new BlockIsotropicHomogeneousTKE3D<T,DESCRIPTOR>(this->_sLattice.getBlock(iC), this-> _sVelocity.getBlockF(iC)));
  }
}

////////////////////////BlockLatticePhysEnstrophyFD3D//////////////////////////////////
template <typename T, typename DESCRIPTOR>
BlockLatticePhysEnstrophyFD3D<T,DESCRIPTOR>::BlockLatticePhysEnstrophyFD3D
(BlockLattice<T,DESCRIPTOR>& blockLattice, BlockF3D<T>& blockVeloGrad, const  UnitConverter<T,DESCRIPTOR>& converter)
  : BlockLatticeF3D<T,DESCRIPTOR>(blockLattice, 1), _blockVeloGrad(blockVeloGrad), _converter(converter)
{
  this->getName() = "PhysEnstrophyFD";
}

template <typename T, typename DESCRIPTOR>
bool BlockLatticePhysEnstrophyFD3D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  T velograd[9];
  _blockVeloGrad(velograd,input);
  output[0] = 0.5 * ( util::pow(velograd[7] - velograd[5], 2) + util::pow(velograd[2] - velograd[6], 2) + util::pow(velograd[3] - velograd[1], 2) );
  return true;
}

////////////////////////SuperLatticePhysEnstrophyFD3D//////////////////////////////////
template <typename T, typename DESCRIPTOR>
SuperLatticePhysEnstrophyFD3D<T,DESCRIPTOR>::SuperLatticePhysEnstrophyFD3D
(SuperGeometry<T,3>& sGeometry, SuperLattice<T,DESCRIPTOR>& sLattice, std::list<int>& matNumber, const  UnitConverter<T,DESCRIPTOR>& converter) : SuperLatticeF3D<T,DESCRIPTOR>(sLattice, 1),
  _sVeloGrad(sGeometry, sLattice, matNumber, converter), _converter(converter)
{
  this->getName() = "PhysEnstrophyFD";
  int maxC = this->_superStructure.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back(new BlockLatticePhysEnstrophyFD3D<T,DESCRIPTOR> (this->_sLattice.getBlock(iC), this->_sVeloGrad.getBlockF(iC), this->_converter));
  }
}


} // end namespace olb
#endif
