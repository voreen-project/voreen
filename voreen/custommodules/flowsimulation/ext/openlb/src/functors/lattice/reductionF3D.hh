/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012-2017 Lukas Baron, Tim Dornieden, Mathias J. Krause,
 *  Albert Mink, Fabian Klemens, Benjamin Förster, Marie-Luise Maier,
 *  Adrian Kummerländer
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

#ifndef REDUCTION_F_3D_HH
#define REDUCTION_F_3D_HH

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include <algorithm>
#include "functors/lattice/reductionF3D.h"
#include "dynamics/lbHelpers.h"  // for computation of lattice rho and velocity

namespace olb {

template<typename T, template<typename U> class DESCRIPTOR>
SuperLatticeFfromAnalyticalF3D<T, DESCRIPTOR>::SuperLatticeFfromAnalyticalF3D(
  AnalyticalF3D<T, T>& f, SuperLattice3D<T, DESCRIPTOR>& sLattice)
  : SuperLatticeF3D<T, DESCRIPTOR>(sLattice, f.getTargetDim()),
    _f(f)
{
  this->getName() = "fromAnalyticalF(" + _f.getName() + ")";

  LoadBalancer<T>&     load   = sLattice.getLoadBalancer();
  CuboidGeometry3D<T>& cuboid = sLattice.getCuboidGeometry();

  for (int iC = 0; iC < load.size(); ++iC) {
    this->_blockF.emplace_back(
      new BlockLatticeFfromAnalyticalF3D<T,DESCRIPTOR>(
        _f,
        sLattice.getExtendedBlockLattice(iC),
        cuboid.get(load.glob(iC)),
        sLattice.getOverlap())
    );
  }
}

template<typename T, template<typename U> class DESCRIPTOR>
bool SuperLatticeFfromAnalyticalF3D<T, DESCRIPTOR>::operator()(
  T output[], const int input[])
{
  T physR[3] = {};
  this->_sLattice.getCuboidGeometry().getPhysR(physR,input);
  return _f(output,physR);
}


template<typename T, template<typename U> class DESCRIPTOR>
BlockLatticeFfromAnalyticalF3D<T, DESCRIPTOR>::BlockLatticeFfromAnalyticalF3D(
  AnalyticalF3D<T, T>&                    f,
  BlockLatticeStructure3D<T, DESCRIPTOR>& lattice,
  Cuboid3D<T>&                            cuboid,
  int                                     overlap)
  : BlockLatticeF3D<T, DESCRIPTOR>(lattice, f.getTargetDim()),
    _f(f),
    _cuboid(cuboid, overlap),
    _overlap(overlap)
{
  this->getName() = "blockFfromAnalyticalF(" + _f.getName() + ")";
}

template<typename T, template<typename U> class DESCRIPTOR>
bool BlockLatticeFfromAnalyticalF3D<T, DESCRIPTOR>::operator()(
  T output[], const int input[])
{
  const int blockInput[3] = {
    input[0] + _overlap,
    input[1] + _overlap,
    input[2] + _overlap,
  };

  OLB_PRECONDITION(blockInput[0] < _cuboid.getNx());
  OLB_PRECONDITION(blockInput[1] < _cuboid.getNy());
  OLB_PRECONDITION(blockInput[2] < _cuboid.getNz());

  T physR[3] = {};
  _cuboid.getPhysR(physR,blockInput);
  return _f(output,physR);
}


//////////// not yet working // symbolically ///////////////////
////////////////////////////////////////////////
template<typename T, template<typename U> class DESCRIPTOR>
SmoothBlockIndicator3D<T, DESCRIPTOR>::SmoothBlockIndicator3D(
  IndicatorF3D<T>& f, T h )
  : BlockDataF3D<T, T>((int)((f.getMax()[0] - f.getMin()[0] + 4*h) / h),
                       (int)((f.getMax()[1] - f.getMin()[1] + 4*h) / h),
                       (int)((f.getMax()[2] - f.getMin()[2] + 4*h) / h)),
    _f(f),
    _h(h)
{
  this->getName() = "SmoothBlockIndicator3D";

  for (int iX=0; iX<this->getBlockData().getNx(); iX++) {
    for (int iY=0; iY<this->getBlockData().getNy(); iY++) {
      for (int iZ=0; iZ<this->getBlockData().getNz(); iZ++) {
        T value = T();
        bool output[1];
        //     for (int iPop=0; iPop<DESCRIPTOR<T>::q; iPop++) {
        //       T input[] = {_f.getMin()[0] - 2*_h + (iX + DESCRIPTOR<T>::c[iPop][0])
        //                    *_h, _f.getMin()[1] - 2*_h + (iY + DESCRIPTOR<T>::c[iPop][1])*_h,
        //                    _f.getMin()[2] - 2*_h + (iZ + DESCRIPTOR<T>::c[iPop][2])*_h
        //                   };
        T input[] = {_f.getMin()[0] + (iX)
                     *_h-2.*h, _f.getMin()[1] + (iY)*_h-2.*h,
                     _f.getMin()[2] + (iZ)*_h-2.*h
                    };
        _f(output,input);
        //        value += (int)output[0]*DESCRIPTOR<T>::t[iPop];
        //      }
        value = output[0];
        //std::cout << iX << "/" << iY << "/" << iZ  << std::endl;
        this->getBlockData().get(iX,iY,iZ,0) = value;
      }
    }
  }
}
/*
template<typename T, template<typename U> class DESCRIPTOR>
bool SmoothBlockIndicator3D<T, DESCRIPTOR>::operator()(
  T output[], const int input[])
{
  T physR[3] = {};
  _superGeometry.getPhysR(physR,input[0],input[1],input[2] );
  _f(output,physR);
  return true;
}*/


template<typename T, template<typename U> class DESCRIPTOR>
SuperLatticeInterpPhysVelocity3Degree3D<T, DESCRIPTOR>::
SuperLatticeInterpPhysVelocity3Degree3D(
  SuperLattice3D<T, DESCRIPTOR>& sLattice, UnitConverter<T,DESCRIPTOR>& conv, int range)
  : SuperLatticeF3D<T, DESCRIPTOR>(sLattice, 3)
{
  this->getName() = "Interp3DegreeVelocity";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    BlockLatticeInterpPhysVelocity3Degree3D<T, DESCRIPTOR>* foo =
      new BlockLatticeInterpPhysVelocity3Degree3D<T, DESCRIPTOR>(
      sLattice.getExtendedBlockLattice(iC),
      conv,
      &sLattice.getCuboidGeometry().get(this->_sLattice.getLoadBalancer().
                                        glob(iC)),
      sLattice.getOverlap(),
      range);
    _bLattices.push_back(foo);
  }
}

template<typename T, template<typename U> class DESCRIPTOR>
void SuperLatticeInterpPhysVelocity3Degree3D<T, DESCRIPTOR>::operator()(
  T output[], const T input[], const int iC)
{
  _bLattices[this->_sLattice.getLoadBalancer().loc(iC)]->operator()(output,
      input);
}

template<typename T, template<typename U> class DESCRIPTOR>
BlockLatticeInterpPhysVelocity3Degree3D<T, DESCRIPTOR>::
BlockLatticeInterpPhysVelocity3Degree3D(
  BlockLatticeStructure3D<T, DESCRIPTOR>& blockLattice, UnitConverter<T,DESCRIPTOR>& conv,
  Cuboid3D<T>* c, int overlap, int range)
  : BlockLatticeF3D<T, DESCRIPTOR>(blockLattice, 3),
    _conv(conv),
    _cuboid(c),
    _overlap(overlap),
    _range(range)
{
  this->getName() = "BlockLatticeInterpVelocity3Degree3D";
}

template<typename T, template<typename U> class DESCRIPTOR>
BlockLatticeInterpPhysVelocity3Degree3D<T, DESCRIPTOR>::
BlockLatticeInterpPhysVelocity3Degree3D(
  const BlockLatticeInterpPhysVelocity3Degree3D<T, DESCRIPTOR>& rhs) :
  BlockLatticeF3D<T, DESCRIPTOR>(rhs._blockLattice, 3),
  _conv(rhs._conv),
  _cuboid(rhs._cuboid),
  _overlap(rhs._overlap),
  _range(rhs._range)
{
}

template<typename T, template<typename U> class DESCRIPTOR>
void BlockLatticeInterpPhysVelocity3Degree3D<T, DESCRIPTOR>::operator()(
  T output[3], const T input[3])
{
  T u[3], rho, volume;
  int latIntPos[3] = {0};
  T latPhysPos[3] = {T()};
  _cuboid->getFloorLatticeR(latIntPos, &input[0]);
  _cuboid->getPhysR(latPhysPos, latIntPos);

  latIntPos[0] += _overlap;
  latIntPos[1] += _overlap;
  latIntPos[2] += _overlap;

  volume=T(1);
  for (int i = -_range; i <= _range+1; ++i) {
    for (int j = -_range; j <= _range+1; ++j) {
      for (int k = -_range; k <= _range+1; ++k) {

        this->_blockLattice.get(latIntPos[0]+i, latIntPos[1]+j,
                                latIntPos[2]+k).computeRhoU(rho, u);
        for (int l = -_range; l <= _range+1; ++l) {
          if (l != i) {
            volume *= (input[0] - (latPhysPos[0]+ l *_cuboid->getDeltaR()))
                      / (latPhysPos[0] + i *_cuboid->getDeltaR()
                         - (latPhysPos[0] + l *_cuboid->getDeltaR()));
          }
        }
        for (int m = -_range; m <= _range+1; ++m) {
          if (m != j) {
            volume *= (input[1]
                       - (latPhysPos[1] + m *_cuboid->getDeltaR()))
                      / (latPhysPos[1] + j * _cuboid->getDeltaR()
                         - (latPhysPos[1] + m * _cuboid->getDeltaR()));
          }
        }
        for (int n = -_range; n <= _range+1; ++n) {
          if (n != k) {
            volume *= (input[2]
                       - (latPhysPos[2] + n * _cuboid->getDeltaR()))
                      / (latPhysPos[2] + k * _cuboid->getDeltaR()
                         - (latPhysPos[2] + n * _cuboid->getDeltaR()));
          }
        }
        output[0] += u[0] * volume;
        output[1] += u[1] * volume;
        output[2] += u[2] * volume;
        volume=T(1);
      }
    }
  }

  output[0] = _conv.getPhysVelocity(output[0]);
  output[1] = _conv.getPhysVelocity(output[1]);
  output[2] = _conv.getPhysVelocity(output[2]);
}


template<typename T, template<typename U> class DESCRIPTOR>
SuperLatticeInterpDensity3Degree3D<T, DESCRIPTOR>::SuperLatticeInterpDensity3Degree3D(
  SuperLattice3D<T, DESCRIPTOR>& sLattice, SuperGeometry3D<T>& sGeometry,
  UnitConverter<T,DESCRIPTOR>& conv, int range) :
  SuperLatticeF3D<T, DESCRIPTOR>(sLattice, 3)
{
  this->getName() = "Interp3DegreeDensity";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int lociC = 0; lociC < maxC; lociC++) {
    int globiC = this->_sLattice.getLoadBalancer().glob(lociC);

    BlockLatticeInterpDensity3Degree3D<T, DESCRIPTOR>* foo =
      new BlockLatticeInterpDensity3Degree3D<T, DESCRIPTOR>(
      sLattice.getExtendedBlockLattice(lociC),
      sGeometry.getExtendedBlockGeometry(lociC),
      conv,
      &sLattice.getCuboidGeometry().get(globiC),
      sLattice.getOverlap(), range);
    _bLattices.push_back(foo);

    if (sLattice.getOverlap() <= range + 1)
      std::cout << "lattice overlap has to be larger than (range + 1)"
                << std::endl;
  }
}

template<typename T, template<typename U> class DESCRIPTOR>
SuperLatticeInterpDensity3Degree3D<T, DESCRIPTOR>::~SuperLatticeInterpDensity3Degree3D()
{
  // first deconstruct vector elements
  for ( auto it : _bLattices) {
    delete it;
  }
  // then delete std::vector
  _bLattices.clear();
}

template<typename T, template<typename U> class DESCRIPTOR>
void SuperLatticeInterpDensity3Degree3D<T, DESCRIPTOR>::operator()(T output[],
    const T input[], const int globiC)
{
  if (this->_sLattice.getLoadBalancer().rank(globiC) == singleton::mpi().getRank()) {
    _bLattices[this->_sLattice.getLoadBalancer().loc(globiC)]->operator()(output,
        input);
  }
}

template<typename T, template<typename U> class DESCRIPTOR>
BlockLatticeInterpDensity3Degree3D<T, DESCRIPTOR>::BlockLatticeInterpDensity3Degree3D(
  BlockLatticeStructure3D<T, DESCRIPTOR>& blockLattice,
  BlockGeometryStructure3D<T>& blockGeometry, UnitConverter<T,DESCRIPTOR>& conv,
  Cuboid3D<T>* c, int overlap, int range) :
  BlockLatticeF3D<T, DESCRIPTOR>(blockLattice, 3), _blockGeometry(blockGeometry),
  _conv(conv), _cuboid(c), _overlap(overlap), _range(range)
{
  this->getName() = "BlockLatticeInterpDensity3Degree3D";
}

template<typename T, template<typename U> class DESCRIPTOR>
BlockLatticeInterpDensity3Degree3D<T, DESCRIPTOR>::BlockLatticeInterpDensity3Degree3D(
  const BlockLatticeInterpDensity3Degree3D<T, DESCRIPTOR>& rhs) :
  BlockLatticeF3D<T, DESCRIPTOR>(rhs._blockLattice, 3),
  _blockGeometry(rhs._blockGeometry),_conv(rhs._conv), _cuboid(
    rhs._cuboid), _overlap(rhs._overlap), _range(rhs._range)
{
}

template<typename T, template<typename U> class DESCRIPTOR>
void BlockLatticeInterpDensity3Degree3D<T, DESCRIPTOR>::operator()(
  T output[DESCRIPTOR<T>::q], const T input[3])
{
  T volume = T(1);
  T f_iPop = 0.;
  /** neighbor position on grid of input value in lattice units
   *referred to local cuboid
   */
  int latIntPos[3] = { 0 };
  // neighbor position on grid of input value in physical units
  T latPhysPos[3] = { T() };
  // input is physical position on grid
  _cuboid->getFloorLatticeR(latIntPos, input);
  // latPhysPos is global physical position on geometry
  _cuboid->getPhysR(latPhysPos, latIntPos);

  // point on cuboid equals cell on blocklattice (extended) shifted by _overlap
  latIntPos[0] += _overlap;
  latIntPos[1] += _overlap;
  latIntPos[2] += _overlap;

  for (unsigned iPop = 0; iPop < DESCRIPTOR<T>::q; ++iPop) {
    output[iPop] = T(0);
    for (int i = -_range; i <= _range + 1; ++i) {
      for (int j = -_range; j <= _range + 1; ++j) {
        for (int k = -_range; k <= _range + 1; ++k) {
          f_iPop = 0.;
          // just if material of cell != 1 there may be information of fluid density
          if (_blockGeometry.getMaterial(latIntPos[0] + i, latIntPos[1] + j,
                                         latIntPos[2] + k) != 0) {
            // because of communication it is possible to get density information
            // from neighboring cuboid
            f_iPop = this->_blockLattice.get(latIntPos[0] + i, latIntPos[1] + j,
                                             latIntPos[2] + k)[iPop];
          }
          for (int l = -_range; l <= _range + 1; ++l) {
            if (l != i) {
              volume *= (input[0] - (latPhysPos[0] + l * _cuboid->getDeltaR()))
                        / (latPhysPos[0] + i * _cuboid->getDeltaR()
                           - (latPhysPos[0] + l * _cuboid->getDeltaR()));
            }
          }
          for (int m = -_range; m <= _range + 1; ++m) {
            if (m != j) {
              volume *= (input[1] - (latPhysPos[1] + m * _cuboid->getDeltaR()))
                        / (latPhysPos[1] + j * _cuboid->getDeltaR()
                           - (latPhysPos[1] + m * _cuboid->getDeltaR()));
            }
          }
          for (int n = -_range; n <= _range + 1; ++n) {
            if (n != k) {
              volume *= (input[2] - (latPhysPos[2] + n * _cuboid->getDeltaR()))
                        / (latPhysPos[2] + k * _cuboid->getDeltaR()
                           - (latPhysPos[2] + n * _cuboid->getDeltaR()));
            }
          }
          output[iPop] += f_iPop * volume;
          volume = T(1);
        }
      }
    }
  }
}

template<typename T, template<typename U> class DESCRIPTOR>
SuperLatticeSmoothDiracDelta3D<T, DESCRIPTOR>::SuperLatticeSmoothDiracDelta3D(
  SuperLattice3D<T, DESCRIPTOR>& sLattice,
  UnitConverter<T,DESCRIPTOR>& conv, SuperGeometry3D<T>& sGeometry) :
  SuperLatticeF3D<T, DESCRIPTOR>(sLattice, 3)
{
  this->getName() = "SuperLatticeSmoothDiracDelta3D";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int lociC = 0; lociC < maxC; lociC++) {
    int globiC = this->_sLattice.getLoadBalancer().glob(lociC);

    BlockLatticeSmoothDiracDelta3D<T, DESCRIPTOR>* foo =
      new BlockLatticeSmoothDiracDelta3D<T, DESCRIPTOR>(
      sLattice.getExtendedBlockLattice(lociC),
      conv, &sLattice.getCuboidGeometry().get(globiC)
    );
    _bLattices.push_back(foo);
  }
}

template<typename T, template<typename U> class DESCRIPTOR>
SuperLatticeSmoothDiracDelta3D<T, DESCRIPTOR>::~SuperLatticeSmoothDiracDelta3D()
{
  for ( auto it : _bLattices) {
    delete it;
  }
  _bLattices.clear();
}

template<typename T, template<typename U> class DESCRIPTOR>
void SuperLatticeSmoothDiracDelta3D<T, DESCRIPTOR>::operator()(T delta[4][4][4],
    const T physPos[3], const int globiC)
{
  if (this->_sLattice.getLoadBalancer().rank(globiC) == singleton::mpi().getRank()) {
    _bLattices[this->_sLattice.getLoadBalancer().loc(globiC)]->operator()(delta,
        physPos);
  }
}

template<typename T, template<typename U> class DESCRIPTOR>
BlockLatticeSmoothDiracDelta3D<T, DESCRIPTOR>::BlockLatticeSmoothDiracDelta3D(
  BlockLattice3D<T, DESCRIPTOR>& blockLattice, UnitConverter<T,DESCRIPTOR>& conv, Cuboid3D<T>* cuboid)
  : BlockLatticeF3D<T, DESCRIPTOR>(blockLattice, 3), _conv(conv), _cuboid(cuboid)
{
  this->getName() = "BlockLatticeSmoothDiracDelta3D";
}

template<typename T, template<typename U> class DESCRIPTOR>
BlockLatticeSmoothDiracDelta3D<T, DESCRIPTOR>::BlockLatticeSmoothDiracDelta3D(
  const BlockLatticeSmoothDiracDelta3D<T, DESCRIPTOR>& rhs)
  :
  BlockLatticeF3D<T, DESCRIPTOR>(rhs._blockLattice, 3), _conv(rhs._conv), _cuboid(rhs._cuboid)
{
}

template<typename T, template<typename U> class DESCRIPTOR>
void BlockLatticeSmoothDiracDelta3D<T, DESCRIPTOR>::operator()(
  T delta[4][4][4], const T physPos[])
{
  int range = 1;
  T a, b, c = T();
  int latticeRoundedPosP[3] = { 0 };
  T physRoundedPosP[3] = { T() };
  T physLatticeL = _conv.getConversionFactorLength();

  T counter = 0.;

  _cuboid->getLatticeR(latticeRoundedPosP, physPos);
  _cuboid->getPhysR(physRoundedPosP, latticeRoundedPosP);

  for (int i = -range; i <= range + 1; ++i) {
    for (int j = -range; j <= range + 1; ++j) {
      for (int k = -range; k <= range + 1; ++k) {
        delta[i+range][j+range][k+range] = T(1);
        // a, b, c in lattice units cause physical ones get cancelled
        a = (physRoundedPosP[0] + i * physLatticeL - physPos[0])
            / physLatticeL;
        b =  (physRoundedPosP[1] + j * physLatticeL - physPos[1])
             / physLatticeL;
        c = (physRoundedPosP[2] + k * physLatticeL - physPos[2])
            / physLatticeL;

        // the for loops already define that a, b, c are smaller than 2
        delta[i+range][j+range][k+range] *= 1. / 4 * (1 + cos(M_PI * a / 2.));
        delta[i+range][j+range][k+range] *= 1. / 4 * (1 + cos(M_PI * b / 2.));
        delta[i+range][j+range][k+range] *= 1. / 4 * (1 + cos(M_PI * c / 2.));

        counter += delta[i+range][j+range][k+range];
      }
    }
  }

  //  if (!util::nearZero(counter - T(1))){
  //    // sum of delta has to be one
  //    std::cout << "[" << this->getName() << "] " <<
  //        "Delta summed up does not equal 1 but = " <<
  //        counter << std::endl;
  //  }

}


}  // end namespace olb

#endif
