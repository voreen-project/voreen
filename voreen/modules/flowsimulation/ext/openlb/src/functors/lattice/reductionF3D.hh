/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012-2017 Lukas Baron, Tim Dornieden, Mathias J. Krause,
 *  Albert Mink, Fabian Klemens, Benjamin Förster, Marie-Luise Maier,
 *  Adrian Kummerlaender
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
#include "dynamics/lbm.h"  // for computation of lattice rho and velocity

namespace olb {

template<typename T, typename DESCRIPTOR>
SuperLatticeFfromAnalyticalF3D<T, DESCRIPTOR>::SuperLatticeFfromAnalyticalF3D(
  FunctorPtr<AnalyticalF3D<T,T>>&& f,
  SuperLattice<T, DESCRIPTOR>&   sLattice)
  : SuperLatticeF3D<T, DESCRIPTOR>(sLattice, f->getTargetDim()),
    _f(std::move(f))
{
  this->getName() = "fromAnalyticalF(" + _f->getName() + ")";

  LoadBalancer<T>&     load   = sLattice.getLoadBalancer();
  CuboidGeometry3D<T>& cuboid = sLattice.getCuboidGeometry();

  for (int iC = 0; iC < load.size(); ++iC) {
    this->_blockF.emplace_back(
      new BlockLatticeFfromAnalyticalF3D<T,DESCRIPTOR>(
        *_f,
        sLattice.getBlock(iC),
        cuboid.get(load.glob(iC)))
    );
  }
}

template<typename T, typename DESCRIPTOR>
bool SuperLatticeFfromAnalyticalF3D<T, DESCRIPTOR>::operator()(
  T output[], const int input[])
{
  T physR[3] = {};
  this->_sLattice.getCuboidGeometry().getPhysR(physR,input);
  return _f(output,physR);
}


template<typename T, typename DESCRIPTOR>
BlockLatticeFfromAnalyticalF3D<T, DESCRIPTOR>::BlockLatticeFfromAnalyticalF3D(
  AnalyticalF3D<T, T>&                    f,
  BlockLattice<T, DESCRIPTOR>& lattice,
  Cuboid3D<T>&                            cuboid)
  : BlockLatticeF3D<T, DESCRIPTOR>(lattice, f.getTargetDim()),
    _f(f),
    _cuboid(cuboid)
{
  this->getName() = "blockFfromAnalyticalF(" + _f.getName() + ")";
}

template<typename T, typename DESCRIPTOR>
bool BlockLatticeFfromAnalyticalF3D<T, DESCRIPTOR>::operator()(
  T output[], const int input[])
{
  T physR[3] = {};
  _cuboid.getPhysR(physR,input);
  return _f(output,physR);
}


//////////// not yet working // symbolically ///////////////////
////////////////////////////////////////////////
template<typename T, typename DESCRIPTOR>
SmoothBlockIndicator3D<T, DESCRIPTOR>::SmoothBlockIndicator3D(
  IndicatorF3D<T>& f, T h, T eps, T sigma)
  : BlockDataF3D<T, T>((int)((f.getMax()[0] - f.getMin()[0]) / h + ( util::round(eps*0.5)*2+2 )),
                       (int)((f.getMax()[1] - f.getMin()[1]) / h + ( util::round(eps*0.5)*2+2 )),
                       (int)((f.getMax()[2] - f.getMin()[2]) / h + ( util::round(eps*0.5)*2+2 ))),
    _h(h),
    _sigma(sigma),
    _eps(util::round(eps*0.5)*2),
    _wa(_eps+1),
    _f(f)
{
  this->getName() = "SmoothBlockIndicator3D";
  T value, dx, dy, dz;
  T weights[this->_wa][this->_wa][this->_wa];
  T sum = 0;
  const int iStart = util::floor(this->_wa*0.5);
  const int iEnd = util::ceil(this->_wa*0.5);

  // calculate weights: they are constants, but calculation here is less error-prone than hardcoding these parameters
  for (int x = -iStart; x < iEnd; ++x) {
    for (int y = -iStart; y < iEnd; ++y) {
      for (int z = -iStart; z < iEnd; ++z) {
        weights[x+iStart][y+iStart][z+iStart] = util::exp(-(x*x+y*y+z*z)/(2*this->_sigma*this->_sigma)) / (util::pow(this->_sigma,3)*util::sqrt(util::pow(2,3)*util::pow(M_PI,3)));
        // important because sum of all weigths only equals 1 for this->_wa -> infinity
        sum += weights[x+iStart][y+iStart][z+iStart];
      }
    }
  }
  const T invSum = 1./sum;

  for (int iX=0; iX<this->getBlockData().getNx(); ++iX) {
    for (int iY=0; iY<this->getBlockData().getNy(); ++iY) {
      for (int iZ=0; iZ<this->getBlockData().getNz(); ++iZ) {
        bool output[1];
        value = 0;

        // input: regarded point (center)
        T input[] = {
          _f.getMin()[0] + iX*_h,
          _f.getMin()[1] + iY*_h,
          _f.getMin()[2] + iZ*_h
        };

        /*
         * three loops to look at every point (which weight is not 0) around the regarded point
         * sum all weighted porosities
         */
        for (int x = -iStart; x < iEnd; ++x) {
          for (int y = -iStart; y < iEnd; ++y) {
            for (int z = -iStart; z < iEnd; ++z) {
              dx = x*_h;
              dy = y*_h;
              dz = z*_h;

              // move from regarded point to point in neighborhood
              input[0] += dx;
              input[1] += dy;
              input[2] += dz;

              // get porosity
              _f(output,input);

              // sum porosity
              value += output[0] * weights[x+iStart][y+iStart][z+iStart];

              // move back to center
              input[0] -= dx;
              input[1] -= dy;
              input[2] -= dz;
            }
          }
        }
        /*
         * Round to 3 decimals
         * See above sum != 1.0, that's the reason for devision, otherwise porosity will never reach 0
         */
        this->getBlockData().get({iX,iY,iZ},0) = value*invSum;//nearbyint(1000*value/sum)/1000.0;
      }
    }
  }
}
/*
template<typename T, typename DESCRIPTOR>
bool SmoothBlockIndicator3D<T, DESCRIPTOR>::operator()(
  T output[], const int input[])
{
  T physR[3] = {};
  _superGeometry.getPhysR(physR,{input[0],input[1],input[2]} );
  _f(output,physR);
  return true;
}*/


template<typename T, typename DESCRIPTOR>
SuperLatticeInterpPhysVelocity3Degree3D<T, DESCRIPTOR>::
SuperLatticeInterpPhysVelocity3Degree3D(
  SuperLattice<T, DESCRIPTOR>& sLattice, UnitConverter<T,DESCRIPTOR>& conv, int range)
  : SuperLatticeF3D<T, DESCRIPTOR>(sLattice, 3)
{
  this->getName() = "Interp3DegreeVelocity";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    BlockLatticeInterpPhysVelocity3Degree3D<T, DESCRIPTOR>* foo =
      new BlockLatticeInterpPhysVelocity3Degree3D<T, DESCRIPTOR>(
      sLattice.getBlock(iC),
      conv,
      &sLattice.getCuboidGeometry().get(this->_sLattice.getLoadBalancer().
                                        glob(iC)),
      range);
    _bLattices.push_back(foo);
  }
}

template<typename T, typename DESCRIPTOR>
void SuperLatticeInterpPhysVelocity3Degree3D<T, DESCRIPTOR>::operator()(
  T output[], const T input[], const int iC)
{
  _bLattices[this->_sLattice.getLoadBalancer().loc(iC)]->operator()(output,
      input);
}

template<typename T, typename DESCRIPTOR>
BlockLatticeInterpPhysVelocity3Degree3D<T, DESCRIPTOR>::
BlockLatticeInterpPhysVelocity3Degree3D(
  BlockLattice<T, DESCRIPTOR>& blockLattice, UnitConverter<T,DESCRIPTOR>& conv,
  Cuboid3D<T>* c, int range)
  : BlockLatticeF3D<T, DESCRIPTOR>(blockLattice, 3),
    _conv(conv),
    _cuboid(c),
    _range(range)
{
  this->getName() = "BlockLatticeInterpVelocity3Degree3D";
}

template<typename T, typename DESCRIPTOR>
BlockLatticeInterpPhysVelocity3Degree3D<T, DESCRIPTOR>::
BlockLatticeInterpPhysVelocity3Degree3D(
  const BlockLatticeInterpPhysVelocity3Degree3D<T, DESCRIPTOR>& rhs) :
  BlockLatticeF3D<T, DESCRIPTOR>(rhs._blockLattice, 3),
  _conv(rhs._conv),
  _cuboid(rhs._cuboid),
  _range(rhs._range)
{
}

template<typename T, typename DESCRIPTOR>
void BlockLatticeInterpPhysVelocity3Degree3D<T, DESCRIPTOR>::operator()(
  T output[3], const T input[3])
{
  T u[3], rho, volume;
  int latIntPos[3] = {0};
  T latPhysPos[3] = {T()};
  _cuboid->getFloorLatticeR(latIntPos, &input[0]);
  _cuboid->getPhysR(latPhysPos, latIntPos);

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


template<typename T, typename DESCRIPTOR>
SuperLatticeInterpDensity3Degree3D<T, DESCRIPTOR>::SuperLatticeInterpDensity3Degree3D(
  SuperLattice<T, DESCRIPTOR>& sLattice, SuperGeometry<T,3>& sGeometry,
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
      sLattice.getBlock(lociC),
      sGeometry.getBlockGeometry(lociC),
      conv,
      &sLattice.getCuboidGeometry().get(globiC),
      range);
    _bLattices.push_back(foo);

    if (sLattice.getOverlap() <= range + 1)
      std::cout << "lattice overlap has to be larger than (range + 1)"
                << std::endl;
  }
}

template<typename T, typename DESCRIPTOR>
SuperLatticeInterpDensity3Degree3D<T, DESCRIPTOR>::~SuperLatticeInterpDensity3Degree3D()
{
  // first deconstruct vector elements
  for ( auto it : _bLattices) {
    delete it;
  }
  // then delete std::vector
  _bLattices.clear();
}

template<typename T, typename DESCRIPTOR>
void SuperLatticeInterpDensity3Degree3D<T, DESCRIPTOR>::operator()(T output[],
    const T input[], const int globiC)
{
  if (this->_sLattice.getLoadBalancer().rank(globiC) == singleton::mpi().getRank()) {
    _bLattices[this->_sLattice.getLoadBalancer().loc(globiC)]->operator()(output,
        input);
  }
}

template<typename T, typename DESCRIPTOR>
BlockLatticeInterpDensity3Degree3D<T, DESCRIPTOR>::BlockLatticeInterpDensity3Degree3D(
  BlockLattice<T, DESCRIPTOR>& blockLattice,
  BlockGeometry<T,3>& blockGeometry, UnitConverter<T,DESCRIPTOR>& conv,
  Cuboid3D<T>* c, int range) :
  BlockLatticeF3D<T, DESCRIPTOR>(blockLattice, 3), _blockGeometry(blockGeometry),
  _conv(conv), _cuboid(c), _range(range)
{
  this->getName() = "BlockLatticeInterpDensity3Degree3D";
}

template<typename T, typename DESCRIPTOR>
BlockLatticeInterpDensity3Degree3D<T, DESCRIPTOR>::BlockLatticeInterpDensity3Degree3D(
  const BlockLatticeInterpDensity3Degree3D<T, DESCRIPTOR>& rhs) :
  BlockLatticeF3D<T, DESCRIPTOR>(rhs._blockLattice, 3),
  _blockGeometry(rhs._blockGeometry),_conv(rhs._conv), _cuboid(
    rhs._cuboid), _range(rhs._range)
{
}

template<typename T, typename DESCRIPTOR>
void BlockLatticeInterpDensity3Degree3D<T, DESCRIPTOR>::operator()(
  T output[DESCRIPTOR::q], const T input[3])
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

  for (unsigned iPop = 0; iPop < DESCRIPTOR::q; ++iPop) {
    output[iPop] = T(0);
    for (int i = -_range; i <= _range + 1; ++i) {
      for (int j = -_range; j <= _range + 1; ++j) {
        for (int k = -_range; k <= _range + 1; ++k) {
          f_iPop = 0.;
          // just if material of cell != 1 there may be information of fluid density
          if (_blockGeometry.getMaterial({latIntPos[0] + i, latIntPos[1] + j,
            latIntPos[2] + k}) != 0) {
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

template<typename T, typename DESCRIPTOR>
SuperLatticeSmoothDiracDelta3D<T, DESCRIPTOR>::SuperLatticeSmoothDiracDelta3D(
  SuperLattice<T, DESCRIPTOR>& sLattice,
  UnitConverter<T,DESCRIPTOR>& conv, SuperGeometry<T,3>& sGeometry) :
  SuperLatticeF3D<T, DESCRIPTOR>(sLattice, 3)
{
  this->getName() = "SuperLatticeSmoothDiracDelta3D";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int lociC = 0; lociC < maxC; lociC++) {
    int globiC = this->_sLattice.getLoadBalancer().glob(lociC);

    BlockLatticeSmoothDiracDelta3D<T, DESCRIPTOR>* foo =
      new BlockLatticeSmoothDiracDelta3D<T, DESCRIPTOR>(
      sLattice.getBlock(lociC),
      conv, &sLattice.getCuboidGeometry().get(globiC)
    );
    _bLattices.push_back(foo);
  }
}

template<typename T, typename DESCRIPTOR>
SuperLatticeSmoothDiracDelta3D<T, DESCRIPTOR>::~SuperLatticeSmoothDiracDelta3D()
{
  for ( auto it : _bLattices) {
    delete it;
  }
  _bLattices.clear();
}

template<typename T, typename DESCRIPTOR>
void SuperLatticeSmoothDiracDelta3D<T, DESCRIPTOR>::operator()(T delta[4][4][4],
    const T physPos[3], const int globiC)
{
  if (this->_sLattice.getLoadBalancer().rank(globiC) == singleton::mpi().getRank()) {
    _bLattices[this->_sLattice.getLoadBalancer().loc(globiC)]->operator()(delta,
        physPos);
  }
}

template<typename T, typename DESCRIPTOR>
BlockLatticeSmoothDiracDelta3D<T, DESCRIPTOR>::BlockLatticeSmoothDiracDelta3D(
  BlockLattice<T, DESCRIPTOR>& blockLattice, UnitConverter<T,DESCRIPTOR>& conv, Cuboid3D<T>* cuboid)
  : BlockLatticeF3D<T, DESCRIPTOR>(blockLattice, 3), _conv(conv), _cuboid(cuboid)
{
  this->getName() = "BlockLatticeSmoothDiracDelta3D";
}

template<typename T, typename DESCRIPTOR>
BlockLatticeSmoothDiracDelta3D<T, DESCRIPTOR>::BlockLatticeSmoothDiracDelta3D(
  const BlockLatticeSmoothDiracDelta3D<T, DESCRIPTOR>& rhs)
  :
  BlockLatticeF3D<T, DESCRIPTOR>(rhs._blockLattice, 3), _conv(rhs._conv), _cuboid(rhs._cuboid)
{
}

template<typename T, typename DESCRIPTOR>
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
        delta[i+range][j+range][k+range] *= 1. / 4 * (1 + util::cos(M_PI * a / 2.));
        delta[i+range][j+range][k+range] *= 1. / 4 * (1 + util::cos(M_PI * b / 2.));
        delta[i+range][j+range][k+range] *= 1. / 4 * (1 + util::cos(M_PI * c / 2.));

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
