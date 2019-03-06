/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012-2017 Lukas Baron, Tim Dornieden, Mathias J. Krause,
 *  Albert Mink, Adrian Kummerländer
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

#ifndef SUPER_CALC_F_3D_HH
#define SUPER_CALC_F_3D_HH

#include "superCalcF3D.h"
#include "blockCalcF3D.h"
#include "core/olbDebug.h"

namespace olb {


template <typename T, typename W, template<typename> class F>
SuperCalc3D<T,W,F>::SuperCalc3D(SuperF3D<T,W>& f, SuperF3D<T,W>& g)
  : SuperF3D<T,W>(
      f.getSuperStructure(),
      f.getTargetDim() > g.getTargetDim() ? f.getTargetDim() : g.getTargetDim()),
  _f(f), _g(g)
{
  OLB_ASSERT(
    f.getTargetDim() == g.getTargetDim() || f.getTargetDim() == 1 || g.getTargetDim() == 1,
    "Componentwise operation must be well defined.");

  this->getName() = "(" + f.getName() + F<T>::symbol + g.getName() + ")";

  std::swap(f._ptrCalcC, this->_ptrCalcC);

  LoadBalancer<T>& load = f.getSuperStructure().getLoadBalancer();
  if ( f.getBlockFSize() == load.size() ) {
    if ( g.getBlockFSize() == load.size() ) {
      // both functors expose the correct count of block level functors
      for (int iC = 0; iC < load.size(); ++iC) {
        this->_blockF.emplace_back(
          new BlockCalc3D<W,F>(f.getBlockF(iC), g.getBlockF(iC))
        );
      }
    } else {
      // operate on super functor `g` and block level functors provided by `f`
      for (int iC = 0; iC < load.size(); ++iC) {
        this->_blockF.emplace_back(
          new BlockCalc3D<W,F>(f.getBlockF(iC), g, load.glob(iC))
        );
      }
    }
  } else if ( g.getBlockFSize() == load.size() ) {
    // operate on block level functors provided by `f` and super functor `g`
    for (int iC = 0; iC < load.size(); ++iC) {
      this->_blockF.emplace_back(
        new BlockCalc3D<W,F>(f, load.glob(iC), g.getBlockF(iC))
      );
    }
  }
}

template <typename T, typename W, template<typename> class F>
bool SuperCalc3D<T,W,F>::operator()(W output[], const int input[])
{
  if ( _f.getTargetDim() == 1 || _g.getTargetDim() == 1 ) {
    // scalar operation
    W scalar;
    if ( _f.getTargetDim() == 1 ) {
      // apply the scalar f to possibly multidimensional g
      _g(output, input);
      _f(&scalar, input);
    } else {
      // apply scalar g to possibly multidimensional f
      _f(output, input);
      _g(&scalar, input);
    }

    for (int i = 0; i < this->getTargetDim(); i++) {
      output[i] = F<T>()(output[i], scalar);
    }
  } else {
    // componentwise operation on equidimensional functors
    W* outputF = output;
    W outputG[this->getTargetDim()];

    _f(outputF, input);
    _g(outputG, input);

    for (int i = 0; i < this->getTargetDim(); i++) {
      output[i] = F<T>()(outputF[i], outputG[i]);
    }
  }
  return true;
}


/////////////////////////////////operator()///////////////////////////////////

template <typename T, typename W>
SuperF3D<T,W>& SuperF3D<T,W>::operator+(SuperF3D<T,W>& rhs)
{
  auto tmp = std::make_shared< SuperPlus3D<T,W> >(*this,rhs);
  this->_ptrCalcC = tmp;
  return *tmp;
}

template <typename T, typename W>
SuperF3D<T,W>& SuperF3D<T,W>::operator-(SuperF3D<T,W>& rhs)
{
  auto tmp = std::make_shared< SuperMinus3D<T,W> >(*this,rhs);
  this->_ptrCalcC = tmp;
  return *tmp;
}

template <typename T, typename W>
SuperF3D<T,W>& SuperF3D<T,W>::operator*(SuperF3D<T,W>& rhs)
{
  auto tmp = std::make_shared< SuperMultiplication3D<T,W> >(*this,rhs);
  this->_ptrCalcC = tmp;
  return *tmp;
}

template <typename T, typename W>
SuperF3D<T,W>& SuperF3D<T,W>::operator/(SuperF3D<T,W>& rhs)
{
  auto tmp = std::make_shared< SuperDivision3D<T,W> >(*this,rhs);
  this->_ptrCalcC = tmp;
  return *tmp;
}

} // end namespace olb

#endif
