/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012-2017 Lukas Baron, Tim Dornieden, Mathias J. Krause,
 *                          Albert Mink, Benjamin Förster, Adrian Kummerländer
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

#ifndef SUPER_BASE_F_3D_HH
#define SUPER_BASE_F_3D_HH

#include "superBaseF3D.h"

namespace olb {

template <typename T, typename W>
SuperF3D<T,W>::SuperF3D(SuperStructure3D<T>& superStructure, int targetDim)
  : GenericF<W,int>(targetDim,4), _superStructure(superStructure) { }

template <typename T, typename W>
SuperStructure3D<T>& SuperF3D<T,W>::getSuperStructure()
{
  return _superStructure;
}

template <typename T, typename W>
int SuperF3D<T,W>::getBlockFSize() const
{
  OLB_ASSERT(_blockF.size() < INT32_MAX,
             "it is safe to cast std::size_t to int");
  return _blockF.size();
}

template <typename T, typename W>
BlockF3D<W>& SuperF3D<T,W>::getBlockF(int iCloc)
{
  OLB_ASSERT(iCloc < int(_blockF.size()) && iCloc >= 0,
             "block functor index within bounds");
  return *(_blockF[iCloc]);
}


template <typename T,typename BaseType>
SuperDataF3D<T,BaseType>::SuperDataF3D(SuperData3D<T,BaseType>& superData)
  : SuperF3D<T,BaseType>( superData, superData.getDataSize() ), _superData(superData)
{
}

template <typename T,typename BaseType>
bool SuperDataF3D<T,BaseType>::operator() (BaseType output[], const int input[])
{
  int iC_loc = _superData.getLoadBalancer().loc(input[0]);
  for (int i=0; i < _superData.getDataSize(); i++) {
    output[i] = _superData.get(iC_loc).get(input[1] + _superData.getOverlap(),
                                           input[2] + _superData.getOverlap(),
                                           input[3] + _superData.getOverlap(), i);
  }
  return true;
}

template <typename T,typename BaseType>
SuperData3D<T,BaseType>& SuperDataF3D<T,BaseType>::getSuperData()
{
  return _superData;
}



template <typename T, typename W>
SuperIdentity3D<T,W>::SuperIdentity3D(SuperF3D<T,W>& f)
  : SuperF3D<T,W>(f.getSuperStructure(), f.getTargetDim() ), _f(f)
{
  this->getName() = _f.getName();
  std::swap( _f._ptrCalcC, this->_ptrCalcC );
}

template <typename T, typename W>
bool SuperIdentity3D<T,W>::operator()(W output[], const int input[])
{
  _f(output, input);
  return true;
}



template <typename T, typename W>
SuperExtractComponentF3D<T,W>::SuperExtractComponentF3D(SuperF3D<T,W>& f, int extractDim)
  : SuperF3D<T,W>(f.getSuperStructure(), 1 ), _f(f), _extractDim(extractDim)
{
  this->getName() = _f.getName();
}

template <typename T, typename W>
int SuperExtractComponentF3D<T,W>::getExtractDim()
{
  return _extractDim;
}

template <typename T, typename W>
bool SuperExtractComponentF3D<T,W>::operator()(W output[], const int input[])
{
  T outTmp[3];
  _f(outTmp, input);
  output[0] = outTmp[_extractDim];
  return true;
}



template <typename T, typename W>
SuperExtractComponentIndicatorF3D<T,W>::SuperExtractComponentIndicatorF3D(SuperF3D<T,W>& f, int extractDim, SuperIndicatorF3D<T>& indicatorF)
  : SuperExtractComponentF3D<T,W>(f , extractDim), _f(f), _extractDim(extractDim), _indicatorF(indicatorF)
{
  this->getName() = _f.getName();
}

template <typename T, typename W>
bool SuperExtractComponentIndicatorF3D<T,W>::operator()(W output[], const int input[])
{
  T outTmp[3];
  bool out;
  _indicatorF(&out, input);
  if (out) {
    _f(outTmp, input);
    output[0] = outTmp[_extractDim];
  }
  return true;
}


template <typename T, typename W>
SuperExtractIndicatorF3D<T,W>::SuperExtractIndicatorF3D(SuperF3D<T,W>& f, SuperIndicatorF3D<T>& indicatorF)
  : SuperF3D<T,W>(f.getSuperStructure(), f.getTargetDim()), _f(f), _indicatorF(indicatorF)
{
  this->getName() = _f.getName();
}

template <typename T, typename W>
bool SuperExtractIndicatorF3D<T,W>::operator()(W output[], const int input[])
{
  bool out;
  _indicatorF(&out, input);
  if (out) {
    _f(output, input);
  }
  return true;
}



template <typename T, typename W>
SuperDotProductF3D<T,W>::SuperDotProductF3D(SuperF3D<T,W>& f, T vector[])
  : SuperF3D<T,W>(f.getSuperStructure() ,1 ), _f(f), _vector(vector)
{
  this->getName() = _f.getName();
  /*if ( (sizeof(_vector)/sizeof(T)) != _f.getTargetDim() ) {
    std::cout << "WARNING: dimension of vectors do not match!" << std::endl;
    exit(-1);
  }*/
}

template <typename T, typename W>
bool SuperDotProductF3D<T,W>::operator()(W output[], const int input[])
{
  T outTmp[3];
  _f(outTmp, input);
  output[0] = T();
  for (int iDim=0; iDim<_f.getTargetDim(); iDim++) {
    output[0] += outTmp[iDim]*_vector[iDim];
  }
  return true;
}



template <typename T, typename W>
SuperIdentityOnSuperIndicatorF3D<T,W>::SuperIdentityOnSuperIndicatorF3D(SuperF3D<T,W>& f,
    SuperIndicatorF3D<T>& indicatorF,
    W defaultValue)
  : SuperF3D<T,W>(f.getSuperStructure() ,f.getTargetDim() ),
    _f(f),
    _indicatorF(indicatorF),
    _defaultValue(defaultValue)
{
  this->getName() = _f.getName() + std::string("_on_") + _indicatorF.getName();
  std::swap( _f._ptrCalcC, this->_ptrCalcC );
}

template <typename T, typename W>
bool SuperIdentityOnSuperIndicatorF3D<T,W>::operator()(W output[], const int input[])
{
  bool indic;
  _indicatorF(&indic, input);
  if (indic) {
    _f(output, input);
  } else {
    for (int i=0; i<_f.getTargetDim(); i++) {
      output[i] = _defaultValue;
    }
  }
  return true;
}


template <typename T, template <typename U> class Lattice>
SuperLatticeF3D<T,Lattice>::SuperLatticeF3D(SuperLattice3D<T,Lattice>& superLattice,
    int targetDim)
  : SuperF3D<T,T>(superLattice, targetDim), _sLattice(superLattice) { }

template <typename T, template <typename U> class Lattice>
SuperLattice3D<T,Lattice>& SuperLatticeF3D<T,Lattice>::getSuperLattice()
{
  return _sLattice;
}

template <typename T, template <typename U> class DESCRIPTOR>
SuperLatticePhysF3D<T,DESCRIPTOR>::SuperLatticePhysF3D
(SuperLattice3D<T,DESCRIPTOR>& sLattice, const UnitConverter<T,DESCRIPTOR>& converter,
 int targetDim)
  : SuperLatticeF3D<T,DESCRIPTOR>(sLattice, targetDim), _converter(converter) { }

template <typename T, template <typename U> class DESCRIPTOR>
UnitConverter<T,DESCRIPTOR> const& SuperLatticePhysF3D<T,DESCRIPTOR>::getConverter() const
{
  return this->_converter;
}
template <typename T, template <typename U> class DESCRIPTOR, template <typename V> class ThermalDESCRIPTOR>
SuperLatticeThermalPhysF3D<T,DESCRIPTOR,ThermalDESCRIPTOR>::SuperLatticeThermalPhysF3D
(SuperLattice3D<T,ThermalDESCRIPTOR>& sLattice, const ThermalUnitConverter<T,DESCRIPTOR,ThermalDESCRIPTOR>& converter,
 int targetDim)
  : SuperLatticeF3D<T,ThermalDESCRIPTOR>(sLattice, targetDim), _converter(converter) { }

template <typename T, template <typename U> class DESCRIPTOR, template <typename V> class ThermalDESCRIPTOR>
ThermalUnitConverter<T,DESCRIPTOR,ThermalDESCRIPTOR> const& SuperLatticeThermalPhysF3D<T,DESCRIPTOR,ThermalDESCRIPTOR>::getConverter() const
{
  return this->_converter;
}

template <typename T, template <typename U> class Lattice>
ComposedSuperLatticeF3D<T,Lattice>::ComposedSuperLatticeF3D
(SuperLatticeF3D<T,Lattice>& f0, SuperLatticeF3D<T,Lattice>& f1,
 SuperLatticeF3D<T,Lattice>& f2)
  : SuperLatticeF3D<T,Lattice>(f0.getSuperLattice(), 3), _f0(f0), _f1(f1), _f2(f2)
{
  this->getName() = "composedSuperLatticeF3D";
}

template <typename T, template <typename U> class Lattice>
bool ComposedSuperLatticeF3D<T,Lattice>::operator() (T output[], const int input[])
{
  T tmp[3] = {};
  SuperIdentity3D<T,T> ff0(_f0), ff1(_f1), ff2(_f2);
  _f0(tmp,input);
  output[0]=tmp[0];
  _f1(tmp,input);
  output[1]=tmp[0];
  _f2(tmp,input);
  output[2]=tmp[0];
  return true;
}


} // end namespace olb

#endif
