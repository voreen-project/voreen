/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012 Lukas Baron, Tim Dornieden, Mathias J. Krause, Albert Mink
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

#ifndef ANALYTICAL_CALC_F_HH
#define ANALYTICAL_CALC_F_HH


#include "analyticCalcF.h"
#include "core/olbDebug.h"

namespace olb {



//////////////////////////////// AnalyticCalc1D ////////////////////////////////
template <typename T, typename S>
AnalyticCalc1D<T,S>::AnalyticCalc1D(AnalyticalF1D<T,S>& f, AnalyticalF1D<T,S>& g)
  : AnalyticalF1D<T,S>(f.getTargetDim()), _f(f), _g(g)
{
  std::swap(f._ptrCalcC, this->_ptrCalcC);
}

template <typename T, typename S>
AnalyticPlus1D<T,S>::AnalyticPlus1D(AnalyticalF1D<T,S>& f, AnalyticalF1D<T,S>& g)
  : AnalyticCalc1D<T,S>(f,g)
{
  OLB_ASSERT(g.getTargetDim() == f.getTargetDim(), "the dimensions of both functors need to be equal");
  this->getName() = "(" + f.getName() + "+" + g.getName() + ")";

}

template <typename T, typename S>
bool AnalyticPlus1D<T,S>::operator()(T output[], const S input[])
{
  T outputTmp[this->_g.getTargetDim()];
  this->_g(outputTmp,input);
  this->_f(output,input);
  for (int i = 0; i < this->_f.getTargetDim(); i++) {
    output[i]+=outputTmp[i];
  }
  return true;
}

template <typename T, typename S>
AnalyticMinus1D<T,S>::AnalyticMinus1D(AnalyticalF1D<T,S>& f, AnalyticalF1D<T,S>& g)
  : AnalyticCalc1D<T,S>(f,g)
{
  OLB_ASSERT(g.getTargetDim() == f.getTargetDim(), "the dimensions of both functors need to be equal");
  this->getName() = "(" + f.getName() + "-" + g.getName() + ")";
}

template <typename T, typename S>
bool AnalyticMinus1D<T,S>::operator()(T output[], const S input[])
{
  T outputTmp[this->_g.getTargetDim()];
  this->_g(outputTmp,input);
  this->_f(output,input);
  for (int i = 0; i < this->_f.getTargetDim(); i++) {
    output[i]-=outputTmp[i];
  }
  return true;
}


template <typename T, typename S>
AnalyticMultiplication1D<T,S>::AnalyticMultiplication1D(AnalyticalF1D<T,S>& f,
    AnalyticalF1D<T,S>& g) : AnalyticCalc1D<T,S>(f,g)
{
  OLB_ASSERT(g.getTargetDim() == f.getTargetDim(), "the dimensions of both functors need to be equal");
  this->getName() = "(" + f.getName() + "*" + g.getName() + ")";
}

template <typename T, typename S>
bool AnalyticMultiplication1D<T,S>::operator()(T output[], const S input[])
{
  T outputTmp[this->_g.getTargetDim()];
  this->_g(outputTmp,input);
  this->_f(output,input);
  for (int i = 0; i < this->_f.getTargetDim(); i++) {
    output[i]*=outputTmp[i];
  }
  return true;
}


template <typename T, typename S>
AnalyticDivision1D<T,S>::AnalyticDivision1D(AnalyticalF1D<T,S>& f,
    AnalyticalF1D<T,S>& g) : AnalyticCalc1D<T,S>(f,g)
{
  OLB_ASSERT(g.getTargetDim() == f.getTargetDim(), "the dimensions of both functors need to be equal");
  this->getName() = "(" + f.getName() + "/" + g.getName() + ")";
}

template <typename T, typename S>
bool AnalyticDivision1D<T,S>::operator()(T output[], const S input[])
{
  T outputTmp[this->_g.getTargetDim()];
  this->_g(outputTmp,input);
  this->_f(output,input);
  for (int i = 0; i < this->_f.getTargetDim(); i++) {
    output[i]/=outputTmp[i];
  }
  return true;
}


/////////////////////////////////operator()/// ////////////////////////////////
template <typename T, typename S>
AnalyticalF1D<T,S>& AnalyticalF1D<T,S>::operator+(AnalyticalF1D<T,S>& rhs)
{
  auto tmp = std::make_shared< AnalyticPlus1D<T,S> >(*this,rhs);
  this->_ptrCalcC = tmp;
  return *tmp;
}

template <typename T, typename S>
AnalyticalF1D<T,S>& AnalyticalF1D<T,S>::operator-(AnalyticalF1D<T,S>& rhs)
{
  auto tmp = std::make_shared< AnalyticMinus1D<T,S> >(*this,rhs);
  this->_ptrCalcC = tmp;
  return *tmp;
}

template <typename T, typename S>
AnalyticalF1D<T,S>& AnalyticalF1D<T,S>::operator*(AnalyticalF1D<T,S>& rhs)
{
  auto tmp = std::make_shared< AnalyticMultiplication1D<T,S> >(*this,rhs);
  this->_ptrCalcC = tmp;
  return *tmp;
}

template <typename T, typename S>
AnalyticalF1D<T,S>& AnalyticalF1D<T,S>::operator/(AnalyticalF1D<T,S>& rhs)
{
  auto tmp = std::make_shared< AnalyticDivision1D<T,S> >(*this,rhs);
  this->_ptrCalcC = tmp;
  return *tmp;
}




//////////////////////////////// AnalyticCalc2D ////////////////////////////////
template <typename T, typename S>
AnalyticCalc2D<T,S>::AnalyticCalc2D(AnalyticalF2D<T,S>& f, AnalyticalF2D<T,S>& g)
  : AnalyticalF2D<T,S>(f.getTargetDim()), _f(f), _g(g)
{
  // pass through the shared_ptr from the first argument f to the arithmetic class itself.
  // used by secsessive calls: e.g. (functorA + functor B) followed by (functorA + functorC)
  // the result of the first operation is overwritten by the second.

  // equivalent operations
  //  std::swap(f._ptrCalcC, this->_ptrCalcC);
  //  this->_ptrCalcC = f._ptrCalcC;
  this->_ptrCalcC.swap(f._ptrCalcC);
}


template <typename T, typename S>
AnalyticPlus2D<T,S>::AnalyticPlus2D(AnalyticalF2D<T,S>& f, AnalyticalF2D<T,S>& g)
  : AnalyticCalc2D<T,S>(f,g)
{
  OLB_ASSERT(g.getTargetDim() == f.getTargetDim(), "the dimensions of both functors need to be equal");
  this->getName() = "(" + f.getName() + "+" + g.getName() + ")";
}

template <typename T, typename S>
bool AnalyticPlus2D<T,S>::operator()(T output[], const S input[])
{
  T outputTmp[this->_g.getTargetDim()];
  this->_g(outputTmp,input);
  this->_f(output,input);
  for (int i = 0; i < this->_f.getTargetDim(); i++) {
    output[i]+=outputTmp[i];
  }
  return true;
}


template <typename T, typename S>
AnalyticMinus2D<T,S>::AnalyticMinus2D(AnalyticalF2D<T,S>& f, AnalyticalF2D<T,S>& g)
  : AnalyticCalc2D<T,S>(f,g)
{
  OLB_ASSERT(g.getTargetDim() == f.getTargetDim(), "the dimensions of both functors need to be equal");
  this->getName() = "(" + f.getName() + "-" + g.getName() + ")";
}

template <typename T, typename S>
bool AnalyticMinus2D<T,S>::operator()(T output[], const S input[])
{
  T outputTmp[this->_g.getTargetDim()];
  this->_g(outputTmp,input);
  this->_f(output,input);
  for (int i = 0; i < this->_f.getTargetDim(); i++) {
    output[i]-=outputTmp[i];
  }
  return true;
}


template <typename T, typename S>
AnalyticMultiplication2D<T,S>::AnalyticMultiplication2D(AnalyticalF2D<T,S>& f,
    AnalyticalF2D<T,S>& g) : AnalyticCalc2D<T,S>(f,g)
{
  OLB_ASSERT(g.getTargetDim() == f.getTargetDim(), "the dimensions of both functors need to be equal");
  this->getName() = "(" + f.getName() + "*" + g.getName() + ")";
}

template <typename T, typename S>
bool AnalyticMultiplication2D<T,S>::operator()(T output[], const S input[])
{
  T outputTmp[this->_g.getTargetDim()];
  this->_g(outputTmp,input);
  this->_f(output,input);
  for (int i = 0; i < this->_f.getTargetDim(); i++) {
    output[i]*=outputTmp[i];
  }
  return true;
}


template <typename T, typename S>
AnalyticDivision2D<T,S>::AnalyticDivision2D(AnalyticalF2D<T,S>& f,
    AnalyticalF2D<T,S>& g) : AnalyticCalc2D<T,S>(f,g)
{
  OLB_ASSERT(g.getTargetDim() == f.getTargetDim(), "the dimensions of both functors need to be equal");
  this->getName() = "(" + f.getName() + "/" + g.getName() + ")";
}

template <typename T, typename S>
bool AnalyticDivision2D<T,S>::operator()(T output[], const S input[])
{
  T outputTmp[this->_g.getTargetDim()];
  this->_g(outputTmp,input);
  this->_f(output,input);
  for (int i = 0; i < this->_f.getTargetDim(); i++) {
    output[i]/=outputTmp[i];
  }
  return true;
}


/////////////////////////////////operator()////////////////////////////////////
template <typename T, typename S>
AnalyticalF2D<T,S>& AnalyticalF2D<T,S>::operator+(AnalyticalF2D<T,S>& rhs)
{
  // version 1
  //    AnalyticalF2D<T,S>* tmp = new AnalyticPlus2D<T,S>(*this,rhs);
  //    std::shared_ptr< GenericF<T,S> > ptr( tmp );
  //    this->_ptrCalcC = ptr;

  // version 2
  //  std::shared_ptr< AnalyticalF2D<T,S> > tmp = std::make_shared< AnalyticPlus2D<T,S> >(*this,rhs);

  // version 2.5
  //  std::shared_ptr< AnalyticPlus2D<T,S> > tmp( new AnalyticPlus2D<T,S>(*this,rhs) );

  // version 3
  auto tmp = std::make_shared< AnalyticPlus2D<T,S> >(*this,rhs);

  this->_ptrCalcC = tmp;

  //  std::cout << "operator+(): " << this->_ptrCalcC.get()->getName() << std::endl;
  return *tmp;
}

template <typename T, typename S>
AnalyticalF2D<T,S>& AnalyticalF2D<T,S>::operator-(AnalyticalF2D<T,S>& rhs)
{
  auto tmp = std::make_shared< AnalyticMinus2D<T,S> >(*this,rhs);
  this->_ptrCalcC = tmp;
  return *tmp;
}

template <typename T, typename S>
AnalyticalF2D<T,S>& AnalyticalF2D<T,S>::operator*(AnalyticalF2D<T,S>& rhs)
{
  auto tmp = std::make_shared< AnalyticMultiplication2D<T,S> >(*this,rhs);
  this->_ptrCalcC = tmp;
  return *tmp;
}

template <typename T, typename S>
AnalyticalF2D<T,S>& AnalyticalF2D<T,S>::operator/(AnalyticalF2D<T,S>& rhs)
{
  auto tmp = std::make_shared< AnalyticDivision2D<T,S> >(*this,rhs);
  this->_ptrCalcC = tmp;
  return *tmp;
}




//////////////////////////////// AnalyticCalc3D ////////////////////////////////
template <typename T, typename S>
AnalyticCalc3D<T,S>::AnalyticCalc3D(AnalyticalF3D<T,S>& f, AnalyticalF3D<T,S>& g)
  : AnalyticalF3D<T,S>(f.getTargetDim()), _f(f), _g(g)
{}


template <typename T, typename S>
AnalyticPlus3D<T,S>::AnalyticPlus3D(AnalyticalF3D<T,S>& f, AnalyticalF3D<T,S>& g)
  : AnalyticCalc3D<T,S>(f,g)
{
  OLB_ASSERT(g.getTargetDim() == f.getTargetDim(), "the dimensions of both functors need to be equal");
  this->getName() = "(" + f.getName() + "+" + g.getName() + ")";
  std::swap(f._ptrCalcC, this->_ptrCalcC);
}

template <typename T, typename S>
bool AnalyticPlus3D<T,S>::operator()(T output[], const S input[])
{
  T outputTmp[this->_g.getTargetDim()];
  this->_g(outputTmp,input);
  this->_f(output,input);
  for (int i = 0; i < this->_f.getTargetDim(); i++) {
    output[i]+=outputTmp[i];
  }
  return true;
}


template <typename T, typename S>
AnalyticMinus3D<T,S>::AnalyticMinus3D(AnalyticalF3D<T,S>& f, AnalyticalF3D<T,S>& g)
  : AnalyticCalc3D<T,S>(f,g)
{
  OLB_ASSERT(g.getTargetDim() == f.getTargetDim(), "the dimensions of both functors need to be equal");
  this->getName() = "(" + f.getName() + "-" + g.getName() + ")";
  std::swap(f._ptrCalcC, this->_ptrCalcC);
}

template <typename T, typename S>
bool AnalyticMinus3D<T,S>::operator()(T output[], const S input[])
{
  T outputTmp[this->_g.getTargetDim()];
  this->_g(outputTmp,input);
  this->_f(output,input);
  for (int i = 0; i < this->_f.getTargetDim(); i++) {
    output[i]-=outputTmp[i];
  }
  return true;
}


template <typename T, typename S>
AnalyticMultiplication3D<T,S>::AnalyticMultiplication3D(AnalyticalF3D<T,S>& f,
    AnalyticalF3D<T,S>& g) : AnalyticCalc3D<T,S>(f,g)
{
  OLB_ASSERT(g.getTargetDim() == f.getTargetDim(), "the dimensions of both functors need to be equal");
  this->getName() = "(" + f.getName() + "*" + g.getName() + ")";
  std::swap(f._ptrCalcC, this->_ptrCalcC);
}

template <typename T, typename S>
bool AnalyticMultiplication3D<T,S>::operator()(T output[], const S input[])
{
  T outputTmp[this->_g.getTargetDim()];
  this->_g(outputTmp,input);
  this->_f(output,input);
  for (int i = 0; i < this->_f.getTargetDim(); i++) {
    output[i]*=outputTmp[i];
  }
  return true;
}


template <typename T, typename S>
AnalyticDivision3D<T,S>::AnalyticDivision3D(AnalyticalF3D<T,S>& f,
    AnalyticalF3D<T,S>& g) : AnalyticCalc3D<T,S>(f,g)
{
  OLB_ASSERT(g.getTargetDim() == f.getTargetDim(), "the dimensions of both functors need to be equal");
  this->getName() = "(" + f.getName() + "/" + g.getName() + ")";
  std::swap(f._ptrCalcC, this->_ptrCalcC);
}

template <typename T, typename S>
bool AnalyticDivision3D<T,S>::operator()(T output[], const S input[])
{
  T outputTmp[this->_g.getTargetDim()];
  this->_g(outputTmp,input);
  this->_f(output,input);
  for (int i = 0; i < this->_f.getTargetDim(); i++) {
    output[i]/=outputTmp[i];
  }
  return true;
}

/////////////////////////////////operator()/// ////////////////////////////////
template <typename T, typename S>
AnalyticalF3D<T,S>& AnalyticalF3D<T,S>::operator+(AnalyticalF3D<T,S>& rhs)
{
  auto tmp = std::make_shared< AnalyticPlus3D<T,S> >(*this,rhs);
  this->_ptrCalcC = tmp;
  return *tmp;
}

template <typename T, typename S>
AnalyticalF3D<T,S>& AnalyticalF3D<T,S>::operator-(AnalyticalF3D<T,S>& rhs)
{
  auto tmp = std::make_shared< AnalyticMinus3D<T,S> >(*this,rhs);
  this->_ptrCalcC = tmp;
  return *tmp;
}

template <typename T, typename S>
AnalyticalF3D<T,S>& AnalyticalF3D<T,S>::operator*(AnalyticalF3D<T,S>& rhs)
{
  auto tmp = std::make_shared< AnalyticMultiplication3D<T,S> >(*this,rhs);
  this->_ptrCalcC = tmp;
  return *tmp;
}

template <typename T, typename S>
AnalyticalF3D<T,S>& AnalyticalF3D<T,S>::operator/(AnalyticalF3D<T,S>& rhs)
{
  auto tmp = std::make_shared< AnalyticDivision3D<T,S> >(*this,rhs);
  this->_ptrCalcC = tmp;
  return *tmp;
}

} // end namespace olb

#endif
