/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2020 Davide Dapelo
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
 * Generic functions for the finite-difference differentiation schemes.
 *  -- header file
 */
#ifndef FD_SCHEME_BASE_H
#define FD_SCHEME_BASE_H

#include "fdTags.h"

namespace olb {

namespace fd {

/////////////////////////////////////////////////////////////////////////////////////////////////////
template <unsigned D, typename T>
struct AdvectionSchemeBase {
  AdvectionSchemeBase() {}
  virtual T operator()(T& f0, std::vector<std::vector<T>>& f, std::vector<std::vector<T>>& F, std::vector<T>& u) =0;
};

template <unsigned D, typename T>
struct DiffusionSchemeBase {
  DiffusionSchemeBase() {}
  virtual T operator()(T& f0, std::vector<std::vector<T>>& f, std::vector<std::vector<T>>& F, T& diff, std::vector<T> u=std::vector<T>(D,T())) =0;
};

template <unsigned D, typename T>
struct AdBoundarySchemeBase {
  AdBoundarySchemeBase() {}
  virtual int getExtraExtent()=0;
  virtual void operator()(std::vector<std::vector<T>>& fOut, T& f0, std::vector<std::vector<T>>& fIn, std::vector<int>& normal, std::vector<T>& u) =0;
};

/////////////////////////////////////////////////////////////////////////////////////////////////////
template <unsigned D, typename T, typename TAG>
struct AdvectionScheme final : AdvectionSchemeBase<D,T> {
  AdvectionScheme() {}
  virtual T operator()(T& f0, std::vector<std::vector<T>>& f, std::vector<std::vector<T>>& F, std::vector<T>& u) override
  {
    throw std::invalid_argument("Wrong advection scheme tag.");
  }
};

template <unsigned D, typename T, typename TAG>
struct DiffusionScheme final : DiffusionSchemeBase<D,T> {
  DiffusionScheme() {}
  virtual T operator()(T& f0, std::vector<std::vector<T>>& f, std::vector<std::vector<T>>& F, T& diff, std::vector<T> u=std::vector<T>(D,T())) override
  {
    throw std::invalid_argument("Wrong diffusion scheme tag.");
  }
};

template <unsigned D, typename T, typename TAG>
struct AdNeumannZeroBoundaryScheme final : AdBoundarySchemeBase<D,T> {
  AdNeumannZeroBoundaryScheme() {}
  virtual int getExtraExtent() override
  {
    return 0;
  }
  virtual void operator() (
    std::vector<std::vector<T>>& fOut, T& f0, std::vector<std::vector<T>>& fIn, std::vector<int>& normal, std::vector<T>& u ) override
  {
    throw std::invalid_argument("Wrong advection scheme tag.");
  }
};


}  // namespace fd

}  // namespace olb

#endif
