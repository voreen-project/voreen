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
 * Specializations for the finite-difference differentiation schemes.
 *  -- header file
 */
#ifndef FD_SCHEMES_H
#define FD_SCHEMES_H

#include "fdTags.h"
#include "fdSchemeBase.h"

namespace olb {

namespace fd {

namespace tag {

struct CENTRAL final : FD_TAG {
  CENTRAL() = delete;
  static constexpr int extent = 1;
};

constexpr int CENTRAL::extent;

struct UPWIND final : FD_TAG {
  UPWIND() = delete;
  static constexpr int extent = 1;
};

constexpr int UPWIND::extent;

struct UPWIND_2_ORDER final : FD_TAG {
  UPWIND_2_ORDER() = delete;
  static constexpr int extent = 2;
};

constexpr int UPWIND_2_ORDER::extent;

} // namespace tag

//////////////////////////////////////////// CENTRAL /////////////////////////////////////////////////////////
template <unsigned D, typename T>
struct AdvectionScheme<D,T,tag::CENTRAL> final : AdvectionSchemeBase<D,T> {
  AdvectionScheme() {}
  virtual T operator()(T& f0, std::vector<std::vector<T>>& f, std::vector<std::vector<T>>& F, std::vector<T>& u) override;
};

template <unsigned D, typename T>
struct DiffusionScheme<D,T,tag::CENTRAL> final : DiffusionSchemeBase<D,T> {
  DiffusionScheme() {}
  virtual T operator()(T& f0, std::vector<std::vector<T>>& f, std::vector<std::vector<T>>& F, T& diff, std::vector<T> u=std::vector<T>(D,T())) override;
};

template <unsigned D, typename T>
struct AdNeumannZeroBoundaryScheme<D,T,tag::CENTRAL> final : AdBoundarySchemeBase<D,T> {
  AdNeumannZeroBoundaryScheme() {}
  virtual int getExtraExtent() override;
  virtual void operator() (
    std::vector<std::vector<T>>& fOut, T& f0, std::vector<std::vector<T>>& fIn, std::vector<int>& normal, std::vector<T>& u ) override;
};

//////////////////////////////////////////// UPWIND /////////////////////////////////////////////////////////
template <unsigned D, typename T>
struct AdvectionScheme<D,T,tag::UPWIND> final : AdvectionSchemeBase<D,T> {
  AdvectionScheme() {}
  virtual T operator()(T& f0, std::vector<std::vector<T>>& f, std::vector<std::vector<T>>& F, std::vector<T>& u) override;
};

template <unsigned D, typename T>
struct AdNeumannZeroBoundaryScheme<D,T,tag::UPWIND> final : AdBoundarySchemeBase<D,T> {
  AdNeumannZeroBoundaryScheme() {}
  virtual int getExtraExtent() override;
  virtual void operator() (
    std::vector<std::vector<T>>& fOut, T& f0, std::vector<std::vector<T>>& fIn, std::vector<int>& normal, std::vector<T>& u ) override;
};

//////////////////////////////////////////// UPWIND_2_ORDER /////////////////////////////////////////////////////////
template <unsigned D, typename T>
struct AdvectionScheme<D,T,tag::UPWIND_2_ORDER> final : AdvectionSchemeBase<D,T> {
  AdvectionScheme() {}
  virtual T operator()(T& f0, std::vector<std::vector<T>>& f, std::vector<std::vector<T>>& F, std::vector<T>& u) override;
};

}  // namespace fd

}  // namespace olb

#endif
