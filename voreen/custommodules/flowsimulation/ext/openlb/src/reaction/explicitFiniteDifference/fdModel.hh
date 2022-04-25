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
 * Finite-difference models.
 *  -- generic implementation
 */
#ifndef FD_MODEL_HH
#define FD_MODEL_HH

#include "fdModel.h"

namespace olb {


/////////////////////// class FdModel ///////////////////////

template<typename T, typename DESCRIPTOR>
FdModel<T,DESCRIPTOR>::FdModel(T diffusivity)
  : _diffusivity(diffusivity)
{ }


/////////////////////// class FdAdvectionDiffusionModel ///////////////////////

template<typename T, typename DESCRIPTOR, typename SCHEME_ADV, typename SCHEME_DIFF>
FdAdvectionDiffusionModel<T,DESCRIPTOR,SCHEME_ADV,SCHEME_DIFF>::FdAdvectionDiffusionModel(T diffusivity)
  : FdModel<T,DESCRIPTOR>(diffusivity),
    _advectionScheme(std::make_shared<fd::AdvectionScheme<DESCRIPTOR::d,T,SCHEME_ADV>>()),
    _diffusionScheme(std::make_shared<fd::DiffusionScheme<DESCRIPTOR::d,T,SCHEME_DIFF>>())
{ }

template<typename T, typename DESCRIPTOR, typename SCHEME_ADV, typename SCHEME_DIFF>
int FdAdvectionDiffusionModel<T,DESCRIPTOR,SCHEME_ADV,SCHEME_DIFF>::extent()
{
  return util::max<int>(SCHEME_ADV::extent, SCHEME_DIFF::extent);
}

template<typename T, typename DESCRIPTOR, typename SCHEME_ADV, typename SCHEME_DIFF>
void FdAdvectionDiffusionModel<T,DESCRIPTOR,SCHEME_ADV,SCHEME_DIFF>::
operator()(T* fNew, T* f0, std::vector<std::vector<T>>& f, std::vector<std::vector<T>>& F, Cell<T,DESCRIPTOR>& cell)
{
  T uArr[DESCRIPTOR::d];
  cell.computeU(uArr);
  std::vector<T> u(DESCRIPTOR::d, T());
  for (int iD=0; iD<DESCRIPTOR::d; ++iD) {
    u[iD] = uArr[iD];
  }
  *fNew = *f0 - _advectionScheme->operator()(*f0, f, F, u) + _diffusionScheme->operator()(*f0, f, F, this->_diffusivity);
}


/////////////////////// class FdAdvectionDiffusionModel ///////////////////////

template<typename T, typename DESCRIPTOR, typename SCHEME_ADV, typename SCHEME_DIFF>
FdAdvectionDiffusionModelWithAntiDiffusion<T,DESCRIPTOR,SCHEME_ADV,SCHEME_DIFF>::FdAdvectionDiffusionModelWithAntiDiffusion (
  T diffusivity, T antiDiffusionTunig)
  : FdAdvectionDiffusionModel<T,DESCRIPTOR,SCHEME_ADV,SCHEME_DIFF>(diffusivity),
    _antiDiffusionTuning(antiDiffusionTunig)
{ }

template<typename T, typename DESCRIPTOR, typename SCHEME_ADV, typename SCHEME_DIFF>
void FdAdvectionDiffusionModelWithAntiDiffusion<T,DESCRIPTOR,SCHEME_ADV,SCHEME_DIFF>::
operator()(T* fNew, T* f0, std::vector<std::vector<T>>& f, std::vector<std::vector<T>>& F, Cell<T,DESCRIPTOR>& cell)
{
  T uArr[DESCRIPTOR::d];
  cell.computeU(uArr);
  std::vector<T> u(DESCRIPTOR::d, T()),  uCorr(DESCRIPTOR::d, T());
  for (int iD=0; iD<DESCRIPTOR::d; ++iD) {
    u[iD] = uArr[iD];
    uCorr[iD] = uArr[iD] * _antiDiffusionTuning;
  }
  *fNew = *f0 - this->_advectionScheme->operator()(*f0, f, F, u) + this->_diffusionScheme->operator()(*f0, f, F, this->_diffusivity, uCorr);
}


}  // namespace olb

#endif
