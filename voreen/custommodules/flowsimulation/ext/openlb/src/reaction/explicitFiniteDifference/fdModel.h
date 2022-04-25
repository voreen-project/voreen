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
 *  -- header file
 */
#ifndef FD_MODEL_H
#define FD_MODEL_H

#include "fdDescriptorField.h"
#include "fdSchemes.h"

namespace olb {

/*
 * Basic virtual class for all the finite-difference models.
 */
template<typename T, typename DESCRIPTOR>
class FdModel {
public:
  FdModel(T diffusivity);
  virtual int extent() =0;
  virtual void operator()(T* fNew, T* f0, std::vector<std::vector<T>>& f, std::vector<std::vector<T>>& F, Cell<T,DESCRIPTOR>& cell) =0;
protected:
  T _diffusivity;
};

/*
 * Class for finite-difference advection-diffusion model
 */
template<typename T, typename DESCRIPTOR, typename SCHEME_ADV, typename SCHEME_DIFF>
class FdAdvectionDiffusionModel : public FdModel<T,DESCRIPTOR> {
public:
  FdAdvectionDiffusionModel(T diffusivity);
  virtual int extent() override;
  virtual void operator()(T* fNew, T* f0, std::vector<std::vector<T>>& f, std::vector<std::vector<T>>& F, Cell<T,DESCRIPTOR>& cell) override;
protected:
  std::shared_ptr<fd::AdvectionSchemeBase<DESCRIPTOR::d,T>> _advectionScheme;
  std::shared_ptr<fd::DiffusionSchemeBase<DESCRIPTOR::d,T>> _diffusionScheme;
};

/*
 * Class for finite-difference advection-diffusion model,
 * with an anti-diffusivity term to counterbalance upwind numerical diffusion
 */
template<typename T, typename DESCRIPTOR, typename SCHEME_ADV, typename SCHEME_DIFF>
class FdAdvectionDiffusionModelWithAntiDiffusion : public FdAdvectionDiffusionModel<T,DESCRIPTOR,SCHEME_ADV,SCHEME_DIFF> {
public:
  FdAdvectionDiffusionModelWithAntiDiffusion(T diffusivity, T antiDiffusionTunig=1.);
  virtual void operator()(T* fNew, T* f0, std::vector<std::vector<T>>& f, std::vector<std::vector<T>>& F, Cell<T,DESCRIPTOR>& cell) override;
protected:
  T _antiDiffusionTuning;
};

}  // namespace olb

#endif
