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
 * 2D postprocessor to locally update the finite-diffusion advection-diffusion external field.
 *  -- header file
 */
#ifndef AD_POST_PROCESSOR_2D_H
#define AD_POST_PROCESSOR_2D_H

#include "fdModel.h"
#include "utilities/aliases.h"
#include <type_traits>

namespace olb {

////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
 * Virtual parent to all the finite-difference postprocessors.
 * Its raison d'etre is to enforce all the finite-difference postprocessors
 * to contain a reference to the simulation's timestep.
 */
template<typename T, typename DESCRIPTOR, typename FIELD=descriptors::AD_FIELD, typename SOURCE=void>
class FdBasePostProcessor2D : public LocalPostProcessor2D<T,DESCRIPTOR> {
protected:
  FdBasePostProcessor2D(size_t& iT);
  int getIT() const;
  std::vector<std::vector<T>> initMatrix(int size0);
  template<typename SOURCE_CHK=SOURCE>
  std::enable_if_t<  std::is_void<SOURCE_CHK>::value, void> applySourceTerm(T* fNew, Cell<T,DESCRIPTOR>& cell);
  template<typename SOURCE_CHK=SOURCE>
  std::enable_if_t<! std::is_void<SOURCE_CHK>::value, void> applySourceTerm(T* fNew, Cell<T,DESCRIPTOR>& cell);
private:
  // Reference to the simulation's timestep in order to take track of odd / even timesteps
  std::size_t& _iT;
};

/*
 * Virtual parent to all the finite-difference postprocessor generators.
 * Its raison d'etre is to enforce all the finite-difference postprocessors
 * to contain a reference to the simulation's timestep.
 */
template<typename T, typename DESCRIPTOR, typename FIELD=descriptors::AD_FIELD, typename SOURCE=void>
class FdBasePostProcessorGenerator2D : public PostProcessorGenerator2D<T,DESCRIPTOR> {
public:
  FdBasePostProcessorGenerator2D(int x0_, int x1_, int y0_, int y1_, std::size_t& iT);
protected:
  std::size_t& _iT;
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
 * Postprocessor to update the finite-difference advection-diffusion external field.
 */
template<typename T, typename DESCRIPTOR, typename FIELD=descriptors::AD_FIELD, typename SOURCE=void>
class FdPostProcessor2D : public FdBasePostProcessor2D<T,DESCRIPTOR,FIELD,SOURCE> {
public:
  FdPostProcessor2D ( int x0, int x1, int y0, int y1, std::size_t& iT, std::shared_ptr<FdModel<T,DESCRIPTOR>> model);
  FdPostProcessor2D ( std::size_t& iT, std::shared_ptr<FdModel<T,DESCRIPTOR>> model);
  int extent() const override
  {
    return 1;
  }
  int extent(int whichDirection) const override
  {
    return 1;
  }
  void process(BlockLattice<T,DESCRIPTOR>& blockLattice) override;
  void processSubDomain ( BlockLattice<T,DESCRIPTOR>& blockLattice,
                          int x0, int x1, int y0, int y1 ) override;
private:
  int _x0, _x1, _y0, _y1;
  std::shared_ptr<FdModel<T,DESCRIPTOR>> _model;
  std::vector<std::vector<T>> f, F;
};

template<typename T, typename DESCRIPTOR, typename FIELD=descriptors::AD_FIELD, typename SOURCE=void>
class FdPostProcessorGenerator2D final : public FdBasePostProcessorGenerator2D<T,DESCRIPTOR,FIELD,SOURCE> {
public:
  FdPostProcessorGenerator2D(int x0_, int x1_, int y0_, int y1_, std::size_t& iT, std::shared_ptr<FdModel<T,DESCRIPTOR>> model);
  FdPostProcessorGenerator2D(size_t& iT, std::shared_ptr<FdModel<T,DESCRIPTOR>> model);
  PostProcessor2D<T,DESCRIPTOR>* generate() const override;
  PostProcessorGenerator2D<T,DESCRIPTOR>* clone() const override;
private:
  std::shared_ptr<FdModel<T,DESCRIPTOR>> _model;
};


}  // namespace olb

#endif
