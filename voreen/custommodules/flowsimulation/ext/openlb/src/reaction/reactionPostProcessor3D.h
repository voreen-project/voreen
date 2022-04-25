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
 * 3D postprocessor to locally perform a generic chemical reactions
 *  -- header file
 */
#ifndef REACTION_POST_PROCESSOR_3D_H
#define REACTION_POST_PROCESSOR_3D_H

#include "rate.h"
#include "reactingSpecies3D.h"

namespace olb {

////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
 * Postprocessor to perform a generic chemical reactions
 */
template<typename T, typename DESCRIPTOR>
class ReactionPostProcessor3D : public LocalPostProcessor3D<T,DESCRIPTOR> {
public:
  ReactionPostProcessor3D ( int x0, int x1, int y0, int y1, int z0, int z1,
                            std::vector<Rate<T>*> rate, std::vector<std::vector<ReactingSpeciesBase3D<T>*>> species );
  ReactionPostProcessor3D ( std::vector<Rate<T>*> rate, std::vector<std::vector<ReactingSpeciesBase3D<T>*>> species );
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
                          int x0, int x1, int y0, int y1, int z0, int z1 ) override;
private:
  int _x0, _x1, _y0, _y1, _z0, _z1;
  std::vector<long unsigned> _size;
  std::vector<Rate<T>*> _rate;
  std::vector<std::vector<ReactingSpeciesBase3D<T>*>> _species;
};

template<typename T, typename DESCRIPTOR>
class ReactionGenerator3D final : public PostProcessorGenerator3D<T,DESCRIPTOR> {
public:
  ReactionGenerator3D ( int x0_, int x1_, int y0_, int y1_, int z0_, int z1_,
                        std::vector<Rate<T>*> rate, std::vector<std::vector<ReactingSpeciesBase3D<T>*>> species );
  ReactionGenerator3D ( std::vector<Rate<T>*> rate, std::vector<std::vector<ReactingSpeciesBase3D<T>*>> species );
  PostProcessor3D<T,DESCRIPTOR>* generate() const override;
  PostProcessorGenerator3D<T,DESCRIPTOR>* clone() const override;
private:
  std::vector<Rate<T>*> _rate;
  std::vector<std::vector<ReactingSpeciesBase3D<T>*>> _species;
};


}  // namespace olb

#endif
