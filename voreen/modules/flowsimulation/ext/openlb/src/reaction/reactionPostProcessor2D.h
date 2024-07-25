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
 * 2D postprocessor to locally perform a generic chemical reactions
 *  -- header file
 */
#ifndef REACTION_POST_PROCESSOR_2D_H
#define REACTION_POST_PROCESSOR_2D_H

#include "rate.h"
#include "reactingSpecies2D.h"

namespace olb {

////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
 * Postprocessor to perform a generic chemical reactions
 */
template<typename T, typename DESCRIPTOR>
class ReactionPostProcessor2D : public LocalPostProcessor2D<T,DESCRIPTOR> {
public:
  ReactionPostProcessor2D ( int x0, int x1, int y0, int y1,
                            std::vector<Rate<T>*> rate, std::vector<std::vector<ReactingSpeciesBase2D<T>*>> species );
  ReactionPostProcessor2D ( std::vector<Rate<T>*> rate, std::vector<std::vector<ReactingSpeciesBase2D<T>*>> species );
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
  std::vector<long unsigned> _size;
  std::vector<Rate<T>*> _rate;
  std::vector<std::vector<ReactingSpeciesBase2D<T>*>> _species;
};

template<typename T, typename DESCRIPTOR>
class ReactionGenerator2D final : public PostProcessorGenerator2D<T,DESCRIPTOR> {
public:
  ReactionGenerator2D ( int x0_, int x1_, int y0_, int y1_,
                        std::vector<Rate<T>*> rate, std::vector<std::vector<ReactingSpeciesBase2D<T>*>> species );
  ReactionGenerator2D ( std::vector<Rate<T>*> rate, std::vector<std::vector<ReactingSpeciesBase2D<T>*>> species );
  PostProcessor2D<T,DESCRIPTOR>* generate() const override;
  PostProcessorGenerator2D<T,DESCRIPTOR>* clone() const override;
private:
  std::vector<Rate<T>*> _rate;
  std::vector<std::vector<ReactingSpeciesBase2D<T>*>> _species;
};


}  // namespace olb

#endif
