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
 * 2D postprocessor to locally perform a generic chemical reaction
 *  -- generic implementation
 */
#ifndef REACTION_POST_PROCESSOR_2D_HH
#define REACTION_POST_PROCESSOR_2D_HH

#include <stdexcept>
namespace olb {


///////////////////////////////////// class ReactionPostProcessor2D /////////////////////////////////////

template <typename T, typename DESCRIPTOR>
ReactionPostProcessor2D<T,DESCRIPTOR>::ReactionPostProcessor2D ( int x0, int x1, int y0, int y1,
    std::vector<Rate<T>*> rate, std::vector<std::vector<ReactingSpeciesBase2D<T>*>> species )
  : _x0(x0), _x1(x1), _y0(y0), _y1(y1),
    _rate(rate), _species(species)
{
  for (auto&& specieOnReaction : species) {
    _size.push_back(specieOnReaction.size());
  }
}

template <typename T, typename DESCRIPTOR>
ReactionPostProcessor2D<T,DESCRIPTOR>::ReactionPostProcessor2D (
  std::vector<Rate<T>*> rate, std::vector<std::vector<ReactingSpeciesBase2D<T>*>> species )
  : ReactionPostProcessor2D<T,DESCRIPTOR>(0,0,0,0,rate,species)
{}

template <typename T, typename DESCRIPTOR>
void ReactionPostProcessor2D<T,DESCRIPTOR>::process(BlockLattice<T,DESCRIPTOR>& blockLattice)
{
  processSubDomain(blockLattice, _x0, _x1, _y0, _y1);
}

template <typename T, typename DESCRIPTOR>
void ReactionPostProcessor2D<T,DESCRIPTOR>::processSubDomain ( BlockLattice<T,DESCRIPTOR>& blockLattice,
    int x0, int x1, int y0, int y1 )
{
  int newX0, newX1, newY0, newY1;
  if ( util::intersect ( _x0, _x1, _y0, _y1,
                         x0, x1, y0, y1,
                         newX0, newX1, newY0, newY1 ) ) {

    for (int iX=newX0-1; iX<=newX1+1; ++iX) {
      for (int iY=newY0-1; iY<=newY1+1; ++iY) {

        // working out reaction rates
        std::vector<T> rate;
        for (long unsigned i=0; i<_size.size(); ++i) {
          // storing the old values of the local field from the partner lattices
          std::vector<T> fields;
          for (auto&& specie : _species[i]) {
            fields.push_back ( specie->getField(iX, iY) );
          }
          // computing the rates from the old values of the local fields as per rate's own law
          rate.push_back ( _rate[i]->compute(fields) );
        }

        // resetting the sources
        for (auto&& specieOnReaction : _species) {
          for (auto&& specie : specieOnReaction) {
            specie->setSource ( T(), iX, iY );
          }
        }

        // updating local field in the partner lattice: source = sum_i rate_i*coeff
        T oldSource = T();
        for (long unsigned i=0; i<_size.size(); ++i) {
          for (long unsigned j=0; j<_size[i]; ++j) {
            oldSource = _species[i][j]->getSource(iX, iY);
            _species[i][j]->setSource ( oldSource + rate[i] * _species[i][j]->getStoichioCoeff(), iX, iY );
          }
        }
      }
    }
  }
}


///////////////////////////////////// class ReactionGenerator2D /////////////////////////////////////

template <typename T, typename DESCRIPTOR>
ReactionGenerator2D<T,DESCRIPTOR>::ReactionGenerator2D ( int x0_, int x1_, int y0_, int y1_,
    std::vector<Rate<T>*> rate, std::vector<std::vector<ReactingSpeciesBase2D<T>*>> species )
  : PostProcessorGenerator2D<T,DESCRIPTOR>(x0_, x1_, y0_, y1_),
    _rate(rate), _species(species)
{}

template <typename T, typename DESCRIPTOR>
ReactionGenerator2D<T,DESCRIPTOR>::ReactionGenerator2D (
  std::vector<Rate<T>*> rate, std::vector<std::vector<ReactingSpeciesBase2D<T>*>> species )
  : ReactionGenerator2D<T,DESCRIPTOR>(0,0,0,0,rate,species)
{}

template<typename T, typename DESCRIPTOR>
PostProcessor2D<T,DESCRIPTOR>* ReactionGenerator2D<T,DESCRIPTOR>::generate () const
{
  return new ReactionPostProcessor2D<T,DESCRIPTOR>(
           this->x0,this->x1,this->y0,this->y1, _rate, _species);
}

template<typename T, typename DESCRIPTOR>
PostProcessorGenerator2D<T,DESCRIPTOR>* ReactionGenerator2D<T,DESCRIPTOR>::clone() const
{
  return new ReactionGenerator2D<T,DESCRIPTOR>(*this);
}

}  // namespace olb

#endif
