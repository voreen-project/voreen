/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2021 Claudius Holeksa
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

namespace olb {

namespace FreeSurface {

template<typename T, size_t S>
std::array<T,S> solvePivotedLU(std::array<std::array<T,S>,S>& matrix, const std::array<T,S>& b, size_t N) {
  std::array<T,S> x;
  std::array<T,S> pivots;
  for(size_t i = 0; i < S; ++i){
    pivots[i] = i;
    x[i] = 0.;
  }

  N = std::min(N,S);

  for(size_t i = 0; i < N; ++i){

    T max = 0.;
    size_t max_index = i;

    for(size_t j = i; j < N; ++j){
      T abs = std::abs(matrix[pivots[j]][i]);
      if(abs > max){
        max_index = j;
        max = abs;
      }
    }

    //Normally we would have a tolerance, but we clamp the result instead. So even if it is NaN it does not matter
    //if(max < 1e-16){
    //  return x;
    //}

    if(max_index != i){
      size_t tmp_index = pivots[i];
      pivots[i] = pivots[max_index];
      pivots[max_index] = tmp_index;
    }

    for(size_t j = i + 1; j < N; ++j){
      matrix[pivots[j]][i] /= matrix[pivots[i]][i];

      for(size_t k = i + 1; k < N; ++k){

        matrix[pivots[j]][k] -= matrix[pivots[j]][i] * matrix[pivots[i]][k];
      }
    }
  }

  for(size_t i = 0; i  < N; ++i){
    x[i] = b[pivots[i]];

    for(size_t j = 0; j < i; ++j){
      x[i] -= matrix[pivots[i]][j] * x[j];
    }
  }

  for(size_t i = N; i > 0; --i){
    for(size_t j = i; j < N; ++j){
      x[i-1] -= matrix[pivots[i-1]][j] * x[j];
    }

    x[i-1] /= matrix[pivots[i-1]][i-1];
  }

  return x;
}

template<typename T, typename DESCRIPTOR>
void initialize(SuperLattice<T,DESCRIPTOR>& lattice) {
  lattice.executePostProcessors(FreeSurface::Stage0());
  lattice.executePostProcessors(FreeSurface::Stage1());
  lattice.executePostProcessors(FreeSurface::Stage2());
  lattice.executePostProcessors(FreeSurface::Stage3());
  lattice.executePostProcessors(FreeSurface::Stage4());
}

}

}
