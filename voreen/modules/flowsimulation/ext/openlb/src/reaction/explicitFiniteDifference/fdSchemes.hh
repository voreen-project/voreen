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
 *  -- generic implementation
 */
#ifndef FD_SCHEMES_HH
#define FD_SCHEMES_HH

namespace olb {

namespace fd {

//////////////////////////////////////////// CENTRAL /////////////////////////////////////////////////////////
template <unsigned D, typename T>
T AdvectionScheme<D,T,tag::CENTRAL>::
operator()(T& f0, std::vector<std::vector<T>>& f, std::vector<std::vector<T>>& F, std::vector<T>& u)
{
  T fNew = 0.;
  for (unsigned iD=0; iD<D; ++iD) {
    fNew += u[iD]*(F[0][iD] - f[0][iD]);
  }
  return fNew / 2.;
}

template <unsigned D, typename T>
T DiffusionScheme<D,T,tag::CENTRAL>::
operator()(T& f0, std::vector<std::vector<T>>& f, std::vector<std::vector<T>>& F, T& diff, std::vector<T> u)
{
  T fNew  = 0.;
  T fCorr = 0.; // Negative-diffusivity correction to counterbalance upwind schemes' numerical diffusivity
  for (unsigned iD=0; iD<D; ++iD) {
    fNew  +=  F[0][iD] + f[0][iD];
    fCorr += (F[0][iD] + f[0][iD] - 2.*f0) * u[iD];
  }
  fNew -= 2.*D*f0;
  return fNew * diff - 0.5*fCorr;
}

template <unsigned D, typename T>
int AdNeumannZeroBoundaryScheme<D,T,tag::CENTRAL>::getExtraExtent()
{
  return 0;
}

template <unsigned D, typename T>
void AdNeumannZeroBoundaryScheme<D,T,tag::CENTRAL>::
operator()(std::vector<std::vector<T>>& fOut, T& f0, std::vector<std::vector<T>>& fIn, std::vector<int>& normal, std::vector<T>& u)
{
  for (unsigned iD=0; iD<D; ++iD) {
    fOut[0][iD] = fIn[0][iD];
  }
}

//////////////////////////////////////////// UPWIND /////////////////////////////////////////////////////////
template <unsigned D, typename T>
T AdvectionScheme<D,T,tag::UPWIND>::
operator()(T& f0, std::vector<std::vector<T>>& f, std::vector<std::vector<T>>& F, std::vector<T>& u)
{
  T fNew = 0.;
  for (unsigned iD=0; iD<D; ++iD) {
    fNew += u[iD] * ( u[iD]==0. ? 0.
                      : ( u[iD] >0. ? f0-f[0][iD]
                          : F[0][iD]-f0 ) );
  }
  return fNew;
}

template <unsigned D, typename T>
int AdNeumannZeroBoundaryScheme<D,T,tag::UPWIND>::getExtraExtent()
{
  return 0;
}

template <unsigned D, typename T>
void AdNeumannZeroBoundaryScheme<D,T,tag::UPWIND>::
operator()(std::vector<std::vector<T>>& fOut, T& f0, std::vector<std::vector<T>>& fIn, std::vector<int>& normal, std::vector<T>& u)
{
  for (unsigned iD=0; iD<D; ++iD) {
    fOut[0][iD] = ( u[iD]*normal[iD]==0 ? fIn[0][iD] : f0 );
  }
}

//////////////////////////////////////////// UPWIND_2_ORDER /////////////////////////////////////////////////////////
template <unsigned D, typename T>
T AdvectionScheme<D,T,tag::UPWIND_2_ORDER>::
operator()(T& f0, std::vector<std::vector<T>>& f, std::vector<std::vector<T>>& F, std::vector<T>& u)
{
  T fNew = 0.;
  for (unsigned iD=0; iD<D; ++iD) {
    fNew += u[iD] * ( u[iD]==0. ? 0.
                      : ( u[iD] >0. ? 3.*f0 - 4.*f[0][iD] + f[1][iD]
                          : -F[1][iD] + 4.*f[0][iD] - 3.*f0 ) );
  }
  return 0.5 * fNew;
}

}  // namespace fd

}  // namespace olb

#endif
