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
 * Classes to implement a generic reaction process
 *  -- generic implementation
 */
#ifndef RATE_HH
#define RATE_HH

namespace olb {

///////////////////////////////////////// class Rate /////////////////////////////////////////
template<typename T>
Rate<T>::Rate(std::vector<T> params)
  : _params(params)
{}


///////////////////////////////////////// class ConstantRate /////////////////////////////////////////
template<typename T>
ConstantRate<T>::ConstantRate()
  : Rate<T>(std::vector<T>{})
{}

template<typename T>
T ConstantRate<T>::compute(std::vector<T> localFieldValues)
{
  return 1.0;
}


///////////////////////////////////////// class ExpOn1stSpecieRate  /////////////////////////////////////////
template<typename T>
ExpOn1stSpecieRate<T>::ExpOn1stSpecieRate(int t0)
  : Rate<T>(std::vector<T>{(T)(t0)})
{}

template<typename T>
T ExpOn1stSpecieRate<T>::compute(std::vector<T> localFieldValues)
{
  return localFieldValues[0] / this->_params[0];
}


///////////////////////////////////////// class MonodRate  /////////////////////////////////////////
template<typename T>
MonodRate<T>::MonodRate(T muMax, T Ks)
  : Rate<T>(std::vector<T>{muMax, Ks})
{}

template<typename T>
T MonodRate<T>::compute(std::vector<T> localFieldValues)
{
  // _params[0] = muMax
  // _params[1] = Ks
  // localFieldValues[0] = [S]
  // localFieldValues[1] = [X]
  return this->_params[0] * localFieldValues[0] * localFieldValues[1] / (localFieldValues[0] + this->_params[1]);
}


///////////////////////////////////////// class HaldaneRate  /////////////////////////////////////////
template<typename T>
HaldaneRate<T>::HaldaneRate(T muMax, T Ks, T KI)
  : Rate<T>(std::vector<T>{muMax, Ks, KI})
{}

template<typename T>
T HaldaneRate<T>::compute(std::vector<T> localFieldValues)
{
  // _params[0] = muMax
  // _params[1] = Ks
  // _params[2] = KI
  // localFieldValues[0] = [S]
  // localFieldValues[1] = [X]
  return this->_params[0] * localFieldValues[0] * localFieldValues[1]
              / (localFieldValues[0] + this->_params[1] + localFieldValues[0]*localFieldValues[0]/this->_params[2]) ;
}


}  // namespace olb

#endif
