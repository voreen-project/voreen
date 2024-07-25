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
 * Base class to implement a generic reaction rate - to be derived to get the specific rate
 *  -- header file
 */
#ifndef RATE_H
#define RATE_H


namespace olb {

// Base abstract class for reaction rates
template<typename T>
class Rate {
public:
  Rate(std::vector<T> params);
  virtual T compute(std::vector<T> localFieldValues) =0;
protected:
  std::vector<T> _params;
};

// Class implementing costant reaction rate
template<typename T>
class ConstantRate final : public Rate<T> {
public:
  // params must be passed as a mock parameter to ensure consistency with parent class's constructor
  ConstantRate();
  virtual T compute(std::vector<T> localFieldValues) override;
};

/** Class implementing exponentially-decreasing reaction rate on the 1st reacting species,
 * that is: nu = [A]/t0, with t0 being the time constant in lattice units.
 * */
template<typename T>
class ExpOn1stSpecieRate final : public Rate<T> {
public:
  ExpOn1stSpecieRate(int t0);
  T compute(std::vector<T> localFieldValues) override;
};

/** Class implementing Monod kinetics, with 1st field being substrate concentration [S], 2nd being bacteria concentration [X]:
 * nu = mu * [X];
 * mu = muMax * [S] / ([S] + Ks)
 * */
template<typename T>
class MonodRate final : public Rate<T> {
public:
  MonodRate(T muMax, T Ks);
  T compute(std::vector<T> localFieldValues) override;
};

/** Class implementing Haldane kinetics, with 1st field being substrate concentration [S], 2nd being bacteria concentration [X]:
 * nu = mu * [X];
 * mu = muMax * [S] / ([S] + Ks + [S]^2/KI)
 * */
template<typename T>
class HaldaneRate final : public Rate<T> {
public:
  HaldaneRate(T muMax, T Ks, T KI);
  T compute(std::vector<T> localFieldValues) override;
};

}  // namespace olb

#endif
