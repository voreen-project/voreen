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
 * Updater processors for finite-difference external fields: they copy the FIELD[1] to FIELD[0] entries.
 * Functional for generic field access.
 *  -- header file
 */
#ifndef FD_UPDATER_H
#define FD_UPDATER_H

namespace olb {

template<typename T, typename DESCRIPTOR>
class FdUpdaterBase {
protected:
  FdUpdaterBase() { }
public:
  virtual void operator()(Cell<T,DESCRIPTOR>& cell) =0;
};

template<typename T, typename DESCRIPTOR, typename FIELD>
class FdUpdater : public FdUpdaterBase<T,DESCRIPTOR> {
public:
  FdUpdater() : FdUpdaterBase<T,DESCRIPTOR>() { }
  virtual void operator()(Cell<T,DESCRIPTOR>& cell) override
  {
    T* f0 = &cell.template getFieldPointer<FIELD>()[0];
    T  f1 =  cell.template getFieldPointer<FIELD>()[1];
    *f0 = f1;
  }
};

}  // namespace olb

#endif
