/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2019 Adrian Kummerlaender
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

#ifndef CELL_D_H
#define CELL_D_H

#include "cell.h"

namespace olb {

template<typename T, typename DESCRIPTOR>
class SingleCellBlockD {
protected:
  BlockStaticPopulationD<T,DESCRIPTOR> _localStaticPopulationD;
  BlockStaticFieldsD<T,DESCRIPTOR>     _localStaticFieldsD;
  BlockDynamicFieldsD<T,DESCRIPTOR>    _localDynamicFieldsD;
  BlockDynamicsMap<T,DESCRIPTOR>       _localDynamicsMap;

public:
  SingleCellBlockD() = default;

};

template<typename T, typename DESCRIPTOR>
class CellD : private SingleCellBlockD<T,DESCRIPTOR>, public Cell<T,DESCRIPTOR> {
public:
  CellD():
    SingleCellBlockD<T,DESCRIPTOR>(),
    Cell<T,DESCRIPTOR>(
      this->_localStaticPopulationD,
      this->_localStaticFieldsD,
      this->_localDynamicFieldsD,
      this->_localDynamicsMap,
      0)
  { }

  using Cell<T,DESCRIPTOR>::operator=;

  CellD(ConstCell<T,DESCRIPTOR> rhs): CellD()
  {
    this->operator=(rhs);
  }

};

}

#endif
