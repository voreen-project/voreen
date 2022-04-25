/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2020 Jan E. Marquardt, Anna Husfeldt
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

#ifndef INDIC_ALT_3D_H
#define INDIC_ALT_3D_H

#include "indicatorBaseF3D.h"
#include "sdf.h"
#include "utilities/functorPtr.h"

namespace olb {

template <typename S>
class IndicElongation3D : public IndicatorF3D<S> {
protected:
  FunctorPtr<IndicatorF3D<S>> _f;
  Vector<S,3> _h;
  Vector<S,3> getEstimatedCenter();
public:
  IndicElongation3D( FunctorPtr<IndicatorF3D<S>> f, Vector<S,3> h );
  S signedDistance( const Vector<S,3>& input );
  Vector<S,3> getElongation();
  void setElongation(Vector<S,3> h);
};

}
#endif