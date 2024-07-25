/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2021 Jan E. Marquardt
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

#ifndef INDIC_MOD_H
#define INDIC_MOD_H

#include "utilities/aliases.h"
#include "utilities/functorPtr.h"
#include "sdf.h"

namespace olb {

template <typename S, unsigned D>
class IndicInverse : public IndicatorF<S,D> {
protected:
  FunctorPtr<IndicatorF<S,D>> _f;
public:
  IndicInverse( FunctorPtr<IndicatorF<S,D>> f, PhysR<S,D> min, PhysR<S,D> max );
  S signedDistance( const Vector<S,D>& input );
  //bool operator() (bool output[], const S input[]) override;
};

}
#endif