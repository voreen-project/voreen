/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2014-2016 Cyril Masquelier, Mathias J. Krause, Albert Mink
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

#include "io/stlReader.h"
#include "indicatorF3D.h"
#include "indicatorF3D.hh"

namespace olb {

// stl indicator functors
template class STLreader<double>;

// indicator functors 3D
template class IndicatorCircle3D<double>;
template class IndicatorSphere3D<double>;
template class IndicatorLayer3D<double>;
template class IndicatorCylinder3D<double>;
template class IndicatorCone3D<double>;
template class IndicatorCuboid3D<double>;

// creator functions for primitives
template IndicatorCircle3D<double>* createIndicatorCircle3D(XMLreader const& params, bool verbose);
template IndicatorSphere3D<double>* createIndicatorSphere3D(XMLreader const& params, bool verbose);
template IndicatorCylinder3D<double>* createIndicatorCylinder3D(XMLreader const& params, bool verbose);
template IndicatorCone3D<double>* createIndicatorCone3D(XMLreader const& params, bool verbose);
template IndicatorCuboid3D<double>* createIndicatorCuboid3D(XMLreader const& params, bool verbose);

// creator functions for arithmetic operations
template IndicatorF3D<double>* createIndicatorUnion3D( XMLreader const& params, bool verbose );
template IndicatorF3D<double>* createIndicatorWithout3D( XMLreader const& params, bool verbose );
template IndicatorF3D<double>* createIndicatorIntersection3D( XMLreader const& params, bool verbose );

// root of all creator functions
template IndicatorF3D<double>* createIndicatorF3D(XMLreader const& params, bool verbose);

}
