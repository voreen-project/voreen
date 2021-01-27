/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2021 University of Muenster, Germany,                        *
 * Department of Computer Science.                                                 *
 * For a list of authors please refer to the file "CREDITS.txt".                   *
 *                                                                                 *
 * This file is part of the Voreen software package. Voreen is free software:      *
 * you can redistribute it and/or modify it under the terms of the GNU General     *
 * Public License version 2 as published by the Free Software Foundation.          *
 *                                                                                 *
 * Voreen is distributed in the hope that it will be useful, but WITHOUT ANY       *
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR   *
 * A PARTICULAR PURPOSE. See the GNU General Public License for more details.      *
 *                                                                                 *
 * You should have received a copy of the GNU General Public License in the file   *
 * "LICENSE.txt" along with this file. If not, see <http://www.gnu.org/licenses/>. *
 *                                                                                 *
 * For non-commercial academic use see the license exception specified in the file *
 * "LICENSE-academic.txt". To get information about commercial licensing please    *
 * contact the authors.                                                            *
 *                                                                                 *
 ***********************************************************************************/

#include "seeds.h"

#include "voreen/core/datastructures/geometry/pointsegmentlistgeometry.h"
#include "voreen/core/datastructures/geometry/geometrysequence.h"

namespace voreen {

void getSeedListsFromPorts(std::vector<PortDataPointer<Geometry>>& geom, PointSegmentListGeometryVec3& seeds) {
    auto mat = tgt::mat4::createIdentity();
    for (size_t i=0; i<geom.size(); i++) {
        if(!seeds.collectSegmentsFromGeometry(*geom.at(i).get())) {
            LWARNINGC("voreen.RandomWalker.util", "Invalid geometry. PointSegmentListGeometry<vec3> or GeometrySequence expected.");
        }
    }
}

}
