/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2019 University of Muenster, Germany,                        *
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

#include "voreen/core/datastructures/roi/roiunion.h"
#include "voreen/core/datastructures/geometry/geometrysequence.h"

namespace voreen {

ROIUnion::ROIUnion(Grid grid) : ROIAggregation(grid) {
}

bool ROIUnion::inROI(tgt::vec3 p) const {
    for(size_t i=0; i<children_.size(); i++)
        if(children_[i]->inROI(p))
            return true;

    return false;
}

Geometry* ROIUnion::generateMesh() const {
    return new GeometrySequence();
}

Geometry* ROIUnion::generateMesh(tgt::plane /*pl*/) const {
    return new GeometrySequence();
}

tgt::Bounds ROIUnion::getBoundingBox() const {
    tgt::Bounds b;
    for(size_t i=0; i<children_.size(); i++)
        b.addVolume(children_[i]->getBoundingBox());
    return b;
}

std::vector<const ControlPoint*> ROIUnion::getControlPoints() const {
    //TODO
    return std::vector<const ControlPoint*>();
}

std::vector<const ControlPoint*> ROIUnion::getControlPoints(tgt::plane /*pl*/) const {
    //TODO
    return std::vector<const ControlPoint*>();
}

} // namespace
