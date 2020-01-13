/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2020 University of Muenster, Germany,                        *
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

#include "voreen/core/datastructures/roi/roisubtract.h"
#include "voreen/core/datastructures/geometry/meshlistgeometry.h"

namespace voreen {

ROISubtract::ROISubtract(Grid grid) : ROIAggregation(grid) {
}

bool ROISubtract::inROI(tgt::vec3 p) const {
    bool in = false;
    if(children_.size() >= 1)
        in = children_[0]->inROI(p);

    if(in) {
        for(size_t i=1; i<children_.size(); i++)
            if(children_[i]->inROI(p))
                return false;
    }

    return in;
}

Geometry* ROISubtract::generateMesh() const {
    return 0;
}

Geometry* ROISubtract::generateMesh(tgt::plane /*pl*/) const {
    return 0;
}

tgt::Bounds ROISubtract::getBoundingBox() const {
    tgt::Bounds b;
    for(size_t i=0; i<children_.size(); i++)
        b.addVolume(children_[i]->getBoundingBox());
    return b;
}

std::vector<const ControlPoint*> ROISubtract::getControlPoints() const {
    //TODO
    return std::vector<const ControlPoint*>();
}

std::vector<const ControlPoint*> ROISubtract::getControlPoints(tgt::plane /*pl*/) const {
    //TODO
    return std::vector<const ControlPoint*>();
}

} // namespace
