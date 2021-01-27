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

#include "voreen/core/datastructures/roi/roiaggregation.h"
#include "voreen/core/datastructures/geometry/meshlistgeometry.h"

namespace voreen {

const std::string ROIAggregation::loggerCat_ = "voreen.ROIAggregation";

ROIAggregation::ROIAggregation(Grid grid) : ROIBase(grid) {
}

ROIAggregation::~ROIAggregation() {
    while(!children_.empty()) {
        delete children_.back();
        children_.pop_back();
    }
}

std::vector<ROIBase*> ROIAggregation::getChildren() {
    return children_;
}

std::vector<const ROIBase*> ROIAggregation::getChildren() const {
    std::vector<const ROIBase*> ret;

    for(size_t i=0; i<children_.size(); i++) {
        ret.push_back(children_[i]);
    }
    return ret;
}

void ROIAggregation::addChild(ROIBase* c) {
    children_.push_back(c);
}

void ROIAggregation::removeChild(ROIBase* c) {
    for(size_t i=0; i<children_.size(); i++) {
        if(children_[i] == c) {
            children_.erase(children_.begin()+i);
            delete c;
            return;
        }
    }
}

void ROIAggregation::deregisterChild(ROIBase* c) {
    for(size_t i=0; i<children_.size(); i++) {
        if(children_[i] == c) {
            children_.erase(children_.begin()+i);
            return;
        }
    }
}

void ROIAggregation::serialize(Serializer& s) const {
    s.serialize("Children", children_, "ROI");
}

void ROIAggregation::deserialize(Deserializer& s) {
    s.deserialize("Children", children_, "ROI");
}

} // namespace
