/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2024 University of Muenster, Germany,                        *
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

#ifndef VRN_ROIUNION_H
#define VRN_ROIUNION_H

#include "roiaggregation.h"

namespace voreen {

class VRN_CORE_API ROIUnion : public ROIAggregation {
public:
    ROIUnion(Grid grid);
    ROIUnion() : ROIAggregation() {}
    virtual std::string getClassName() const { return "ROIUnion"; }
    virtual ROIBase* create() const { return new ROIUnion(); }

    virtual bool inROI(tgt::vec3 p) const;

    virtual Geometry* generateMesh() const;
    virtual Geometry* generateMesh(tgt::plane pl) const;

    virtual tgt::Bounds getBoundingBox() const;

    virtual std::vector<const ControlPoint*> getControlPoints() const;
    virtual std::vector<const ControlPoint*> getControlPoints(tgt::plane pl) const;

private:
    static const std::string loggerCat_;
};

} //namespace

#endif
