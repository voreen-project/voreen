/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany.                        *
 * Visualization and Computer Graphics Group <http://viscg.uni-muenster.de>        *
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

#ifndef VRN_ROINORMALIZEDGEOMETRY_H
#define VRN_ROINORMALIZEDGEOMETRY_H

#include "voreen/core/datastructures/roi/roisingle.h"
#include "voreen/core/properties/optionproperty.h"

namespace voreen {

class VRN_CORE_API ROINormalizedGeometry : public ROISingle {
public:
    ROINormalizedGeometry(Grid grid, tgt::vec3 center, tgt::vec3 dimensions);
    ROINormalizedGeometry();
    virtual bool inROI(tgt::vec3 p) const;

    tgt::mat4 getPhysicalToNormalizedMatrix() const;
    tgt::mat4 getNormalizedToPhysicalMatrix() const;

    virtual Geometry* generateMesh() const;
    virtual Geometry* generateMesh(tgt::plane pl) const;

    virtual tgt::Bounds getBoundingBox() const;

    virtual std::vector<const ControlPoint*> getControlPoints() const;
    virtual std::vector<const ControlPoint*> getControlPoints(tgt::plane pl) const;

    virtual bool moveControlPoint(const ControlPoint* cp, tgt::vec3 to);

    tgt::vec3 getCenter() const { return center_.get(); }
    void setCenter(tgt::vec3 c) { center_.set(c); }
    tgt::vec3 getDimensions() const { return dimensions_.get(); }
    void setDimensions(tgt::vec3 d) { dimensions_.set(d); }
    int getPrimaryDir() const { return primaryDir_.getValue(); }
    void setPrimaryDir(int pd) { primaryDir_.selectByValue(pd); }
    int getResizingDOF() const { return resizingDOF_.getValue(); }
    void setResizingDOF(int dof) { return resizingDOF_.selectByValue(dof); }

    tgt::vec3 getLLF() const { return center_.get() - dimensions_.get() * 0.5f; }
    tgt::vec3 getURB() const { return center_.get() + dimensions_.get() * 0.5f; }
    TriangleMeshGeometrySimple* getPlanePolygon(tgt::plane pl) const;

    // Methods to implement:
    virtual bool inROINormalized(tgt::vec3 p) const = 0;
    virtual Geometry* generateNormalizedMesh() const = 0;
    /**
     * Creates a (intersection) mesh in the given plane, usually for display in slice renderings.
     * @param pl Plane in normalized coordinates.
     */
    virtual Geometry* generateNormalizedMesh(tgt::plane pl) const = 0;
private:
    FloatVec3Property center_;
    FloatVec3Property dimensions_;
    IntOptionProperty primaryDir_;
    IntOptionProperty resizingDOF_;

    static const std::string loggerCat_;
};

} //namespace

#endif
