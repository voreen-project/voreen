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

#include "voreen/core/datastructures/roi/roinormalizedgeometry.h"
#include "voreen/core/datastructures/geometry/trianglemeshgeometry.h"

#include "voreen/core/io/serialization/serialization.h"

namespace voreen {

using tgt::vec3;
using tgt::vec4;
using tgt::mat4;

const std::string ROINormalizedGeometry::loggerCat_ = "voreen.ROINormalizedGeometry";

ROINormalizedGeometry::ROINormalizedGeometry() : ROISingle()
    , center_("center", "Center", vec3(0.0f), vec3(-999.9f), vec3(999.9f), ROI_CHANGE)
    , dimensions_("dimensions", "Dimensions", vec3(1.0f), vec3(-999.9f), vec3(999.9f), ROI_CHANGE)
    , primaryDir_("primaryDir", "Primary Direction")
    , resizingDOF_("resizingBehavior", "Resizing Behavior")
{
    addProperty(center_);
    addProperty(dimensions_);

    primaryDir_.addOption("none", "None", -1);
    primaryDir_.addOption("x", "X", 0);
    primaryDir_.addOption("y", "Y", 1);
    primaryDir_.addOption("z", "Z", 2);
    primaryDir_.select("none");
    addProperty(primaryDir_);

    resizingDOF_.addOption("3d", "No Restrictions", 3);
    resizingDOF_.addOption("2d", "1+2 Locked", 2);
    resizingDOF_.addOption("1d", "All Locked", 1);
    resizingDOF_.select("3d");
    addProperty(resizingDOF_);
}

ROINormalizedGeometry::ROINormalizedGeometry(Grid grid, tgt::vec3 center, tgt::vec3 dimensions) : ROISingle(grid)
    , center_("center", "Center", vec3(center), vec3(-999.9f), vec3(999.9f), ROI_CHANGE)
    , dimensions_("dimensions", "Dimensions", vec3(dimensions), vec3(-999.9f), vec3(999.9f), ROI_CHANGE)
    , primaryDir_("primaryDir", "Primary Direction")
    , resizingDOF_("resizingDOF", "Resizing Behavior")
{
    addProperty(center_);
    addProperty(dimensions_);

    primaryDir_.addOption("none", "None", -1);
    primaryDir_.addOption("x", "X", 0);
    primaryDir_.addOption("y", "Y", 1);
    primaryDir_.addOption("z", "Z", 2);
    primaryDir_.select("none");
    addProperty(primaryDir_);

    resizingDOF_.addOption("3d", "No Restrictions", 3);
    resizingDOF_.addOption("2d", "1+2 Locked", 2);
    resizingDOF_.addOption("1d", "All Locked", 1);
    resizingDOF_.select("3d");
    addProperty(resizingDOF_);
}

bool ROINormalizedGeometry::inROI(tgt::vec3 p) const {
    vec3 pNorm = getPhysicalToNormalizedMatrix() * p;
    if(inROINormalized(pNorm))
        return true;
    else
        return false;
}

tgt::mat4 ROINormalizedGeometry::getPhysicalToNormalizedMatrix() const {
    return mat4::createScale(vec3(2.0f/dimensions_.get().x, 2.0f/dimensions_.get().y, 2.0f/dimensions_.get().z)) * mat4::createTranslation(-center_.get());
}

tgt::mat4 ROINormalizedGeometry::getNormalizedToPhysicalMatrix() const {
    return mat4::createTranslation(center_.get()) * mat4::createScale(dimensions_.get() * 0.5f);
}

Geometry* ROINormalizedGeometry::generateMesh() const {
    Geometry* geometry = generateNormalizedMesh();
    geometry->transform(getNormalizedToPhysicalMatrix());
    geometry->transform(getGrid().getPhysicalToWorldMatrix());
    return geometry;
}

TriangleMeshGeometrySimple* ROINormalizedGeometry::getPlanePolygon(tgt::plane pl) const {
    // Transform plane to grid coordinates:
    mat4 m = getGrid().getWorldToPhysicalMatrix();
    const tgt::plane physicalPlane = pl.transform(m);

    TriangleMeshGeometrySimple* planeFace = createQuad(physicalPlane);

    tgt::plane xm(-1.0f, 0.0f, 0.0f, -getLLF().x);
    tgt::plane xp(1.0f, 0.0f, 0.0f, getURB().x);
    tgt::plane ym(0.0f, -1.0f, 0.0f, -getLLF().y);
    tgt::plane yp(0.0f, 1.0f, 0.0f, getURB().y);
    tgt::plane zm(0.0f, 0.0f, -1.0f, -getLLF().z);
    tgt::plane zp(0.0f, 0.0f, 1.0f, getURB().z);

    planeFace->clip(xm);
    planeFace->clip(xp);
    planeFace->clip(ym);
    planeFace->clip(yp);
    planeFace->clip(zm);
    planeFace->clip(zp);

    return planeFace;
}

Geometry* ROINormalizedGeometry::generateMesh(tgt::plane pl) const {
    const tgt::plane physicalPlane = pl.transform(getGrid().getWorldToPhysicalMatrix());
    const tgt::plane normalizedPlane = physicalPlane.transform(getPhysicalToNormalizedMatrix());

    Geometry* mlg = generateNormalizedMesh(normalizedPlane);

    mlg->transform(getNormalizedToPhysicalMatrix());
    mlg->transform(getGrid().getPhysicalToWorldMatrix());

    return mlg;
}

tgt::Bounds ROINormalizedGeometry::getBoundingBox() const {
    return tgt::Bounds(getLLF(), getURB());
}

std::vector<const ControlPoint*> ROINormalizedGeometry::getControlPoints() const {
    std::vector<const ControlPoint*> cps;
    cps.push_back(new ControlPoint(this, ControlPoint::CENTER, getGrid().getPhysicalToWorldMatrix() * center_.get()));
    return cps;
}

std::vector<const ControlPoint*> ROINormalizedGeometry::getControlPoints(tgt::plane pl) const {
    std::vector<const ControlPoint*> cps;

    TriangleMeshGeometrySimple* planeFace = getPlanePolygon(pl);
    if(!planeFace->isEmpty()) {
        std::vector<vec4> planes;

        planes.push_back(vec4(1.0f, 0.0f, 0.0f, getURB().x));
        planes.push_back(vec4(-1.0f, 0.0f, 0.0f, -getLLF().x));
        planes.push_back(vec4(0.0f, 1.0f, 0.0f, getURB().y));
        planes.push_back(vec4(0.0f, -1.0f, 0.0f, -getLLF().y));
        planes.push_back(vec4(0.0f, 0.0f, 1.0f, getURB().z));
        planes.push_back(vec4(0.0f, 0.0f, -1.0f, -getLLF().z));

        vec3 com(0.0f); //Center of mass
        //TODO
        //for(size_t i=0; i<planeFace.getVertexCount(); i++) {
            //com += planeFace.getVertex(i).getCoords();

            //size_t prev;
            //if(i == 0)
                //prev = planeFace.getVertexCount() - 1;
            //else
                //prev = i - 1;

            //vec3 a = planeFace.getVertex(prev).getCoords();
            //vec3 b = planeFace.getVertex(i).getCoords();
            //tgt::vec3 p = 0.5f * (a + b);

            //// Find plane id:
            //int closest = 0;
            //for(size_t j=1; j<planes.size(); j++) {
                //float d = fabs(dot(p, planes[j].xyz()) - planes[j].w);
                //if(d < fabs(dot(p, planes[closest].xyz()) - planes[closest].w))
                    //closest = static_cast<int>(j);
            //}

            //cps.push_back(new ControlPoint(this, closest+1, getGrid().getPhysicalToWorldMatrix() * p));
        //}

        //com /= (float) planeFace.getVertexCount();
        ////Center:
        //cps.push_back(new ControlPoint(this, ControlPoint::CENTER, getGrid().getPhysicalToWorldMatrix() * com));
    }

    return cps;
}

bool ROINormalizedGeometry::moveControlPoint(const ControlPoint* cp, tgt::vec3 offset) {
    if(cp->getROI() != this)
        return false;

    offset = getGrid().getWorldToPhysicalMatrix().getRotationalPart() * offset;

    vec3 dir(0.0f);
    switch(cp->getID()) {
        case ControlPoint::CENTER:
            center_.set(center_.get() + offset);
            return true;
        case ControlPoint::X_PLUS:
            dir = vec3(1.0f, 0.0f, 0.0f);
            break;
        case ControlPoint::X_MINUS:
            dir = vec3(-1.0f, 0.0f, 0.0f);
            break;
        case ControlPoint::Y_PLUS:
            dir = vec3(0.0f, 1.0f, 0.0f);
            break;
        case ControlPoint::Y_MINUS:
            dir = vec3(0.0f, -1.0f, 0.0f);
            break;
        case ControlPoint::Z_PLUS:
            dir = vec3(0.0f, 0.0f, 1.0f);
            break;
        case ControlPoint::Z_MINUS:
            dir = vec3(0.0f, 0.0f, -1.0f);
            break;
    }
    if(dir != vec3(0.0f)) {
        bool neg = hadd(dir) < 0.0f;

        dir *= dot(dir, offset);
        if(!neg)
            dimensions_.set(dimensions_.get() + dir);
        else
            dimensions_.set(dimensions_.get() - dir);

        dir *= 0.5f;
        center_.set(center_.get() + dir);
        return true;
    }

    return false;
}

} // namespace
