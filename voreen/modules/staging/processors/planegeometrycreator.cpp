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

#include "planegeometrycreator.h"

#include "voreen/core/datastructures/geometry/trianglemeshgeometry.h"
#include "tgt/glmath.h"

namespace voreen {

const std::string PlaneGeometryCreator::loggerCat_("voreen.PlaneGeometryCreator");

PlaneGeometryCreator::PlaneGeometryCreator()
    : Processor()
    //properties
    , adaptFromVolume_("adaptFromVolume", "Adapt Range from Inport", true)
        //settings
    , normal_("planeNormal", "Plane Normal", tgt::vec3(0.f, 1.f, 0.f), tgt::vec3(-1.f), tgt::vec3(1.f))
    , position_("planePosition", "Plane Position", tgt::vec3(0.0f), tgt::vec3(-1000000.f), tgt::vec3(1000000.f))
    , size_("size", "Plane Size", 10.0f, 0.01f, 1000000.f)
    , reset_("reset","Reset Values")
    , volInport_(Port::INPORT,"volInport","Volume Inport")
    , geomOutport_(Port::OUTPORT, "geometry", "Geometry Output")
{
    addProperty(adaptFromVolume_);
    adaptFromVolume_.onChange(MemberFunctionCallback<PlaneGeometryCreator>(this, &PlaneGeometryCreator::onVolumePortChangeCallback));

    addProperty(normal_);
        normal_.setGroupID("plane");
    addProperty(position_);
        position_.setGroupID("plane");
        position_.setNumDecimals(4);
    addProperty(size_);
        size_.setGroupID("plane");
        size_.setNumDecimals(4);
    addProperty(reset_);
        reset_.setGroupID("plane");
        reset_.onChange(MemberFunctionCallback<PlaneGeometryCreator>(this, &PlaneGeometryCreator::onResetClickedCallback));
    setPropertyGroupGuiName("plane","Plane Settings");

    addPort(volInport_);
        volInport_.onNewData(MemberFunctionCallback<PlaneGeometryCreator>(this, &PlaneGeometryCreator::onVolumePortChangeCallback));
    addPort(geomOutport_);
}

PlaneGeometryCreator::~PlaneGeometryCreator() {
}

bool PlaneGeometryCreator::isReady() const {
    return geomOutport_.isReady();
}

void PlaneGeometryCreator::process() {

    if(normal_.get() == tgt::vec3::zero) {
        LERROR ("Plane normal of (0,0,0) is not allowed! No output gererated!");
        geomOutport_.setData(nullptr);
        return;
    }

    tgt::vec3 n = normalize(normal_.get());

    tgt::vec3 temp = tgt::vec3(1.0, 0.0, 0.0);
    if(std::abs(dot(temp, n)) > 0.9)
        temp = tgt::vec3(0.0, 1.0, 0.0);

    tgt::vec3 inPlaneA = normalize(cross(n, temp)) * 0.5f * size_.get();
    tgt::vec3 inPlaneB = normalize(cross(n, inPlaneA)) * 0.5f * size_.get();

    tgt::vec3 base = position_.get();

    TriangleMeshGeometryNormal* plane = new TriangleMeshGeometryNormal();
    plane->addQuad(VertexNormal(base + inPlaneA + inPlaneB, base + inPlaneA + inPlaneB),
                   VertexNormal(base - inPlaneA + inPlaneB, base - inPlaneA + inPlaneB),
                   VertexNormal(base - inPlaneA - inPlaneB, base - inPlaneA - inPlaneB),
                   VertexNormal(base + inPlaneA - inPlaneB, base + inPlaneA - inPlaneB));

    geomOutport_.setData(plane);
}

void PlaneGeometryCreator::onVolumePortChangeCallback() {
    //should not be needed
    if(volInport_.hasData() && adaptFromVolume_.get()) {
        const VolumeBase* volume = volInport_.getData();

        float diagonal = tgt::length(volume->getCubeSize());
        tgt::vec3 center = volume->getOffset() + volume->getCubeSize() / 2.f;

        //set to 4 decimals
        position_.setMaxValue(((center + volume->getCubeSize())*10000.f) / 10000.f);
        position_.setMinValue(((center - volume->getCubeSize())*10000.f) / 10000.f);
        size_.setMaxValue(std::floor((2*diagonal)*10000.f) / 10000.f + 0.5f);
        size_.set(size_.getMaxValue());
    } else {
        position_.setMaxValue(tgt::vec3( 1000000.f));
        position_.setMinValue(tgt::vec3(-1000000.f));
        size_.setMaxValue(1000000.f);
    }
}

void PlaneGeometryCreator::onResetClickedCallback() {
    normal_.set(tgt::vec3(0.f,1.f,0.f));

    if(tgt::hand(tgt::greaterThanEqual(position_.getMaxValue(),tgt::vec3(1000000.f - 1.f)))) //default (-"epsilon")
        position_.set(tgt::vec3::zero);
    else
        position_.set((position_.getMinValue() + position_.getMaxValue()) / 2.f);

    if(size_.getMaxValue() >= 1000000.f - 1.f) //default ( -"epsilon")
        size_.set(10.f);
    else
        size_.set(size_.getMaxValue() / 2.f);
}

} // namespace voreen
