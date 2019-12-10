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

#include "meshfrustumclipping.h"

namespace voreen {

const std::string MeshFrustumClipping::loggerCat_("voreen.MeshFrustumClipping");

MeshFrustumClipping::MeshFrustumClipping()
    : Processor()
    , inport_(Port::INPORT, "geometry.geometry")
    , outport_(Port::OUTPORT, "geometry.clippedgeometry")
    , clipLeft_("clipLeft", "Clip against left frustum plane")
    , offsetLeft_("offsetLeft", "Left plane offset", 0.0f, -1.0f, 1.0f)
    , clipRight_("clipRight", "Clip against right frustum plane")
    , offsetRight_("offsetRight", "Right plane offset", 0.0f, -1.0f, 1.0f)
    , clipTop_("clipTop", "Clip against top frustum plane")
    , offsetTop_("offsetTop", "Top plane offset", 0.0f, -1.0f, 1.0f)
    , clipBottom_("clipBottom", "Clip against bottom frustum plane")
    , offsetBottom_("offsetBottom", "Bottom plane offset", 0.0f, -1.0f, 1.0f)
    , clipNear_("clipNear", "Clip against near frustum plane")
    , offsetNear_("offsetNear", "Near plane offset", 0.0f, -1.0f, 1.0f)
    , clipFar_("clipFar", "Clip against far frustum plane")
    , offsetFar_("offsetFar", "Far plane offset", 0.0f, -10.0f, 10.0f)
    , camera_("camera", "Camera", tgt::Camera(tgt::vec3(0.f, 0.f, 3.5f), tgt::vec3(0.f, 0.f, 0.f), tgt::vec3(0.f, 1.f, 0.f)))
{
    addPort(inport_);
    addPort(outport_);

    addProperty(clipLeft_);
    addProperty(offsetLeft_);

    addProperty(clipRight_);
    addProperty(offsetRight_);

    addProperty(clipTop_);
    addProperty(offsetTop_);

    addProperty(clipBottom_);
    addProperty(offsetBottom_);

    addProperty(clipNear_);
    addProperty(offsetNear_);

    addProperty(clipFar_);
    addProperty(offsetFar_);

    addProperty(camera_);
    camera_.setVisibleFlag(false);
}

MeshFrustumClipping::~MeshFrustumClipping()
{}

Processor* MeshFrustumClipping::create() const {
    return new MeshFrustumClipping();
}

void MeshFrustumClipping::process() {
    tgtAssert(inport_.getData(), "no geometry");

    const MeshListGeometry* inportGeometry = dynamic_cast<const MeshListGeometry*>(inport_.getData());
    if (!inportGeometry) {
        LWARNING("Input geometry type not supported, expecting MeshListGeometry.");
        outport_.setData(0);
        return;
    }

    geometry_ = *inportGeometry;

    tgt::Camera cam = camera_.get();
    cam.updateFrustum();
    const tgt::Frustum frust = cam.getFrustum();
    float epsilon = 0.0001f;

    if(clipLeft_.get())
        geometry_.clip(tgt::plane(frust.leftn(), offsetLeft_.get()+tgt::dot(frust.leftn(), frust.campos())), epsilon);

    if(clipRight_.get())
        geometry_.clip(tgt::plane(frust.rightn(), offsetRight_.get()+tgt::dot(frust.rightn(), frust.campos())), epsilon);

    if(clipTop_.get())
        geometry_.clip(tgt::plane(frust.topn(), offsetTop_.get()+tgt::dot(frust.topn(), frust.campos())), epsilon);

    if(clipBottom_.get())
        geometry_.clip(tgt::plane(frust.bottomn(), offsetBottom_.get()+tgt::dot(frust.bottomn(), frust.campos())), epsilon);

    if(clipNear_.get())
        geometry_.clip(tgt::plane(frust.nearn(), offsetNear_.get()+tgt::dot(frust.nearn(), frust.nearp())), epsilon);

    if(clipFar_.get())
        geometry_.clip(tgt::plane(frust.farn(), offsetFar_.get()+tgt::dot(frust.farn(), frust.farp())), epsilon);

    outport_.setData(&geometry_, false);
}

}  //namespace
