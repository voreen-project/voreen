/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany,                        *
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

#include "flowindicatorrenderer.h"

#include "voreen/core/datastructures/callback/memberfunctioncallback.h"

#include "tgt/immediatemode/immediatemode.h"

namespace {

tgt::mat4 createTransformationMatrix(const tgt::vec3& position, const tgt::vec3& velocity) {

    tgt::vec3 tangent(tgt::normalize(velocity));

    tgt::vec3 temp(0.0f, 0.0f, 1.0f);
    if(1.0f - std::abs(tgt::dot(temp, tangent)) <= std::numeric_limits<float>::epsilon())
        temp = tgt::vec3(0.0f, 1.0f, 0.0f);

    tgt::vec3 binormal(tgt::normalize(tgt::cross(temp, tangent)));
    tgt::vec3 normal(tgt::normalize(tgt::cross(tangent, binormal)));

    return tgt::mat4(normal.x, binormal.x, tangent.x, position.x,
                     normal.y, binormal.y, tangent.y, position.y,
                     normal.z, binormal.z, tangent.z, position.z,
                     0.0f, 0.0f, 0.0f, 1.0f);
}

}

namespace voreen {

FlowIndicatorRenderer::FlowIndicatorRenderer()
    : GeometryRendererBase()
    , inport_(Port::INPORT, "parametrization", "Parametrization Input")
    , enable_("enable", "Enable", true)
    , geometry_(nullptr)
{
    addPort(inport_);

    addProperty(enable_);
}

void FlowIndicatorRenderer::initialize() {
    GeometryRendererBase::initialize();
    geometry_ = new GlMeshGeometryUInt16Simple();
    geometry_->setDiskGeometry(0.0f, 1.0f, 16);
}

void FlowIndicatorRenderer::deinitialize() {
    delete geometry_;
    geometry_ = nullptr;
    GeometryRendererBase::deinitialize();
}

tgt::Bounds FlowIndicatorRenderer::getBoundingBox() const {
    if (geometry_)
        return geometry_->getBoundingBox();

    return GeometryRendererBase::getBoundingBox();
}

void FlowIndicatorRenderer::process() {}

void FlowIndicatorRenderer::render() {
    if (!inport_.isReady() || !enable_.get())
        return;

    // TODO: set color dependend on in-/ outflow

    for(const FlowIndicator& indicator : inport_.getData()->getFlowIndicators()) {

        MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
        MatStack.pushMatrix();
        MatStack.multMatrix(createTransformationMatrix(indicator.center_, indicator.normal_));
        MatStack.scale(tgt::vec3(indicator.radius_));

        geometry_->render();

        MatStack.popMatrix();
    }
}

}

