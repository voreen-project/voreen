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

#include "flowindicatorrenderer.h"

#include "tgt/immediatemode/immediatemode.h"

#include "../../utils/utils.h"

namespace voreen {

FlowIndicatorRenderer::FlowIndicatorRenderer()
    : GeometryRendererBase()
    , inport_(Port::INPORT, "parametrization", "Parametrization Input")
    , enable_("enable", "Enable", true)
    , colorMode_("colorMode", "Color Mode", false)
    , velocityBoundaryColor_("velocityBoundaryColor", "Velocity Boundary Color", tgt::vec4(0.0f, 1.0f, 0.0f, 0.8f))
    , pressureBoundaryColor_("pressureBoundaryColor", "Pressure Boundary Color", tgt::vec4(1.0f, 0.0f, 0.0f, 0.8f))
    , measureFluxColor_("measureFluxColor", "Flux Measure Color", tgt::vec4(0.0f, 0.0f, 1.0f, 0.8f))
{
    addPort(inport_);

    addProperty(enable_);
    addProperty(colorMode_);
    colorMode_.addOption("type", "Indicator Type");
    colorMode_.addOption("color", "Indicator Color");

    addProperty(velocityBoundaryColor_);
    addProperty(pressureBoundaryColor_);
    addProperty(measureFluxColor_);
}

void FlowIndicatorRenderer::initialize() {
    GeometryRendererBase::initialize();

    const int numSlices = 16;

    // Disk.
    diskGeometry_.reset(new GlMeshGeometryUInt16Simple());
    diskGeometry_->setCylinderGeometry(tgt::vec4::zero, 1.0f, 1.0f, 1.0f, numSlices, 2, true, true);

    // Cone.
    coneGeometry_.reset(new GlMeshGeometryUInt16Simple());
    coneGeometry_->setPrimitiveType(GL_TRIANGLE_FAN);
    coneGeometry_->addVertex(VertexBase(tgt::vec3(0.0f, 0.0f, 1.0f)));
    for(size_t i=0; i <= numSlices; ++i) {
        float s = std::sin(tgt::PIf*2*i/numSlices);
        float c = std::cos(tgt::PIf*2*i/numSlices);
        tgt::vec3 vertPos(s, c, 0.0f);
        coneGeometry_->addVertex(vertPos);
    }
}

void FlowIndicatorRenderer::deinitialize() {
    diskGeometry_.reset();
    coneGeometry_.reset();
    GeometryRendererBase::deinitialize();
}

tgt::Bounds FlowIndicatorRenderer::getBoundingBox() const {

    tgt::Bounds bounds;
    if (inport_.hasData()) {
        for (const FlowIndicator& indicator : inport_.getData()->getFlowIndicators()) {
            bounds.addPoint(indicator.center_ - indicator.radius_);
            bounds.addPoint(indicator.center_ + indicator.radius_);
        }
    }

    return bounds;
}

void FlowIndicatorRenderer::renderDisk(const FlowIndicator& indicator) const {
    MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
    MatStack.pushMatrix();
    MatStack.multMatrix(utils::createTransformationMatrix(indicator.center_ - indicator.normal_ * indicator.length_ * 0.5f , indicator.normal_));
    MatStack.scale(indicator.radius_, indicator.radius_, indicator.length_);
    diskGeometry_->render();
    MatStack.popMatrix();
}

void FlowIndicatorRenderer::renderCone(const FlowIndicator& indicator) const {
    MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
    MatStack.pushMatrix();
    MatStack.multMatrix(utils::createTransformationMatrix(indicator.center_ + indicator.normal_ * indicator.length_ * 0.5f, indicator.normal_));
    MatStack.scale(indicator.radius_, indicator.radius_, indicator.radius_);
    coneGeometry_->render();
    MatStack.popMatrix();
}

tgt::vec4 FlowIndicatorRenderer::getColor(const FlowIndicator& indicator) const {
    if(colorMode_.get() == "type") {
        switch (indicator.type_) {
        case FIT_VELOCITY:
            return velocityBoundaryColor_.get();
        case FIT_PRESSURE:
            return pressureBoundaryColor_.get();
        case FIT_MEASURE:
            return measureFluxColor_.get();
        default:
            return tgt::vec4::zero;
        }
    }
    else if(colorMode_.get() == "color") {
        return indicator.color_;
    }
    else {
        LERROR("Unknown color mode");
        return tgt::vec4::zero;
    }
}

void FlowIndicatorRenderer::render() {
    if (!inport_.isReady() || !enable_.get())
        return;

    for(const FlowIndicator& indicator : inport_.getData()->getFlowIndicators()) {
        IMode.color(getColor(indicator));
        if(indicator.type_ == FIT_VELOCITY) {
            renderDisk(indicator);
            renderCone(indicator);
        }
        else if(indicator.type_ == FIT_PRESSURE) {
            renderDisk(indicator);
            renderCone(indicator);
        }
        else if(indicator.type_ == FIT_MEASURE) {
            renderDisk(indicator);
        }
    }

    IMode.color(tgt::vec4::one);
}

}

