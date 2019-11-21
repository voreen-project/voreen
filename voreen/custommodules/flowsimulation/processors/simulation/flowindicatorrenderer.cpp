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

#include "tgt/immediatemode/immediatemode.h"

#include "../../utils/utils.h"

namespace voreen {

FlowIndicatorRenderer::FlowIndicatorRenderer()
    : GeometryRendererBase()
    , inport_(Port::INPORT, "parametrization", "Parametrization Input")
    , enable_("enable", "Enable", true)
    , flowGeneratorColor_("flowGeneratorColor", "Flow Generator Color", tgt::vec4(0.0f, 1.0f, 0.0f, 0.8f))
    , pressureBoundaryColor_("pressureBoundaryColor", "Pressure Boundary Color", tgt::vec4(1.0f, 0.0f, 0.0f, 0.8f))
    , measureFluxColor_("measureFluxColor", "Flux Measure Color", tgt::vec4(0.0f, 0.0f, 1.0f, 0.8f))
{
    addPort(inport_);

    addProperty(enable_);
    addProperty(flowGeneratorColor_);
    addProperty(pressureBoundaryColor_);
    addProperty(measureFluxColor_);
}

void FlowIndicatorRenderer::initialize() {
    GeometryRendererBase::initialize();

    const int numSlices = 16;

    // Disk.
    diskGeometry_.reset(new GlMeshGeometryUInt16Simple());
    diskGeometry_->setDiskGeometry(0.0f, 1.0f, numSlices);

    // Cone.
    coneGeometry_.reset(new GlMeshGeometryUInt16Simple());
    coneGeometry_->setPrimitiveType(GL_TRIANGLE_FAN);
    coneGeometry_->addVertex(VertexBase(tgt::vec3(0.0f, 0.0f, 1.0f)));
    for(size_t i=0; i <= numSlices; ++i) {
        float s = std::sin(tgt::PIf*2*i/numSlices);
        float c = std::cos(tgt::PIf*2*i/numSlices);
        tgt::vec3 vertPos(s, c, 0);
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

void FlowIndicatorRenderer::process() {}

void FlowIndicatorRenderer::render() {
    if (!inport_.isReady() || !enable_.get())
        return;

    for(const FlowIndicator& indicator : inport_.getData()->getFlowIndicators()) {

        MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
        MatStack.pushMatrix();
        MatStack.multMatrix(utils::createTransformationMatrix(indicator.center_, indicator.normal_));
        MatStack.scale(tgt::vec3(indicator.radius_));

        if(indicator.type_ == FIT_GENERATOR) {
            IMode.color(flowGeneratorColor_.get());
            coneGeometry_->render();
        }
        else if(indicator.type_ == FIT_PRESSURE) {
            IMode.color(pressureBoundaryColor_.get());
            diskGeometry_->render();
        }
        else if(indicator.type_ == FIT_MEASURE) {
            IMode.color(measureFluxColor_.get());
            diskGeometry_->render();
        }

        MatStack.popMatrix();
    }

    IMode.color(tgt::vec4::one);
}

}

