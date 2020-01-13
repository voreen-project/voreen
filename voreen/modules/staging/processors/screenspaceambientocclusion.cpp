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

#include "screenspaceambientocclusion.h"

#include "voreen/core/utils/voreenqualitymode.h"
#include "tgt/textureunit.h"
#include <time.h>

using tgt::TextureUnit;

namespace voreen {

ScreenSpaceAmbientOcclusion::ScreenSpaceAmbientOcclusion()
    : ImageProcessorBypassable("ssao", true)
    , inportOrig_(Port::INPORT, "image.orig", "Image Input")
    , outport_(Port::OUTPORT, "image.outport", "Image Output", true)
    , camera_("camera", "Camera", tgt::Camera(tgt::vec3(0.f, 0.f, 3.5f), tgt::vec3(0.f, 0.f, 0.f), tgt::vec3(0.f, 1.f, 0.f)))
    , radius_("radius", "Radius", 3, 1, 20)
    , occFac_("aofac", "Occlusion multiplicator", 1.f, 0.f, 10.f)
    , colFac_("gradfac", "Color bleeding multiplicator", 3.f, 0.f, 10.f)
{
    addPort(inportOrig_);
    addPort(outport_);

    addProperty(camera_);
    addProperty(radius_);
    addProperty(occFac_);
    addProperty(colFac_);
}

Processor* ScreenSpaceAmbientOcclusion::create() const {
    return new ScreenSpaceAmbientOcclusion();
}

bool ScreenSpaceAmbientOcclusion::isReady() const {
    return inportOrig_.isReady() && outport_.isReady();
}

void ScreenSpaceAmbientOcclusion::process() {

    if (!enableSwitch_.get()){
        bypass(&inportOrig_, &outport_);
        return;
    }

    // compile program if needed
    if (getInvalidationLevel() >= Processor::INVALID_PROGRAM)
        compile();
    LGL_ERROR;

    if(!program_ || !program_->isLinked())
        return;

    outport_.activateTarget();
    outport_.clearTarget();

    TextureUnit colorUnit, depthUnit, normUnit;
    inportOrig_.bindTextures(colorUnit.getEnum(), depthUnit.getEnum());

    // initialize shader
    program_->activate();
    setGlobalShaderParameters(program_, &camera_.get());

    program_->setUniform("colorTex_", colorUnit.getUnitNumber());
    program_->setUniform("depthTex_", depthUnit.getUnitNumber());

    int kernelRadius = radius_.get();
    if(QualityMode.isInteractionMode() && interactionAdapt_.get())
        kernelRadius = std::max(1, (int)(kernelRadius * (1.f / interactionFactor_.get())));
    program_->setUniform("numIts_", kernelRadius);
    program_->setUniform("occFac_", occFac_.get());
    program_->setUniform("colFac_", colFac_.get());

    inportOrig_.setTextureParameters(program_, "colorParams_");

    glDepthFunc(GL_ALWAYS);
    renderQuad();
    glDepthFunc(GL_LESS);

    program_->deactivate();
    outport_.deactivateTarget();
    TextureUnit::setZeroUnit();
    LGL_ERROR;
}

std::string ScreenSpaceAmbientOcclusion::generateHeader(const tgt::GpuCapabilities::GlVersion* version = 0) {
    std::string header = ImageProcessor::generateHeader(version);
    return header;
}

} // voreen namespace
