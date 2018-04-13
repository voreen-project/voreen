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

#include "depthdarkening.h"

#include "tgt/textureunit.h"

using tgt::TextureUnit;

namespace voreen {

DepthDarkening::DepthDarkening()
    : ImageProcessorBypassable("image/depthdarkening")
    , sigma_("sigma", "Sigma", 2.0f, 0.1f, 10.0f)
    , lambda_("lambda", "Lambda", 10.0f, 0.0f, 20.0f)
    , inport_(Port::INPORT, "image.inport", "Image Input")
    , outport_(Port::OUTPORT, "image.outport", "Image Output")
    , privatePort_(Port::OUTPORT, "image.privateport")
{
    addProperty(sigma_);
    addProperty(lambda_);

    addPort(inport_);
    addPort(outport_);
    addPrivateRenderPort(privatePort_);
}

Processor* DepthDarkening::create() const {
    return new DepthDarkening();
}

void DepthDarkening::process() {

    if (!enableSwitch_.get()){
        bypass(&inport_, &outport_);
        return;
    }

    privatePort_.activateTarget();
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    //get min and max depth values
    tgt::vec2 depthRange = computeDepthRange(&inport_);

    // since the blurring is implemented as a separable filter,
    // two rendering passes are needed
    // first horizontal pass
    TextureUnit colorUnit, depthUnit0, depthUnit1;
    inport_.bindTextures(colorUnit.getEnum(), depthUnit0.getEnum());
    inport_.bindDepthTexture(depthUnit1.getEnum());
    // initialize shader
    program_->activate();
    setGlobalShaderParameters(program_);
    program_->setUniform("colorTex_", colorUnit.getUnitNumber());
    program_->setUniform("depthTex0_", depthUnit0.getUnitNumber());
    program_->setUniform("depthTex1_", depthUnit1.getUnitNumber());
    inport_.setTextureParameters(program_, "texParams_");
    program_->setUniform("sigma_", sigma_.get());
    program_->setUniform("lambda_", lambda_.get());
    program_->setUniform("minDepth_", depthRange.x);
    program_->setUniform("maxDepth_", depthRange.y);
    program_->setUniform("dir_", tgt::vec2(1.0,0.0));

    renderQuad();

    program_->deactivate();
    LGL_ERROR;

    // second vertical pass
    outport_.activateTarget();
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    inport_.bindColorTexture(colorUnit.getEnum());
    privatePort_.bindDepthTexture(depthUnit0.getEnum());

    // initialize shader
    program_->activate();
    setGlobalShaderParameters(program_);
    program_->setUniform("colorTex_", colorUnit.getUnitNumber());
    program_->setUniform("depthTex0_", depthUnit0.getUnitNumber());
    program_->setUniform("depthTex1_", depthUnit1.getUnitNumber());
    inport_.setTextureParameters(program_, "texParams_");
    program_->setUniform("sigma_", sigma_.get());
    program_->setUniform("lambda_", lambda_.get());
    program_->setUniform("minDepth_", depthRange.x);
    program_->setUniform("maxDepth_", depthRange.y);
    program_->setUniform("dir_", tgt::vec2(0.0,1.0));

    renderQuad();

    program_->deactivate();
    outport_.deactivateTarget();
    LGL_ERROR;
}

} // voreen namespace
