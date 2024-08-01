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

#include "samplingpositiontransformation.h"

#include "tgt/textureunit.h"

using tgt::TextureUnit;

namespace voreen {

SamplingPositionTransformation::SamplingPositionTransformation()
    : ImageProcessor("samplingpositiontransformation"),
      mode_("mode_", "Mode"),
      inport_(Port::INPORT, "inport", "FHP Input"),
      outport_(Port::OUTPORT, "outport", "FHP Transformed", true, INVALID_RESULT, RenderPort::RENDERSIZE_DEFAULT, GL_RGBA16F),
      volPort_(Port::INPORT, "volport", "Volume Input")
{
    addPort(inport_);
    addPort(outport_);
    addPort(volPort_);

    mode_.addOption("textureToWorld", "Texture to World");
    mode_.addOption("worldToTexture", "World to Texture");
    mode_.select("textureToWorld");
    addProperty(mode_);
}

SamplingPositionTransformation::~SamplingPositionTransformation() {}

void SamplingPositionTransformation::process() {
    outport_.activateTarget();
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    TextureUnit colorUnit, depthUnit;
    inport_.bindTextures(colorUnit.getEnum(), depthUnit.getEnum());
    LGL_ERROR;

    // initialize shader
    program_->activate();
    setGlobalShaderParameters(program_);
    inport_.setTextureParameters(program_, "texParams_");
    program_->setUniform("colorTex_", colorUnit.getUnitNumber());
    program_->setUniform("depthTex_", depthUnit.getUnitNumber());

    if(mode_.get() == "worldToTexture")
        program_->setUniform("transformation_", volPort_.getData()->getWorldToTextureMatrix());
    else if(mode_.get() == "textureToWorld")
        program_->setUniform("transformation_", volPort_.getData()->getTextureToWorldMatrix());

    renderQuad();

    program_->deactivate();
    LGL_ERROR;
    outport_.deactivateTarget();
}

} // voreen namespace
