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

#include "fxaa.h"

#include "tgt/textureunit.h"

using tgt::TextureUnit;

namespace voreen {

FXAA::FXAA()
    : ImageProcessorBypassable("image/fxaa"), // loads fragment shader fxaa.frag
      inport_(Port::INPORT, "inport", "Image Input"),
      outport_(Port::OUTPORT, "outport", "Image Output"),
      spanMax_("spanMax", "Maximum Edge Detection Pixel Span", 8.f, 1.f, 16.f),
      reduceMul_("reduceMul", "Reduction Multiplier: 1.0 / ", 8.f, 1.f, 128.f)
{
    // register ports and properties
    addPort(inport_);
    addPort(outport_);

    addProperty(spanMax_);
    addProperty(reduceMul_);
}

Processor* FXAA::create() const {
    return new FXAA();
}

void FXAA::process() {

    if (!enableSwitch_.get()) {
        bypass(&inport_, &outport_);
        return;
    }

    // activate and clear output render target
    outport_.activateTarget();
    outport_.clearTarget();

    // bind input rendering to texture units
    TextureUnit colorUnit, depthUnit;
    inport_.bindTextures(colorUnit.getEnum(), depthUnit.getEnum());

    // activate shader and set uniforms
    program_->activate();
    setGlobalShaderParameters(program_);
    inport_.setTextureParameters(program_, "textureParameters_");
    program_->setUniform("colorTex_", colorUnit.getUnitNumber());
    program_->setUniform("depthTex_", depthUnit.getUnitNumber());
    program_->setUniform("spanMax_", spanMax_.get());
    program_->setUniform("reduceMul_", reduceMul_.get());

    // render screen aligned quad
    renderQuad();

    // cleanup
    program_->deactivate();
    outport_.deactivateTarget();
    TextureUnit::setZeroUnit();

    // check for OpenGL errors
    LGL_ERROR;
}

} // namespace
