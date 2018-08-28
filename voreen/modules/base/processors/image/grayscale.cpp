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

#include "grayscale.h"

#include "tgt/textureunit.h"

using tgt::TextureUnit;

namespace voreen {

Grayscale::Grayscale()
    : ImageProcessorBypassable("image/grayscale"), // loads fragment shader pp_grayscale.frag
      inport_(Port::INPORT, "inport", "Image Input"),
      outport_(Port::OUTPORT, "outport", "Image Output"),
      saturation_("saturation", "Saturation", 0.f)
{
    // register ports and properties
    addPort(inport_);
    addPort(outport_);

    addProperty(saturation_);
}

Processor* Grayscale::create() const {
    return new Grayscale();
}

void Grayscale::process() {

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
    program_->setUniform("saturation_", saturation_.get());

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
