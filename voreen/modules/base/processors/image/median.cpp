/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2021 University of Muenster, Germany,                        *
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

#include "median.h"

#include "tgt/textureunit.h"

using tgt::TextureUnit;

namespace voreen {

Median::Median()
    : ImageProcessorBypassable("image/median"),
      halfKernelDim_("halfKernelSize", "Half kernel size", 1, 1, 3),
      inport_(Port::INPORT, "image.inport", "Image Input"),
      outport_(Port::OUTPORT, "image.outport", "Image Output")
{
    addProperty(halfKernelDim_);

    addPort(inport_);
    addPort(outport_);
}

void Median::process() {
    if (!enableSwitch_.get()){
        bypass(&inport_, &outport_);
        return;
    }

    outport_.activateTarget();
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    TextureUnit colorUnit, depthUnit;
    inport_.bindTextures(colorUnit.getEnum(), depthUnit.getEnum());

    // initialize shader
    program_->activate();
    setGlobalShaderParameters(program_);
    program_->setUniform("colorTex_", colorUnit.getUnitNumber());
    program_->setUniform("depthTex_", depthUnit.getUnitNumber());
    inport_.setTextureParameters(program_, "textureParameters_");
    program_->setUniform("halfKernelDim_", halfKernelDim_.get());

    renderQuad();

    program_->deactivate();
    outport_.deactivateTarget();
    LGL_ERROR;
}

} // voreen namespace
