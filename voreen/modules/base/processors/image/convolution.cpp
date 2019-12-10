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

#include "convolution.h"

namespace voreen {

Convolution::Convolution()
    : ImageProcessorBypassable("image/convolution"),
      filterSize_("filtersize", "Filter Size (NxN)", 7),
      inport_(Port::INPORT, "inport", "Image Input"),
      filterPort_(Port::INPORT, "filterport", "Filter Input"),
      outport_(Port::OUTPORT, "outport", "Image Output")
{
    // register properties and ports:
    addProperty(filterSize_);

    addPort(inport_);
    addPort(filterPort_);
    addPort(outport_);
}

Processor* Convolution::create() const {
    return new Convolution();
}

void Convolution::process() {
    if (!enableSwitch_.get()) {
        bypass(&inport_, &outport_);
        return;
    }

    // activate and clear output render target
    outport_.activateTarget();
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // bind result from previous processor
    inport_.bindTextures(GL_TEXTURE0, GL_TEXTURE1);
    filterPort_.bindColorTexture(GL_TEXTURE2);

    // activate shader and set uniforms:
    program_->activate();
    setGlobalShaderParameters(program_);
    LGL_ERROR;
    program_->setUniform("colorTex_", 0);
    program_->setUniform("depthTex_", 1);
    inport_.setTextureParameters(program_, "textureParameters_");
    LGL_ERROR;

    program_->setUniform("filter_", 2);
    LGL_ERROR;
    filterPort_.setTextureParameters(program_, "filterParameters_");
    LGL_ERROR;

    program_->setUniform("filterSize_", filterSize_.get());

    LGL_ERROR;
    // render screen aligned quad:
    glDepthFunc(GL_ALWAYS);
    renderQuad();
    glDepthFunc(GL_LESS);

    program_->deactivate();
    outport_.deactivateTarget();
    glActiveTexture(GL_TEXTURE0);
    LGL_ERROR;
}

} // namespace
