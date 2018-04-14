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

#include "gabor.h"

namespace voreen {

Gabor::Gabor()
    : ImageProcessor("gabor"), // loads fragment shader pp_greyscale.frag
      angle_("angle", "Angle", 0.0f, 0.0f, 180.0f),
      wavelength_("wavelength", "Wavelength", 0.85f, 0.05f, 5.0f),
      offset_("offset", "Offset", 0.0f, 0.0f, 5.0f),
      sigma_("sigma", "Sigma", 0.45f, 0.05f, 5.0f),
      aspectRatio_("aspectratio", "Aspect ratio", 0.9f, 0.0f, 50.0f),
      resolution_("resolution_", "Resolution", 7, 3, 512),
      outport_(Port::OUTPORT, "outport")
{
    // register properties and ports:
    addProperty(angle_);
    addProperty(wavelength_);
    addProperty(offset_);
    addProperty(sigma_);
    addProperty(aspectRatio_);
    addProperty(resolution_);

    addPort(outport_);
}

Processor* Gabor::create() const {
    return new Gabor();
}

void Gabor::process() {
    outport_.resize(tgt::ivec2(resolution_.get(), resolution_.get()));
    // activate and clear output render target
    outport_.activateTarget();
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // activate shader and set uniforms:
    program_->activate();
    setGlobalShaderParameters(program_);
    program_->setUniform("orientation_", (angle_.get()/360.0f)*2.0f*tgt::PIf);
    program_->setUniform("wavelength_", wavelength_.get());
    program_->setUniform("offset_", offset_.get());
    program_->setUniform("sigma_", sigma_.get());
    program_->setUniform("aRatio_", aspectRatio_.get());
    //program_->setUniform("aspectRatio_", aspectRatio_.get());

    // render screen aligned quad:
    glDepthFunc(GL_ALWAYS);
    renderQuad();
    glDepthFunc(GL_LESS);

    program_->deactivate();
    outport_.deactivateTarget();
    LGL_ERROR;
}

} // namespace
