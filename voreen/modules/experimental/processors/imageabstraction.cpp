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

#include "imageabstraction.h"

#include "tgt/textureunit.h"

using tgt::TextureUnit;

namespace voreen {

ImageAbstraction::ImageAbstraction()
    : ImageProcessor("pp_imageabstraction"),
      minSigma_("minSigma", "Minimum sigma", 0.1f, 0.1f, 5.0f),
      maxSigma_("maxSigma", "Maximum sigma", 4.0f, 0.1f, 10.0f),
      mappingFunc_("mappingFunction", "Mapping function", Processor::INVALID_RESULT),
      inport_(Port::INPORT, "image.inport"),
      outport_(Port::OUTPORT, "image.outport")
{
    addProperty(minSigma_);
    addProperty(maxSigma_);
    addProperty(mappingFunc_);

    addPort(inport_);
    addPort(outport_);
}

void ImageAbstraction::process() {
    outport_.activateTarget();
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    tgt::vec2 depthRange = computeDepthRange(&inport_);

    TextureUnit colorUnit, depthUnit;
    inport_.bindTextures(colorUnit.getEnum(), depthUnit.getEnum());

    // bind mapping function
    TextureUnit mappingUnit;
    mappingUnit.activate();

    if (mappingFunc_.get())
        mappingFunc_.get()->getTexture()->bind();

    // initialize shader
    program_->activate();
    setGlobalShaderParameters(program_);
    program_->setUniform("colorTex_", colorUnit.getUnitNumber());
    program_->setUniform("depthTex_", depthUnit.getUnitNumber());
    program_->setUniform("minDepth_", depthRange.x);
    program_->setUniform("maxDepth_", depthRange.y);
    inport_.setTextureParameters(program_, "textureParameters_");
    program_->setUniform("minSigma_", minSigma_.get());
    program_->setUniform("maxSigma_", maxSigma_.get());
    program_->setUniform("mappingFunc_", mappingUnit.getUnitNumber());

    renderQuad();

    program_->deactivate();
    outport_.deactivateTarget();
    TextureUnit::setZeroUnit();
    LGL_ERROR;
}

} // voreen namespace
