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

#include "binaryimageprocessor.h"

#include "tgt/textureunit.h"

using tgt::TextureUnit;

namespace voreen {

BinaryImageProcessor::BinaryImageProcessor()
    : ImageProcessor("image/binary")
    , inport0_(Port::INPORT, "image.inport0", "Image0 Input")
    , inport1_(Port::INPORT, "image.inport1", "Image1 Input")
    , outport_(Port::OUTPORT, "image.outport", "Image Output", true)
{
    addPort(inport0_);
    addPort(inport1_);
    addPort(outport_);
}

Processor* BinaryImageProcessor::create() const {
    return new BinaryImageProcessor();
}

void BinaryImageProcessor::process() {
    // compile program if needed
    if (getInvalidationLevel() >= Processor::INVALID_PROGRAM)
        compile();
    LGL_ERROR;

    if(!program_ || !program_->isLinked())
        return;

    outport_.activateTarget();
    outport_.clearTarget();

    TextureUnit colorUnit0, depthUnit0, colorUnit1, depthUnit1;
    inport0_.bindTextures(colorUnit0.getEnum(), depthUnit0.getEnum());
    inport1_.bindTextures(colorUnit1.getEnum(), depthUnit1.getEnum());

    // initialize shader
    program_->activate();
    setGlobalShaderParameters(program_);

    program_->setIgnoreUniformLocationError(true);
    program_->setUniform("colorTex0_", colorUnit0.getUnitNumber());
    program_->setUniform("depthTex0_", depthUnit0.getUnitNumber());
    program_->setUniform("colorTex1_", colorUnit1.getUnitNumber());
    program_->setUniform("depthTex1_", depthUnit1.getUnitNumber());
    inport0_.setTextureParameters(program_, "texParams0_");
    inport1_.setTextureParameters(program_, "texParams1_");
    program_->setIgnoreUniformLocationError(false);

    glDepthFunc(GL_ALWAYS);
    renderQuad();
    glDepthFunc(GL_LESS);

    program_->deactivate();
    outport_.deactivateTarget();
    TextureUnit::setZeroUnit();
    LGL_ERROR;
}

} // voreen namespace
