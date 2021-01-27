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

#include "nonminmaxsuppression.h"

namespace voreen {

NonMinMaxSuppression::NonMinMaxSuppression()
    : ImageProcessorBypassable("image/nonminmax"),
      inportImage_(Port::INPORT, "image.inport", "Image Input"),
      inportSobel_(Port::INPORT, "image.sobel", "SobelImage Input"),
      outport_(Port::OUTPORT, "image.outport", "Image Output"),
      mode_("mode", "Mode", Processor::INVALID_PROGRAM),
      isotropic_("isotropic", "Isotropic (ignore edges)", true, Processor::INVALID_PROGRAM),
      kernelRadius_("kernelRadius", "Kernel radius", tgt::ivec2(1))
{
    addPort(inportImage_);
    addPort(inportSobel_);
    addPort(outport_);

    mode_.addOption("nonmax", "Non-Maximum", "NON_MAX_SUPPRESSION");
    mode_.addOption("nonmin", "Non-Minimum", "NON_MIN_SUPPRESSION");
    addProperty(mode_);
    addProperty(isotropic_);
    addProperty(kernelRadius_);
}

Processor* NonMinMaxSuppression::create() const {
    return new NonMinMaxSuppression();
}

void NonMinMaxSuppression::process() {
    if (!enableSwitch_.get()) {
        bypass(&inportImage_, &outport_);
        return;
    }

    if (getInvalidationLevel() >= Processor::INVALID_PROGRAM)
        compile();

    outport_.activateTarget();
    outport_.clearTarget();

    inportImage_.bindTextures(GL_TEXTURE0, GL_TEXTURE1);
    if (!isotropic_.get()) {
        tgtAssert(inportSobel_.isReady(), "No sobel image");
        inportSobel_.bindColorTexture(GL_TEXTURE2);
    }

    // initialize shader
    program_->activate();
    setGlobalShaderParameters(program_);
    program_->setUniform("colorTex_", 0);
    program_->setUniform("depthTex_", 1);
    inportImage_.setTextureParameters(program_, "textureParameters_");
    if (!isotropic_.get()) {
        program_->setUniform("sobelTex_", 2);
        inportSobel_.setTextureParameters(program_, "textureParametersSobel_");
    }
    else {
        program_->setUniform("kernelRadius_", kernelRadius_.get());
    }

    renderQuad();

    program_->deactivate();
    outport_.deactivateTarget();
    glActiveTexture(GL_TEXTURE0);
    LGL_ERROR;
}

bool NonMinMaxSuppression::isReady() const {
    if (!isInitialized())
        return false;

    if (!inportImage_.isReady() || !outport_.isReady())
        return false;

    if (!isotropic_.get() && !inportSobel_.isReady())
        return false;

    return true;
}

std::string NonMinMaxSuppression::generateHeader(const tgt::GpuCapabilities::GlVersion* version) {
    std::string header = ImageProcessor::generateHeader(version);
    header += "#define " + mode_.getValue() + "\n";
    if (isotropic_.get())
        header += "#define NON_MIN_MAX_ISOTROPIC \n";
    return header;
}

void NonMinMaxSuppression::compile() {
    if (program_)
        program_->setHeaders(generateHeader());
    ImageProcessor::compile();
}

} // voreen namespace
