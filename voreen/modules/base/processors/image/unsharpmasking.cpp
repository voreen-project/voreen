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

#include "unsharpmasking.h"

#include "voreen/core/utils/stringutils.h"
#include "tgt/textureunit.h"

using tgt::TextureUnit;

namespace voreen {

UnsharpMasking::UnsharpMasking()
    : ImageProcessorBypassable("image/unsharpmaskingblur"),
      radius_("radius", "Radius", 5, 0, 25),
      sigma_("sigma", "Sigma", 2.0f, 0.1f, 10.0f),
      amount_("amount", "Amount", 1.0f, 0.f, 10.0f),
      threshold_("threshold", "Threshhold", 0.f, 0.f, 0.5f),
      alpha_("alpha", "Process Alpha-Channel", true),
      inport_(Port::INPORT, "image.inport", "Image Input"),
      outport_(Port::OUTPORT, "image.outport", "Image Output"),
      privatePort_(Port::OUTPORT, "image.privateport")
{
    addProperty(sigma_);
    addProperty(radius_);
    addProperty(amount_);
    addProperty(threshold_);
    addProperty(alpha_);

    sigma_.onChange(MemberFunctionCallback<UnsharpMasking>(this, &UnsharpMasking::sigmaChanged));
    radius_.onChange(MemberFunctionCallback<UnsharpMasking>(this, &UnsharpMasking::radiusChanged));

    addPort(inport_);
    addPort(outport_);

    addPrivateRenderPort(&privatePort_);

    sigmaFlag_ = false;
    radiusFlag_ = false;
}

void UnsharpMasking::initialize() {
    ImageProcessorBypassable::initialize(); //call parent's initialize method, initializes first shader program (blur)

    program_->deactivate();

    //initialize second shader program (masking)
    secondProgram_ = ShdrMgr.loadSeparate("passthrough.vert", "image/unsharpmaskingmask.frag", generateHeader(), false);
    invalidate(Processor::INVALID_PROGRAM);
    if (getInvalidationLevel() >= Processor::INVALID_PROGRAM)
        secondProgram_->rebuild();

    secondProgram_->deactivate();
}

void UnsharpMasking::deinitialize() {
    ShdrMgr.dispose(secondProgram_);
    secondProgram_ = 0;
    LGL_ERROR;

    ImageProcessorBypassable::deinitialize();
}

void UnsharpMasking::process() {

    if (!enableSwitch_.get()) {
        bypass(&inport_, &outport_);
        return;
    }

    TextureUnit firstColorUnit, secondColorUnit, depthUnit;

    //Compute the Gauss kernel and norm
    float gaussKernel[26];
    float sigma = sigma_.get();
    int kernelRadius = static_cast<int>(sigma*2.5f);
    tgtAssert(kernelRadius >= 0 && kernelRadius <= 25, std::string("invalid kernel radius: " + itos(kernelRadius)).c_str());

    for (int i=0; i<=kernelRadius; i++)
        gaussKernel[i] = exp(-float(i*i)/(2.f*sigma*sigma));

    // compute norm
    float norm = 0.0;
    for (int i=1; i<=kernelRadius; i++)
        norm += gaussKernel[i];

    // so far we have just computed norm for one half
    norm = 2.f * norm + gaussKernel[0];

    //create a vector for the channels the masking should be applied to
    tgt::vec4 channelWeights;
    if (alpha_.get())
        channelWeights = tgt::vec4(1.f,1.f,1.f,1.f);
    else
        channelWeights = tgt::vec4(1.f,1.f,1.f,0.f);

    //Gaussian blur: first pass
    outport_.activateTarget();
    outport_.clearTarget();

    inport_.bindTextures(firstColorUnit.getEnum(), depthUnit.getEnum());

    program_->activate();
    setGlobalShaderParameters(program_);
    program_->setUniform("colorTex_", firstColorUnit.getUnitNumber());
    program_->setUniform("depthTex_", depthUnit.getUnitNumber());
    inport_.setTextureParameters(program_, "textureParameters_");
    program_->setUniform("dir_", tgt::vec2(1.f,0.f));
    program_->setUniform("gaussKernel_", gaussKernel, 26);
    program_->setUniform("norm_", norm);
    program_->setUniform("kernelRadius_", kernelRadius);
    program_->setUniform("channelWeights_", channelWeights);

    renderQuad();

    program_->deactivate();
    outport_.deactivateTarget();

    //Gaussian blur: second pass
    privatePort_.activateTarget();
    privatePort_.clearTarget();

    outport_.bindColorTexture(firstColorUnit.getEnum());

    program_->activate();
    setGlobalShaderParameters(program_);
    program_->setUniform("colorTex_", firstColorUnit.getUnitNumber());
    program_->setUniform("depthTex_", depthUnit.getUnitNumber());
    inport_.setTextureParameters(program_, "textureParameters_");
    program_->setUniform("dir_", tgt::vec2(0.f,1.f));
    program_->setUniform("gaussKernel_", gaussKernel, 26);
    program_->setUniform("norm_", norm);
    program_->setUniform("kernelRadius_", kernelRadius);
    program_->setUniform("channelWeights_", channelWeights);

    renderQuad();

    program_->deactivate();
    privatePort_.deactivateTarget();

    //use the masking shader with the original image and the blurred one to do the actual masking
    outport_.activateTarget();
    outport_.clearTarget();

    inport_.bindTextures(firstColorUnit.getEnum(), depthUnit.getEnum());
    privatePort_.bindColorTexture(secondColorUnit.getEnum());

    secondProgram_->activate();
    setGlobalShaderParameters(secondProgram_);
    secondProgram_->setUniform("originalColorTex_", firstColorUnit.getUnitNumber());
    secondProgram_->setUniform("blurredColorTex_", secondColorUnit.getUnitNumber());
    secondProgram_->setUniform("depthTex_", depthUnit.getUnitNumber());
    secondProgram_->setUniform("amount_", amount_.get()+1.f);
    secondProgram_->setUniform("threshold_", threshold_.get());
    inport_.setTextureParameters(secondProgram_, "textureParameters_");

    renderQuad();

    secondProgram_->deactivate();
    outport_.deactivateTarget();

    TextureUnit::setZeroUnit();

    LGL_ERROR;
}

/*
  Set Radius according to new sigma value
*/
void UnsharpMasking::sigmaChanged(){
    sigmaFlag_ = true;
    if (!radiusFlag_)
        radius_.set((int) (sigma_.get()*2.5f));
    else {
        sigmaFlag_ = false;
        radiusFlag_ = false;
    }
}

/*
  Set Sigma according to new radius
*/
void UnsharpMasking::radiusChanged(){
    radiusFlag_ = true;
    if (!sigmaFlag_){
       (radius_.get() > 0) ? sigma_.set((float) radius_.get()/2.5f) : sigma_.set(0.1f);
    }
    else {
        sigmaFlag_ = false;
        radiusFlag_ = false;
    }
}

} // voreen namespace
