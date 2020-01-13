/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2020 University of Muenster, Germany,                        *
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

#include "scale.h"

namespace voreen {

using tgt::vec2;
using tgt::ivec2;

ScalingProcessor::ScalingProcessor()
    : ImageProcessor("image/scale"),
      distributeEvents_("distributeEvents", "Distribute Events", true, Processor::VALID)
{
    addProperty(distributeEvents_);
}

void ScalingProcessor::applyScalingMatrix(int scalingMode, RenderPort* inport, RenderPort* outport) {
    float aspectRatioIn = (float)inport->getSize().x / (float)inport->getSize().y;
    float aspectRatioOut = (float)outport->getSize().x / (float)outport->getSize().y;
    MatStack.loadIdentity();
    switch(scalingMode) {
        case 0:
            break;
        case 1:
            //map pixels 1:1
            MatStack.scale((float)inport->getSize().x / (float)outport->getSize().x, (float)inport->getSize().y / (float)outport->getSize().y, 1.0f);
            break;
        case 2:
            //scale respecting aspect ratio and view all of the inport:
            if (aspectRatioOut < aspectRatioIn) {
                MatStack.scale(1.0f, aspectRatioOut/aspectRatioIn, 1.0f);
            }
            else {
                MatStack.scale(aspectRatioIn/aspectRatioOut, 1.0f, 1.0f);
            }
            break;
        case 3:
            //scale respecting aspect ratio and fill outport:
            if (aspectRatioOut < aspectRatioIn) {
                MatStack.scale(aspectRatioIn/aspectRatioOut, 1.0f, 1.0f);
            }
            else {
                MatStack.scale(1.0f, aspectRatioOut/aspectRatioIn, 1.0f);
            }
            break;
        case 4:
            //scale to height:
            MatStack.scale(aspectRatioIn/aspectRatioOut, 1.0f, 1.0f);
            break;
        case 5:
            //scale to width:
            MatStack.scale(1.0f, aspectRatioIn*aspectRatioOut, 1.0f);
            break;
        default:
            break;
    }
    setGlobalShaderParameters(program_);
}

void ScalingProcessor::onEvent(tgt::Event* e) {
    e->ignore();
    if (distributeEvents_.get())
        ImageProcessor::onEvent(e);
}

tgt::MouseEvent* ScalingProcessor::transformMouseCoordinates(tgt::MouseEvent* e, int scalingMode, RenderPort* inport, RenderPort* outport ) const {
    float aspectRatioIn = (float)inport->getSize().x / (float)inport->getSize().y;
    float aspectRatioOut = (float)outport->getSize().x / (float)outport->getSize().y;
    vec2 scalingFactor(1.f);
    switch(scalingMode) {
        case 0:
            break;
        case 1:
            //map pixels 1:1
            scalingFactor = vec2((float)inport->getSize().x / (float)outport->getSize().x, (float)inport->getSize().y / (float)outport->getSize().y);
            break;
        case 2:
            //scale respecting aspect ratio and view all of the inport:
            if (aspectRatioOut < aspectRatioIn) {
                scalingFactor = vec2(1.0f, aspectRatioOut/aspectRatioIn);
            }
            else {
                scalingFactor = vec2(aspectRatioIn/aspectRatioOut, 1.0f);
            }
            break;
        case 3:
            //scale respecting aspect ratio and fill outport:
            if (aspectRatioOut < aspectRatioIn) {
                scalingFactor = vec2(aspectRatioIn/aspectRatioOut, 1.0f);
            }
            else {
                scalingFactor = vec2(1.0f, aspectRatioOut/aspectRatioIn);
            }
            break;
        case 4:
            //scale to height:
            scalingFactor = vec2(aspectRatioIn/aspectRatioOut, 1.0f);
            break;
        case 5:
            //scale to width:
            scalingFactor = vec2(1.0f, aspectRatioIn*aspectRatioOut);
            break;
        default:
            break;
    }

    // compute transformation
    vec2 pixelOffset = ((1.f - scalingFactor) * vec2(outport->getSize())) / 2.f;
    vec2 pixelScale = (vec2(outport->getSize())*scalingFactor) / vec2(inport->getSize());
    vec2 trafoCoords = (vec2(e->coord()) - pixelOffset)/pixelScale;

    // clone event and assign transformed coords
    tgt::MouseEvent* trafoEvent = new tgt::MouseEvent(*e);
    trafoEvent->setCoord(ivec2(trafoCoords));
    trafoEvent->setViewport(inport->getSize());
    return trafoEvent;
}

//-----------------------------------------------------------------------------

SingleScale::SingleScale()
    : ScalingProcessor(),
    scalingMode_("scaling", "Scaling Mode", Processor::INVALID_RESULT),
    inport_(Port::INPORT, "image.inport", "Image Input"),
    outport_(Port::OUTPORT, "image.outport", "Image Output", true, Processor::INVALID_RESULT, RenderPort::RENDERSIZE_RECEIVER)
{

    scalingMode_.addOption("brute-force",     "Brute force",            0);
    scalingMode_.addOption("1:1",             "Pixels 1:1",             1);
    scalingMode_.addOption("keep-ratio-all",  "Keep ratio, view all",   2);
    scalingMode_.addOption("keep-ratio-fill", "Keep ratio, fill",       3);
    scalingMode_.addOption("scale-to-height", "Scale to height",        4);
    scalingMode_.addOption("scale-to-width",  "Scale to width",         5);
    scalingMode_.select("keep-ratio-all");
    addProperty(scalingMode_);

    addPort(inport_);
    addPort(outport_);
}

void SingleScale::process() {
    if (!inport_.isReady() || !outport_.isReady())
        return;

    inport_.bindTextures(GL_TEXTURE0, GL_TEXTURE1);
    LGL_ERROR;

    // initialize shader
    program_->activate();
    program_->setUniform("colorTex_", 0);
    program_->setUniform("depthTex_", 1);
    inport_.setTextureParameters(program_, "colorTexParameters_");
    inport_.setTextureParameters(program_, "depthTexParameters_");

    MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
    applyScalingMatrix(scalingMode_.getValue(), &inport_, &outport_);

    outport_.activateTarget();
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glDepthFunc(GL_ALWAYS);
    renderQuad();
    glDepthFunc(GL_LESS);
    outport_.deactivateTarget();

    MatStack.loadIdentity();

    program_->deactivate();
    glActiveTexture(GL_TEXTURE0);
    LGL_ERROR;
}

void SingleScale::onEvent(tgt::Event* e) {
    e->ignore();
    if (distributeEvents_.get()) {
        tgt::MouseEvent* me = dynamic_cast<tgt::MouseEvent*>(e);
        if (me) {
            me = transformMouseCoordinates(me, scalingMode_.getValue(), &inport_, &outport_);
            ScalingProcessor::onEvent(me);
            delete me;  // has been cloned by transformMouseCoordinates
        }
        else {
            ScalingProcessor::onEvent(e);
        }
    }
}
} // voreen namespace
