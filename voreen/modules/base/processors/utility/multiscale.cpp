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

#include "multiscale.h"

namespace voreen {

MultiScale::MultiScale()
    : ScalingProcessor(),
      scalingMode1_("scalingMode1", "Scaling Mode 1", Processor::INVALID_RESULT),
      scalingMode2_("scalingMode2", "Scaling Mode 2", Processor::INVALID_RESULT),
      scalingMode3_("scalingMode3", "Scaling Mode 3", Processor::INVALID_RESULT),
      scalingMode4_("scalingMode4", "Scaling Mode 4", Processor::INVALID_RESULT),
      selectionMode_("selectionMode", "Selection Mode", Processor::INVALID_RESULT),
      inport_(Port::INPORT, "image.inport", "inport", false, Processor::INVALID_RESULT, RenderPort::RENDERSIZE_ORIGIN),
      outport1_(Port::OUTPORT, "image.outport1", "Image1 Output", true, Processor::INVALID_RESULT, RenderPort::RENDERSIZE_RECEIVER),
      outport2_(Port::OUTPORT, "image.outport2", "Image2 Output", true, Processor::INVALID_RESULT, RenderPort::RENDERSIZE_RECEIVER),
      outport3_(Port::OUTPORT, "image.outport3", "Image3 Output", true, Processor::INVALID_RESULT, RenderPort::RENDERSIZE_RECEIVER),
      outport4_(Port::OUTPORT, "image.outport4", "Image4 Output", true, Processor::INVALID_RESULT, RenderPort::RENDERSIZE_RECEIVER)
{

    scalingMode1_.addOption("brute-force",     "Brute force",            0);
    scalingMode1_.addOption("1:1",             "Pixels 1:1",             1);
    scalingMode1_.addOption("keep-ratio-all",  "Keep ratio, view all",   2);
    scalingMode1_.addOption("keep-ratio-fill", "Keep ratio, fill",       3);
    scalingMode1_.addOption("scale-to-height", "Scale to height",        4);
    scalingMode1_.addOption("scale-to-width",  "Scale to width",        5);
    scalingMode1_.select("keep-ratio-all");
    addProperty(scalingMode1_);

    scalingMode2_.addOption("brute-force",     "Brute force",            0);
    scalingMode2_.addOption("1:1",             "Pixels 1:1",             1);
    scalingMode2_.addOption("keep-ratio-all",  "Keep ratio, view all",   2);
    scalingMode2_.addOption("keep-ratio-fill", "Keep ratio, fill",       3);
    scalingMode2_.addOption("scale-to-height", "Scale to height",        4);
    scalingMode2_.addOption("scale-to-width",  "Scale to width",        5);
    scalingMode2_.select("keep-ratio-all");
    addProperty(scalingMode2_);

    scalingMode3_.addOption("brute-force",     "Brute force",            0);
    scalingMode3_.addOption("1:1",             "Pixels 1:1",             1);
    scalingMode3_.addOption("keep-ratio-all",  "Keep ratio, view all",   2);
    scalingMode3_.addOption("keep-ratio-fill", "Keep ratio, fill",       3);
    scalingMode3_.addOption("scale-to-height", "Scale to height",        4);
    scalingMode3_.addOption("scale-to-width",  "Scale to width",        5);
    scalingMode3_.select("keep-ratio-all");
    addProperty(scalingMode3_);

    scalingMode4_.addOption("brute-force",     "Brute force",            0);
    scalingMode4_.addOption("1:1",             "Pixels 1:1",             1);
    scalingMode4_.addOption("keep-ratio-all",  "Keep ratio, view all",   2);
    scalingMode4_.addOption("keep-ratio-fill", "Keep ratio, fill",       3);
    scalingMode4_.addOption("scale-to-height", "Scale to height",        4);
    scalingMode4_.addOption("scale-to-width",  "Scale to width",        5);
    scalingMode4_.select("keep-ratio-all");
    addProperty(scalingMode4_);

    selectionMode_.addOption("size",           "Size");
    selectionMode_.addOption("aspect-ratio",  "Aspect ratio");
    selectionMode_.addOption("width",          "Width");
    selectionMode_.addOption("height",         "Height");
    addProperty(selectionMode_);

    addPort(inport_);

    addPort(outport1_);
    addPort(outport2_);
    addPort(outport3_);
    addPort(outport4_);

    outport1_.onSizeReceiveChange(this, &MultiScale::portSizeReceiveChanged);
    outport2_.onSizeReceiveChange(this, &MultiScale::portSizeReceiveChanged);
    outport3_.onSizeReceiveChange(this, &MultiScale::portSizeReceiveChanged);
    outport4_.onSizeReceiveChange(this, &MultiScale::portSizeReceiveChanged);
}

MultiScale::~MultiScale() {}

void MultiScale::initialize() {
    ScalingProcessor::initialize();
    inport_.requestSize(selectBest());
}

bool MultiScale::isReady() const {
    if (!inport_.isReady() || (!outport1_.isReady() && !outport2_.isReady() && !outport3_.isReady() && !outport4_.isReady()))
        return false;

    return true;
}

void MultiScale::process() {
    inport_.bindTextures(GL_TEXTURE0, GL_TEXTURE1);
    LGL_ERROR;

    // initialize shader
    program_->activate();
    program_->setUniform("colorTex_", 0);
    program_->setUniform("depthTex_", 1);
    inport_.setTextureParameters(program_, "colorTexParameters_");
    inport_.setTextureParameters(program_, "depthTexParameters_");

    MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);

    glDepthFunc(GL_ALWAYS);
    if (outport1_.isReady()) {
        outport1_.activateTarget();
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        applyScalingMatrix(scalingMode1_.getValue(), &inport_, &outport1_);
        renderQuad();
        outport1_.deactivateTarget();
    }
    if (outport2_.isReady()) {
        outport2_.activateTarget();
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        applyScalingMatrix(scalingMode2_.getValue(), &inport_, &outport2_);
        renderQuad();
        outport2_.deactivateTarget();
    }
    if (outport3_.isReady()) {
        outport3_.activateTarget();
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        applyScalingMatrix(scalingMode3_.getValue(), &inport_, &outport3_);
        renderQuad();
        outport3_.deactivateTarget();
    }
    if (outport4_.isReady()) {
        outport4_.activateTarget();
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        applyScalingMatrix(scalingMode4_.getValue(), &inport_, &outport4_);
        renderQuad();
        outport4_.deactivateTarget();
    }
    glDepthFunc(GL_LESS);

    MatStack.loadIdentity();

    program_->deactivate();
    glActiveTexture(GL_TEXTURE0);
    LGL_ERROR;
}

void MultiScale::portSizeReceiveChanged() {
    inport_.requestSize(selectBest());
}


tgt::ivec2 MultiScale::selectBest() {
    if (!outport1_.isConnected() && !outport2_.isConnected() && !outport3_.isConnected() && !outport4_.isConnected())
        return tgt::ivec2(128,128);

    tgt::ivec2 best(0,0);

    if (selectionMode_.get() == "size") {
        if (outport1_.isConnected() && (hmul(best) < hmul(outport1_.getReceivedSize())))
            best = outport1_.getReceivedSize();

        if (outport2_.isConnected() && (hmul(best) < hmul(outport2_.getReceivedSize())))
            best = outport2_.getReceivedSize();

        if (outport3_.isConnected() && (hmul(best) < hmul(outport3_.getReceivedSize())))
            best = outport3_.getReceivedSize();

        if (outport4_.isConnected() && (hmul(best) < hmul(outport4_.getReceivedSize())))
            best = outport4_.getReceivedSize();
    }
    else if (selectionMode_.get() == "aspect-ratio") {

        float bestAR = 999999.0;
        if (outport1_.isConnected()) {
            float ar = outport1_.getReceivedSize().x / (float)outport1_.getReceivedSize().y;
            ar = std::max(ar, 1.0f/ar);
               if (ar < bestAR) {
                best = outport1_.getReceivedSize();
                bestAR = best.x / (float)best.y;
                bestAR = std::max(bestAR, 1.0f/bestAR);
            }
        }

        if (outport2_.isConnected()) {
            float ar = outport2_.getReceivedSize().x / (float)outport2_.getReceivedSize().y;
            ar = std::max(ar, 1.0f/ar);
               if (ar < bestAR) {
                best = outport2_.getReceivedSize();
                bestAR = best.x / (float)best.y;
                bestAR = std::max(bestAR, 1.0f/bestAR);
            }
        }

        if(outport3_.isConnected()) {
            float ar = outport3_.getReceivedSize().x / (float)outport3_.getReceivedSize().y;
            ar = std::max(ar, 1.0f/ar);
               if (ar < bestAR) {
                best = outport3_.getReceivedSize();
                bestAR = best.x / (float)best.y;
                bestAR = std::max(bestAR, 1.0f/bestAR);
            }
        }

        if(outport4_.isConnected()) {
            float ar = outport4_.getReceivedSize().x / (float)outport4_.getReceivedSize().y;
            ar = std::max(ar, 1.0f/ar);
               if (ar < bestAR) {
                best = outport4_.getReceivedSize();
                bestAR = best.x / (float)best.y;
                bestAR = std::max(bestAR, 1.0f/bestAR);
            }
        }
    }
    else if (selectionMode_.get() == "width") {
        if (outport1_.isConnected() && (best.x < outport1_.getReceivedSize().x))
            best = outport1_.getReceivedSize();

        if (outport2_.isConnected() && (best.x < outport2_.getReceivedSize().x))
            best = outport2_.getReceivedSize();

        if (outport3_.isConnected() && (best.x < outport3_.getReceivedSize().x))
            best = outport3_.getReceivedSize();

        if (outport4_.isConnected() && (best.x < outport4_.getReceivedSize().x))
            best = outport4_.getReceivedSize();
    }
    else if (selectionMode_.get() == "height") {
        if (outport1_.isConnected() && (best.y < outport1_.getReceivedSize().y))
            best = outport1_.getReceivedSize();

        if (outport2_.isConnected() && (best.y < outport2_.getReceivedSize().y))
            best = outport2_.getReceivedSize();

        if (outport3_.isConnected() && (best.y < outport3_.getReceivedSize().y))
            best = outport3_.getReceivedSize();

        if (outport4_.isConnected() && (best.y < outport4_.getReceivedSize().y))
            best = outport4_.getReceivedSize();
    }
    return best;
}

} // voreen namespace
