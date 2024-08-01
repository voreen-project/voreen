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

#include "renderstore.h"

#include "tgt/textureunit.h"

using tgt::TextureUnit;

namespace voreen {

const std::string RenderStore::loggerCat_("voreen.RenderStore");

RenderStore::RenderStore()
    : RenderProcessor()
    , inport_(Port::INPORT, "image.input", "Image Input")
    , privatePort_(Port::OUTPORT, "image.stored", "image.stored", true)
    , outport_(Port::OUTPORT, "image.output", "Image Output")
    , shaderPrg_(0)
{

    addPort(inport_);
    addPrivateRenderPort(&privatePort_);
    addPort(outport_);
}

Processor* RenderStore::create() const {
    return new RenderStore();
}

bool RenderStore::isReady() const {
    return (inport_.isReady() && privatePort_.hasRenderTarget());
}

void RenderStore::process() {

    tgtAssert(inport_.isReady(), "Inport not ready");
    tgtAssert(privatePort_.hasRenderTarget(), "Private port has no data");

    //
    // render private port to outport, if connected
    //
    if (outport_.isConnected()) {

        // put out rendering currently stored in private port, if available
        // otherwise, pass-through input rendering
        RenderPort& storePort = privatePort_.hasValidResult() ? privatePort_ : inport_;

        // adjust size
        if (outport_.getSize() != storePort.getSize()) {
            outport_.resize(storePort.getSize());
            LGL_ERROR;
        }

        outport_.activateTarget();
        glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
        LGL_ERROR;

        TextureUnit storeUnit, storeUnitDepth;
        storePort.bindTextures(storeUnit.getEnum(), storeUnitDepth.getEnum());
        LGL_ERROR;

        shaderPrg_->activate();
        setGlobalShaderParameters(shaderPrg_, 0, storePort.getSize());
        storePort.setTextureParameters(shaderPrg_, "texParams_");
        shaderPrg_->setUniform("colorTex_", storeUnit.getUnitNumber());
        shaderPrg_->setUniform("depthTex_", storeUnitDepth.getUnitNumber());
        LGL_ERROR;

        glDepthFunc(GL_ALWAYS);
        renderQuad();
        glDepthFunc(GL_LESS);
        shaderPrg_->deactivate();
        outport_.deactivateTarget();
        LGL_ERROR;
    }


    //
    // render inport to private port
    //
    if (inport_.getSize() != privatePort_.getSize()) {
        privatePort_.resize(inport_.getSize());
        LGL_ERROR;
    }

    privatePort_.activateTarget("RenderStore::process() image.private");
    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
    LGL_ERROR;

    TextureUnit storeUnit, storeUnitDepth;
    inport_.bindTextures(storeUnit.getEnum(), storeUnitDepth.getEnum());
    LGL_ERROR;

    shaderPrg_->activate();
    setGlobalShaderParameters(shaderPrg_, 0, inport_.getSize());
    inport_.setTextureParameters(shaderPrg_, "texParams_");
    shaderPrg_->setUniform("colorTex_", storeUnit.getUnitNumber());
    shaderPrg_->setUniform("depthTex_", storeUnitDepth.getUnitNumber());
    LGL_ERROR;

    glDepthFunc(GL_ALWAYS);
    renderQuad();
    glDepthFunc(GL_LESS);
    shaderPrg_->deactivate();
    privatePort_.deactivateTarget();
    LGL_ERROR;
}

void RenderStore::initialize() {
    RenderProcessor::initialize();

    shaderPrg_ = ShdrMgr.loadSeparate("passthrough.vert", "copyimage.frag",
        generateHeader(), false);
}

void RenderStore::deinitialize() {
    ShdrMgr.dispose(shaderPrg_);
    shaderPrg_ = 0;

    RenderProcessor::deinitialize();
}

tgt::vec4 RenderStore::getStoredTargetPixel(const tgt::ivec2 &pos) {

    if (!privatePort_.hasValidResult()) {
        LWARNING("Private port has no valid result");
        return tgt::vec4(-1.f);
    }
    else {
        return privatePort_.getRenderTarget()->getColorAtPos(pos);
    }
}

} // namespace voreen
