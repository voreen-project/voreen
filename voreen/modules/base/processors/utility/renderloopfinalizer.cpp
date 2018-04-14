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

#include "renderloopfinalizer.h"

#include "tgt/textureunit.h"

using tgt::TextureUnit;

namespace voreen {

RenderLoopFinalizer::RenderLoopFinalizer()
    : RenderProcessor(),
      shaderPrg_(0),
      inport_(Port::INPORT, "image.inport", "Image Input"),
      outport_(Port::OUTPORT, "image.outport", "Image Output"),
      loopOutport_(Port::OUTPORT, "loop.outport", "Loop Outport")
{
    loopOutport_.setLoopPort(true);

    addPort(inport_);
    addPort(outport_);
    addPort(loopOutport_);
}

Processor* RenderLoopFinalizer::create() const {
    return new RenderLoopFinalizer();
}

void RenderLoopFinalizer::initialize() {
    RenderProcessor::initialize();

    shaderPrg_ = ShdrMgr.loadSeparate("passthrough.vert", "copyimage.frag",
        generateHeader(), false);
}

void RenderLoopFinalizer::deinitialize() {
    ShdrMgr.dispose(shaderPrg_);

    RenderProcessor::deinitialize();
}

bool RenderLoopFinalizer::isReady() const {
    return (inport_.isReady() && outport_.isReady());
}

void RenderLoopFinalizer::process() {
    if (loopOutport_.getLoopIteration() == loopOutport_.getNumLoopIterations()-1) {
        // last iteration: write input image to outport
        outport_.activateTarget();
    }
    else {
        // loop iteration: write input image to loop port
        loopOutport_.activateTarget();
    }

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    TextureUnit colorUnit, depthUnit;
    inport_.bindTextures(colorUnit.getEnum(), depthUnit.getEnum());
    LGL_ERROR;

    // initialize shader
    shaderPrg_->activate();

    // set common uniforms used by all shaders
    setGlobalShaderParameters(shaderPrg_);

    // pass the remaining uniforms to the shader
    shaderPrg_->setUniform("useTexcoords", false);
    shaderPrg_->setUniform("colorTex_", colorUnit.getUnitNumber());
    shaderPrg_->setUniform("depthTex_", depthUnit.getUnitNumber());
    inport_.setTextureParameters(shaderPrg_, "texParams_");

    renderQuad();

    shaderPrg_->deactivate();

    if (loopOutport_.getLoopIteration() == loopOutport_.getNumLoopIterations()-1) {
        outport_.deactivateTarget();
    } else {
        loopOutport_.deactivateTarget();
    }
    TextureUnit::setZeroUnit();

    LGL_ERROR;
}

} // voreen namespace
