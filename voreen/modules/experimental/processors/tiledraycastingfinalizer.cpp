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

#include "tiledraycastingfinalizer.h"

#include "tgt/textureunit.h"

using tgt::TextureUnit;

namespace voreen {

const std::string TiledRaycastingFinalizer::loggerCat_("voreen.experimental.TiledRaycastingFinalizer");

TiledRaycastingFinalizer::TiledRaycastingFinalizer()
    : RenderProcessor(),
      shaderPrg_(0),
      inportTile_(Port::INPORT, "inportTile", "Rendered Tile Input", false,
        Processor::INVALID_RESULT, RenderPort::RENDERSIZE_ORIGIN),
      outportRendering_(Port::OUTPORT, "outportRendering", "Rendering Output", true,
        Processor::INVALID_RESULT, RenderPort::RENDERSIZE_RECEIVER),
      loopOutport_(Port::OUTPORT, "loopOutport", "Loop Outport"),
      subdivisions_("subdivisions", "Subdivisions Per Dim", 2, 1, 6)
{
    loopOutport_.setLoopPort(true);

    addPort(inportTile_);
    addPort(outportRendering_);
    addPort(loopOutport_);

    addProperty(subdivisions_);
    ON_PROPERTY_CHANGE(subdivisions_, TiledRaycastingFinalizer, updateInputSizeRequest);

    outportRendering_.onSizeReceiveChange<TiledRaycastingFinalizer>(this, &TiledRaycastingFinalizer::updateInputSizeRequest);
}

Processor* TiledRaycastingFinalizer::create() const {
    return new TiledRaycastingFinalizer();
}

void TiledRaycastingFinalizer::initialize() {
    RenderProcessor::initialize();

    updateInputSizeRequest();

    shaderPrg_ = ShdrMgr.loadSeparate("passthrough.vert", "copytosubimage.frag",
        generateHeader(), false);
}

void TiledRaycastingFinalizer::deinitialize() {
    ShdrMgr.dispose(shaderPrg_);

    RenderProcessor::deinitialize();
}

bool TiledRaycastingFinalizer::isReady() const {
    return (inportTile_.isReady() && outportRendering_.isReady());
}

void TiledRaycastingFinalizer::process() {
    tgt::ivec2 tileDim = inportTile_.getSize();
    int iteration = loopOutport_.getLoopIteration();

    tgt::ivec2 tile;
    tile.y = iteration / subdivisions_.get();
    tile.x = iteration % subdivisions_.get();
    tgt::ivec2 tileLL = tile*tileDim;
    //tgt::ivec2 tileUR = tileLL + tileDim;
    //tgt::ivec2 outputDim = tileDim*subdivisions_.get();
    //LINFO("Tile: " << tile << ", LL: " << tileLL << ", UR: " << tileUR << ", outputDim: " << outputDim);

    if (iteration == 0) {
        //outportRendering_.resize(outputDim);
        outportRendering_.activateTarget();
        outportRendering_.clearTarget();
    }
    else {
        outportRendering_.activateTarget();
    }

    TextureUnit colorUnit, depthUnit;
    inportTile_.bindTextures(colorUnit.getEnum(), depthUnit.getEnum());
    LGL_ERROR;

    // initialize shader
    shaderPrg_->activate();

    // set common uniforms used by all shaders
    setGlobalShaderParameters(shaderPrg_);

    // pass the remaining uniforms to the shader
    shaderPrg_->setUniform("colorTex_", colorUnit.getUnitNumber());
    shaderPrg_->setUniform("depthTex_", depthUnit.getUnitNumber());
    inportTile_.setTextureParameters(shaderPrg_, "texParams_");

    shaderPrg_->setUniform("tileLL_", tileLL);
    shaderPrg_->setUniform("tileDim_", tileDim);
    renderQuad();

    shaderPrg_->deactivate();
    outportRendering_.deactivateTarget();
    TextureUnit::setZeroUnit();

    if (loopOutport_.getLoopIteration() < loopOutport_.getNumLoopIterations()-1) {
        loopOutport_.invalidatePort();
    }
}

void TiledRaycastingFinalizer::updateInputSizeRequest() {
    tgt::ivec2 receivedSize = outportRendering_.getReceivedSize();

    tgt::ivec2 sizeRequest = receivedSize / subdivisions_.get();
    inportTile_.requestSize(sizeRequest);
}

} // voreen namespace
