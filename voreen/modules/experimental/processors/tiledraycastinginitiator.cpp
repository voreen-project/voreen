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

#include "tiledraycastinginitiator.h"

#include "tgt/textureunit.h"
#include "voreen/core/datastructures/callback/memberfunctioncallback.h"

using tgt::TextureUnit;

namespace voreen {

const std::string TiledRaycastingInitiator::loggerCat_("voreen.experimental.TiledRaycastingInitiator");

TiledRaycastingInitiator::TiledRaycastingInitiator()
    : RenderProcessor(),
      shader_(0),
      inportEntryPoints_(Port::INPORT, "inportEntryPoints", "Entry Points Input"),
      inportExitPoints_(Port::INPORT, "inportExitPoints", "Exit Points Input"),
      outportEntryPoints_(Port::OUTPORT, "outportEntryPoints", "Entry Points Output", true,
        Processor::INVALID_RESULT, RenderPort::RENDERSIZE_STATIC, GL_RGBA16F_ARB),
      outportExitPoints_(Port::OUTPORT, "outportExitPoints", "Exit Points Output", true,
        Processor::INVALID_RESULT, RenderPort::RENDERSIZE_STATIC, GL_RGBA16F_ARB),
      loopInport_(Port::INPORT, "loopInport", "Loop Inport"),
      subdivisions_("subdivisions", "Subdivisions Per Dim", 2, 1, 6)
{
    loopInport_.setLoopPort(true);
    loopInport_.setNumLoopIterations(subdivisions_.get()*subdivisions_.get());

    subdivisions_.onChange(MemberFunctionCallback<TiledRaycastingInitiator>(this, &TiledRaycastingInitiator::updateSubdivision));
    addProperty(subdivisions_);

    addPort(inportEntryPoints_);
    addPort(inportExitPoints_);
    addPort(outportEntryPoints_);
    addPort(outportExitPoints_);
    addPort(loopInport_);
}

Processor* TiledRaycastingInitiator::create() const {
    return new TiledRaycastingInitiator();
}

void TiledRaycastingInitiator::initialize() {
    RenderProcessor::initialize();

    shader_ = ShdrMgr.loadSeparate("passthrough.vert", "copyfromsubimage.frag",
        generateHeader(), false);

    updateSubdivision();
}

void TiledRaycastingInitiator::deinitialize() {
    ShdrMgr.dispose(shader_);

    RenderProcessor::deinitialize();
}

bool TiledRaycastingInitiator::isReady() const {
    return (inportEntryPoints_.isReady() && outportEntryPoints_.isConnected()) ||
           (inportExitPoints_.isReady()  && outportExitPoints_.isConnected());
}

void TiledRaycastingInitiator::process() {
    tgt::ivec2 tileDim;
    if (inportEntryPoints_.isReady())
        tileDim = inportEntryPoints_.getSize() / subdivisions_.get();
    else if (inportExitPoints_.isReady())
        tileDim = inportExitPoints_.getSize() / subdivisions_.get();
    else {
        LERROR("No inport is ready");
        return;
    }

    if (!tgt::hand(tgt::greaterThanEqual(tileDim, tgt::ivec2(2)))) {
        LERROR("Tile size too small: " << tileDim);
        return;
    }

    int iteration = loopInport_.getLoopIteration();
    tgt::ivec2 tile;
    tile.y = iteration / subdivisions_.get();
    tile.x = iteration % subdivisions_.get();
    tgt::ivec2 tileLL = tile*tileDim;
    //tgt::ivec2 tileUR = tileLL + tileDim;
    //LINFO("Tile: " << tile << ", LL: " << tileLL << ", UR: " << tileUR);

    if (iteration == 0)
        setProgress(0.f);

    shader_->activate();

    // 1. copy entry points texture tile
    if (inportEntryPoints_.isReady() && outportEntryPoints_.isConnected()) {
        outportEntryPoints_.resize(tileDim);
        outportEntryPoints_.activateTarget();
        outportEntryPoints_.clearTarget();

        TextureUnit colorUnit, depthUnit;
        inportEntryPoints_.bindTextures(colorUnit.getEnum(), depthUnit.getEnum());
        LGL_ERROR;

        setGlobalShaderParameters(shader_);
        shader_->setUniform("colorTex_", colorUnit.getUnitNumber());
        shader_->setUniform("depthTex_", depthUnit.getUnitNumber());
        inportEntryPoints_.setTextureParameters(shader_, "texParams_");
        shader_->setUniform("subImageOffset_", tileLL);
        LGL_ERROR;

        renderQuad();
        outportEntryPoints_.deactivateTarget();
    }

    // 2. copy exit points texture tile
    if (inportExitPoints_.isReady() && outportExitPoints_.isReady()) {
        outportExitPoints_.resize(tileDim);
        outportExitPoints_.activateTarget();
        outportExitPoints_.clearTarget();

        TextureUnit colorUnit, depthUnit;
        inportExitPoints_.bindTextures(colorUnit.getEnum(), depthUnit.getEnum());
        LGL_ERROR;

        setGlobalShaderParameters(shader_);
        shader_->setUniform("colorTex_", colorUnit.getUnitNumber());
        shader_->setUniform("depthTex_", depthUnit.getUnitNumber());
        inportExitPoints_.setTextureParameters(shader_, "texParams_");
        shader_->setUniform("subImageOffset_", tileLL);
        LGL_ERROR;

        renderQuad();
        outportExitPoints_.deactivateTarget();
    }

    shader_->deactivate();

    float progress;
    if (loopInport_.getNumLoopIterations() > 0)
        progress = (float)(iteration) / (loopInport_.getNumLoopIterations()-1.f);
    else
        progress = 1.f;
    setProgress(progress);
}

void TiledRaycastingInitiator::updateSubdivision() {
    loopInport_.setNumLoopIterations(subdivisions_.get()*subdivisions_.get());
}

} // voreen namespace
