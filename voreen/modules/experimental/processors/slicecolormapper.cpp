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

#include "slicecolormapper.h"
#include "voreen/core/datastructures/transfunc/transfuncbase.h"

using tgt::vec3;

namespace voreen {

SliceColorMapper::SliceColorMapper()
    : RenderProcessor()
    , transferFunc_("transferFunction", "Transfer function", Processor::INVALID_RESULT)
    , inport_(Port::INPORT, "inport")
    , outport_(Port::OUTPORT, "output")
    , shader_(0)
{
    addProperty(transferFunc_);

    addPort(inport_);
    addPort(outport_);
}

void SliceColorMapper::deinitialize() {
    ShdrMgr.dispose(shader_);
    RenderProcessor::deinitialize();
}

SliceColorMapper::~SliceColorMapper() {
}

std::string SliceColorMapper::getProcessorInfo() const {
    return "Applies a transfer function to a slice.";
}

Processor* SliceColorMapper::create() const {
    return new SliceColorMapper();
}

void SliceColorMapper::initialize() {
    RenderProcessor::initialize();

    shader_ = ShdrMgr.loadSeparate("passthrough.vert", "colormapper.frag", generateHeader() + transferFunc_.get()->getShaderDefines(), false);

    processorState_ = PROCESSOR_STATE_NOT_READY;
}

void SliceColorMapper::compile() {
    if (!shader_)
        return;

    shader_->setHeaders(generateHeader() + transferFunc_.get()->getShaderDefines());
    shader_->rebuild();
}

void SliceColorMapper::process() {
    outport_.activateTarget();
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    if (getInvalidationLevel() >= Processor::INVALID_PROGRAM)
        compile();
    LGL_ERROR;

    inport_.bindColorTexture(GL_TEXTURE0);

    glActiveTexture(GL_TEXTURE1);
    if (transferFunc_.get())
        transferFunc_.get()->getTexture()->bind();

    shader_->activate();

    shader_->setUniform("input_", 0);
    inport_.setTextureParameters(shader_, "inputParameters_");

    transferFunc_.get()->setUniform(shader_, "transferFunc_", "transferFuncTex_", 1);

    renderQuad();

    shader_->deactivate();
    glActiveTexture(GL_TEXTURE0);

    outport_.deactivateTarget();
    LGL_ERROR;
}

} // namespace voreen
