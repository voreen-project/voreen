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

//header file
#include "samplerenderprocessor.h"
//needed headers (used in process())
#include "tgt/textureunit.h"

//we are in namespace voreen
namespace voreen {

SampleRenderProcessor::SampleRenderProcessor()
    : RenderProcessor()
    , inport_(Port::INPORT, "inport","Unmodified Image")
    , outport_(Port::OUTPORT, "outport","Modified Image")
    , saturationProp_("saturation", "Saturation"    ///< property ID and GUI-label
                     ,0.5f                          ///< default value
                     ,0.f, 1.f)                     ///< min and max value
{
    //register ports
    addPort(inport_);
    addPort(outport_);
    //register properties
    addProperty(saturationProp_);
}

void SampleRenderProcessor::initialize() {
    // call superclass function first
    RenderProcessor::initialize();

    // load fragment shader 'sample.frag'
    shader_ = ShdrMgr.loadSeparate("passthrough.vert", "sample.frag", generateHeader(), false);
}

void SampleRenderProcessor::deinitialize() {
    // free shader
    ShdrMgr.dispose(shader_);
    shader_ = 0;

    // call superclass function last
    RenderProcessor::deinitialize();
}

void SampleRenderProcessor::process() {
    // activate and clear output render target
    outport_.activateTarget();
    outport_.clearTarget();

    // bind input image to texture units
    tgt::TextureUnit colorUnit, depthUnit;
    inport_.bindTextures(colorUnit.getEnum(), depthUnit.getEnum());

    // activate shader and pass data
    shader_->activate();
    setGlobalShaderParameters(shader_); //(see RenderProcessor)

    // pass input image to shader
    inport_.setTextureParameters(shader_, "textureParameters_");
    shader_->setUniform("colorTex_", colorUnit.getUnitNumber());
    shader_->setUniform("depthTex_", depthUnit.getUnitNumber());

    // pass property value to shader
    shader_->setUniform("saturation_", saturationProp_.get());

    // render screen aligned quad to run the fragment shader
    renderQuad();//(see RenderProcessor)

    // cleanup
    shader_->deactivate();
    outport_.deactivateTarget();
    tgt::TextureUnit::setZeroUnit();

    // check for OpenGL errors
    LGL_ERROR;
}

} // namespace
