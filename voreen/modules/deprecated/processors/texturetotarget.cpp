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

#include "texturetotarget.h"

using tgt::TextureUnit;

namespace voreen {

TextureToTarget::TextureToTarget()
    : ImageProcessor("image/texturetotarget")
    , inport_(Port::INPORT, "inport")
    , outport_(Port::OUTPORT, "outport")
{
    // register ports and properties
    addPort(inport_);
    addPort(outport_);
}

Processor* TextureToTarget::create() const {
    return new TextureToTarget();
}

void TextureToTarget::process() {
    // activate and clear output render target
    outport_.activateTarget();
    outport_.clearTarget();

    // bind input rendering to texture units
    TextureUnit colorUnit;
    const tgt::Texture* texture = inport_.getData();
    colorUnit.activate();
    texture->bind();

    // activate shader and set uniforms
    program_->activate();
    setGlobalShaderParameters(program_);
    //inport_.setTextureParameters(program_, "textureParameters_");
    program_->setUniform("colorTex_", colorUnit.getUnitNumber());

    // render screen aligned quad
    renderQuad();

    // cleanup
    program_->deactivate();
    outport_.deactivateTarget();
    TextureUnit::setZeroUnit();

    // check for OpenGL errors
    LGL_ERROR;
}

} // namespace
