/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2019 University of Muenster, Germany,                        *
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

#include "imageselector.h"

#include "voreen/core/datastructures/imagesequence.h"
#include "voreen/core/interaction/mwheelnumpropinteractionhandler.h"

namespace voreen {

const std::string ImageSelector::loggerCat_("voreen.core.ImageSelector");

ImageSelector::ImageSelector()
    : RenderProcessor(),
      inport_(Port::INPORT, "imagesequence.in", "ImageSequence Input", false),
      outport_(Port::OUTPORT, "image.out", "image.out", false),
      imageID_("imageID", "Selected Image", -1, -1, std::numeric_limits<int>::max() - 1),
      imageSize_("imageSize", "Image Size", tgt::ivec2(0), tgt::ivec2(0), tgt::ivec2(1 << 12), VALID),
      wheelHandler_("wheelHandler.imageCycling", "Image Cycling", &imageID_),
      shader_(0)
{
    addPort(inport_);
    addPort(outport_);

    addProperty(imageID_);
    addProperty(imageSize_);
    imageSize_.setReadOnlyFlag(true);
    addInteractionHandler(&wheelHandler_);
}

Processor* ImageSelector::create() const {
    return new ImageSelector();
}

void ImageSelector::process() {

    if (imageID_.get() == -1) {
        outport_.clear();
        outport_.invalidatePort();
        return;
    }

    // get texture to render
    tgt::Texture* tex = inport_.getData()->at(imageID_.get());
    tgtAssert(tex, "Texture is null");

    // adjust outport size to texture dimensions
    outport_.resize(tex->getDimensions().xy());
    imageSize_.set(tex->getDimensions().xy());

    // clear outport
    outport_.activateTarget();
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    LGL_ERROR;

    // activate shader
    shader_->activate();

    // set common uniforms
    setGlobalShaderParameters(shader_);

    // bind image texture
    glActiveTexture(GL_TEXTURE0);
    tex->bind();

    // pass texture parameters to the shader
    shader_->setUniform("colorTex_", 0);
    shader_->setIgnoreUniformLocationError(true);
    shader_->setUniform("texParams_.dimensions_", tgt::vec2(tex->getDimensions().xy()));
    shader_->setUniform("texParams_.dimensionsRCP_", tgt::vec2(1.f) / tgt::vec2(tex->getDimensions().xy()));
    shader_->setUniform("texParams_.matrix_", tgt::mat4::identity);
    shader_->setIgnoreUniformLocationError(false);

    // execute the shader
    renderQuad();
    LGL_ERROR;

    // clean up
    shader_->deactivate();
    outport_.deactivateTarget();
    LGL_ERROR;
}

void ImageSelector::initialize() {
    RenderProcessor::initialize();

    shader_ = ShdrMgr.loadSeparate("passthrough.vert", "copyimage.frag",
        generateHeader() + "#define NO_DEPTH_TEX\n", false);

    adjustPropertiesToInput();
}

void ImageSelector::deinitialize() {
    ShdrMgr.dispose(shader_);
    shader_ = 0;

    RenderProcessor::deinitialize();
}

void ImageSelector::adjustPropertiesToInput() {

    const ImageSequence* sequence = inport_.getData();

    if(!sequence)
        return;

    //if we have a list, adapt min and max values
    imageID_.setMinValue(std::min(0, static_cast<int>(sequence->size()) - 1));
    imageID_.setMaxValue(static_cast<int>(sequence->size()) - 1);
    // set to first image if no image was present earlier
    if (!sequence->empty() && imageID_.get() == -1)
        imageID_.set(0);
}

} // namespace
