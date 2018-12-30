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

#include "imagesource.h"

#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/voreenapplication.h"
#include "voreen/core/utils/stringutils.h"

#include "tgt/texturemanager.h"
#include "tgt/filesystem.h"
#include "tgt/gpucapabilities.h"

#ifdef VRN_MODULE_DEVIL
#include "modules/devil/devilmodule.h"
#endif

namespace voreen {

const std::string ImageSource::loggerCat_("voreen.core.ImageSource");

ImageSource::ImageSource()
    : RenderProcessor()
    , outport_(Port::OUTPORT, "image.out", "Image Output", false)
    , imageFile_("imageFile", "Image File", "Image File", VoreenApplication::app()->getUserDataPath())
    , reloadImage_("reloadImage", "Reload Image")
    , clearImage_("clearButton", "Clear Image")
    , imageSize_("imageSize", "Image Size", tgt::ivec2(0), tgt::ivec2(0), tgt::ivec2(1 << 12), Processor::VALID)
    , texture_(0)
    , shader_(0)
    , textureOwner_(false)
    , forceReload_(false)
{
    addPort(outport_);

    addProperty(imageFile_);
    ON_CHANGE(imageFile_, ImageSource, forceReload);
    addProperty(reloadImage_);
    ON_CHANGE(reloadImage_, ImageSource, forceReload);
    addProperty(clearImage_);
    ON_CHANGE(clearImage_, ImageSource, clearImage);
    addProperty(imageSize_);
    imageSize_.setReadOnlyFlag(true);
}

ImageSource::~ImageSource() {
}

Processor* ImageSource::create() const {
    return new ImageSource();
}

void ImageSource::beforeProcess() {

    if (forceReload_) {
        loadImage(imageFile_.get());
        forceReload_ = false;
        LGL_ERROR;
    }

    // resize port to image dimensions, if image is present
    if (texture_) {
        outport_.resize(texture_->getDimensions().xy());
        LGL_ERROR;
    }
}

void ImageSource::process() {

    // clear output
    outport_.activateTarget();
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    LGL_ERROR;

    // no image to render
    if (!texture_) {
        outport_.invalidateResult();
        outport_.deactivateTarget();
        LGL_ERROR;
        return;
    }

    // activate shader
    shader_->activate();

    // set common uniforms
    setGlobalShaderParameters(shader_);

    // bind image texture
    glActiveTexture(GL_TEXTURE0);
    texture_->bind();
    LGL_ERROR;

    // pass texture parameters to the shader
    shader_->setUniform("colorTex_", 0);
    shader_->setUniform("useTexcoords", false);
    shader_->setIgnoreUniformLocationError(true);
    shader_->setUniform("texParams_.dimensions_", tgt::vec2(texture_->getDimensions().xy()));
    shader_->setUniform("texParams_.dimensionsRCP_", tgt::vec2(1.f) / tgt::vec2(texture_->getDimensions().xy()));
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

void ImageSource::initialize() {
    RenderProcessor::initialize();

#ifdef VRN_MODULE_DEVIL
    std::string filterStr = "Image files (*." + strJoin(DevILModule::getSupportedReadExtensions(), " *.") + ")";
    imageFile_.setFileFilter(filterStr);
#endif

    shader_ = ShdrMgr.loadSeparate("passthrough.vert", "copyimage.frag",
        generateHeader() + "#define NO_DEPTH_TEX\n", false);

    loadImage(imageFile_.get());
}

void ImageSource::deinitialize() {
    ShdrMgr.dispose(shader_);
    shader_ = 0;
    clearImage();

    RenderProcessor::deinitialize();
}

void ImageSource::loadImage(const std::string& fname) {

    if (!tgt::TextureManager::isInited()) {
        return;
    }

    // necessary since the passed string reference might be changed during clearImage/invalidate,
    std::string filename = fname;

    // clear image and check for empty filename
    if (texture_) {
        clearImage();
        LGL_ERROR;
    }
    if (filename.empty())
        return;

    // load image as texture
    texture_ = TexMgr.load(filename, tgt::Texture::LINEAR,
        false, false, true, false);
    LGL_ERROR;
    if (texture_) {
        imageFile_.set(filename);
        imageSize_.set(texture_->getDimensions().xy());
        textureOwner_ = true;
        LINFO("Loaded image with dimensions " << texture_->getDimensions().xy() <<
            " from file " << tgt::FileSystem::cleanupPath(filename));
    }
    else {
        LWARNING("Failed to load image: " << tgt::FileSystem::cleanupPath(filename));
        imageFile_.set("");
        imageSize_.set(tgt::ivec2(0));
    }

    LGL_ERROR;
    invalidate();
}

void ImageSource::clearImage() {

    LGL_ERROR;
    if (texture_) {
        if (textureOwner_)
            delete texture_;
        texture_ = 0;
        LGL_ERROR;
    }
    textureOwner_ = false;
    imageFile_.set("");
    imageSize_.set(tgt::ivec2(0));
}

void ImageSource::forceReload() {
    forceReload_ = true;
    invalidate();
}

void ImageSource::setImage(tgt::Texture* texture) {

    if (!isInitialized()) {
        LWARNING("loadImage(): not initialized");
        return;
    }

    clearImage();
    texture_ = texture;
    textureOwner_ = false;
    if (texture_) {
        imageFile_.set(texture_->getOptionalName());
        imageSize_.set(texture_->getDimensions().xy());
    }
    invalidate();
}

const tgt::Texture* ImageSource::getImage() const {
    return texture_;
}

} // namespace
