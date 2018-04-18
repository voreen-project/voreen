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

#include "rawtexturesource.h"

#include "tgt/texturemanager.h"
#include "tgt/gpucapabilities.h"
#include "voreen/core/datastructures/callback/memberfunctioncallback.h"

namespace voreen {

const std::string RawTextureSource::loggerCat_("voreen.RawTextureSource");

RawTextureSource::RawTextureSource()
    : Processor(),
      texturePort_(Port::OUTPORT, "texture.out"),
      textureFile_("textureFile", "Texture File", "Texture File", VoreenApplication::app()->getUserDataPath()),
      loadTexButton_("loadButton", "Load Texture"),
      textureSize_("textureSize", "Texture Size", tgt::ivec2(0), tgt::ivec2(0), tgt::ivec2(1 << 12), VALID),
      textureFormat_("textureFormat", "Texture Format"),
      textureDataType_("textureDataType", "Texture Data Type"),
      texture_(0),
      textureOwner_(false),
      forceReload_(false)
{
    addPort(texturePort_);

    loadTexButton_.onChange(MemberFunctionCallback<RawTextureSource>(this, &RawTextureSource::loadTexture));
    addProperty(textureFile_);
    addProperty(loadTexButton_);
    addProperty(textureSize_);
    addProperty(textureFormat_);
    textureFormat_.addOption("la", "LUMINANCE_ALPHA", GL_LUMINANCE_ALPHA);
    textureFormat_.addOption("rgb", "RGB", GL_RGB);
    textureFormat_.addOption("rgba", "RGBA", GL_RGBA);
    addProperty(textureDataType_);
    textureDataType_.addOption("ubyte", "UBYTE", GL_UNSIGNED_BYTE);
    textureDataType_.addOption("ushort", "USHORT", GL_UNSIGNED_SHORT);
    textureDataType_.addOption("float", "FLOAT", GL_FLOAT);
}

RawTextureSource::~RawTextureSource() {
}

Processor* RawTextureSource::create() const {
    return new RawTextureSource();
}

void RawTextureSource::beforeProcess() {
    // perform enforced reload
    if (forceReload_) {
        loadTexture();
        LGL_ERROR;
    }
}

void RawTextureSource::process() {
    texturePort_.setData(texture_, false);
}

void RawTextureSource::initialize() {
    Processor::initialize();
    loadTexture();
}

void RawTextureSource::deinitialize() {
    clearTexture();
    Processor::deinitialize();
}

void RawTextureSource::loadTexture() {
    if (!tgt::TextureManager::isInited())
        return;

    // necessary since the passed string reference might be changed during clearTexture/invalidate,
    std::string filename = textureFile_.get();

    // clear image and check for empty filename
    if (texture_) {
        clearTexture();
        LGL_ERROR;
    }
    if (filename.empty())
        return;

    unsigned int numChannels = 0;
    switch (textureFormat_.getValue()) {
        case GL_LUMINANCE_ALPHA: numChannels = 2; break;
        case GL_RGB: numChannels = 3; break;
        case GL_RGBA: numChannels = 4; break;
        default:
            LWARNING("Unknown texture format selected!");
    }
    unsigned int dataSize = 0;
    switch (textureDataType_.getValue()) {
        case GL_UNSIGNED_BYTE: dataSize = 1; break;
        case GL_UNSIGNED_SHORT: dataSize = 2; break;
        case GL_FLOAT: dataSize = 4; break;
        default:
            LWARNING("Unknown texture data size selected!");
    }

    // load texture
    texture_ = new tgt::Texture(tgt::ivec3(textureSize_.get().x, textureSize_.get().y, 1), textureFormat_.getValue(), textureDataType_.getValue(), tgt::Texture::NEAREST);
    std::ifstream dataFile(filename.c_str(), std::ios::in | std::ios::binary);
    unsigned int numBytesToRead = textureSize_.get().x*textureSize_.get().y*numChannels*dataSize;
    char* texData = new char[numBytesToRead];
    if (!dataFile.read(texData, numBytesToRead))
        LERROR("Could not read texture, please check parameters.");
    dataFile.close();

    LGL_ERROR;
    if (texture_) {
        textureOwner_ = true;
        LINFO("Loaded texture with dimensions " << texture_->getDimensions().xy() << " from file " << filename);
    } else {
        LWARNING("Failed to load texture: " << filename);
    }

    LGL_ERROR;
    texture_->setCpuTextureData((GLubyte*)texData,true);
    texture_->uploadTexture();
    invalidate();
}

void RawTextureSource::clearTexture() {
    LGL_ERROR;
    if (texture_) {
        if (textureOwner_)
            delete texture_;
        texture_ = 0;
        LGL_ERROR;
    }
    textureOwner_ = false;
}

} // namespace
