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

#include "rawtexturesave.h"
#include "voreen/core/datastructures/callback/memberfunctioncallback.h"

namespace voreen {

const std::string RawTextureSave::loggerCat_("voreen.RawTextureSave");

RawTextureSave::RawTextureSave()
    : RenderProcessor(),
      texturePort_(Port::INPORT, "texture.in"),
      textureFile_("textureFile", "texture File", "Texture File", VoreenApplication::app()->getUserDataPath(), "Raw texture file (*.raw)", FileDialogProperty::SAVE_FILE),
      saveTexButton_("saveButton", "Save Texture")
{
    addPort(texturePort_);

    saveTexButton_.onChange(MemberFunctionCallback<RawTextureSave>(this, &RawTextureSave::saveTexture));
    addProperty(textureFile_);
    addProperty(saveTexButton_);
}

RawTextureSave::~RawTextureSave() {
}

Processor* RawTextureSave::create() const {
    return new RawTextureSave();
}

void RawTextureSave::process() {
}

void RawTextureSave::initialize() {
    Processor::initialize();
}

void RawTextureSave::deinitialize() {
    Processor::deinitialize();
}

void RawTextureSave::saveTexture() {
    const Texture* texture = texturePort_.getData();
    char* texData = (char*)texture->getCpuTextureData();
    std::ofstream dataFile(textureFile_.get().c_str(), std::ios::out | std::ios::binary);
    unsigned int numBytes = texture->getWidth()*texture->getHeight();
    unsigned int numChannels = 0;
    switch (texture->getGLFormat()) {
        case GL_LUMINANCE_ALPHA: numChannels = 2; break;
        case GL_RGB: numChannels = 3; break;
        case GL_RGBA: numChannels = 4; break;
        default:
            LWARNING("Unknown texture format selected!");
    }
    numBytes *= numChannels;
    numBytes *= texture->getBpt();
    if (!dataFile.write(texData, numBytes))
        LERROR("Could not write to texture.");
    dataFile.close();
}

} // namespace
