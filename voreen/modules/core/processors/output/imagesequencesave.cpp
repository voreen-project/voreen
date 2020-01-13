/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2020 University of Muenster, Germany,                        *
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

#include "imagesequencesave.h"
#include "voreen/core/utils/stringutils.h"
#include "tgt/filesystem.h"

#ifdef VRN_MODULE_DEVIL
#include <IL/il.h>
#include <IL/ilu.h>
#include "modules/devil/devilmodule.h"
#endif

namespace voreen {

const std::string ImageSequenceSave::loggerCat_("voreen.core.ImageSequenceSave");

ImageSequenceSave::ImageSequenceSave()
    : RenderProcessor()
    , inport_(Port::INPORT, "imagesequence.input", "ImageSequence Input", false)
    , outputDirectory_("outputDirectory", "Directory", "Select directory...",
        "", "", FileDialogProperty::DIRECTORY, Processor::INVALID_RESULT, Property::LOD_DEFAULT, VoreenFileWatchListener::ALWAYS_OFF)
    , baseName_("basename", "Basename")
    , saveButton_("save", "Save")
    , continousSave_("continousSave", "Save continuously", false)
{
    addPort(inport_);

    addProperty(outputDirectory_);
    addProperty(baseName_);

    saveButton_.onChange(MemberFunctionCallback<ImageSequenceSave>(this, &ImageSequenceSave::saveSequence));
    addProperty(saveButton_);
    addProperty(continousSave_);
}

Processor* ImageSequenceSave::create() const {
    return new ImageSequenceSave();
}

void ImageSequenceSave::initialize() {
    RenderProcessor::initialize();
}

void ImageSequenceSave::process() {
    if (inport_.hasChanged() && continousSave_.get())
        saveSequence();
}

void ImageSequenceSave::saveSequence() {
    if (!isInitialized())
        return;

    if (!inport_.hasData() || inport_.getData()->empty())
        return;
    if (outputDirectory_.get() == "")
        return;

    std::string directory = outputDirectory_.get();
    const ImageSequence* inputSequence = inport_.getData();
    tgtAssert(inputSequence, "no collection");
    const std::string maxNr = itos(inputSequence->size() - 1);
    for (size_t i=0; i<inputSequence->size(); i++) {
        const tgt::Texture* texture = inputSequence->at(i);
        tgtAssert(texture, "imagesequence contains null pointer");
        std::string imageFilename;

        std::string nrSuffix = itos(i);
        nrSuffix = std::string(maxNr.length()-nrSuffix.length(), '0') + nrSuffix;
        if (baseName_.get() != "")
            imageFilename += baseName_.get() + nrSuffix;
        else if (texture->getOptionalName()!= "")
            imageFilename = tgt::FileSystem::fileName(inputSequence->at(i)->getOptionalName());
        else
            imageFilename = "image" + nrSuffix;

        tgtAssert(imageFilename != "", "no imageFilename");
        if (tgt::FileSystem::fileExtension(imageFilename) == "")
            imageFilename += ".png";

        std::string outputFilename = tgt::FileSystem::cleanupPath(directory + "/" + imageFilename);
        try {
            LINFO("Saving image to file: " << outputFilename);
            saveImage(outputFilename, inputSequence->at(i));
        }
        catch (VoreenException& e) {
            LERROR("Failed to save image to file '" << outputFilename << "': " << e.what());
        }

    }
}

#ifdef VRN_MODULE_DEVIL
void ImageSequenceSave::saveImage(const std::string& filename, tgt::Texture* image) {
    tgtAssert(filename != "", "filename is empty");
    tgtAssert(tgt::FileSystem::fileExtension(filename) != "", "filename has no extension");
    tgtAssert(image, "no texture");
    if (image->getDepth() > 1)
        throw VoreenException("Passed image is a 3D texture");

    // get color buffer content
    tgt::ivec2 dim = image->getDimensions().xy();

    bool luminance = (image->getGLFormat() != GL_RGBA && image->getGLFormat() != GL_RGB);
    GLushort* colorBuffer = 0;

    if (luminance) {
        colorBuffer = reinterpret_cast<GLushort*>(image->downloadTextureToBuffer(GL_RED, GL_UNSIGNED_SHORT));
    }
    else {
        colorBuffer = reinterpret_cast<GLushort*>(image->downloadTextureToBuffer(GL_RGBA, GL_UNSIGNED_SHORT));
    }

    // create Devil image from image data and write it to file
    ILuint img;
    ilGenImages(1, &img);
    ilBindImage(img);
    // put pixels into IL-Image
    ilTexImage(dim.x, dim.y, 1, (luminance ? 1 : 4), (luminance ? IL_LUMINANCE : IL_RGBA), IL_UNSIGNED_SHORT, colorBuffer);
    ilEnable(IL_FILE_OVERWRITE);
    ilResetWrite();
    ILboolean success = ilSaveImage(const_cast<char*>(filename.c_str()));
    ilDeleteImages(1, &img);

    delete[] colorBuffer;

    if (!success) {
        throw VoreenException(DevILModule::getDevILError());
    }
}
#else
void ImageSequenceSave::saveImage(const std::string& /*filename*/, tgt::Texture* /*image*/) {
    throw VoreenException("DevIL module required");
}
#endif

}   // namespace
