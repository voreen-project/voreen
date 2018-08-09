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

#include "imagesequencesource.h"

#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/voreenapplication.h"
#include "voreen/core/datastructures/imagesequence.h"
#include "voreen/core/io/progressbar.h"

#include "tgt/filesystem.h"
#include "tgt/texturemanager.h"
#include "tgt/gpucapabilities.h"

namespace voreen {

const std::string ImageSequenceSource::loggerCat_("voreen.core.ImageSequenceSource");

ImageSequenceSource::ImageSequenceSource()
    : RenderProcessor(),
      outport_(Port::OUTPORT, "imagesequence.out", "ImageSequence Output"),
      imageDirectory_("imageDirectory","Image Directory", "Select Image Directory",
          "", "", FileDialogProperty::DIRECTORY),
      textureFiltering_("textureFiltering", "Enable Texture Filtering", true),
      showProgressBar_("showProgressBar", "Show Progress Bar", true),
      reloadSequence_("reloadSequence", "Reload"),
      clearSequence_("clearSequence", "Clear Sequence"),
      numImages_("numImages", "Num Images", 0, 0, 10000, VALID),
      imageSequence_(0),
      currentDir_(""),
      sequenceOwner_(false),
      forceReload_(false)
{
    addPort(outport_);

    textureFiltering_.onChange(MemberFunctionCallback<ImageSequenceSource>(this, &ImageSequenceSource::forceReload));
    reloadSequence_.onClick(MemberFunctionCallback<ImageSequenceSource>(this, &ImageSequenceSource::forceReload));
    clearSequence_.onClick(MemberFunctionCallback<ImageSequenceSource>(this, &ImageSequenceSource::unsetDirectoryName));
    numImages_.setReadOnlyFlag(true);

    addProperty(imageDirectory_);
    addProperty(textureFiltering_);
    addProperty(showProgressBar_);
    addProperty(reloadSequence_);
    addProperty(clearSequence_);
    addProperty(numImages_);
}

ImageSequenceSource::~ImageSequenceSource() {
}

Processor* ImageSequenceSource::create() const {
    return new ImageSequenceSource();
}

void ImageSequenceSource::process() {
    if (forceReload_) {
        loadImageSequence(imageDirectory_.get());
    }
}

void ImageSequenceSource::initialize() {

    RenderProcessor::initialize();

    imageSequence_ = new ImageSequence();
    sequenceOwner_ = true;

    if (!imageDirectory_.get().empty()) {
        //currentDir_ = imageDirectory_.get();
        try {
            loadImageSequence(imageDirectory_.get());
        }
        catch (std::exception& e) {
            LERROR("Failed to load image sequence: " << e.what());
            numImages_.set(0);
        }
    }
    else {
        numImages_.set(0);
    }
}

void ImageSequenceSource::deinitialize() {
    outport_.setData(0);
    clearSequence();
    delete imageSequence_;
    imageSequence_ = 0;

    RenderProcessor::deinitialize();
}

void ImageSequenceSource::invalidate(int inv) {

    RenderProcessor::invalidate(inv);

    if (inv == Processor::VALID)
        return;

    if (!isInitialized())
        return;

    if (imageDirectory_.get() != currentDir_) {
        forceReload_ = true;
    }
}

void ImageSequenceSource::loadImageSequence(const std::string& d) {

    // important: d might be cleared by clearSequence
    std::string dir(d);

    forceReload_ = false;

    if (!imageSequence_) {
        LERROR("No image sequence present");
        return;
    }

    // clear current sequence
    clearSequence();

    if (dir.empty())
        return;

    LINFO("Loading images from directory " << dir << " ...");
    currentDir_ = dir;
    imageDirectory_.set(dir);

    // load images as textures and collect them in an image sequence
    std::vector<std::string> filenames = tgt::FileSystem::readDirectory(dir, true, false);

    // create progress bar
    ProgressBar* progressDialog = 0;
    if (showProgressBar_.get() && !filenames.empty()) {
        progressDialog = VoreenApplication::app()->createProgressDialog();
        if (progressDialog) {
            progressDialog->setTitle("Loading Images");
            progressDialog->show();
            progressDialog->setProgress(0.f);
            progressDialog->forceUpdate();
        }
    }

    tgt::Texture::Filter filterMode = textureFiltering_.get() ? tgt::Texture::LINEAR : tgt::Texture::NEAREST;
    for (size_t i=0; i<filenames.size(); ++i) {
        if (progressDialog) {
            progressDialog->setProgressMessage("Loading " + filenames[i] + " ...");
            progressDialog->setProgress(static_cast<float>(i) / static_cast<float>(filenames.size()));
        }
        LINFO("Loading image " << filenames[i] << " ...");
        tgt::Texture* texture = TexMgr.load(dir + "/" + filenames[i], filterMode,
            false, false, true, false);
        if (texture)
            imageSequence_->add(texture);
        else
            LWARNING("Failed to load image: " << filenames[i]);
    }
    LGL_ERROR;

    // clear progress
    if (progressDialog) {
        progressDialog->hide();
        delete progressDialog;
        progressDialog = 0;
    }

    // output sequence
    outport_.setData(imageSequence_, false);

    LINFO("Successfully loaded " << imageSequence_->size() << " images.");
    numImages_.set(static_cast<int>(imageSequence_->size()));
    LGL_ERROR;
}


void ImageSequenceSource::reloadImageSequence() {
    if (!currentDir_.empty())
        forceReload();
}

void ImageSequenceSource::clearSequence() {

    if (sequenceOwner_ && imageSequence_) {
        if (!imageSequence_->empty()) {
            LDEBUG("Clearing sequence");
            for (size_t i=0; i<imageSequence_->size(); ++i) {
                delete imageSequence_->at(i);
            }
            imageSequence_->clear();
            LGL_ERROR;
        }
    }
    else {
        imageSequence_ = new ImageSequence();
        sequenceOwner_ = true;
    }
    currentDir_ = "";
    imageDirectory_.set("");
    numImages_.set(0);
}

void ImageSequenceSource::setImageSequence(ImageSequence* sequence) {

    clearSequence();  // now owner of the sequence -> delete it before assigning
    delete imageSequence_;

    imageSequence_ = sequence;
    sequenceOwner_ = false;

    outport_.setData(imageSequence_, false);
    outport_.invalidatePort();
    invalidate();

    numImages_.set(sequence ? static_cast<int>(sequence->size()) : 0);
}

const ImageSequence* ImageSequenceSource::getImageSequence() const {
    return imageSequence_;
}

void ImageSequenceSource::forceReload() {
    forceReload_ = true;
    invalidate();
}

void ImageSequenceSource::unsetDirectoryName() {
    imageDirectory_.set("");
}

} // namespace
