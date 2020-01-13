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

#include "voreen/core/io/volumeserializerpopulator.h"

#include "voreen/core/voreenapplication.h"
#include "voreen/core/voreenmodule.h"
#include "voreen/core/io/volumeserializer.h"
#include "voreen/core/io/volumereader.h"
#include "voreen/core/io/volumewriter.h"

#include "voreen/core/io/progressbar.h"

namespace voreen {

VolumeSerializerPopulator::VolumeSerializerPopulator(ProgressBar* progressBar)
    : vs_(new VolumeSerializer()),
      progressBar_(progressBar)
{

    if (VoreenApplication::app()) {
        // retrieve volume readers/writers from modules
        std::vector<VoreenModule*> modules = VoreenApplication::app()->getModules();
        for (size_t i=0; i<modules.size(); i++) {
            for (size_t j=0; j<modules.at(i)->getRegisteredVolumeReaders().size(); j++) {
                readers_.push_back(modules.at(i)->getRegisteredVolumeReaders().at(j)->create(progressBar_));
            }
            for (size_t j=0; j<modules.at(i)->getRegisteredVolumeWriters().size(); j++) {
                writers_.push_back(modules.at(i)->getRegisteredVolumeWriters().at(j)->create(progressBar_));
            }
        }
    }
    else {
        LWARNINGC("voreen.VolumeSerializerPopulator", "VoreenApplication not instantiated");
    }

    // populate array with all known VolumeReaders
    for (size_t i = 0; i < readers_.size(); ++i) {
        vs_->registerReader(readers_[i]);
    }

    // populate array with all known VolumeWriters
    for (size_t i = 0; i < writers_.size(); ++i) {
        vs_->registerWriter(writers_[i]);
    }
}

VolumeSerializerPopulator::~VolumeSerializerPopulator() {
    delete vs_;

    for (size_t i = 0; i < readers_.size(); ++i)
        if (readers_[i])
            delete readers_[i];

    for (size_t i = 0; i < writers_.size(); ++i)
        delete writers_[i];
}

VolumeSerializer* VolumeSerializerPopulator::getVolumeSerializer() const {
    return vs_;
}

std::vector<std::string> VolumeSerializerPopulator::getSupportedReadExtensions() const {
    std::vector<std::string> extensions;
    for (size_t i=0; i<readers_.size(); ++i) {
        for (size_t j=0; j<readers_[i]->getSupportedExtensions().size(); ++j) {
            if (std::find(extensions.begin(), extensions.end(), readers_[i]->getSupportedExtensions()[j]) == extensions.end())
                extensions.push_back(readers_[i]->getSupportedExtensions()[j]);
        }
    }
    return extensions;
}

std::vector<std::string> VolumeSerializerPopulator::getSupportedReadFilenames() const {
    std::vector<std::string> filenames;
    for (size_t i=0; i<readers_.size(); ++i) {
        for (size_t j=0; j<readers_[i]->getSupportedFilenames().size(); ++j) {
            if (std::find(filenames.begin(), filenames.end(), readers_[i]->getSupportedFilenames()[j]) == filenames.end())
                filenames.push_back(readers_[i]->getSupportedFilenames()[j]);
        }
    }
    return filenames;
}

std::vector<std::string> VolumeSerializerPopulator::getSupportedWriteExtensions() const {
    std::vector<std::string> extensions;
    for (size_t i=0; i<writers_.size(); ++i) {
        for (size_t j=0; j<writers_[i]->getSupportedExtensions().size(); ++j) {
            if (std::find(extensions.begin(), extensions.end(), writers_[i]->getSupportedExtensions()[j]) == extensions.end())
                extensions.push_back(writers_[i]->getSupportedExtensions()[j]);
        }
    }
    return extensions;
}

std::vector<std::string> VolumeSerializerPopulator::getSupportedWriteFilenames() const {
    std::vector<std::string> filenames;
    for (size_t i=0; i<writers_.size(); ++i) {
        for (size_t j=0; j<writers_[i]->getSupportedFilenames().size(); ++j) {
            if (std::find(filenames.begin(), filenames.end(), writers_[i]->getSupportedFilenames()[j]) == filenames.end())
                filenames.push_back(writers_[i]->getSupportedFilenames()[j]);
        }
    }
    return filenames;
}

} // namespace voreen
