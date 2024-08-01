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

#include "volumebricksave.h"

#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/datastructures/volume/volumeram.h"
#include "voreen/core/datastructures/volume/volumedisk.h"
#include "voreen/core/datastructures/octree/volumeoctree.h"
#include "voreen/core/voreenapplication.h"
#include "voreen/core/datastructures/callback/lambdacallback.h"

#include "modules/hdf5/io/hdf5volumewriter.h"

#include "tgt/filesystem.h"

#include <memory>

namespace voreen {

const std::string VolumeBrickSave::loggerCat_("voreen.bigdataimageprocessingextra.volumebricksave");

VolumeBrickSave::VolumeBrickSave()
    : VolumeProcessor()
    , inport_(Port::INPORT, "connectedcomponentanalysis.inport", "Volume Input")
    , volumeInfo_("volumeInfo","Volume info")
    , volumeFilePath_("volumeFilePath", "Volume file output", "Path", "", "HDF5 (*.h5)", FileDialogProperty::SAVE_FILE, Processor::INVALID_RESULT, Property::LOD_DEFAULT)
    , allowTruncateFile_("allowTruncateFile_", "Allow truncation", false)
    , volumeDimensions_("volumeDimensions", "Volume dimensions", tgt::ivec3(1), tgt::ivec3(1), tgt::ivec3(std::numeric_limits<int>::max()))
    , brickDimensions_("brickDimensions", "Brick dimensions", tgt::ivec3(1), tgt::ivec3(1), tgt::ivec3(std::numeric_limits<int>::max()))
    , brickOffset_("brickOffset", "Brick offset", tgt::ivec3(0), tgt::ivec3(0), tgt::ivec3(std::numeric_limits<int>::max()))
    , saveButton_("saveButton", "Save brick", Processor::INVALID_RESULT, Property::LOD_DEFAULT)
{
    addPort(inport_);

    addProperty(volumeInfo_);
    addProperty(volumeFilePath_);

    addProperty(allowTruncateFile_);
    addProperty(volumeDimensions_);
    ON_CHANGE_LAMBDA(volumeDimensions_, [this] () {
            brickOffset_.setMaxValue(volumeDimensions_.get());
            });
    addProperty(brickDimensions_);
    addProperty(brickOffset_);

    addProperty(saveButton_);
        ON_CHANGE(saveButton_, VolumeBrickSave, saveBrick);
}

VolumeBrickSave::~VolumeBrickSave() {
}
VoreenSerializableObject* VolumeBrickSave::create() const {
    return new VolumeBrickSave();
}

void VolumeBrickSave::initialize() {
}
void VolumeBrickSave::deinitialize() {
}

void VolumeBrickSave::process() {
    if(inport_.hasChanged()) {
        volumeInfo_.setVolume(inport_.getData());
    }
}

void VolumeBrickSave::saveBrick() {
    const VolumeBase* input = inport_.getData();
    if(!input) {
        return;
    }

    std::unique_ptr<HDF5FileVolume> currentVolume(nullptr);

    const std::string volumeLocation = HDF5VolumeWriter::VOLUME_DATASET_NAME;

    // First try to load existing file
    if(tgt::FileSystem::fileExists(volumeFilePath_.get())) {
        try {
            currentVolume = HDF5FileVolume::openVolume(volumeFilePath_.get(), volumeLocation, false /*readonly*/);
        } catch(tgt::IOException& e) {
            if(allowTruncateFile_.get()) {
                currentVolume = nullptr;
            } else {
                const std::string errorMsg("Could not open output file, but it exists. Allow truncation to overwrite.");
                LERROR(errorMsg);
                VoreenApplication::app()->showMessageBox("Cannot write brick", errorMsg);
                return;
            }
        }
    }

    // Do the property states and input brick fit the chunk size of the output volume?
    if(currentVolume && (tgt::svec3(brickDimensions_.get()) % currentVolume->getChunkSize() != tgt::svec3::zero)) {
        // No? Then...
        if(allowTruncateFile_.get()) {
            // ... if we are allowed, create a new volume later
            currentVolume = nullptr;
        } else {
            // ... if not: Issue a warning.
            const std::string errorMsg("Brick size is not a multiple output volume chunk size. Consider allowing truncation.");
            LWARNING(errorMsg);
        }
    }
    if(currentVolume && (tgt::svec3(brickOffset_.get()) % currentVolume->getChunkSize() != tgt::svec3::zero)) {
        // No? Then...
        if(allowTruncateFile_.get()) {
            // ... if we are allowed, create a new volume later
            currentVolume = nullptr;
        } else {
            // ... if not: Issue a warning.
            const std::string errorMsg("Brick offset is not a multiple output volume chunk size. That might impact performance.");
            LWARNING(errorMsg);
        }
    }

    // Is it possible at all to write the chunk into the output volume?
    if(currentVolume && (tgt::ivec3(currentVolume->getDimensions()) != volumeDimensions_.get()
            || currentVolume->getBaseType() != input->getBaseType()))
    {
        // No? Then...
        if(allowTruncateFile_.get()) {
            // ... if we are allowed, create a new volume later
            currentVolume = nullptr;
        } else {
            // ... if not: Put out an error message.
            const std::string errorMsg("Volume dimensions, and/or format of output file do not fit. Consider allowing truncation.");
            LERROR(errorMsg);
            VoreenApplication::app()->showMessageBox("Cannot write brick", errorMsg);
            return;
        }
    }

    // Okay, apparently we need to create a new file
    if(!currentVolume) {
        LINFO("Creating new output file.");
        try {
            const tgt::svec3 chunkSize(brickDimensions_.get());
            currentVolume = HDF5FileVolume::createVolume(volumeFilePath_.get(), volumeLocation, input->getBaseType(), volumeDimensions_.get(), 1, true, 0 /* deflate level */, chunkSize, false /* shuffle */);

            // Write meta data that is available now
            currentVolume->writeSpacing(input->getSpacing());
            if(input->hasMetaData("RealWorldMapping")) {
                currentVolume->writeRealWorldMapping(input->getRealWorldMapping());
            }

        } catch(tgt::IOException& e) {
            VoreenApplication::app()->showMessageBox("Could not create Output Volume",
                    std::string("Could not create Output Volume:\n\n") + e.what(),
                    true);
            currentVolume = nullptr;
            return;
        }
    }

    tgt::svec3 offset = brickOffset_.get();
    if(tgt::hor(tgt::greaterThan(offset + input->getDimensions(), currentVolume->getDimensions()))) {
        std::stringstream errorMsg;
        errorMsg << "Input with dimensions " << input->getDimensions() << " and offset " << offset << "does not fit into volume of dimensions " << currentVolume->getDimensions() << "!";
        LERROR(errorMsg.str());
        VoreenApplication::app()->showMessageBox("Cannot write brick", errorMsg.str());
        return;
    }

    LINFO("Writing volume brick.");
    currentVolume->writeBrick(input->getRepresentation<VolumeRAM>(), offset);

    // Disallow truncation now so that the user does not accidentally overwrite the volume
    allowTruncateFile_.set(false);
}

} // namespace voreen
