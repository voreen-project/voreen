/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2021 University of Muenster, Germany,                        *
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

#include "volumelistsave.h"

#include "voreen/core/io/volumeserializer.h"
#include "voreen/core/io/volumeserializerpopulator.h"
#include "voreen/core/utils/stringutils.h"

#include "tgt/filesystem.h"

#ifdef VRN_MODULE_HDF5
#include "modules/hdf5/io/hdf5volumewriter.h"
#endif

#include <set>

namespace voreen {

const std::string VolumeListSave::loggerCat_("voreen.core.VolumeListSave");

VolumeListSave::VolumeListSave()
    : VolumeProcessor()
    , inport_(Port::INPORT, "volumehandle.input", "VolumeList Input", false)
    , outputFormat_("outputFormat", "Output Format")
    , fileNameHDF5_("fileNameHDF5", "Save as", "Select output file...",
        "", "Hierarchical Data Format 5 (*.h5 *.hdf5)", FileDialogProperty::SAVE_FILE)
    , folderNameAll_("folderNameVVD", "Save as", "Select output file...",
        "", "Voreen Volume Data (*.vvd)", FileDialogProperty::DIRECTORY, Processor::INVALID_RESULT, Property::LOD_DEFAULT, VoreenFileWatchListener::ALWAYS_OFF)
    , saveButton_("save", "Save")
    , progressProp_("progress", "Progress")
    , useOriginFileNames_("useOriginFileNames_", "Use origin volume names", true)
    , baseName_("baseName","Base name","volume_",Processor::INVALID_RESULT,Property::LOD_DEFAULT)
    , continousSave_("continousSave", "Save continuously", false)
    , enableCompression_("enableCompression", "Enable Compression", false, Processor::INVALID_PROGRAM, Property::LOD_ADVANCED)
    , compressionLevel_("compressionLevel", "Deflate Level", 1, 1, 9, Processor::INVALID_RESULT, NumericProperty<int>::STATIC, Property::LOD_ADVANCED)
    , enableShuffling_("enableShuffling", "Enable Shuffling", false, Processor::INVALID_RESULT, Property::LOD_ADVANCED)
{
    addPort(inport_);

    addProperty(outputFormat_);
        outputFormat_.setOptions(constructFormats());
        ON_CHANGE(outputFormat_, VolumeListSave, adaptToOutputFormat);

    addProperty(fileNameHDF5_);
    addProperty(folderNameAll_);

    addProperty(saveButton_);
        ON_CHANGE(saveButton_, VolumeListSave, saveVolumeList);

    addProperty(progressProp_);
        addProgressBar(&progressProp_);

    addProperty(useOriginFileNames_);
        ON_CHANGE(useOriginFileNames_, VolumeListSave, adaptToUseOriginFileNames);

    addProperty(baseName_);

    addProperty(continousSave_);

    addProperty(enableCompression_);
        enableCompression_.setGroupID("compression");
        ON_CHANGE(enableCompression_, VolumeListSave, adjustCompressionProperties);
    addProperty(compressionLevel_);
        compressionLevel_.setGroupID("compression");
    addProperty(enableShuffling_);
        enableShuffling_.setGroupID("compression");
    setPropertyGroupGuiName("compression", "Compression");
    adjustCompressionProperties();
}

VolumeListSave::~VolumeListSave() {
}

Processor* VolumeListSave::create() const {
    return new VolumeListSave();
}

void VolumeListSave::initialize() {
    VolumeProcessor::initialize();

    // Adapt to changed properties after deserialization
    adaptToOutputFormat();
    adaptToUseOriginFileNames();
}

void VolumeListSave::process() {
    if(!isInitialized()) {
        return;
    }
    // New data => no save progress
    progressProp_.setProgress(0);
    if (inport_.hasChanged()) {
        // Adapt UI to changed data
        adaptToInput();

        // Save if continousSave_ is enabled
        if(continousSave_.get()) {
            saveVolumeList();
        }
    }
}

void VolumeListSave::saveVolumeList() {
    if (!isInitialized())
        return;

    const VolumeList* inputList = inport_.getData();
    // If there is nothing to save, we are done.
    if (!inputList || inputList->empty())
        return;

    switch(outputFormat_.getValue()) {

        case VOL_ALL:
            // Save to vvd
            if (folderNameAll_.get().empty()) {
                LWARNING("No output folder specified.");
                return;
            }
            saveVolumes(inputList, std::bind(&VolumeListSave::saveVolume, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3));
            break;

#ifdef VRN_MODULE_HDF5
        case VOL_HDF5:
            //Save to hdf5
            if (fileNameHDF5_.get().empty()) {
                LWARNING("No output file specified.");
                return;
            }
            saveVolumes(inputList, std::bind(&VolumeListSave::saveVolumeHDF5, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3));
            break;
#endif

        default:
            // We should not get here.
            tgtAssert(false, "Invalid output format.");
    }
}
void VolumeListSave::saveVolume(const std::string& volumeName, const VolumeList* volumeList, size_t i) {
    tgtAssert(!volumeName.empty(), "no volumeName");
    const VolumeSerializer* serializer = volumeSerializerPopulator_.getVolumeSerializer();

    // Change extension to vvd
    std::string basename = tgt::FileSystem::fullBaseName(tgt::FileSystem::cleanupPath(volumeName));
    std::string fileName = basename + "." + outputFormat_.getKey();

    try {
        serializer->write(folderNameAll_.get() + "/" + fileName, volumeList->at(i));
    }
    catch (tgt::FileException& e) {
        LERROR("Failed to save volume to file '" << fileName << "': " << e.what());
    }
}

#ifdef VRN_MODULE_HDF5
void VolumeListSave::saveVolumeHDF5(const std::string& volumeName, const VolumeList* volumeList, size_t i) {
    tgtAssert(!volumeName.empty(), "no volumeName");

    const VolumeBase* vol = volumeList->at(i);

    // Set chunksize to a single z slice.
    const tgt::svec3 chunkSize(vol->getDimensions().xy(), 1);

    // Set compressionlevel, if enabled.
    int compressionLevel = enableCompression_.get() ? compressionLevel_.get() : 0;

    try {
        // Write to HDF5 file, but only truncate for the first volume (i == 0)!
        HDF5VolumeWriter().write(fileNameHDF5_.get(), vol, volumeName, (i == 0), compressionLevel, chunkSize, enableShuffling_.get());
    }
    catch (tgt::IOException& e) {
        LERROR("Failed to save volume " << volumeName << " to file '" << fileNameHDF5_.get() << "': " << e.what());
    }
}
#endif

void VolumeListSave::saveVolumes(const VolumeList* inputList, std::function<void (const std::string&, const VolumeList*, size_t)> write) {
    tgtAssert(inputList, "No input");

    const size_t listSize = inputList->size();
    const size_t nrLength = std::to_string(listSize - 1).size();

    for (size_t i=0; i<listSize; i++) {
        // Update progress
        progressProp_.setProgress(static_cast<float>(i)/listSize);

        std::string volumeName;
        if(useOriginFileNames_.get()) {
            // Derive name from volume origin
            volumeName = tgt::FileSystem::fileName(static_cast<Volume*>(inputList->at(i))->getOrigin().getFilename());
            // We checked if all volumes hat an origin on changed inport, so volumeName should now be != ""
        } else  {
            // Enumerate volumes
            std::ostringstream ss;
            ss << std::setw(nrLength) << std::setfill('0') << i;
            std::string nrSuffix(ss.str());
            volumeName = baseName_.get() + nrSuffix;
        }
        tgtAssert(!volumeName.empty(), "No file name.");

        write(volumeName, inputList, i);
    }
    // We are done, set progress to 100%
    progressProp_.setProgress(1);
}

void VolumeListSave::adaptToOutputFormat() {
    // Only show the filedialogproperty suitable for current output format
    switch(outputFormat_.getValue()) {
        case VOL_ALL:
            fileNameHDF5_.setVisibleFlag(false);
            folderNameAll_.setVisibleFlag(true);
            break;
#ifdef VRN_MODULE_HDF5
        case VOL_HDF5:
            fileNameHDF5_.setVisibleFlag(true);
            folderNameAll_.setVisibleFlag(false);
            break;
#endif
        default:
            tgtAssert(false, "Invalid output format.");
    }

    adjustCompressionProperties();
}

void VolumeListSave::adaptToInput() {
    const VolumeList* input = inport_.getData();
    if (!input || input->empty()) {
        return;
    }
    std::set<std::string> volumeNames;

    // Check if there are duplicate file names in input
    for (size_t i=0; i<input->size(); i++) {
        std::string volFilename = tgt::FileSystem::fileName(static_cast<Volume*>(input->at(i))->getOrigin().getFilename());
        if (!volFilename.empty() && volumeNames.find(volFilename) == volumeNames.end()) {
            volumeNames.insert(volFilename);
        } else {
            // Duplicate found
            useOriginFileNames_.set(false);
            useOriginFileNames_.setReadOnlyFlag(true);
            return;
        }
    }
    // No duplicates found
    useOriginFileNames_.setReadOnlyFlag(false);
}
void VolumeListSave::adaptToUseOriginFileNames() {
    // Show the basename property iff useOriginFileNames_ is disabled.
    baseName_.setVisibleFlag(!useOriginFileNames_.get());
}

void VolumeListSave::adjustCompressionProperties() {
    // Adjust commpression properties
    bool visible = enableCompression_.get();
    compressionLevel_.setVisibleFlag(visible);
    enableShuffling_.setVisibleFlag(visible);

    bool enabled = (outputFormat_.getValue() == VOL_HDF5);
    enableCompression_.setReadOnlyFlag(!enabled);
    compressionLevel_.setReadOnlyFlag(!enabled);
    enableShuffling_.setReadOnlyFlag(!enabled);
}

std::deque<Option<VolumeListSave::OutputFormat>> VolumeListSave::constructFormats() const {
    std::deque<Option<OutputFormat>> formats;

    const std::vector<VolumeWriter*> volumeWriters = volumeSerializerPopulator_.getVolumeSerializer()->getWriters();
    std::string vvdWriterFilter;
    for (size_t i=0; i<volumeWriters.size(); i++) {
        VolumeWriter* curWriter = volumeWriters.at(i);

        // extensions
        std::vector<std::string> extensionVec = curWriter->getSupportedExtensions();
        for (size_t j=0; j<extensionVec.size(); j++) {
            std::string extension = extensionVec.at(j);
            std::string filterStr = curWriter->getFormatDescription() + " (*." + extension + ")";
            if (extension == "vvd")
                vvdWriterFilter = filterStr;
            else
#ifdef VRN_MODULE_HDF5
                if (extension != "h5")
#endif
                    formats.push_back({extension, filterStr, VOL_ALL});
        }
    }

    if (vvdWriterFilter != "")
        formats.insert(formats.begin(), {"vvd", "VVD (Voreen Volume Data)", VOL_ALL});

#ifdef VRN_MODULE_HDF5
    formats.insert(formats.begin(), {"hdf5", "HDF5 (Hierarchical Data Format 5)", VOL_HDF5});
#endif

    return formats;
}

} // namespace voreen
