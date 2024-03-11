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

#ifndef VRN_VOLUMELISTSAVE_H
#define VRN_VOLUMELISTSAVE_H

#include "voreen/core/processors/volumeprocessor.h"

#include "voreen/core/io/volumeserializerpopulator.h"
#include "voreen/core/properties/filedialogproperty.h"
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/buttonproperty.h"
#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/stringproperty.h"
#include "voreen/core/properties/progressproperty.h"
#include "voreen/core/ports/genericport.h"

#include <string>
#include <functional>

namespace voreen {

class Volume;
class VolumeSerializer;

class VRN_CORE_API VolumeListSave : public VolumeProcessor {
public:
    VolumeListSave();
    virtual ~VolumeListSave();
    virtual Processor* create() const;

    virtual std::string getClassName() const  { return "VolumeListSave";   }
    virtual std::string getCategory() const   { return "Output";           }
    virtual CodeState getCodeState() const    { return CODE_STATE_STABLE;  }
    virtual bool isEndProcessor() const       { return true;               }

protected:
    virtual void setDescriptions() {
        setDescription("Saves all volumes in the input VolumeList at once.");
        outputFormat_.setDescription("Switch between file formats as output format.");
        useOriginFileNames_.setDescription(
                "Whether or not to use the origin attribute of the volumes as an output name. "
                "This property will be disabled automatically if a volume does not have an origin name or "
                "if multiple volumes share the save origin name.");
        baseName_.setDescription(
                "Base name for volumes to be saved by this processor. The whole name will be the basename concatenated"
                "with the position of the volume in the list (including leading zeros if needed)."
                "Only available if origin volume names are not used."
                );
    }

    virtual void process();
    virtual void initialize();

    /**
     * Save the volume list using the method specified by outputFormat_.
     */
    void saveVolumeList();
    /**
     * Save a volume to a file.
     * @param volumeName identifier for the given volume. Has to be unique in the list of all volumes to be saved.
     * @param volumeList list containing all volumes
     * @param i number of volume to save this call. volumeList->at(i) is the volume to be saved.
     */
    void saveVolume(const std::string& volumeName, const VolumeList* volumeList, size_t i);
#ifdef VRN_MODULE_HDF5
    /**
     * Save a volume to a hdf5 file.
     * @param volumeName identifier for the given volume. Has to be unique in the list of all volumes to be saved.
     * @param volumeList list containing all volumes
     * @param i number of volume to save this call. volumeList->at(i) is the volume to be saved.
     */
    void saveVolumeHDF5(const std::string& volumeName, const VolumeList* volumeList, size_t i);
#endif
    /**
     * Save all volumes of the input list using the provided save function write.
     * @param inputList list of volumes to be saved.
     * @param write Function to save a single volume. Subject to the same restrictions as saveVolumeVVD or saveVolumeHDF5.
     */
    void saveVolumes(const VolumeList* inputList, std::function<void (const std::string&, const VolumeList*, size_t)> write);

    /**
     * Adjust other properties etc. to the current output format.
     */
    void adaptToOutputFormat();
    /**
     * Adjust properties to the current output format.
     */
    void adaptToInput();
    /**
     * Adjust properties to the state of useOriginFileNames_
     */
    void adaptToUseOriginFileNames();

private:

    enum OutputFormat {
        VOL_ALL,
        VOL_HDF5,
    };

    void adjustCompressionProperties();
    std::deque<Option<OutputFormat>> constructFormats() const;

    VolumeListPort inport_;

    // Format to save volumes as.
    OptionProperty<OutputFormat> outputFormat_;
    // Output file for HDF5 format.
    FileDialogProperty fileNameHDF5_;
    // Output folder for VVD format.
    FileDialogProperty folderNameAll_;
    // Whether or not to use the origin attribute of the volumes as an output name
    // This property will be disabled automatically if a volume does not have an origin name
    // or if there are duplicates.
    BoolProperty useOriginFileNames_;
    // Use this as the base name for volumes. The volume names will be baseName_.get() + number.
    StringProperty baseName_;
    // Save all volumes in the list.
    ButtonProperty saveButton_;
    // Progress for saving the whole list. One step per volume.
    ProgressProperty progressProp_;
    // Save whenever the input changes.
    BoolProperty continousSave_;

    //compression
    BoolProperty enableCompression_;
    IntProperty compressionLevel_;
    BoolProperty enableShuffling_;

    // Used to serialize volumes to vvd.
    VolumeSerializerPopulator volumeSerializerPopulator_;

    static const std::string loggerCat_; ///< category used in logging
};

} //namespace

#endif
