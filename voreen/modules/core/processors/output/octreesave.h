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

#ifndef VRN_OCTREESAVE_H
#define VRN_OCTREESAVE_H

#include <string>
#include "voreen/core/processors/volumeprocessor.h"

#include "voreen/core/properties/filedialogproperty.h"
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/buttonproperty.h"
#include "voreen/core/properties/volumeinfoproperty.h"
#include "voreen/core/properties/progressproperty.h"
#include "voreen/core/properties/boundingboxproperty.h"
#include "voreen/core/properties/intproperty.h"
#include "voreen/core/datastructures/volume/volumeram.h"
#include "voreen/core/properties/vectorproperty.h"
#include "voreen/core/properties/stringproperty.h"

#include "modules/core/io/vvodvolumewriter.h"

#ifdef VRN_MODULE_HDF5
#include "modules/hdf5/io/hdf5filevolume.h"
#endif

namespace voreen {

//----------------------------------------------------------------------------------------------------------
// Classes used to write to volume files slice by slice
//----------------------------------------------------------------------------------------------------------

class VolumeSliceWriter {
public:
    /**
     * Write a slice to this volume.
     * Slices are written in z dimension from 0 to dimensions.z-1
     * The caller must make sure to call writeSlice only dimensions.z times per VolumeSliceWriter and channel.
     */
    virtual void writeSlice(const VolumeRAM* slice, size_t channel) = 0;
    virtual ~VolumeSliceWriter() { }
};

class VVDSliceWriter : public VolumeSliceWriter {
public:
    VVDSliceWriter(std::string filename, const VolumeBase* volume, const tgt::svec3& volumeLLF, const tgt::svec3& volumeURB);
    ~VVDSliceWriter() { }
    virtual void writeSlice(const VolumeRAM* slice, size_t channel);

private:
    const std::string filename_;
    const size_t numChannels_;
    const VolumeBase* volume_;
    const tgt::svec3 volumeLLF_;
    const tgt::svec3 volumeURB_;
    std::vector<std::fstream> rawout_;
};

#ifdef VRN_MODULE_HDF5
class HDF5SliceWriter : public VolumeSliceWriter {
public:
    /**
     * Construct a HDF5SliceWriter.
     * @throws tgt::IOException if the file could not be opened/created or the dataset that contains the volume could not be created.
     */
    HDF5SliceWriter(const std::string& filename, const tgt::svec3& volumeLLF, const tgt::svec3& volumeURB, const size_t numberOfChannels, const tgt::vec3& spacing, const tgt::vec3& offset, const tgt::mat4& physicalToWorldTransformation, const RealWorldMapping* rwm, size_t deflateLevel, tgt::svec3 chunkSize, bool shuffle);
    ~HDF5SliceWriter() { }
    virtual void writeSlice(const VolumeRAM* slice, size_t channel);

private:
    const std::unique_ptr<const HDF5FileVolume> fileVolume_;
    std::vector<size_t> slicesWritten_;
};
#endif



//----------------------------------------------------------------------------------------------------------
// The OctreeSave processor:
//----------------------------------------------------------------------------------------------------------


/**
 * Class to store a clipped octree as vvd/hdf5 file and apply channel shifts.
 * Alternatively stores and octree as vvod.
 */
class VRN_CORE_API OctreeSave : public VolumeProcessor {
public:
    OctreeSave();
    virtual Processor* create() const;

    virtual bool usesExpensiveComputation() const { return true; }

    virtual std::string getClassName() const  { return "OctreeSave";      }
    virtual std::string getCategory() const   { return "Output";          }
    virtual CodeState getCodeState() const    { return CODE_STATE_EXPERIMENTAL; }
    virtual bool isEndProcessor() const       { return true;              }

    // Save the volume according to properties currently available at inport_
    virtual void saveVolume();

protected:
    virtual void setDescriptions() {
        setDescription("Provides additional options to save a volume that has an octree representation to disk.");
        inport_.setDescription("The Volume to be saved. <strong>Note:</strong> An OctreeRepresentation has to be available or saving will fail.");
        filename_.setDescription("The name of the file output will be written to.");
        outputFormat_.setDescription("The Format used to save the Volume to disk. Different Formats provide different additional options:"
                "<ul>"
                    "<li><strong>VVD</strong> (Voreen Volume Data): Saves the highest resolution level of the octree, but enables usage of clipping and channel shift.</li>"
                    "<li><strong>VVOD</strong> (Voreen Volume Octree Data): Saves the complete octree, but does not provide an option to specify a clipregion or aplly channel shift. </li>"
#ifdef VRN_MODULE_HDF5
                    "<li><strong>HDF5</strong> (Hierarchical Data Format 5): Saves the highest resolution level of the octree, as a volume in a HDF5 file compatible with many other tools. "
                    "It also enables usage of clipping and channel shift, but additionally provides the option to use compression.</li>"
#endif
                "</ul>"
                );
        outputProp_.setDescription("Gives an estimate for the size of the volume on disk:"
                "<ul>"
                    "<li><strong>VVD</strong>: The value is the actual size of the .raw file(s) (of all channels combined) saved beside the .vvd files.</li>"
                    "<li><strong>VVOD</strong>: Not available."
#ifdef VRN_MODULE_HDF5
                    "<li><strong>HDF5</strong>: The value is only an estimate for the size of the .hdf5 file (as it also contains metadata), but can be "
                "significantly lower if compression is used.</li>"
#endif
                "</ul>"
                );
        enableClipping_.setDescription("Turns clipping on or off. Only available for VVD "
#ifdef VRN_MODULE_HDF5
                "and HDF5 "
#endif
                "output formats.");
        clipRegion_.setDescription("Can be used to specify the Subregion of the volume that will be saved to disk.");
        applyChannelShift_.setDescription("Turns channel shift on or off. Only available for VVD "
#ifdef VRN_MODULE_HDF5
                "and HDF5 "
#endif
                "output formats.");
        channelShift0_.setDescription("Channel shift (in voxel coordinates) to be applied to channel 1. (Also available for other channels.)");
        defaultValue_.setDescription("Value to fill the regions of the output volume that lie outside the input volume (due to channel shift). "
                "As the volume will be saved as uint16 values range from 0 (minimum) to 2<sup>16</sup>=65535 (maximum).");
        enableCompression_.setDescription("Turns compression on or off. Only available for HDF5 format."
#ifndef VRN_MODULE_HDF5
                "<i>HDF5 support has to be enabled at compile time. Enable VRN_MODULE_HDF5 during cmake configuration.</i>"
#endif
                );
#ifdef VRN_MODULE_HDF5
        compressionLevel_.setDescription("The Level of compression. Values range from 1 to 9. The Compression uses the DEFLATE algorithm after optionally applying prefilters.");
        enableShuffling_.setDescription("Turns the shuffle prefilter on or off. The filter groups bits of same significance together so that in many cases the DEFLATE algorithm will yield better results.");
        chunkSize_.setDescription("The size of the chunks the volume will be saved as. The compression pipeline operates on a single chunk of the volume at a time. Choosing larger chunk sizes yields generally in better"
                "compression rates but in longer compression times. As the volume is saved to disk slice by slice, it is generally advisable to pick smaller values for z, but larger values for x and y. By default x and "
                "y will cover the whole clip region, but z will be set to 1.");
#endif
    }

    virtual void process();
    virtual void initialize();
    virtual void deinitialize();
    virtual void invalidate(int inv = 1);

    /*
     * Update filters of filename_ according to the current selection of output format
     */
    void updateFilenameFilters();
    /*
     * Will be called whenever outputFormat_ changes
     */
    void outputFormatChanged();
    /**
     * Used to clean up the file name (e.g. remove __channel and put the correct suffix according to outputFormat_)
     * after filename_ changes.
     */
    void cleanupFilename();

    /**
     * Used to save vvd or hdf5 files.
     */
    void saveToVolume(std::string filename, const VolumeBase* volume);
    /**
     * Used to store vvod files.
     */
    void saveToOctree(std::string filename, const VolumeBase* volume);
    /**
     * Create a slice writer according to the current selection of outputFormat_.
     * @throws tgt::IOException if the writer could not be created at filename.
     */
    VolumeSliceWriter* createSliceWriter(const std::string& filename, const VolumeBase* volume, const tgt::svec3& volumeLLF, const tgt::svec3& volumeURB) const;

    /*
     * Adjust properites in groups to changed input/parameter conditions
     */
    virtual void adjustClipProperties();
    virtual void adjustShiftProperties();
    virtual void adjustCompressionProperties();
    virtual void adjustOutputProperty();
    virtual void adjustPropertiesToInput();
private:
    enum OutputFormat {
        VOL_VVOD,
        VOL_VVD,
        VOL_HDF5,
    };
    VolumePort inport_; ///< volume inport. must contain an octree representation.

    //saving
    FileDialogProperty filename_;
    VolumeInfoProperty volumeInfo_;
    StringProperty outputProp_;
    ProgressProperty progressProp_;
    ButtonProperty saveButton_;
    OptionProperty<OutputFormat> outputFormat_;

    //clipping
    BoolProperty enableClipping_;
    IntBoundingBoxProperty clipRegion_;

    //channel shift
    BoolProperty applyChannelShift_;
    IntVec3Property channelShift0_;
    IntVec3Property channelShift1_;
    IntVec3Property channelShift2_;
    IntVec3Property channelShift3_;
    IntProperty defaultValue_;

    //compression
    BoolProperty enableCompression_;
    IntProperty compressionLevel_;
    BoolProperty enableShuffling_;
    IntVec3Property chunkSize_;

    VvodVolumeWriter* vvodWriter_;

    static const std::string loggerCat_;
};

}   //namespace

#endif
