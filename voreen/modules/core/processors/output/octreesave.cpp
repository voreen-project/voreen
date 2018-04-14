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

#include "octreesave.h"

#include "voreen/core/io/volumeserializer.h"
#include "voreen/core/io/volumeserializerpopulator.h"
#include "voreen/core/io/volumewriter.h"
#include "voreen/core/datastructures/octree/volumeoctree.h"
#include "voreen/core/datastructures/volume/volumefactory.h"
#include "voreen/core/datastructures/volume/volumederiveddata.h"
#include "voreen/core/datastructures/volume/volumehash.h"
#include "voreen/core/datastructures/volume/volumepreview.h"
#include "voreen/core/utils/stringutils.h"

#include "modules/core/io/vvdformat.h"

#include "tgt/filesystem.h"

namespace { // anonymous namespace for helper functions
tgt::svec3 getDimensions(const tgt::svec3& llf, const tgt::svec3 urb) {
    return urb - llf + tgt::svec3::one;
}
} // anonymous namespace for helper functions

namespace voreen {

//----------------------------------------------------------------------------------------------------------
// Classes used to write to volume files slice by slice
//----------------------------------------------------------------------------------------------------------

VVDSliceWriter::VVDSliceWriter(std::string filename, const VolumeBase* volume, const tgt::svec3& volumeLLF, const tgt::svec3& volumeURB)
    : filename_(filename)
    , numChannels_(volume->getNumChannels())
    , volume_(volume)
    , volumeLLF_(volumeLLF)
    , volumeURB_(volumeURB)
    , rawout_(volume->getNumChannels())
{
    tgt::svec3 dim = getDimensions(volumeLLF_, volumeURB_);
    tgtAssert(dim.x > 0 && dim.y > 0 && dim.z > 0, "Dimensions must all be > 0");
    tgtAssert(volume_, "Volume must not be 0.");
    tgtAssert(volume_->getNumChannels() >= 1 && volume_->getNumChannels() <= 4, "Invalid number of channels");

    std::string basefilename =  tgt::FileSystem::fullBaseName(tgt::FileSystem::cleanupPath(filename));
    //init filename arrays
    std::vector<std::string> vvdname(volume->getNumChannels());
    std::vector<std::string> rawname(volume->getNumChannels());

    //
    // VVD: ---------------------------
    //
    for(size_t channel = 0; channel < volume_->getNumChannels(); channel++) {
        //create names
        vvdname[channel] = basefilename + "__Channel" + std::to_string(channel+1) + ".vvd";
        rawname[channel] = basefilename + "__Channel" + std::to_string(channel+1) + ".raw";

        XmlSerializer s(vvdname[channel]);
        s.setUseAttributes(true);

        // TODO: removing the tgt::FileSystem::fileName fixed a crash, however, this should be checked again
        //VvdObject o = VvdObject(volume, tgt::FileSystem::fileName(rawname[channel]));
        VvdObject o(volume_, rawname[channel]);

        o.representationData_->dimensions_ = tgt::ivec3(dim);
        o.representationData_->format_ = "uint16";
        //remove derived data since they are wrong for the clipped volume
        o.derivedData_.erase(volume_->getDerivedData<VolumeHash>());
        o.derivedData_.erase(volume_->getDerivedData<VolumeMinMax>());
        o.derivedData_.erase(volume_->getDerivedData<VolumePreview>());
        o.derivedData_.erase(volume_->getDerivedData<VolumeHistogramIntensity>());
        //update offset
        tgt::vec3 tmpSpacing = (o.metaData_.hasMetaData(VolumeBase::META_DATA_NAME_SPACING) ? static_cast<Vec3MetaData*>(o.metaData_.getMetaData(VolumeBase::META_DATA_NAME_SPACING))->getValue() : tgt::vec3::one);
        tgt::vec3 tmpOffset = (o.metaData_.hasMetaData(VolumeBase::META_DATA_NAME_OFFSET) ? static_cast<Vec3MetaData*>(o.metaData_.getMetaData(VolumeBase::META_DATA_NAME_OFFSET))->getValue() : tgt::vec3::zero);
        o.metaData_.addMetaData(VolumeBase::META_DATA_NAME_OFFSET,new Vec3MetaData(tmpOffset+tgt::vec3(volumeLLF_)*tmpSpacing));

        std::vector<VvdObject> vec;
        vec.push_back(o);

        Serializer serializer(s);
        serializer.serialize("Volumes", vec, "Volume");

        // write serialization data to temporary string stream
        std::ostringstream textStream;
        try {
            s.write(textStream);
            if (textStream.fail())
                throw tgt::IOException("Failed to write serialization data to string stream.");
        }
        catch (std::exception& e) {
            throw tgt::IOException("Failed to write serialization data to string stream: " + std::string(e.what()));
        }
        catch (...) {
            throw tgt::IOException("Failed to write serialization data to string stream (unknown exception).");
        }

        // Now we have a valid StringStream containing the serialization data.
        // => Open output file and write it to the file.
        std::fstream fileStream(vvdname[channel].c_str(), std::ios_base::out);
        if (fileStream.fail())
            throw tgt::IOException("Failed to open file '" + vvdname[channel] + "' for writing.");

        try {
            fileStream << textStream.str();
        }
        catch (std::exception& e) {
            throw tgt::IOException("Failed to write serialization data stream to file '"
                                         + vvdname[channel] + "': " + std::string(e.what()));
        }
        catch (...) {
            throw tgt::IOException("Failed to write serialization data stream to file '"
                                         + vvdname[channel] + "' (unknown exception).");
        }
        fileStream.close();
    }
    //open streams
    for(size_t channel = 0; channel < volume_->getNumChannels(); channel++) {
        rawout_[channel].open(rawname[channel].c_str(), std::ios::out | std::ios::binary);
        if (!rawout_[channel].is_open() || rawout_[channel].bad()) {
            throw tgt::IOException("Error opening rawout");
        }
    }
}
void VVDSliceWriter::writeSlice(const VolumeRAM* slice, size_t channel) {
    tgtAssert(slice, "Slice must not be 0.");
    tgtAssert(channel >= 0 && channel < volume_->getNumChannels(), "Slice must not be 0.");
    tgt::svec3 volDim = getDimensions(volumeLLF_, volumeURB_);
    tgt::svec3 sliceDim = slice->getDimensions();
    tgtAssert(volDim.x == sliceDim.x && volDim.y == sliceDim.y && sliceDim.z == 1, "Invalid slice dimensions.");

    rawout_[channel].write(reinterpret_cast<const char*>(slice->getData()), slice->getNumBytes());

    //check if everthing went right
    if (rawout_[channel].bad()) {
        throw tgt::IOException("Failed to write volume data to file (bad stream)");
    }
}

#ifdef VRN_MODULE_HDF5
HDF5SliceWriter::HDF5SliceWriter(const std::string& filename, const tgt::svec3& volumeLLF, const tgt::svec3& volumeURB, const size_t numberOfChannels, const tgt::vec3& spacing, const tgt::vec3& offset, const tgt::mat4& physicalToWorldTransformation, const RealWorldMapping* rwm, size_t deflateLevel, tgt::svec3 chunkSize, bool shuffle) 
    : VolumeSliceWriter()
    , fileVolume_(HDF5FileVolume::createVolume(filename, "/volume", "uint16", getDimensions(volumeLLF, volumeURB), numberOfChannels, true, deflateLevel, chunkSize, shuffle))
    , slicesWritten_(numberOfChannels)
{
    fileVolume_->writeSpacing(spacing);
    fileVolume_->writeOffset(offset);
    fileVolume_->writePhysicalToWorldTransformation(physicalToWorldTransformation);
    if(rwm) {
        fileVolume_->writeRealWorldMapping(*rwm);
    }
}
void HDF5SliceWriter::writeSlice(const VolumeRAM* slice, size_t channel) {
    fileVolume_->writeSlices(slice, slicesWritten_.at(channel)++, channel);
}
#endif

//----------------------------------------------------------------------------------------------------------
// The OctreeSave processor:
//----------------------------------------------------------------------------------------------------------


const std::string OctreeSave::loggerCat_("voreen.core.OctreeSave");

OctreeSave::OctreeSave()
    : VolumeProcessor()
    , inport_(Port::INPORT, "volumehandle.input", "Volume Input", false)
    , filename_("outputFilename", "File", "Select file...", "",
            "PROPERTY NOT INITIALIZED", FileDialogProperty::SAVE_FILE, Processor::INVALID_PATH)
    , saveButton_("save", "Save")
    , progressProp_("progress", "Progress")
    , volumeInfo_("volumeInfo","info")
    , outputFormat_("outputFormat", "Output Format")
    , enableClipping_("enableClipping", "Enable Clipping", false, Processor::INVALID_PROGRAM, Property::LOD_ADVANCED)
    , clipRegion_("clipRegion", "Clip Region", tgt::IntBounds(tgt::ivec3(0), tgt::ivec3(std::numeric_limits<int>::max()-1)), tgt::ivec3(0), tgt::ivec3(std::numeric_limits<int>::max()-1 /* -1 so that the total range value fits into an int */), tgt::ivec3(static_cast<int>(0)), tgt::ivec3(std::numeric_limits<int>::max()), Processor::INVALID_RESULT, Property::LOD_ADVANCED)
    , applyChannelShift_("applyChannelShift", "Apply Channel Shift", false, Processor::INVALID_PROGRAM, Property::LOD_ADVANCED)
    , channelShift0_("channelShift0", "Channel Shift 1", tgt::ivec3(0), tgt::ivec3(-50), tgt::ivec3(50), Processor::INVALID_RESULT, NumericProperty<tgt::ivec3>::STATIC, Property::LOD_ADVANCED)
    , channelShift1_("channelShift1", "Channel Shift 2", tgt::ivec3(0), tgt::ivec3(-50), tgt::ivec3(50), Processor::INVALID_RESULT, NumericProperty<tgt::ivec3>::STATIC, Property::LOD_ADVANCED)
    , channelShift2_("channelShift2", "Channel Shift 3", tgt::ivec3(0), tgt::ivec3(-50), tgt::ivec3(50), Processor::INVALID_RESULT, NumericProperty<tgt::ivec3>::STATIC, Property::LOD_ADVANCED)
    , channelShift3_("channelShift3", "Channel Shift 4", tgt::ivec3(0), tgt::ivec3(-50), tgt::ivec3(50), Processor::INVALID_RESULT, NumericProperty<tgt::ivec3>::STATIC, Property::LOD_ADVANCED)
    , defaultValue_("defaultValue","Outside Volume Value",0,0,65535,Processor::INVALID_RESULT,NumericProperty<int>::STATIC,Property::LOD_ADVANCED)
    , outputProp_("outputProp","Required Data Memory","no volume",Processor::INVALID_RESULT,Property::LOD_DEBUG)
    , enableCompression_("enableCompression", "Enable Compression", false, Processor::INVALID_PROGRAM, Property::LOD_ADVANCED)
    , compressionLevel_("compressionLevel_","Deflate Level",1,1,9,Processor::INVALID_RESULT,NumericProperty<int>::STATIC,Property::LOD_ADVANCED)
    , enableShuffling_("enableShuffling", "Enable Shuffling", false, Processor::INVALID_RESULT, Property::LOD_ADVANCED)
    , chunkSize_("chunkSize", "Chunk Size", tgt::ivec3(std::numeric_limits<int>::max()), tgt::ivec3(1), tgt::ivec3(std::numeric_limits<int>::max()), Processor::INVALID_RESULT, NumericProperty<tgt::ivec3>::STATIC, Property::LOD_ADVANCED)
    , vvodWriter_(0)
{
    addPort(inport_);

    // saving
        addProperty(filename_);
            filename_.setGroupID("saving");
            ON_CHANGE(filename_, OctreeSave, cleanupFilename);
        addProperty(saveButton_);
            saveButton_.setGroupID("saving");
            ON_CHANGE(saveButton_, OctreeSave, saveVolume);
        addProperty(outputFormat_);
            outputFormat_.setGroupID("saving");
#ifdef VRN_MODULE_HDF5
            outputFormat_.addOption("hdf5", "HDF5 (Hierarchical Data Format 5)", VOL_HDF5);
#endif
            outputFormat_.addOption("vvd", "VVD (Voreen Volume Data)", VOL_VVD);
            outputFormat_.addOption("vvod", "VVOD (Voreen Volume Octree Data)", VOL_VVOD);
            ON_CHANGE(outputFormat_, OctreeSave, outputFormatChanged);
        addProperty(outputProp_);
            outputProp_.setGroupID("saving");
            outputProp_.setReadOnly(true);
            adjustOutputProperty();
        addProperty(progressProp_);
            progressProp_.setGroupID("saving");
            addProgressBar(&progressProp_);
    setPropertyGroupGuiName("saving", "Saving");


    // clipping
        addProperty(enableClipping_);
            enableClipping_.setGroupID("clipping");
            ON_CHANGE(enableClipping_, OctreeSave, adjustClipProperties);
            ON_CHANGE(enableClipping_, OctreeSave, adjustOutputProperty);
            ON_CHANGE(enableClipping_, OctreeSave, adjustCompressionProperties);
        addProperty(clipRegion_);
            clipRegion_.setGroupID("clipping");
            ON_CHANGE(clipRegion_, OctreeSave, adjustOutputProperty); // Update output size
            ON_CHANGE(clipRegion_, OctreeSave, adjustCompressionProperties); // Update max chunk size
    setPropertyGroupGuiName("clipping", "Octree Clipping");
    adjustClipProperties();

    // channel shift
        addProperty(applyChannelShift_);
            applyChannelShift_.setGroupID("shift");
            ON_CHANGE(applyChannelShift_, OctreeSave, adjustShiftProperties);
        addProperty(channelShift0_);
            channelShift0_.setGroupID("shift");
        addProperty(channelShift1_);
            channelShift1_.setGroupID("shift");
        addProperty(channelShift2_);
            channelShift2_.setGroupID("shift");
        addProperty(channelShift3_);
            channelShift3_.setGroupID("shift");
        addProperty(defaultValue_);
            defaultValue_.setGroupID("shift");
    setPropertyGroupGuiName("shift", "Channel Shift Correction (in voxels)");
    adjustShiftProperties();

    // compression
        addProperty(enableCompression_);
            enableCompression_.setGroupID("compression");
            ON_CHANGE(enableCompression_, OctreeSave, adjustCompressionProperties);
        addProperty(compressionLevel_);
            compressionLevel_.setGroupID("compression");
        addProperty(enableShuffling_);
            enableShuffling_.setGroupID("compression");
        addProperty(chunkSize_);
            chunkSize_.setGroupID("compression");
    setPropertyGroupGuiName("compression", "Compression");
    adjustCompressionProperties();

    // Information about the volume at inport_.
    addProperty(volumeInfo_);
}

Processor* OctreeSave::create() const {
    return new OctreeSave();
}
void OctreeSave::updateFilenameFilters() {
    std::vector<std::string> filters;
    switch(outputFormat_.getValue()) {
        case VOL_VVOD:
            {
                // file extensions
                std::vector<std::string> extensionVec = vvodWriter_->getSupportedExtensions();
                for (size_t j=0; j<extensionVec.size(); j++) {
                    std::string extension = extensionVec.at(j);
                    std::string filterStr = vvodWriter_->getFormatDescription() + " (*." + extension + ")";
                    filters.push_back(filterStr);
                }

                // filename
                std::vector<std::string> filenamesVec = vvodWriter_->getSupportedFilenames();
                for (size_t j=0; j<filenamesVec.size(); j++) {
                    std::string filename = filenamesVec.at(j);
                    std::string filterStr = vvodWriter_->getFormatDescription() + " (" + filename + ")";
                    filters.push_back(filterStr);
                }
            }
            break;
        case VOL_VVD:
            filters.push_back("Voreen Volume Data (*.vvd)");
            break;
        case VOL_HDF5:
#ifdef VRN_MODULE_HDF5
            filters.push_back("Hierachical Data Format 5 (*.h5)");
            filters.push_back("Hierachical Data Format 5 (*.hdf5)");
            break;
#endif
        default:
            tgtAssert(false, "Invalid output format.");
    }

    filename_.setFileFilter(strJoin(filters, ";;"));
}
void OctreeSave::outputFormatChanged() {
    if (!isInitialized())
        return;

    // Update UI
    updateFilenameFilters();
    adjustClipProperties();
    adjustShiftProperties();
    adjustCompressionProperties();
    adjustOutputProperty();
    cleanupFilename();
}

void OctreeSave::cleanupFilename() {
    if (filename_.get() != "") {
        //get base filename and remove __ChannelX
        std::string basefilename =  tgt::FileSystem::fullBaseName(tgt::FileSystem::cleanupPath(filename_.get()));
        if(basefilename.substr(basefilename.length()-10,9).compare("__Channel") == 0)
            basefilename = basefilename.substr(0,basefilename.length()-10);

        std::string extension;
        switch(outputFormat_.getValue()) {
            case VOL_VVOD:
                extension = vvodWriter_->getSupportedExtensions().at(0);
                break;
            case VOL_VVD:
                extension = "vvd";
                break;
#ifdef VRN_MODULE_HDF5
            case VOL_HDF5:
                extension = "h5";
                break;
#endif
            default:
                tgtAssert(false, "Invalid output format.");
        }

        // Block callbacks, because otherwise this callback would be called again.
        filename_.set(basefilename + "." + extension);
    }

}

void OctreeSave::initialize() {
    VolumeProcessor::initialize();

    vvodWriter_ = new VvodVolumeWriter(progressBars_.empty() ? 0 : progressBars_.front());

    updateFilenameFilters();
    // Select correct extension after deserialization
    cleanupFilename();
}

void OctreeSave::deinitialize() {
    delete vvodWriter_;
    vvodWriter_ = 0;

    volumeInfo_.setVolume(0);
    VolumeProcessor::deinitialize();
}

void OctreeSave::invalidate(int inv) {
    Processor::invalidate(inv);

    if (inport_.hasChanged()) {
        volumeInfo_.setVolume(inport_.getData());

        if (inport_.getData() && inport_.getData()->hasRepresentation<VolumeOctree>())
           saveButton_.setReadOnlyFlag(false);
        else
           saveButton_.setReadOnlyFlag(true);

        adjustPropertiesToInput();
    }
}

void OctreeSave::process() {

    if (isInitialized())
        progressProp_.setProgress(0);

    // do nothing
}

void OctreeSave::saveVolume() {
    if (!isInitialized())
        return;
    if (!inport_.getData()) {
        LWARNING("no input volume");
        return;
    }
    if (!inport_.getData()->hasRepresentation<VolumeOctree>()) {
        LWARNING("Input volume does not have an octree representation");
        return;
    }
    if (filename_.get() == "") {
        LWARNING("no filename specified");
        return;
    }

    if (inport_.getData()->getOrigin().getPath() == filename_.get()) {
        // TODO: this should display a message box in application mode...
        LERROR("Cannot overwrite octree volume by itself!");
    }

    switch(outputFormat_.getValue()) {
        case VOL_VVOD:
            try {
                saveToOctree(filename_.get(), inport_.getData());
            }
            catch(tgt::FileException e) {
                LERROR(e.what());
                filename_.set("");
            }
            break;
#ifdef VRN_MODULE_HDF5
        case VOL_HDF5: // Fall through
#endif
        case VOL_VVD:
            try {
                saveToVolume(filename_.get(), inport_.getData());
            }
            catch(tgt::Exception e) {
                LERROR(e.what());
                filename_.set("");
            }
            catch(std::invalid_argument e) {
                LERROR(e.what());
                filename_.set("");
            }
            break;
        default:
            tgtAssert(false, "Invalid output format.");
    }
}
void OctreeSave::saveToOctree(std::string filename, const VolumeBase* volume) {
    tgtAssert(volume, "No volume");
    tgtAssert(vvodWriter_, "No VVOD writer instantiated");
    vvodWriter_->write(filename, volume);
}
void OctreeSave::saveToVolume(std::string filename, const VolumeBase* volume) {
    //get and check for octree
    const VolumeOctree* octree = volume->getRepresentation<VolumeOctree>();
    tgtAssert(octree, "No octree");
    size_t numChannels = octree->getNumChannels();

    //calculate dimensions of new volumes
    tgt::svec3 volumeLLF, volumeURB, volumeDim;
    if(enableClipping_.get()) {
        volumeLLF = tgt::svec3(clipRegion_.get().getLLF());
        volumeURB = tgt::svec3(clipRegion_.get().getURB());
    } else {
        volumeLLF = tgt::svec3::zero;
        volumeURB = volume->getDimensions() - tgt::svec3::one;
    }
    volumeDim = getDimensions(volumeLLF, volumeURB);

    //calculate dimensions of loaded slices
    tgt::svec3 sliceLLF = volumeLLF, sliceURB = volumeURB;
    tgt::ivec3 minShift = tgt::ivec3::zero, maxShift = tgt::ivec3::zero;
    tgt::svec3 offsetInSlice;
    std::vector<tgt::ivec3> shiftArray(numChannels);
    size_t numSliceVoxels = (volumeDim.x) * (volumeDim.y);
    uint16_t defaultValue = defaultValue_.get(); //magic number
    if(applyChannelShift_.get()) {
        switch(numChannels) {
        case 4:
            minShift = tgt::min(minShift,channelShift3_.get());
            maxShift = tgt::max(maxShift,channelShift3_.get());
            shiftArray[3] = channelShift3_.get();
            //break not needed
        case 3:
            minShift = tgt::min(minShift,channelShift2_.get());
            maxShift = tgt::max(maxShift,channelShift2_.get());
            shiftArray[2] = channelShift2_.get();
            //break not needed
        case 2:
            minShift = tgt::min(minShift,channelShift1_.get());
            maxShift = tgt::max(maxShift,channelShift1_.get());
            shiftArray[1] = channelShift1_.get();
            //break not needed
        case 1:
            minShift = tgt::min(minShift,channelShift0_.get());
            maxShift = tgt::max(maxShift,channelShift0_.get());
            shiftArray[0] = channelShift0_.get();
            break;
        default:
            tgtAssert(false,"should not get here");
            break;
        }
        sliceLLF = tgt::svec3((std::max(0l, static_cast<long>(sliceLLF.x)+minShift.x)),
                              (std::max(0l, static_cast<long>(sliceLLF.y)+minShift.y)),
                              (std::max(0l, static_cast<long>(sliceLLF.z)+minShift.z)));
        sliceURB = tgt::svec3((sliceURB.x+maxShift.x < volume->getDimensions().x ? sliceURB.x+maxShift.x : volume->getDimensions().x-1),
                              (sliceURB.y+maxShift.y < volume->getDimensions().y ? sliceURB.y+maxShift.y : volume->getDimensions().y-1),
                              (sliceURB.z+maxShift.z < volume->getDimensions().z ? sliceURB.z+maxShift.z : volume->getDimensions().z-1));
    } else {
        //set default values
        for(size_t c = 0; c < numChannels; c++)
            shiftArray[c] = tgt::ivec3::zero;
    }
    tgt::svec3 sliceDim = getDimensions(sliceLLF, sliceURB);
    offsetInSlice = volumeLLF - sliceLLF;

    std::unique_ptr<VolumeSliceWriter> output(createSliceWriter(filename, volume, volumeLLF, volumeURB));


    //fill up files with default values at front
    std::unique_ptr<VolumeRAM> outputSlice(VolumeFactory().create("uint16", tgt::svec3(volumeDim.xy(), 1)));
    //TODO numslicevoxels ersetzen?
    uint16_t* channelData = static_cast<uint16_t*>(outputSlice->getData());
    std::fill_n(channelData,numSliceVoxels,defaultValue);
    for(size_t channel = 0; channel < numChannels; channel++) {
        int shiftedLLF_z = (int)(volumeLLF.z) + (int)(shiftArray[channel].z);
        //While the slice we want is outside the volume (z<0) we write default values
        for(int slicePosZ = shiftedLLF_z; slicePosZ < 0; slicePosZ++) {
            output->writeSlice(outputSlice.get(), channel);
        }
    }

    //copy volume slices
    for(size_t slice = sliceLLF.z; slice <= sliceURB.z; slice++) {
        progressProp_.setProgress((float)(slice-sliceLLF.z)/(float)(sliceURB.z));

        std::unique_ptr<VolumeRAM> tmpSlice(octree->createSlice(XY_PLANE,slice,0,0,0,tgt::svec3(sliceLLF),tgt::svec3(sliceURB)));

        const uint16_t* data = reinterpret_cast<const uint16_t*>(tmpSlice->getData());

        for(size_t channel = 0; channel < numChannels; channel++) {
            if((int)slice - shiftArray[channel].z < (int)volumeLLF.z || (int)slice - shiftArray[channel].z > (int)volumeURB.z) {
                continue; //slice not needed for this volume
            } else {
                for (int x = 0; x < (int)(volumeDim.x); x++) {
                    for (int y = 0; y < (int)(volumeDim.y); y++) {
                        const size_t w = y*volumeDim.x + x;
                        //tgtAssert(w<outputSlice->getNumBytes(), "Invalid write");
                        if((int)(volumeLLF.x + x) + shiftArray[channel].x < 0 || x + shiftArray[channel].x >= (int)(sliceDim.x) ||
                           (int)(volumeLLF.y + y) + shiftArray[channel].y < 0 || y + shiftArray[channel].y >= (int)(sliceDim.y)) {
                            // Sampleposition is outside the volume. We set the defaultValue
                            channelData[w] = defaultValue;
                        } else {
                            // We are in the volume! Read data from the slice previously read from the octree.
                            const size_t r = ((offsetInSlice.y+y+shiftArray[channel].y)*sliceDim.x + (offsetInSlice.x+x+shiftArray[channel].x))*numChannels+channel;
                            //tgtAssert(r<tmpSlice->getNumBytes(), "Invalid read");
                            channelData[w] = data[r];
                        }
                    }
                }
            }

            //write data to file
            output->writeSlice(outputSlice.get(), channel);
        }
    }

    //fill up files with default values at back
    std::fill_n(channelData,numSliceVoxels,defaultValue);
    for(size_t channel = 0; channel < numChannels; channel++) {
        int shiftedURB_z = (int)(volumeURB.z) + (int)(shiftArray[channel].z);
        //While the slice we want is outside the volume (z>=dim.z) we write default values
        for(int slicePosZ = shiftedURB_z; slicePosZ >= static_cast<int>(volume->getDimensions().z); slicePosZ--) {
            output->writeSlice(outputSlice.get(), channel);
        }
    }

    //clean up
    progressProp_.setProgress(1.f);
}
VolumeSliceWriter* OctreeSave::createSliceWriter(const std::string& filename, const VolumeBase* volume, const tgt::svec3& volumeLLF, const tgt::svec3& volumeURB) const {
    tgtAssert(volume, "No volume.");
    switch(outputFormat_.getValue()) {
        case VOL_VVD:
            return new VVDSliceWriter(filename, volume, volumeLLF, volumeURB);
        case VOL_HDF5:
#ifdef VRN_MODULE_HDF5
            {
                tgt::vec3 spacing = volume->getSpacing();
                tgt::vec3 offset = volume->getOffset()+tgt::vec3(volumeLLF)*spacing;
                tgt::mat4 physicalToWorldTransformation = volume->getPhysicalToWorldMatrix();
                RealWorldMapping rwm;
                RealWorldMapping* rwmPointer = nullptr;
                if(volume->hasMetaData("RealWorldMapping")) {
                    rwm = volume->getRealWorldMapping();
                    rwmPointer = &rwm;
                }
                if(enableCompression_.get()) {
                    return new HDF5SliceWriter(filename, volumeLLF, volumeURB, volume->getNumChannels(), spacing, offset, physicalToWorldTransformation, rwmPointer, compressionLevel_.get(), chunkSize_.get(), enableShuffling_.get());
                } else {
                    return new HDF5SliceWriter(filename, volumeLLF, volumeURB, volume->getNumChannels(), spacing, offset, physicalToWorldTransformation, rwmPointer, 0, tgt::svec3::zero, false);
                }
            }
#endif
        default:
            tgtAssert(false, "Cannot create VolumeSliceWriter for given output format.");
            return nullptr;
    }
}
//---------------------------------------------------------------------------------------------
//      adjustment functions
//---------------------------------------------------------------------------------------------
void OctreeSave::adjustClipProperties() {
    bool visible = enableClipping_.get();
    clipRegion_.setVisibleFlag(visible);

    bool enabled = (outputFormat_.getValue() != VOL_VVOD);
    enableClipping_.setReadOnlyFlag(!enabled);
    clipRegion_.setReadOnlyFlag(!enabled);
}

void OctreeSave::adjustShiftProperties() {
    bool visible = applyChannelShift_.get();
    channelShift0_.setVisibleFlag(visible);
    channelShift1_.setVisibleFlag(visible);
    channelShift2_.setVisibleFlag(visible);
    channelShift3_.setVisibleFlag(visible);
    defaultValue_.setVisibleFlag(visible);

    bool enabled = (outputFormat_.getValue() != VOL_VVOD);
    applyChannelShift_.setReadOnlyFlag(!enabled);
    channelShift0_.setReadOnlyFlag(!enabled);
    channelShift1_.setReadOnlyFlag(!enabled);
    channelShift2_.setReadOnlyFlag(!enabled);
    channelShift3_.setReadOnlyFlag(!enabled);
    defaultValue_.setReadOnlyFlag(!enabled);
}

void OctreeSave::adjustCompressionProperties() {
    bool visible = enableCompression_.get();
    compressionLevel_.setVisibleFlag(visible);
    enableShuffling_.setVisibleFlag(visible);
    chunkSize_.setVisibleFlag(visible);

    bool enabled = (outputFormat_.getValue() == VOL_HDF5);
    enableCompression_.setReadOnlyFlag(!enabled);
    compressionLevel_.setReadOnlyFlag(!enabled);
    enableShuffling_.setReadOnlyFlag(!enabled);
    chunkSize_.setReadOnlyFlag(!enabled);

    if(!inport_.getData()) {
        return;
    }
    tgt::svec3 dim;
    if(enableClipping_.get()) {
        // Set max chunk size to size of clipped volume
        dim = tgt::svec3(clipRegion_.get().getURB() - clipRegion_.get().getLLF()) + tgt::svec3::one;
    } else {
        // Set max chunk size to size of whole volume
        dim = inport_.getData()->getDimensions();
    }
    chunkSize_.setMaxValue(dim);
    // Set default chunk size to a single slice of the (clipped) volume
    chunkSize_.set(tgt::svec3(dim.xy(), 1));
}
void OctreeSave::adjustOutputProperty() {
    if(outputFormat_.getValue() == VOL_VVOD) {
        outputProp_.setVisibleFlag(false);
        return;
    }
    outputProp_.setVisibleFlag(true);
    if(const VolumeBase* vol = inport_.getData()) {
        std::string numChannels = itos(vol->getNumChannels());
        std::string numVoxels;
        std::string memory;
        if(enableClipping_.get()) {
            /*tgt::svec3 dim = tgt::svec3(clipLeft_.get() - clipRight_.get(),
                                        clipBack_.get() - clipFront_.get(),
                                        clipTop_.get() - clipBottom_.get()) + tgt::svec3::one;*/
            tgt::svec3 dim = tgt::svec3(clipRegion_.get().getURB() - clipRegion_.get().getLLF()) + tgt::svec3::one;
            numVoxels = itos(tgt::hmul(dim));
            memory = formatMemorySize(tgt::hmul(dim)*vol->getBytesPerVoxel());
        } else {
            numVoxels = itos(tgt::hmul(vol->getDimensions()));
            memory = formatMemorySize(tgt::hmul(vol->getDimensions())*vol->getBytesPerVoxel());
        }
        outputProp_.set(numChannels + "C * " + numVoxels + "V = " + memory);
    } else {
        outputProp_.set("no volume");
    }
}

void OctreeSave::adjustPropertiesToInput() {
    //progressProp_.setProgress(0.f); //This creates a segfault apparently?

    if(!inport_.getData()) {
        outputProp_.set("no volume");
        return;
    }

    // clipping
    //if (oldVolumeDimensions_ == tgt::svec3::zero)
    //    oldVolumeDimensions_ = inport_.getData()->getDimensions();
    tgt::svec3 oldVolumeDimensions = tgt::svec3(clipRegion_.getMaxValue()) + tgt::svec3::one;

    tgt::svec3 numSlices = inport_.getData()->getDimensions();

    // reset clip region if the volume dimensions have changed, else keep the value (e.g. for keeping the clipping for multiple time steps)
    if (numSlices != oldVolumeDimensions) {
        // set max values and reset clip region
        clipRegion_.setMaxValue(tgt::ivec3(numSlices) - tgt::ivec3::one);
        clipRegion_.set(tgt::IntBounds(tgt::ivec3(0), tgt::ivec3(numSlices) - tgt::ivec3::one));
    }

    adjustOutputProperty();
    adjustCompressionProperties();
}
}   // namespace voreen
