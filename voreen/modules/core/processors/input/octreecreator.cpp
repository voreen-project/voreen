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

#include "octreecreator.h"

#include "voreen/core/datastructures/octree/volumeoctree.h"
#include "voreen/core/datastructures/octree/octreebrickpoolmanager.h"
#include "voreen/core/datastructures/octree/octreebrickpoolmanagerdisk.h"
#include "voreen/core/datastructures/octree/octreeutils.h"

#include "voreen/core/voreenapplication.h"
#include "voreen/core/datastructures/volume/volumeram.h"
#include "voreen/core/datastructures/volume/volumedisk.h"
#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/datastructures/volume/operators/volumeoperatorresizepoweroftwo.h"
#include "voreen/core/datastructures/volume/operators/volumeoperatorconvert.h"
#include "voreen/core/datastructures/volume/volumedecorator.h"
#include "voreen/core/datastructures/meta/realworldmappingmetadata.h"
#include "voreen/core/datastructures/volume/volumeminmax.h"
#include "voreen/core/datastructures/volume/histogram.h"
#include "voreen/core/datastructures/volume/volumehash.h"
#include "voreen/core/datastructures/volume/volumepreview.h"
#include "voreen/core/datastructures/octree/octreebrickpoolmanagerdisk.h"
#include "voreen/core/datastructures/octree/octreeutils.h"
#include "voreen/core/io/serialization/serialization.h"
#include "voreen/core/utils/hashing.h"
#include "voreen/core/utils/memoryinfo.h"

#include "tgt/vector.h"
#include "tgt/tgt_math.h"
#include "tgt/filesystem.h"

#include "tinyxml/tinyxml.h"

#include <queue>
#include <stack>

using tgt::ivec2;
using tgt::ivec3;
using tgt::ivec4;
using tgt::vec2;
using tgt::vec3;
using tgt::vec4;
using tgt::svec3;

namespace {

const std::string CACHE_SUBDIR =             "OctreeCreator";
const std::string OCTREE_FILENAME =          "octree.xml";
const std::string BRICK_BUFFER_SUBDIR =      "brickBuffer";
const std::string BRICK_BUFFER_FILE_PREFIX = "buffer_";

/// Retrieves the octree cache limit from the VoreenApplication.
uint64_t getCacheLimit() {
    const std::string propertyID = "octreeCacheLimit";
    tgtAssert(voreen::VoreenApplication::app(), "VoreenApplication not instantiated");
    voreen::IntProperty* cacheSizeProp = dynamic_cast<voreen::IntProperty*>(voreen::VoreenApplication::app()->getProperty(propertyID));
    if (cacheSizeProp) {
        return (uint64_t)cacheSizeProp->get() * (1<<30);
    }
    else {
        LWARNINGC("voreen.OctreeCreator", "IntProperty '" << propertyID << "' not found");
        return 0;
    }
}

/// Returns the brick dimension for a given volume dimension and tree depth.
size_t computeBrickDim(const tgt::svec3& volumeDim, size_t numLevels) {
    tgtAssert(numLevels > 0, "numLevels is 0");
    size_t octreeDim = tgt::nextLargerPowerOfTwo(tgt::max(volumeDim));
    size_t brickDim = octreeDim / (size_t)(1 << (numLevels-1));
    return brickDim;
}

/// Returns the octree depth for the passed brick and volume dimensions
size_t computeTreeDepth(const tgt::svec3& brickDim, const tgt::svec3& volumeDim) {
    tgtAssert(isCubicAndPot(brickDim), "brick dim not cubic-pot");
    size_t octreeDim = tgt::nextLargerPowerOfTwo(tgt::max(volumeDim));
    size_t treeDepth = tgt::ilog2((int)octreeDim / (int)brickDim.x) + 1;
    return treeDepth;
}

size_t computeMaxTreeDepth(const tgt::svec3& volumeDim) {
    size_t octreeDim = tgt::nextLargerPowerOfTwo(tgt::max(volumeDim));
    return (size_t)tgt::ilog2((int)octreeDim);
}

} // namespace anonymous

namespace voreen {

const std::string OctreeCreator::loggerCat_("voreen.OctreeCreator");

OctreeCreator::OctreeCreator()
    : AsyncComputeProcessor<ComputeInput, ComputeOutput>()
    , volumeInport_(Port::INPORT, "volumeInput", "Volume Input")
    , volumeInport2_(Port::INPORT, "volumeInput2", "Volume Input 2")
    , volumeInport3_(Port::INPORT, "volumeInput3", "Volume Input 3")
    , volumeInport4_(Port::INPORT, "volumeInput4", "Volume Input 4")
    , volumeOutport_(Port::OUTPORT, "volumeOutport", "Volume Output")
    //, saveOctreeFile_("saveOctreeName","Octree File","Select Octree File", VoreenApplication::app()->getUserDataPath(),
    //"Voreen Volume Octree Data (*.vvod)", FileDialogProperty::SAVE_FILE,Processor::VALID)
    //, saveOctreeButton_("saveOctreeButton","Save Octree")
    , brickDimensions_("brickDimensions", "Brick Dimensions")
    , treeDepth_("treeDepth", "Tree Depth", 6, 1, 8)
    , homogeneityThreshold_("homogeneityThreshold", "Homogeneity Threshold", 0.00f, 0.f, 0.1f)
    , useRelativeThreshold_("useRelativeThreshold", "Relative Threshold", true)
    , brickPoolManager_("brickPoolManager", "Brick Pool Manager")
    , singleBufferMemorySize_("singleBufferMemorySize", "Page File Size (MB)", 32, 1, 256)
    , numThreads_("numThreads", "Num Threads", 8, 1, 16, VALID)
    , clearOctree_("clearOctree", "Clear Octree")
    , statusMessage_("progressMessage", "Progress Message", "", Processor::VALID)
    , currentConfigurationHash_("octreeHash", "Octree Hash", "", INVALID_RESULT, Property::LOD_ADVANCED)
{
    addPort(volumeInport_);
    addPort(volumeInport2_);
    addPort(volumeInport3_);
    addPort(volumeInport4_);
    addPort(volumeOutport_);

    brickDimensions_.addOption("treeDepth", "Derive from Tree Depth",     0);
    //brickDimensions_.addOption("2",       "[2 2 2]",                    2);
    brickDimensions_.addOption("4",         "[4 4 4]",                    4);
    brickDimensions_.addOption("8",         "[8 8 8]",                    8);
    brickDimensions_.addOption("16",        "[16 16 16]",                16);
    brickDimensions_.addOption("32",        "[32 32 32]",                32);
    brickDimensions_.addOption("64",        "[64 64 64]",                64);
    brickDimensions_.addOption("128",       "[128 128 128]",            128);
    brickDimensions_.selectByKey("treeDepth");
    brickDimensions_.onChange(MemberFunctionCallback<OctreeCreator>(this, &OctreeCreator::updatePropertyConfiguration));
    addProperty(brickDimensions_);
    treeDepth_.setReadOnlyFlag(true);
    addProperty(treeDepth_);
    //homogeneityThreshold_.setNumDecimals(2);
    //homogeneityThreshold_.setStepping(0.01f);
    homogeneityThreshold_.setTracking(false);
    addProperty(homogeneityThreshold_);
    addProperty(useRelativeThreshold_);
    treeDepth_.setGroupID("configuration");
    brickDimensions_.setGroupID("configuration");
    homogeneityThreshold_.setGroupID("configuration");
    useRelativeThreshold_.setGroupID("configuration");

    brickPoolManager_.addOption("brickPoolManagerRAM",  "RAM (non-persistent)");
    brickPoolManager_.addOption("brickPoolManagerDisk", "Disk");
    brickPoolManager_.select("brickPoolManagerDisk");
    addProperty(brickPoolManager_);
    brickPoolManager_.onChange(MemberFunctionCallback<OctreeCreator>(this, &OctreeCreator::updatePropertyConfiguration));
    addProperty(singleBufferMemorySize_);
    brickPoolManager_.setGroupID("configuration");
    singleBufferMemorySize_.setGroupID("configuration");
    numThreads_.setGroupID("configuration");
    addProperty(numThreads_);
    setPropertyGroupGuiName("configuration", "Octree Configuration");

    /*addProperty(saveOctreeFile_);
    saveOctreeFile_.setGroupID("save");
    addProperty(saveOctreeButton_);
    saveOctreeButton_.setGroupID("save");
    saveOctreeButton_.onChange(MemberFunctionCallback<OctreeCreator>(this, &OctreeCreator::saveOctreeToVVOD));
    setPropertyGroupGuiName("save", "Save Current Octree");*/

    addProperty(clearOctree_);
    clearOctree_.onChange(MemberFunctionCallback<OctreeCreator>(this, &OctreeCreator::clearOctree));

    addProperty(statusMessage_);
    statusMessage_.setEditable(false);

    addProperty(currentConfigurationHash_);
    currentConfigurationHash_.setEditable(false);
}

OctreeCreator::~OctreeCreator() {
    // Remove old commands.
    VoreenApplication::app()->getCommandQueue()->removeAll(this);
}

Processor* OctreeCreator::create() const {
    return new OctreeCreator();
}

void OctreeCreator::initialize() {
    AsyncComputeProcessor<ComputeInput, ComputeOutput>::initialize();

    updatePropertyConfiguration();

    limitCacheSize(VoreenApplication::app()->getCachePath(), getCacheLimit(), false);
}

void OctreeCreator::deinitialize() {
    volumeOutport_.clear();
    if (!VoreenApplication::app()->useCaching() && currentConfigurationHash_.get() != "")
        clearOctree();

    limitCacheSize(VoreenApplication::app()->getCachePath(), getCacheLimit(), false);

    AsyncComputeProcessor<ComputeInput, ComputeOutput>::deinitialize();
}

void OctreeCreator::adjustPropertiesToInput() {
    updatePropertyConfiguration();
}

bool OctreeCreator::isReady() const {
    return volumeInport_.isReady() && volumeOutport_.isReady();
}

OctreeCreatorInput OctreeCreator::prepareComputeInput() {
    VolumeBase* inputVolume = const_cast<VolumeBase*>(volumeInport_.getData());
    if (!inputVolume) {
        updatePropertyConfiguration();
        throw InvalidInputException("No input", InvalidInputException::S_WARNING);
    }

    // this is necessary for handling octree volumes which have been loaded, e.g. by a VolumeSource processor
    if (inputVolume->hasRepresentation<VolumeOctree>() && !volumeInport2_.hasData() && !volumeInport3_.hasData() && !volumeInport3_.hasData()) {
        // we already have a single octree input -> set it to outport
        volumeOutport_.setData(inputVolume, false);
        throw InvalidInputException("Input already is has an octree representation", InvalidInputException::S_WARNING);
    }
    else if (inputVolume->hasRepresentation<VolumeOctree>()) {
        // input volume 1 has an octree, but other inputs are present
        volumeOutport_.clear();
        throw InvalidInputException("Input volume 1 is octree volume, but other input volumes are present", InvalidInputException::S_ERROR);
    }
    else if ((volumeInport2_.getData() && volumeInport2_.getData()->hasRepresentation<VolumeOctree>())
        || (volumeInport3_.getData() && volumeInport3_.getData()->hasRepresentation<VolumeOctree>())
        || (volumeInport4_.getData() && volumeInport4_.getData()->hasRepresentation<VolumeOctree>()))
    {
        volumeOutport_.clear();
        throw InvalidInputException("Octree volume input is only allowed in first inport", InvalidInputException::S_ERROR);
    }

    //get current hash, if exists
    if (!volumeOutport_.hasData())
        currentConfigurationHash_.set("");
    //get current input hash
    std::string configHash = getConfigurationHash();
    //we are up to date
    if (configHash == currentConfigurationHash_.get()) {
        updatePropertyConfiguration();
        throw InvalidInputException("Up to date!", InvalidInputException::S_IGNORE);
    }

    //remove old octree
    volumeOutport_.clear();
    if (!VoreenApplication::app()->useCaching() && currentConfigurationHash_.get() != "")
        clearOctree();
    currentConfigurationHash_.set("");
    // remove existing octree representation
    inputVolume->removeRepresentation<VolumeOctreeBase>();

    // retrieve input RAM volumes
    std::vector<const VolumeBase*> inputVolumes;
    if (volumeInport_.hasData())
        inputVolumes.push_back(volumeInport_.getData());
    if (volumeInport2_.hasData())
        inputVolumes.push_back(volumeInport2_.getData());
    if (volumeInport3_.hasData())
        inputVolumes.push_back(volumeInport3_.getData());
    if (volumeInport4_.hasData())
        inputVolumes.push_back(volumeInport4_.getData());

    if (inputVolumes.empty())
        throw InvalidInputException("No input volumes", InvalidInputException::S_ERROR);

    const tgt::svec3 volumeDim = inputVolumes.front()->getDimensions();

    // select brick pool manager
    OctreeBrickPoolManagerBase* brickPoolManager = nullptr;
    if (brickPoolManager_.isSelected("brickPoolManagerRAM")) {
        // determine if the octree size stays within the set main memory limit -> otherwise, it can not be constructed using a RAM octree manager
        size_t allowedMemory = VoreenApplication::app()->getCpuRamLimit();
        size_t requiredMemory = 0;
        for (auto it : inputVolumes) {
            requiredMemory += it->getNumVoxels() * it->getNumChannels() * 2; // 2 bytes per voxel, since we store the data as uint16_t
        }
        // required size of the octree has some overhead due to the various resolution levels
        requiredMemory += requiredMemory / 4;

        if (requiredMemory > allowedMemory)
            throw InvalidInputException("CPU memory limit not sufficient for data size required by RAM brick pool manager", InvalidInputException::S_ERROR);

        // if the memory is sufficient: create the brick pool manager
        brickPoolManager = new OctreeBrickPoolManagerRAM(singleBufferMemorySize_.get() << 20);
    }
    else if (brickPoolManager_.isSelected("brickPoolManagerDisk")) {
        std::string brickPoolPath = tgt::FileSystem::cleanupPath(getOctreeStoragePath() + "/" + BRICK_BUFFER_SUBDIR);
        if (!tgt::FileSystem::dirExists(brickPoolPath))
            tgt::FileSystem::createDirectoryRecursive(brickPoolPath);
        brickPoolManager = new OctreeBrickPoolManagerDisk(static_cast<size_t>(singleBufferMemorySize_.get()) << 20,
            VoreenApplication::app()->getCpuRamLimit(), brickPoolPath, BRICK_BUFFER_FILE_PREFIX);
    }
    else {
        throw InvalidInputException("Unknown brick pool manager: " + brickPoolManager_.get(), InvalidInputException::S_ERROR);
    }

    // determine homogeneity threshold
    float homogeneityThreshold = homogeneityThreshold_.get();
    if (useRelativeThreshold_.get()) {
        // use middle slice of first volume to calibrate threshold
        size_t sliceID = volumeDim.z / 2;
        try {
            VolumeRAM* slice = 0;
            if (inputVolumes.front()->hasRepresentation<VolumeRAM>()) {
                slice = inputVolumes.front()->getRepresentation<VolumeRAM>()->getSubVolume(tgt::svec3(volumeDim.x, volumeDim.y, 1), tgt::svec3(0, 0, sliceID));
            }
            else if (inputVolumes.front()->hasRepresentation<VolumeDisk>()) {
                slice = inputVolumes.front()->getRepresentation<VolumeDisk>()->loadSlices(sliceID, sliceID);
            }
            tgtAssert(slice, "no slice");

            float min = slice->minNormalizedValue();
            float max = slice->maxNormalizedValue();
            tgtAssert(min <= max, "invalid min/max values");
            homogeneityThreshold *= (max - min);

            delete slice;

            LDEBUG("Middle slice min/max: " << min << "/" << max);
            LDEBUG("Homogeneity threshold: " << homogeneityThreshold);
        }
        catch (std::exception& e) {
            LWARNING("Failed to extract slice " << sliceID << " from input volume for min/max range computation: " << e.what());
        }
    }

    size_t brickDim;
    if (brickDimensions_.isSelected("treeDepth")) {
        brickDim = computeBrickDim(inputVolumes.front()->getDimensions(), treeDepth_.get());
        brickDim = std::max<size_t>(brickDim, 4);
    }
    else {
        brickDim = brickDimensions_.getValue();
    }
    tgtAssert(brickDim >= 4, "brick dim too small");
    tgtAssert(tgt::isPowerOfTwo((int)brickDim), "brick dim not power-of-two");

    return OctreeCreatorInput(
        VoreenApplication::app()->useCaching(),
        inputVolumes,
        brickPoolManager,
        homogeneityThreshold,
        brickDim,
        numThreads_.get()
    );
}

OctreeCreatorOutput OctreeCreator::compute(OctreeCreatorInput input, ProgressReporter& progressReporter) const {

    std::unique_ptr<VolumeOctreeBase> octree;
    std::string errorMessage;

    // Load Octree from Cache, if available.
    if (input.loadCached) {
        updateStatusMessage("Loading octree from cache...");
        octree.reset(restoreOctreeFromCache());
        if (octree) {
            updateStatusMessage("Octree loaded from cache.");
        }
    }

    // Create Octree, also when cache was not available.
    if(!octree) {
        updateStatusMessage("Generating octree...");
        
        try {
            octree.reset(new VolumeOctree(input.input, (int)input.brickDim, input.homogeneityThreshold, input.brickPoolManager, input.numThreads, &progressReporter));
            updateStatusMessage("Octree generated.");
        }
        catch (VoreenException& e) {
            updateStatusMessage("Failed to generate octree!");
            LINFO(MemoryInfo::getProcessMemoryUsageAsString());
            LINFO(MemoryInfo::getAvailableMemoryAsString());
            errorMessage = e.what();
        }
    }

    // Done.
    return OctreeCreatorOutput(
        std::move(octree),
        errorMessage
    );
}

void OctreeCreator::processComputeOutput(OctreeCreatorOutput output) {

    VolumeBase* inputVolume = const_cast<VolumeBase*>(volumeInport_.getData());
    tgtAssert(inputVolume, "input volume null");

    // Retrieve Octree.
    std::unique_ptr<VolumeOctreeBase> octree = std::move(output.octree);
    if(!octree || !output.errorMessage.empty()) {
        std::string errorMsg = "Failed to generate octree: " + output.errorMessage;
        LERROR(errorMsg);
        VoreenApplication::app()->showMessageBox(getGuiName(), errorMsg, true);
        return;
    }

    LDEBUG("- " << MemoryInfo::getProcessMemoryUsageAsString());
    LDEBUG("- " << MemoryInfo::getAvailableMemoryAsString());

    // store octree to cache, always store, although caching might be disabled to enable save option
    if (!brickPoolManager_.isSelected("brickPoolManagerRAM")) {
        storeOctreeToCache(octree.get());
    }

    // Configuration did not change since the last call of prepareComputeInput,
    // since the calculation would have been aborted otherwise.
    currentConfigurationHash_.set(getConfigurationHash());
    if (!(octree->getOctreeConfigurationHash() == currentConfigurationHash_.get()))
        octree->setOctreeConfigurationHash(currentConfigurationHash_.get());

    // assign RAM limit
    size_t ramLimit = VoreenApplication::app()->getCpuRamLimit();
    if (const OctreeBrickPoolManagerDisk* brickPoolManager =
        dynamic_cast<const OctreeBrickPoolManagerDisk*>(static_cast<VolumeOctree*>(octree.get())->getBrickPoolManager())) {
        const_cast<OctreeBrickPoolManagerDisk*>(brickPoolManager)->setRAMLimit(ramLimit);
    }

    octree->logDescription();

    // min/max values
    std::vector<float> minValues, maxValues, minNormValues, maxNormValues;
    for (size_t i = 0; i < octree->getNumChannels(); i++) {
        float minNorm = octree->getRootNode()->getMinValue(i) / 65535.f;
        float maxNorm = octree->getRootNode()->getMaxValue(i) / 65535.f;
        tgtAssert(minNorm <= maxNorm, "invalid min/max values");
        float min = inputVolume->getRealWorldMapping().normalizedToRealWorld(minNorm);
        float max = inputVolume->getRealWorldMapping().normalizedToRealWorld(maxNorm);
        minValues.push_back(min);
        maxValues.push_back(max);
        minNormValues.push_back(minNorm);
        maxNormValues.push_back(maxNorm);
    }
    VolumeMinMax* volumeMinMax = new VolumeMinMax(minValues, maxValues, minNormValues, maxNormValues);

    // histograms
    std::vector<Histogram1D> histograms;
    for (size_t i = 0; i < octree->getNumChannels(); i++) {
        histograms.push_back(Histogram1D(*(octree->getHistogram(i))));
    }
    VolumeHistogramIntensity* histogramData = new VolumeHistogramIntensity(histograms);

    // single-channel uint16_t octree => add octree to input volume handle
    if (octree->getNumChannels() == 1 && inputVolume->getFormat() == "uint16") {
        // assign minmax derived data to inputVolume volume, if not present
        if (!inputVolume->hasDerivedData<VolumeMinMax>())
            const_cast<VolumeBase*>(inputVolume)->addDerivedData(volumeMinMax);
        else
            delete volumeMinMax;

        // assign histograms, if not present
        if (!inputVolume->hasDerivedData<VolumeHistogramIntensity>())
            const_cast<VolumeBase*>(inputVolume)->addDerivedData(histogramData);
        else
            delete histogramData;

        inputVolume->addRepresentation(octree.release());

        volumeOutport_.setData(inputVolume, false);
    }
    else { // multi-channel volume or differing data type => create new volume handle
        VolumeBase* outputVolume = new Volume(octree.release(), inputVolume);
        outputVolume->setOrigin(inputVolume->getOrigin());

        outputVolume->addDerivedData(volumeMinMax);
        outputVolume->addDerivedData(histogramData);
        outputVolume->getDerivedData<VolumePreview>(); //< prevent creation in background thread (brick pool access)

        // TODO: add hash computation to VolumeOctreeBase interface
        VolumeHash* volumeHash = new VolumeHash();
        volumeHash->setHash(VoreenHash::getHash(getConfigurationHash()));
        outputVolume->addDerivedData(volumeHash);

        volumeOutport_.setData(outputVolume, true);
    }
    tgtAssert(volumeOutport_.hasData(), "outport has no data");

    //updatePropertyConfiguration();
}

void OctreeCreator::storeOctreeToCache(const VolumeOctreeBase* octree) const {
    tgtAssert(octree, "null pointer passed");

    // detect/create base cache path
    const std::string cachePath = getOctreeStoragePath();
    if (!tgt::FileSystem::dirExists(cachePath) && !tgt::FileSystem::createDirectoryRecursive(cachePath)) {
        LWARNING("Failed to create cache directory: " << cachePath);
        return;
    }
    else {
        //tgt::FileSystem::clearDirectory(cachePath);
    }

    // serialize octree
    std::string octreePath = tgt::FileSystem::cleanupPath(cachePath + "/" + OCTREE_FILENAME);
    XmlSerializer serializer(octreePath);
    try {
        serializer.serialize("Octree", octree);
    }
    catch (SerializationException& e) {
        LWARNING("Failed to serialize octree: " << e.what());
        return;
    }

    // write serialization to stream
    std::ostringstream textStream;
    try {
        serializer.write(textStream);
        if (textStream.fail()) {
            LWARNING("Failed to write octree serialization to string stream");
            return;
        }
    }
    catch (std::exception& e) {
        LWARNING("Failed to write octree serialization to string stream: " << e.what());
        return;
    }

    // now we have a valid string stream containing the serialized octree
    // => open output file and write it to the file
    std::fstream fileStream(octreePath.c_str(), std::ios_base::out);
    if (fileStream.fail()) {
        LWARNING("Failed to open file '" << octreePath << "' for writing.");
        return;
    }

    try {
        fileStream << textStream.str();
    }
    catch (std::exception& e) {
        LWARNING("Failed to write serialization data stream to file '" << octreePath << "': " << std::string(e.what()));
    }
    fileStream.close();

    LINFO("Serialized octree to file: " << octreePath);

    limitCacheSize(VoreenApplication::app()->getCachePath(), getCacheLimit(), true);
}

VolumeOctreeBase* OctreeCreator::restoreOctreeFromCache() const {
    const std::string octreeFile = tgt::FileSystem::cleanupPath(getOctreeStoragePath() + "/" + OCTREE_FILENAME);
    if (!tgt::FileSystem::fileExists(octreeFile)) {
        LINFO("No matching octree found in cache.");
        return 0;
    }

    // open file for reading
    std::fstream fileStream(octreeFile.c_str(), std::ios_base::in);
    if (fileStream.fail()) {
        LWARNING("Failed to open cached octree file '" << octreeFile << "' for reading.");
        return 0;
    }

    // read data stream into deserializer
    XmlDeserializer d(octreeFile);
    try {
        d.read(fileStream);
    }
    catch (std::exception& e) {
        LWARNING("Failed to read serialization data stream from cached octree file '" + octreeFile + "': " + e.what());
        fileStream.close();
        return 0;
    }

    fileStream.close();

    // deserialize octree from data stream
    VolumeOctreeBase* octree = 0;
    try {
        LDEBUG("Restoring cached octree from file: " << octreeFile );
        LDEBUG("- Before: " << MemoryInfo::getProcessMemoryUsageAsString());
        LDEBUG("- Before: " << MemoryInfo::getAvailableMemoryAsString());

        d.deserialize("Octree", octree);
        tgtAssert(octree, "null pointer after deserialization");
        tgtAssert(octree->getRootNode(), "deserialized octree has no root node");

        LDEBUG("- After: " << MemoryInfo::getProcessMemoryUsageAsString());
        LDEBUG("- After: " << MemoryInfo::getAvailableMemoryAsString());

        // update file access time
        tgt::FileSystem::updateFileTime(octreeFile);
    }
    catch (std::exception& e) {
        LWARNING("Failed to restore cached octree from file '" + octreeFile + "': " + e.what());
        delete octree;
        octree = 0;
    }

    return octree;
}

void OctreeCreator::clearOctree() {
    // clear octree
    volumeOutport_.clear();
    currentConfigurationHash_.set("");

    // clear disk cache
    std::string octreePath = getOctreeStoragePath();
    if (tgt::FileSystem::dirExists(octreePath)) {
        if (!tgt::FileSystem::deleteDirectoryRecursive(octreePath))
            LWARNING("Failed to delete octree on disk: " + octreePath);
    }

    updatePropertyConfiguration();
}

std::string OctreeCreator::getOctreeStoragePath() const {
    return VoreenApplication::app()->getCachePath(CACHE_SUBDIR + "/" + getConfigurationHash());
}

std::string OctreeCreator::getConfigurationHash() const {
    // get input volume hash
    std::string volumeHash;
    if (volumeInport_.hasData())
        volumeHash += volumeInport_.getData()->getHash();
    if (volumeInport2_.hasData())
        volumeHash += volumeInport2_.getData()->getHash();
    if (volumeInport3_.hasData())
        volumeHash += volumeInport3_.getData()->getHash();
    if (volumeInport4_.hasData())
        volumeHash += volumeInport4_.getData()->getHash();
    volumeHash = VoreenHash::getHash(volumeHash);

    // compute property string
    std::map<std::string, const Property*> propertyMap;
    propertyMap[brickDimensions_.getID()] = &brickDimensions_;
    propertyMap[treeDepth_.getID()] = &treeDepth_;
    propertyMap[homogeneityThreshold_.getID()] = &homogeneityThreshold_;
    propertyMap[useRelativeThreshold_.getID()] = &useRelativeThreshold_;
    propertyMap[brickPoolManager_.getID()] = &brickPoolManager_;
    propertyMap[singleBufferMemorySize_.getID()] = &singleBufferMemorySize_;

    XmlSerializer s;
    Serializer serializer(s);
    const bool usePointerContentSerialization = serializer.getUsePointerContentSerialization();
    serializer.setUsePointerContentSerialization(true);
    try {
        serializer.serialize("Properties", propertyMap, "Property", "name");
    }
    catch (SerializationException& e) {
        LWARNING(e.what());
    }
    serializer.setUsePointerContentSerialization(usePointerContentSerialization);
    std::stringstream propertyConfig;
    s.write(propertyConfig);
    std::string propertyHash = VoreenHash::getHash(propertyConfig.str());

    // concatenate property and volume hash
    return volumeHash + "-" + propertyHash;
}

void OctreeCreator::updatePropertyConfiguration() {
    treeDepth_.setReadOnlyFlag(!brickDimensions_.isSelected("treeDepth"));

    if (volumeInport_.hasData()) {
        tgt::svec3 volumeDim = volumeInport_.getData()->getDimensions();

        /*size_t maxDepth = tgt::clamp(computeMaxTreeDepth(volumeDim), (size_t)1, (size_t)8);
        treeDepth_.setMaxValue((int)maxDepth);

        size_t brickDim = computeBrickDim(volumeDim, treeDepth_.get());
        tgtAssert(tgt::isPowerOfTwo((int)brickDim), "brick dim not power-of-two");
        brickDimensions_.set((int)brickDim); */

        if (!brickDimensions_.isSelected("treeDepth")) {
            tgt::svec3 brickDim(brickDimensions_.getValue());
            size_t treeDepth = computeTreeDepth(brickDim, volumeDim);
            treeDepth_.set((int)treeDepth);
        }

        bool inCache = false;
        if (VoreenApplication::app()->useCaching()) {
            std::string octreeFile = tgt::FileSystem::cleanupPath(getOctreeStoragePath() + "/octree.xml");
            inCache = tgt::FileSystem::fileExists(octreeFile);
        }
        
        if (inCache)
            updateStatusMessage("Octree found in cache!");
        else
            updateStatusMessage("Octree not found in cache!");
    }
    else {
        updateStatusMessage("No input data.");
    }
}



void OctreeCreator::updateStatusMessage(const std::string& message) const {
    // Enqueue new command for an ui update.
    VoreenApplication::app()->getCommandQueue()->enqueue(this, LambdaFunctionCallback([this, message] {
        statusMessage_.set(message);
    }));
}

// statics

void OctreeCreator::limitCacheSize(const std::string& cacheDirectory, const uint64_t maxCacheSize, bool keepLatest) {
    std::string octreeCachePath = tgt::FileSystem::cleanupPath(cacheDirectory + "/" + CACHE_SUBDIR);
    if (!tgt::FileSystem::dirExists(octreeCachePath))
        return;

    if (maxCacheSize == 0)
        return;

    uint64_t cacheSize = tgt::FileSystem::dirSize(octreeCachePath);
    if (cacheSize <= maxCacheSize)
        return;

    LINFO("Octree cache size (" << formatMemorySize(cacheSize) << ") exceeds limit (" << formatMemorySize(maxCacheSize) << "). " << "Cleaning cache...");

    // collect octree directories and their access time in chronologically ordered list
    std::list< std::pair<std::string, time_t> > octreeDirs;
    std::vector<std::string> subDirs = tgt::FileSystem::listSubDirectories(octreeCachePath);
    for (size_t i=0; i<subDirs.size(); i++) {
        std::string octreeDir = tgt::FileSystem::cleanupPath(octreeCachePath + "/" + subDirs.at(i));
        std::string octreeFile = octreeDir + "/" + OCTREE_FILENAME;
        if (tgt::FileSystem::fileExists(octreeFile)) {
            time_t octreeTime = tgt::FileSystem::fileAccessTime(octreeFile);
            for (std::list< std::pair<std::string, time_t> >::iterator it = octreeDirs.begin(); ;++it) {
                if (it == octreeDirs.end() || it->second > octreeTime) {
                    octreeDirs.insert(it, std::make_pair(octreeDir, octreeTime));
                    break;
                }
            }
        }
        else { // invalid directory (no octree.xml) => delete first
            octreeDirs.insert(octreeDirs.begin(), std::make_pair(octreeDir, 0));
        }
    }

    // delete collected octree dir in chronological order until cache does not exceed the limit anymore
    size_t minCacheElems = (keepLatest ? 1 : 0);
    while (cacheSize > maxCacheSize && (octreeDirs.size() > minCacheElems) ) {
        std::string octreeDir = octreeDirs.front().first;
        uint64_t dirSize = tgt::FileSystem::dirSize(octreeDir);
        LINFO("Deleting cached octree: " << octreeDir << " (size: " << formatMemorySize(dirSize) << ")");
        if (!tgt::FileSystem::deleteDirectoryRecursive(octreeDir))
            LWARNING("Failed to delete octree directory: " << octreeDir);
        octreeDirs.pop_front();
        cacheSize -= std::min(dirSize, cacheSize);
    }
}

void OctreeCreator::deleteCache(const std::string& cacheDirectory) {
    std::string octreeCachePath = tgt::FileSystem::cleanupPath(cacheDirectory + "/" + CACHE_SUBDIR);
    if (tgt::FileSystem::dirExists(octreeCachePath)) {
        LINFOC("voreen.OctreeCreator", "Clearing octree cache directory: " << octreeCachePath);
        if (!tgt::FileSystem::deleteDirectoryRecursive(octreeCachePath))
            LWARNINGC("voreen.OctreeCreator", "Failed to delete octree cache directory: " << octreeCachePath);
    }
}

//save function
/*void OctreeCreator::saveOctreeToVVOD() {
    const std::string octreeFile = tgt::FileSystem::cleanupPath(getOctreeStoragePath() + "/" + OCTREE_FILENAME);
    if(!volumeOutport_.hasData() || !tgt::FileSystem::fileExists(octreeFile)) {
        LWARNING("No octree has been created, nothing to save!");
        return;
    }
    if(saveOctreeFile_.get().empty()) {
        LWARNING("No save file name selected!");
        return;
    }

    //determine prefix
    std::string prefix = tgt::FileSystem::baseName(saveOctreeFile_.get());
    std::string oldNodeBufferString, newNodeBufferString;
    std::string oldBrickBufferDir, newBrickBufferString;
    //copy octree.xml
    tgt::FileSystem::copyFile(octreeFile,saveOctreeFile_.get());
        //alter xml
        TiXmlDocument doc(saveOctreeFile_.get());
        if(!doc.LoadFile()) {
            LERROR(saveOctreeFile_.get() + " could not be loaded!");
            return;
        }
        TiXmlHandle curRoot(&doc);
        TiXmlElement* pElem;
        pElem = curRoot.FirstChild("VoreenData").Element();
        if (!pElem) { LERROR("No VoreenData element in " + saveOctreeFile_.get()); return; }
        else { curRoot = TiXmlHandle(pElem); }
        pElem = curRoot.FirstChild("Octree").Element();
        if (!pElem) { LERROR("No Octree element in " + saveOctreeFile_.get()); return; }
        else { curRoot = TiXmlHandle(pElem); }
        //rename nodebuffer
        pElem = curRoot.FirstChild("nodeBufferName").Element();
        if (!pElem) { LERROR("No nodeBufferName element in " + saveOctreeFile_.get()); return; }
        else {
            oldNodeBufferString = std::string(pElem->Attribute("value"));
            newNodeBufferString = std::string(prefix + "_" + oldNodeBufferString);
            pElem->SetAttribute("value",newNodeBufferString);
        }
        //rename brickbuffer
        pElem = curRoot.FirstChild("brickPoolManager").Element();
        if (!pElem) { LERROR("No brickPoolManager element in " + saveOctreeFile_.get()); return; }
        else { curRoot = TiXmlHandle(pElem); }
        pElem = curRoot.FirstChild("brickPoolPath").Element();
        if (!pElem) { LERROR("No brickPoolPath element in " + saveOctreeFile_.get()); return; }
        else {
            oldBrickBufferDir =  tgt::FileSystem::cleanupPath(tgt::FileSystem::parentDir(octreeFile) + "/" + pElem->Attribute("value"), true);
            newBrickBufferString = prefix + "_" + BRICK_BUFFER_SUBDIR;
            pElem->SetAttribute("value","../"+ newBrickBufferString);
        }
        //rename bricks
        pElem = curRoot.FirstChild("bufferFiles").Element();
        if (!pElem) { LERROR("No bufferFiles element in " + saveOctreeFile_.get()); return; }
        else { curRoot = TiXmlHandle(pElem); }
        pElem = curRoot.FirstChild("item").Element();
        size_t pos = 0;
        if(pElem) {
            pos = std::string(pElem->Attribute("value")).find(BRICK_BUFFER_SUBDIR + "/" + BRICK_BUFFER_FILE_PREFIX);
            tgtAssert(pos != std::string::npos,"unknown buffer file name");
        }
        for( pElem; pElem; pElem=pElem->NextSiblingElement()) {
            pElem->SetAttribute("value","../"+ prefix + "_" + std::string(pElem->Attribute("value")).substr(pos));
        }
        //save file
        if(!doc.SaveFile()) {
            LERROR(saveOctreeFile_.get() + " could not be saved!");
            return;
        }
    //copy nodebuffer.raw
    std::string oldNodeBufferPath = tgt::FileSystem::parentDir(octreeFile) + "/" + oldNodeBufferString;
    std::string newNodeBufferPath = tgt::FileSystem::parentDir(saveOctreeFile_.get()) + "/" + newNodeBufferString;
    if(tgt::FileSystem::fileExists(oldNodeBufferPath))
        tgt::FileSystem::copyFile(oldNodeBufferPath,newNodeBufferPath);
    else {
        LERROR("missing nodebuffer.raw");
        return;
    }
    //copy brickBuffer
    std::string newBrickBufferDir = tgt::FileSystem::parentDir(saveOctreeFile_.get()) + "/" + newBrickBufferString;
    tgt::FileSystem::createDirectory(newBrickBufferDir);
    tgt::FileSystem::clearDirectory(newBrickBufferDir);
    if(tgt::FileSystem::dirExists(oldBrickBufferDir)) {
        std::vector<std::string> buffers = tgt::FileSystem::listFiles(oldBrickBufferDir);
        progressProperty_.setProgress(0.f);
        for(size_t i = 0; i < buffers.size(); i++) {
            tgt::FileSystem::copyFile(oldBrickBufferDir + "/" + buffers[i], newBrickBufferDir + "/" + buffers[i]);
            progressProperty_.setProgress(static_cast<float>(i+1)/static_cast<float>(buffers.size()));
        }
    } else {
        LERROR("Missing brick buffer directory: " << oldBrickBufferDir);
        return;
    }
    progressProperty_.setProgress(1.f);
    LINFO("Saving octree complete");
}*/

}   // namespace
