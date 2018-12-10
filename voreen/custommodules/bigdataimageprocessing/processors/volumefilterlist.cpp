/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany,                        *
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

#include "volumefilterlist.h"

#include "modules/hdf5/io/hdf5volumereader.h"
#include "modules/hdf5/io/hdf5volumewriter.h"
#include "modules/hdf5/io/hdf5filevolume.h"

#include "../volumefiltering/binarymedianfilter.h"
#include "../volumefiltering/gaussianfilter.h"
#include "../volumefiltering/medianfilter.h"

#include <algorithm>

namespace voreen {

class FilterProperties : public Serializable {
public:

    static const std::string DEFAULT_SETTINGS;

    virtual const std::vector<Property*> getProperties() const {
        return properties_;
    }

    virtual std::string getVolumeFilterName() const = 0;
    virtual void adjustPropertiesToInput(const VolumeBase& input) = 0;
    virtual VolumeFilter* getVolumeFilter(const VolumeBase& volume, const std::string& name) const = 0;
    virtual void applyInstance(const std::string& name) = 0;
    virtual void storeInstance(const std::string& name) = 0;
    virtual void removeInstance(const std::string& name) = 0;

    virtual ~FilterProperties() {}

protected:

    virtual void addProperties() = 0;

    std::string getId(const std::string& id) {
        return getVolumeFilterName() + "_" + id;
    }

    std::vector<Property*> properties_;

    static const std::string loggerCat_;
};

const std::string FilterProperties::DEFAULT_SETTINGS = "__default__";
const std::string FilterProperties::loggerCat_ = "FilterProperties";

class BinaryMedianFilterProperties : public FilterProperties {
public:
    BinaryMedianFilterProperties()
        : extentX_(getId("extentx"), "Extent X", 1, 1, 10)
        , extentY_(getId("extenty"), "Extent Y", 1, 1, 10)
        , extentZ_(getId("extentz"), "Extent Z", 1, 1, 10)
        , outsideVolumeValue_(getId("outsideVolumeValue"), "Outside Volume Value", 0, 0, 1)
        , binarizationThreshold_(getId("binarizationThreshold"), "Threshold", 0.5f, 0.0f, std::numeric_limits<float>::max(), Processor::INVALID_RESULT, FloatProperty::STATIC, Property::LOD_ADVANCED)
        , samplingStrategyType_(getId("samplingStrategyType"), "Sampling Strategy", SamplingStrategyType::CLAMP_T)
        , objectVoxelThreshold_(getId("objectVoxelThreshold"), "Object Voxel Threshold", 0, 0, std::numeric_limits<int>::max())
    {
        samplingStrategyType_.addOption("clamp", "Clamp", SamplingStrategyType::CLAMP_T);
        samplingStrategyType_.addOption("mirror", "Mirror", SamplingStrategyType::MIRROR_T);
        samplingStrategyType_.addOption("set", "Set", SamplingStrategyType::SET_T);
        ON_CHANGE_LAMBDA(samplingStrategyType_, [this] () {
            outsideVolumeValue_.setVisibleFlag(samplingStrategyType_.getValue() == SamplingStrategyType::SET_T);
        });

        // Store default settings.
        storeInstance(DEFAULT_SETTINGS);

        // Add properties to list.
        addProperties();
    }

    virtual std::string getVolumeFilterName() const {
        return "Binary Median Filter";
    }

    virtual void adjustPropertiesToInput(const VolumeBase& input) {
        if(!input.hasDerivedData<VolumeMinMax>()) {
            LINFOC("voreen.base.VolumeFilter", "Calculating VolumeMinMax. This may take a while...");
        }
        const VolumeMinMax* mm = input.getDerivedData<VolumeMinMax>();

        binarizationThreshold_.setMinValue(mm->getMin());
        binarizationThreshold_.setMaxValue(mm->getMax());
        binarizationThreshold_.adaptDecimalsToRange(2);
    }

    virtual VolumeFilter* getVolumeFilter(const VolumeBase& volume, const std::string& name) const {
        Settings settings = instanceSettings_.at(DEFAULT_SETTINGS);
        if(instanceSettings_.find(name) != instanceSettings_.end()) {
            settings = instanceSettings_.at(name);
        }
        else {
            LWARNING("Filter '" + name + "' has not been configured yet");
        }
        RealWorldMapping rwm;
        if(volume.hasMetaData("RealWorldMapping")) {
            rwm = volume.getRealWorldMapping();
        }
        return new BinaryMedianFilter(
                tgt::ivec3(settings.extentX_, settings.extentY_, settings.extentZ_),
                rwm.realWorldToNormalized(settings.binarizationThreshold_),
                settings.objectVoxelThreshold_,
                SamplingStrategy<float>(settings.samplingStrategyType_, static_cast<float>(settings.outsideVolumeValue_))
        );
    }
    virtual void applyInstance(const std::string& name) {
        auto iter = instanceSettings_.find(name);
        if(iter == instanceSettings_.end()) {
            instanceSettings_[name] = instanceSettings_[DEFAULT_SETTINGS];
        }

        Settings settings = instanceSettings_[name];
        extentX_.set(settings.extentX_);
        extentY_.set(settings.extentY_);
        extentZ_.set(settings.extentZ_);
        outsideVolumeValue_.set(settings.outsideVolumeValue_);
        binarizationThreshold_.set(settings.binarizationThreshold_);
        samplingStrategyType_.selectByValue(settings.samplingStrategyType_);
        objectVoxelThreshold_.set(settings.objectVoxelThreshold_);
    }
    virtual void storeInstance(const std::string& name) {
        Settings& settings = instanceSettings_[name];
        settings.extentX_ = extentX_.get();
        settings.extentY_ = extentY_.get();
        settings.extentZ_ = extentZ_.get();
        settings.outsideVolumeValue_ = outsideVolumeValue_.get();
        settings.binarizationThreshold_ = binarizationThreshold_.get();
        settings.samplingStrategyType_ = samplingStrategyType_.getValue();
        settings.objectVoxelThreshold_ = objectVoxelThreshold_.get();
    }
    virtual void removeInstance(const std::string& name) {
        instanceSettings_.erase(name);
    }
    virtual void addProperties() {
        properties_.push_back(&extentX_);
        properties_.push_back(&extentY_);
        properties_.push_back(&extentZ_);
        properties_.push_back(&outsideVolumeValue_);
        properties_.push_back(&binarizationThreshold_);
        properties_.push_back(&samplingStrategyType_);
        properties_.push_back(&objectVoxelThreshold_);
    }
    virtual void serialize(Serializer& s) const {
        std::vector<std::string> names;
        std::vector<Settings> settings;
        for(const auto& pair : instanceSettings_) {
            names.push_back(pair.first);
            settings.push_back(pair.second);
        }
        s.serializeBinaryBlob("names", names);
        s.serializeBinaryBlob("settings", settings);
    }
    virtual void deserialize(Deserializer& s) {
        std::vector<std::string> names;
        std::vector<Settings> settings;
        s.deserializeBinaryBlob("names", names);
        s.deserializeBinaryBlob("settings", settings);
        tgtAssert(names.size() == settings.size(), "number of keys and values does not match");
        for(size_t i = 0; i < names.size(); i++) {
            instanceSettings_[names[i]] = settings[i];
        }
    }

private:

    struct Settings {
        int extentX_;
        int extentY_;
        int extentZ_;
        int outsideVolumeValue_;
        float binarizationThreshold_;
        SamplingStrategyType samplingStrategyType_;
        int objectVoxelThreshold_;
    };
    std::map<std::string, Settings> instanceSettings_;

    IntProperty extentX_;
    IntProperty extentY_;
    IntProperty extentZ_;
    IntProperty outsideVolumeValue_;
    FloatProperty binarizationThreshold_;
    OptionProperty<SamplingStrategyType> samplingStrategyType_;
    IntProperty objectVoxelThreshold_;
};

class MedianFilterProperties : public FilterProperties {
public:
    MedianFilterProperties()
            : extent_(getId("extent"), "Extent", 1, 1, 10)
            , outsideVolumeValue_(getId("outsideVolumeValue"), "Outside Volume Value", 0, 0, 1)
            , samplingStrategyType_(getId("samplingStrategyType"), "Sampling Strategy", SamplingStrategyType::CLAMP_T)
    {
        samplingStrategyType_.addOption("clamp", "Clamp", SamplingStrategyType::CLAMP_T);
        samplingStrategyType_.addOption("mirror", "Mirror", SamplingStrategyType::MIRROR_T);
        samplingStrategyType_.addOption("set", "Set", SamplingStrategyType::SET_T);
        ON_CHANGE_LAMBDA(samplingStrategyType_, [this] () {
            outsideVolumeValue_.setVisibleFlag(samplingStrategyType_.getValue() == SamplingStrategyType::SET_T);
        });

        // Store default settings.
        storeInstance(DEFAULT_SETTINGS);

        // Add properties to list.
        addProperties();
    }

    virtual std::string getVolumeFilterName() const {
        return "Median Filter";
    }

    virtual void adjustPropertiesToInput(const VolumeBase& input) {
        // unused
    }

    virtual VolumeFilter* getVolumeFilter(const VolumeBase& volume, const std::string& name) const {
        Settings settings = instanceSettings_.at(DEFAULT_SETTINGS);
        if(instanceSettings_.find(name) != instanceSettings_.end()) {
            settings = instanceSettings_.at(name);
        }
        else {
            LWARNING("Filter '" + name + "' has not been configured yet");
        }
        return new MedianFilter(
                settings.extent_,
                SamplingStrategy<float>(settings.samplingStrategyType_, static_cast<float>(settings.outsideVolumeValue_)),
                volume.getBaseType()
        );
    }
    virtual void applyInstance(const std::string& name) {
        auto iter = instanceSettings_.find(name);
        if(iter == instanceSettings_.end()) {
            instanceSettings_[name] = instanceSettings_[DEFAULT_SETTINGS];
        }

        Settings settings = instanceSettings_[name];
        extent_.set(settings.extent_);
        outsideVolumeValue_.set(settings.outsideVolumeValue_);
        samplingStrategyType_.selectByValue(settings.samplingStrategyType_);
    }
    virtual void storeInstance(const std::string& name) {
        Settings& settings = instanceSettings_[name];
        settings.extent_ = extent_.get();
        settings.outsideVolumeValue_ = outsideVolumeValue_.get();
        settings.samplingStrategyType_ = samplingStrategyType_.getValue();
    }
    virtual void removeInstance(const std::string& name) {
        instanceSettings_.erase(name);
    }
    virtual void addProperties() {
        properties_.push_back(&extent_);
        properties_.push_back(&outsideVolumeValue_);
        properties_.push_back(&samplingStrategyType_);
    }
    virtual void serialize(Serializer& s) const {
        std::vector<std::string> names;
        std::vector<Settings> settings;
        for(const auto& pair : instanceSettings_) {
            names.push_back(pair.first);
            settings.push_back(pair.second);
        }
        s.serializeBinaryBlob("names", names);
        s.serializeBinaryBlob("settings", settings);
    }
    virtual void deserialize(Deserializer& s) {
        std::vector<std::string> names;
        std::vector<Settings> settings;
        s.deserializeBinaryBlob("names", names);
        s.deserializeBinaryBlob("settings", settings);
        tgtAssert(names.size() == settings.size(), "number of keys and values does not match");
        for(size_t i = 0; i < names.size(); i++) {
            instanceSettings_[names[i]] = settings[i];
        }
    }

private:

    struct Settings {
        int extent_;
        int outsideVolumeValue_;
        SamplingStrategyType samplingStrategyType_;
    };
    std::map<std::string, Settings> instanceSettings_;

    IntProperty extent_;
    IntProperty outsideVolumeValue_;
    OptionProperty<SamplingStrategyType> samplingStrategyType_;
};

class GaussianFilterProperties : public FilterProperties {
public:
    GaussianFilterProperties()
        : extentX_(getId("extentx"), "Extent X", 1, 1, 10)
        , extentY_(getId("extenty"), "Extent Y", 1, 1, 10)
        , extentZ_(getId("extentz"), "Extent Z", 1, 1, 10)
        , outsideVolumeValue_(getId("outsideVolumeValue"), "Outside Volume Value", 0, 0, 1)
        , samplingStrategyType_(getId("samplingStrategyType"), "Sampling Strategy", SamplingStrategyType::CLAMP_T)
    {
        samplingStrategyType_.addOption("clamp", "Clamp", SamplingStrategyType::CLAMP_T);
        samplingStrategyType_.addOption("mirror", "Mirror", SamplingStrategyType::MIRROR_T);
        samplingStrategyType_.addOption("set", "Set", SamplingStrategyType::SET_T);

        // Store default settings.
        storeInstance(DEFAULT_SETTINGS);

        // Add properties to list.
        addProperties();
    }

    virtual std::string getVolumeFilterName() const {
        return "Gaussian Filter";
    }

    virtual void adjustPropertiesToInput(const VolumeBase&) {
        // unused
    }

    virtual VolumeFilter* getVolumeFilter(const VolumeBase& volume, const std::string& name) const {
        Settings settings = instanceSettings_.at(DEFAULT_SETTINGS);
        if(instanceSettings_.find(name) != instanceSettings_.end()) {
            settings = instanceSettings_.at(name);
        }
        else {
            LWARNING("Filter '" + name + "' has not been configured yet");
        }
        return new GaussianFilter(
                tgt::ivec3(settings.extentX_, settings.extentY_, settings.extentZ_),
                SamplingStrategy<float>(settings.samplingStrategyType_, static_cast<float>(settings.outsideVolumeValue_)),
                volume.getBaseType()
        );
    }
    virtual void applyInstance(const std::string& name) {
        auto iter = instanceSettings_.find(name);
        if(iter == instanceSettings_.end()) {
            instanceSettings_[name] = instanceSettings_[DEFAULT_SETTINGS];
        }

        Settings settings = instanceSettings_[name];
        extentX_.set(settings.extentX_);
        extentY_.set(settings.extentY_);
        extentZ_.set(settings.extentZ_);
        outsideVolumeValue_.set(settings.outsideVolumeValue_);
        samplingStrategyType_.selectByValue(settings.samplingStrategyType_);
    }
    virtual void storeInstance(const std::string& name) {
        Settings& settings = instanceSettings_[name];
        settings.extentX_ = extentX_.get();
        settings.extentY_ = extentY_.get();
        settings.extentZ_ = extentZ_.get();
        settings.outsideVolumeValue_ = outsideVolumeValue_.get();
        settings.samplingStrategyType_ = samplingStrategyType_.getValue();
    }
    virtual void removeInstance(const std::string& name) {
        instanceSettings_.erase(name);
    }
    virtual void addProperties() {
        properties_.push_back(&extentX_);
        properties_.push_back(&extentY_);
        properties_.push_back(&extentZ_);
        properties_.push_back(&outsideVolumeValue_);
        properties_.push_back(&samplingStrategyType_);
    }
    virtual void serialize(Serializer& s) const {
        std::vector<std::string> names;
        std::vector<Settings> settings;
        for(const auto& pair : instanceSettings_) {
            names.push_back(pair.first);
            settings.push_back(pair.second);
        }
        s.serializeBinaryBlob("names", names);
        s.serializeBinaryBlob("settings", settings);
    }
    virtual void deserialize(Deserializer& s) {
        std::vector<std::string> names;
        std::vector<Settings> settings;
        s.deserializeBinaryBlob("names", names);
        s.deserializeBinaryBlob("settings", settings);
        tgtAssert(names.size() == settings.size(), "number of keys and values does not match");
        for(size_t i = 0; i < names.size(); i++) {
            instanceSettings_[names[i]] = settings[i];
        }
    }

private:

    struct Settings {
        int extentX_;
        int extentY_;
        int extentZ_;
        int outsideVolumeValue_;
        SamplingStrategyType samplingStrategyType_;
    };
    std::map<std::string, Settings> instanceSettings_;

    IntProperty extentX_;
    IntProperty extentY_;
    IntProperty extentZ_;
    IntProperty outsideVolumeValue_;
    OptionProperty<SamplingStrategyType> samplingStrategyType_;
};



const std::string VolumeFilterList::loggerCat_("voreen.bigdataimageprocessing.VolumeFilterList");

VolumeFilterList::VolumeFilterList()
    : AsyncComputeProcessor<VolumeFilterListInput, VolumeFilterListOutput>()
    , inport_(Port::INPORT, "volumehandle.input", "Volume Input")
    , outport_(Port::OUTPORT, "volumehandle.output", "Volume Output", false)
    , enabled_("enabled", "Enabled", true)
    , outputVolumeFilePath_("outputVolumeFilePath", "Output Volume", "Path", "", "HDF5 (*.h5)", FileDialogProperty::SAVE_FILE, Processor::INVALID_RESULT, Property::LOD_DEFAULT)
    , outputVolumeDeflateLevel_("outputVolumeDeflateLevel", "Deflate Level", 1, 0, 9, Processor::INVALID_RESULT, IntProperty::STATIC, Property::LOD_DEFAULT)
    , filterList_("filterList", "Filter List", true)
    , numInstances_(0)
{
    addPort(inport_);
        ON_CHANGE(inport_, VolumeFilterList, adjustPropertiesToInput);
    addPort(outport_);

    addProperty(enabled_);
        enabled_.setGroupID("output");
    addProperty(outputVolumeFilePath_);
        outputVolumeFilePath_.setGroupID("output");
    addProperty(outputVolumeDeflateLevel_);
        outputVolumeDeflateLevel_.setGroupID("output");

    addProperty(filterList_);
        filterList_.setGroupID("output");
        filterList_.setDuplicationAllowed(true);
        ON_CHANGE(filterList_, VolumeFilterList, onFilterListChange);

    setPropertyGroupGuiName("output", "Output");

    // Add filters.
    addFilter(new BinaryMedianFilterProperties());
    addFilter(new MedianFilterProperties());
    addFilter(new GaussianFilterProperties());
}

VolumeFilterList::~VolumeFilterList() {
}

bool VolumeFilterList::isReady() const {
    if(!isInitialized()) {
        setNotReadyErrorMessage("Not initialized.");
        return false;
    }
    if(!inport_.isReady()) {
        setNotReadyErrorMessage("Inport not ready.");
        return false;
    }
    return true;
}

Processor* VolumeFilterList::create() const {
    return new VolumeFilterList();
}

void VolumeFilterList::serialize(Serializer& s) const {
    AsyncComputeProcessor<ComputeInput, ComputeOutput>::serialize(s);
    //s.serialize("filterProperties", filterProperties_);
}

void VolumeFilterList::deserialize(Deserializer& s) {
    AsyncComputeProcessor<ComputeInput, ComputeOutput>::deserialize(s);
    //s.deserialize("filterProperties", filterProperties_);
}

VolumeFilterListInput VolumeFilterList::prepareComputeInput() {
    if(!enabled_.get()) {
        return VolumeFilterListInput(
                nullptr,
                nullptr
        );
    }

    if(!inport_.hasData()) {
        throw InvalidInputException("No input", InvalidInputException::S_WARNING);
    }

    if(filterList_.getInstances().empty()) {
        throw InvalidInputException("No filter selected", InvalidInputException::S_ERROR);
    }

    auto inputVolPtr = inport_.getThreadSafeData();
    const VolumeBase& inputVolume = *inputVolPtr;

    if(inputVolume.getNumChannels() != 1) {
        throw InvalidInputException("Input volume has multiple channels, but a single channel volume is expected!", InvalidInputException::S_ERROR);
    }

    // Reset output volume to make sure it (and the hdf5filevolume) are not used any more
    outport_.setData(nullptr);

    const std::string volumeFilePath = outputVolumeFilePath_.get();
    const std::string volumeLocation = HDF5VolumeWriter::VOLUME_DATASET_NAME;
    const std::string baseType = "uint8";
    const tgt::svec3 dim = inputVolume.getDimensions();

    if(volumeFilePath.empty()) {
        throw InvalidInputException("No volume file path specified!", InvalidInputException::S_ERROR);
    }

    std::unique_ptr<HDF5FileVolume> outputVolume = nullptr;
    try {
        outputVolume = std::unique_ptr<HDF5FileVolume>(HDF5FileVolume::createVolume(volumeFilePath, volumeLocation, baseType, dim, 1, true, outputVolumeDeflateLevel_.get(), tgt::svec3(dim.x, dim.y, 1), false));
    } catch(tgt::IOException e) {
        throw InvalidInputException("Could not create output volume.", InvalidInputException::S_ERROR);
    }

    outputVolume->writeSpacing(inputVolume.getSpacing());
    outputVolume->writeOffset(inputVolume.getOffset());
    outputVolume->writeRealWorldMapping(RealWorldMapping(1,0,""));
    // For all zero or all one volumes the following is not correct,
    // and we cannot easily get the real min/max values without iterating
    // through the whole resulting volume.
    //const VolumeMinMax vmm(0, 1, 0, 1);
    //outputVolume->writeVolumeMinMax(&vmm);

    VolumeFilterStackBuilder builder(inputVolume);
    for(const InteractiveListProperty::Instance& instance : filterList_.getInstances()) {
        VolumeFilter* filter = filterProperties_[instance.itemId_]->getVolumeFilter(inputVolume, instance.name_);
        builder.addLayer(std::unique_ptr<VolumeFilter>(filter));
    }

    std::unique_ptr<SliceReader> sliceReader = builder.build(0);

    return VolumeFilterListInput(
            std::move(sliceReader),
            std::move(outputVolume)
    );
}
VolumeFilterListOutput VolumeFilterList::compute(VolumeFilterListInput input, ProgressReporter& progressReporter) const {
    if(!enabled_.get()) {
        return { "" };
    }
    tgtAssert(input.sliceReader, "No sliceReader");
    tgtAssert(input.outputVolume, "No outputVolume");

    writeSlicesToHDF5File(*input.sliceReader, *input.outputVolume, &progressReporter);

    return {
        input.outputVolume->getFileName()
    };
    //outputVolume will be destroyed and thus closed now.
}
void VolumeFilterList::processComputeOutput(VolumeFilterListOutput output) {
    if(!enabled_.get()) {
        outport_.setData(inport_.getData(), false);
    } else {
        // outputVolume has been destroyed and thus closed by now.
        // So we can open it again (and use HDF5VolumeReader's implementation to read all the metadata with the file)
        const VolumeBase* vol = HDF5VolumeReader().read(output.outputVolumeFilePath)->at(0);
        outport_.setData(vol);
    }
}

void VolumeFilterList::adjustPropertiesToInput() {
    const VolumeBase* input = inport_.getData();
    if(!input) {
        return;
    }

    for(auto& filterProperties : filterProperties_) {
        filterProperties->adjustPropertiesToInput(*input);
    }
}

// private methods
//

void VolumeFilterList::onFilterListChange() {

    // Check if instance was deleted.
    bool numInstancesChanged = filterList_.getInstances().size() != numInstances_;
    if(numInstancesChanged) {
        // Handle removal.
        if(numInstances_ > filterList_.getInstances().size() && selectedInstance_) {
            // Assumes that only the selected item can be removed!
            setPropertyGroupVisible(filterList_.getItems()[selectedInstance_->itemId_], false);
            filterProperties_[selectedInstance_->itemId_]->removeInstance(selectedInstance_->name_);
            selectedInstance_.reset();
        }
        numInstances_ = filterList_.getInstances().size();
    }

    // Hide old group.
    if(selectedInstance_) {
        filterProperties_[selectedInstance_->itemId_]->storeInstance(selectedInstance_->name_);
        setPropertyGroupVisible(filterList_.getItems()[selectedInstance_->itemId_], false);
    }

    // Show new group.
    boost::optional<InteractiveListProperty::Instance> currentInstance;
    if(filterList_.getSelectedInstance() != -1) {
        currentInstance = filterList_.getInstances()[filterList_.getSelectedInstance()];
        filterProperties_[currentInstance->itemId_]->applyInstance(currentInstance->name_);
        setPropertyGroupVisible(filterList_.getItems()[currentInstance->itemId_], true);
    }

    selectedInstance_ = currentInstance;
}

void VolumeFilterList::addFilter(FilterProperties* filterProperties) {
    filterList_.addItem(filterProperties->getVolumeFilterName());
    filterProperties_.push_back(std::unique_ptr<FilterProperties>(filterProperties));
    for(Property* property : filterProperties->getProperties()) {
        addProperty(property);
        property->setGroupID(filterProperties->getVolumeFilterName());
    }
    setPropertyGroupGuiName(filterProperties->getVolumeFilterName(), filterProperties->getVolumeFilterName());
    setPropertyGroupVisible(filterProperties->getVolumeFilterName(), false);
}

}   // namespace
