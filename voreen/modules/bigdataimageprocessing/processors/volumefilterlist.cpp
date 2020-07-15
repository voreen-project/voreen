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

#include "volumefilterlist.h"

#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/floatproperty.h"

#include "modules/hdf5/io/hdf5volumereader.h"
#include "modules/hdf5/io/hdf5volumewriter.h"
#include "modules/hdf5/io/hdf5filevolume.h"

// Include actual Volume Filters.
#include "../volumefiltering/binarizationfilter.h"
#include "../volumefiltering/binarymedianfilter.h"
#include "../volumefiltering/gaussianfilter.h"
#include "../volumefiltering/gradientfilter.h"
#include "../volumefiltering/medianfilter.h"
#include "../volumefiltering/morphologyfilter.h"
#include "../volumefiltering/resamplefilter.h"
#include "../volumefiltering/rescalefilter.h"
#include "../volumefiltering/thresholdingfilter.h"
#include "../volumefiltering/valuemapfilter.h"
#include "../volumefiltering/vorticityfilter.h"

// Include their properties.
#include "../volumefilterproperties/binarizationfilterproperties.h"
#include "../volumefilterproperties/binarymedianfilterproperties.h"
#include "../volumefilterproperties/gaussianfilterproperties.h"
#include "../volumefilterproperties/gradientfilterproperties.h"
#include "../volumefilterproperties/medianfilterproperties.h"
#include "../volumefilterproperties/morphologyfilterproperties.h"
#include "../volumefilterproperties/resamplefilterproperties.h"
#include "../volumefilterproperties/rescalefilterproperties.h"
#include "../volumefilterproperties/thresholdingfilterproperties.h"
#include "../volumefilterproperties/valuemapfilterproperties.h"
#include "../volumefilterproperties/vorticityfilterproperties.h"

namespace voreen {

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
    , propertyDisabler_(*this)
    , skipPropertySync_(false)
{
    addPort(inport_);
        ON_CHANGE(inport_, VolumeFilterList, inputOutputChannelCheck);
    addPort(outport_);

    addProperty(filterList_);
        filterList_.setGroupID("filter");
        filterList_.setDuplicationAllowed(true);
        ON_CHANGE(filterList_, VolumeFilterList, onFilterListChange);
    setPropertyGroupGuiName("filter", "Filter");

    // Add filters (this will add their properties!)
    // Note: The items will appear in the order below.
    // Reordering and removal of single items is possible.
    addFilter<BinarizationFilterSettings>();
    addFilter<BinaryMedianFilterSettings>();
    addFilter<GaussianFilterSettings>();
    addFilter<GradientFilterSettings>();
    addFilter<MedianFilterSettings>();
    addFilter<MorphologyFilterSettings>();
    addFilter<ResampleFilterSettings>();
    addFilter<RescaleFilterSettings>();
    addFilter<ThresholdingFilterSettings>();
    addFilter<ValueMapFilterSettings>();
    addFilter<VorticityFilterSettings>();

    // Technical stuff.
    addProperty(enabled_);
        enabled_.setGroupID("output");
        ON_CHANGE_LAMBDA(enabled_, [this] () {
            if(enabled_.get()) {
                this->outport_.setData(nullptr);
                propertyDisabler_.restore();
            } else {
                this->forceComputation();
                propertyDisabler_.saveState([this] (Property* p) { return p == &enabled_; });
                propertyDisabler_.disable();
            }
        });
    addProperty(outputVolumeFilePath_);
        outputVolumeFilePath_.setGroupID("output");
    addProperty(outputVolumeDeflateLevel_);
        outputVolumeDeflateLevel_.setGroupID("output");
    setPropertyGroupGuiName("output", "Output");

    propertyDisabler_.saveState([this] (Property* p) { return p == &enabled_; });
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

void VolumeFilterList::initialize() {
    AsyncComputeProcessor::initialize();

    for(auto& prop : filterProperties_) {
        prop->initialize();
    }
}

void VolumeFilterList::deinitialize() {

    for(auto& prop : filterProperties_) {
        prop->deinitialize();
    }
    AsyncComputeProcessor::deinitialize();
}

Processor* VolumeFilterList::create() const {
    return new VolumeFilterList();
}

void VolumeFilterList::serialize(Serializer& s) const {
    AsyncComputeProcessor<ComputeInput, ComputeOutput>::serialize(s);
    for(size_t i=0; i < filterProperties_.size(); i++) {
        filterProperties_[i]->serialize(s);
    }
}

void VolumeFilterList::deserialize(Deserializer& s) {
    for(size_t i=0; i < filterProperties_.size(); i++) {
        // In case a new filter was added, it won't be able to be deserialized.
        try {
            filterProperties_[i]->deserialize(s);

            // In a previous implementation all filterproperties (at least
            // those used by VolumeFilterList) had a default instance with id
            // -1. The implementation has now changed so that such an instance
            // is no longer required. Old workspaces may have serialized this
            // instance which would now count as an additional one. That can
            // problems, for example, when linking properties (see
            // VolumeFilterList::onFilterPropertyChange). For backwards
            // compatibility we therefore explicitly remove all former default
            // instances here:
            filterProperties_[i]->removeInstance(-1);

        } catch(SerializationException& e) {
            s.removeLastError();
            LWARNING("Failed to deserialize Filterproperty '" << filterProperties_[i]->getVolumeFilterName() << "': " << e.what());
        }
    }

    AsyncComputeProcessor<ComputeInput, ComputeOutput>::deserialize(s);

    // If we deleted an instance of an item the related properties remain set to the instance's values
    // until another instance is selected which restores proper values.
    // Here we enforce valid values from the last instance of each item to not accidentally serialize incorrect values.
    for(auto& instance : filterList_.getInstances()) {
        restoreInstance(instance);
    }

    inputOutputChannelCheck();
}

VolumeFilterListInput VolumeFilterList::prepareComputeInput() {
    if(!enabled_.get() || filterList_.getInstances().empty()) {
        outport_.setData(inport_.getData(), false);
        throw InvalidInputException("", InvalidInputException::S_IGNORE);
    }

    if(!inport_.hasData()) {
        throw InvalidInputException("No input", InvalidInputException::S_WARNING);
    }

    auto inputVolPtr = inport_.getThreadSafeData();
    const VolumeBase& inputVolume = *inputVolPtr;

    VolumeFilterStackBuilder builder(inputVolume);

    SliceReaderMetaData metadata = SliceReaderMetaData::fromVolume(inputVolume);
    for(const InteractiveListProperty::Instance& instance : filterList_.getInstances()) {

        if(!instance.isActive()) {
            LINFO("Filter: '" << instance.getName() << "' is not active. Skipping.");
            continue;
        }

        VolumeFilter* filter = filterProperties_[instance.getItemId()]->getVolumeFilter(metadata, instance.getInstanceId());
        tgtAssert(filter, "filter was null");

        metadata = filter->getMetaData(metadata);
        builder.addLayer(std::unique_ptr<VolumeFilter>(filter));
    }

    std::unique_ptr<SliceReader> sliceReader = builder.build(0);
    std::string baseType = sliceReader->getMetaData().getBaseType();
    size_t numOutputChannels = sliceReader->getNumChannels();

    // Reset output volume to make sure it (and the hdf5filevolume) are not used any more
    outport_.setData(nullptr);

    const std::string volumeFilePath = outputVolumeFilePath_.get();
    const std::string volumeLocation = HDF5VolumeWriter::VOLUME_DATASET_NAME;
    const tgt::svec3 dim = sliceReader->getDimensions();

    if(volumeFilePath.empty()) {
        throw InvalidInputException("No volume file path specified!", InvalidInputException::S_ERROR);
    }

    std::unique_ptr<HDF5FileVolume> outputVolume = nullptr;
    try {
        outputVolume = std::unique_ptr<HDF5FileVolume>(HDF5FileVolume::createVolume(volumeFilePath, volumeLocation, baseType, dim, numOutputChannels, true, outputVolumeDeflateLevel_.get(), tgt::svec3(dim.x, dim.y, 1), false));
    } catch(tgt::IOException& e) {
        throw InvalidInputException("Could not create output volume.", InvalidInputException::S_ERROR);
    }

    tgt::vec3 scale(tgt::vec3(inputVolume.getDimensions()) / tgt::vec3(dim));
    tgt::vec3 additionalOffset(scale * tgt::vec3(0.5) - tgt::vec3(0.5));

    outputVolume->writeSpacing(inputVolume.getSpacing() * scale);
    outputVolume->writeOffset(inputVolume.getOffset() + additionalOffset * inputVolume.getSpacing());
    outputVolume->writePhysicalToWorldTransformation(inputVolume.getPhysicalToWorldMatrix());

    const auto& srmd = sliceReader->getMetaData();
    auto vmm = srmd.getVolumeMinMax();
    if(vmm) {
        outputVolume->writeVolumeMinMax(vmm.get());
    }
    outputVolume->writeRealWorldMapping(srmd.getRealworldMapping());

    return VolumeFilterListInput(
            std::move(sliceReader),
            std::move(outputVolume)
    );
}
VolumeFilterListOutput VolumeFilterList::compute(VolumeFilterListInput input, ProgressReporter& progressReporter) const {
    tgtAssert(input.sliceReader, "No sliceReader");
    tgtAssert(input.outputVolume, "No outputVolume");

    writeSlicesToHDF5File(*input.sliceReader, *input.outputVolume, &progressReporter);

    return {
        input.outputVolume->getFileName()
    };
    //outputVolume will be destroyed and thus closed now.
}
void VolumeFilterList::processComputeOutput(VolumeFilterListOutput output) {
    // outputVolume has been destroyed and thus closed by now.
    // So we can open it again (and use HDF5VolumeReader's implementation to read all the metadata with the file)
    const VolumeBase* vol = HDF5VolumeReaderOriginal().read(output.outputVolumeFilePath)->at(0);
    outport_.setData(vol);
}

void VolumeFilterList::adjustPropertiesToInput() {

    const VolumeBase* input = inport_.getData();
    if(!input) {
        return;
    }

    // As we only adjust list properties, we first must sync towards list properties...
    if(selectedInstance_) {
        storeInstance(*selectedInstance_);
    }

    SliceReaderMetaData metadata = SliceReaderMetaData::fromVolume(*input);

    if(selectedInstance_) {
        for(const InteractiveListProperty::Instance& instance : filterList_.getInstances()) {

            filterProperties_[instance.getItemId()]->adjustPropertiesToInput(metadata, instance.getInstanceId());

            if(!instance.isActive()) {
                continue;
            }

            VolumeFilter* filter = filterProperties_[instance.getItemId()]->getVolumeFilter(metadata, instance.getInstanceId());

            tgtAssert(filter, "filter was null");

            metadata = filter->getMetaData(metadata);
        }
    }
    // ... and later restore the current
    if(selectedInstance_) {
        restoreInstance(*selectedInstance_);
    }
}

void VolumeFilterList::dataWillChange(const Port* source) {
    outport_.clear();
    AsyncComputeProcessor::dataWillChange(source);
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
            tgtAssert(numInstances_ == filterList_.getInstances().size() + 1, "Only single instance removal allowed!");
            setPropertyGroupVisible(filterList_.getItems()[selectedInstance_->getItemId()], false);
            filterProperties_[selectedInstance_->getItemId()]->removeInstance(selectedInstance_->getInstanceId());
            selectedInstance_.reset();
        }
        numInstances_ = filterList_.getInstances().size();
    }

    // Hide old group.
    if(selectedInstance_) {
        filterProperties_[selectedInstance_->getItemId()]->storeVisibility();
        // No need to store the settings here, since it is done on change anyways.
        setPropertyGroupVisible(filterList_.getItems()[selectedInstance_->getItemId()], false);

        // We need to reset here, because otherwise onFilterPropertyChange
        // will be triggered while the current instance is restored.
        selectedInstance_.reset();
    }

    // Show new group.
    boost::optional<InteractiveListProperty::Instance> currentInstance;
    if(filterList_.getSelectedInstance() != -1) {
        currentInstance = filterList_.getInstances()[filterList_.getSelectedInstance()];
        setPropertyGroupVisible(filterList_.getItems()[currentInstance->getItemId()], true);
        filterProperties_[currentInstance->getItemId()]->restoreVisibility();

        restoreInstance(*currentInstance);
    }

    selectedInstance_ = currentInstance;

    // Check, if channels match at filter interfaces.
    inputOutputChannelCheck();

    // Set min/max values etc. for new filters
    adjustPropertiesToInput();
}

void VolumeFilterList::storeInstance(InteractiveListProperty::Instance& instance) {
    filterProperties_[instance.getItemId()]->storeInstance(instance.getInstanceId());
}
void VolumeFilterList::restoreInstance(InteractiveListProperty::Instance& instance) {
    // Restoring the current instance of a kind of filterproperties invalidates
    // the associated properties.  As we _just_ restored them from a concrete
    // instance, we don't need to sync back to the FilterProperties class (see
    // check in onFilterPropertyChange). This avoids (among useless work)
    // tripping the "multiple-instance-linking check".
    skipPropertySync_ = true;
    filterProperties_[instance.getItemId()]->restoreInstance(instance.getInstanceId());
    skipPropertySync_ = false;
}
void VolumeFilterList::onFilterPropertyChange(Property* property) {
    if(skipPropertySync_) {
        return;
    }

    if(!isInitialized()) {
        return;
    }

    // If any filter property was modified, we need to store the settings
    // immediately.
    //
    // However, this function may also be called if (e.g. via linking) a
    // property of a different filter has been changed. So we only use the
    // selected filter if it does have a matching property.
    //
    // FIXME: This can still go wrong/unexpected if:
    // 1. There are more than one instance of a particular Filter
    // 2. One of the two instances is selected
    // 3. The corresponding property is linked in the expectation that the
    //    value of the OTHER instance changes.
    // Unfortunately, we cannot even warn in this case, because here we cannot
    // distinguish property changes through links and changes through user
    // interaction. It would be quite annoying to warn about linking in case of
    // multiple filter instances when the user has not linked anything and just
    // wants to change values manually.
    if(selectedInstance_) {
        FilterProperties& fp = *filterProperties_[selectedInstance_->getItemId()];
        auto& props = fp.getProperties();
        if(std::find(props.begin(), props.end(), property) != props.end()) {
            fp.storeInstance(selectedInstance_->getInstanceId());
            return;
        }
    }


    // Otherwise, search through all filterProperties and see if one matches
    // the changed property
    for(auto& fp_ptr : filterProperties_) {
        FilterProperties& fp = *fp_ptr;
        auto& props = fp.getProperties();
        if(std::find(props.begin(), props.end(), property) != props.end()) {

            // We found the corresponding type of filter.
            auto instances = fp.getStoredInstances();
            switch (instances.size()) {
                case 1: {
                    // If there is exactly one we are lucky and can (trivially)
                    // match the corresponding instance
                    fp.storeInstance(instances[0]);
                    break;
                }
                case 0: {
                    // There is no instance. Either something has gone horribly
                    // wrong or someone has linked a property without a
                    // corresponding filter instance (probably by mistake).
                    LWARNING("Property without corresponding filter instance changed.");
                    break;
                }
                default: {
                    // We have no chance of finding the corresponding filter.
                    // It does not make much sense to change one at random or
                    // both. Instead just warn the user about her/his
                    // configuration mistake.
                    LWARNING("Linking properties with multiple filter instances of the same type is not supported.");
                    break;
                }
            }
            return;
        }
    }
}

void VolumeFilterList::inputOutputChannelCheck() {
    if(inport_.hasData()) {
        const VolumeBase& volume = *inport_.getData();
        size_t numOutputChannels = volume.getNumChannels();
        SliceReaderMetaData metadata = SliceReaderMetaData::fromVolume(volume);

        for (InteractiveListProperty::Instance& instance : filterList_.getInstances()) {
            VolumeFilter* filter = filterProperties_[instance.getItemId()]->getVolumeFilter(metadata, instance.getInstanceId());
            tgtAssert(filter, "filter was null");

            if (numOutputChannels == filter->getNumInputChannels()) {
                instance.setActive(true);
                numOutputChannels = filter->getNumOutputChannels();
            }
            else if(instance.isActive()) {
                instance.setActive(false);
                LERROR("Input channel count of filter '" << instance.getName() << "' is not satisfied. Deactivating.");
            }

            metadata = filter->getMetaData(metadata);
        }
    }
    else {
        // Reset filter active state.
        for(InteractiveListProperty::Instance& instance : filterList_.getInstances()) {
            instance.setActive(true);
        }
    }

    // Don't invalidate here, since this will lead to infinite recursion.
    // We just need to update the widgets only, anyways.
    filterList_.updateWidgets();
}


void VolumeFilterList::disableTracking(Property* property) {
    // In case of numeric properties, we disable tracking to not clash with auto update.
    if (auto numericProperty = dynamic_cast<NumericPropertyBase*>(property)) {
        numericProperty->setTracking(false);
    }
}

}   // namespace
