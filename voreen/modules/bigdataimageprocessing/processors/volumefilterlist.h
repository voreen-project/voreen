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

#ifndef VRN_VOLUMEFILTERLIST_H
#define VRN_VOLUMEFILTERLIST_H

#include "voreen/core/processors/asynccomputeprocessor.h"
#include "voreen/core/properties/temppathproperty.h"
#include "modules/base/properties/interactivelistproperty.h"
#include "../volumefiltering/slicereader.h"
#include "../volumefilterproperties/templatefilterproperties.h"

namespace voreen {

class SliceReader;
class VolumeFilter;
class FilterProperties;

struct VolumeFilterListInput {
    std::unique_ptr<SliceReader> sliceReader;
    std::unique_ptr<HDF5FileVolume> outputVolume;

    VolumeFilterListInput(std::unique_ptr<SliceReader>&& pSliceReader, std::unique_ptr<HDF5FileVolume>&& pOutputVolume)
    : sliceReader(std::move(pSliceReader))
    , outputVolume(std::move(pOutputVolume))
    {
    }

    VolumeFilterListInput(const VolumeFilterListInput&) = delete;
    VolumeFilterListInput(VolumeFilterListInput&& old)
    : sliceReader(old.sliceReader.release())
    , outputVolume(old.outputVolume.release())
    {
    }
};

struct VolumeFilterListOutput {
    std::string outputVolumeFilePath;
};

/**
 * Applies multiple filters onto a volume.
 */
class VRN_CORE_API VolumeFilterList : public AsyncComputeProcessor<VolumeFilterListInput, VolumeFilterListOutput> {
public:
    VolumeFilterList();
    virtual ~VolumeFilterList();

    Processor* create() const;

    std::string getClassName() const { return "VolumeFilterList"; }

    std::string getCategory() const { return "Volume Processing"; }

    CodeState getCodeState() const { return CODE_STATE_TESTING; }

    virtual bool isReady() const;
    virtual bool isEndProcessor() const       { return true; }

    /** @see Property::serialize */
    virtual void serialize(Serializer& s) const;
    /** @see Property::deserialize */
    virtual void deserialize(Deserializer& s);


    virtual ComputeInput prepareComputeInput();
    virtual ComputeOutput compute(ComputeInput input, ProgressReporter& progressReporter) const;
    virtual void processComputeOutput(ComputeOutput output);

protected:
    void setDescriptions() {
        setDescription(
                "This processor allows to apply multiple filters onto an input volume."
                "<br>"
                "The left list provides all available filters. "
                "Use drag and drop or the keyboard to add/move/delete filters. "
                "The right list provides all added filters. "
                "Clicking on one of them allows to configure the filter's settings."
                "<br><br>"
                "Note that the order of application is from top to bottom!"
                "<br><br>"
                "<b>Extent:</b> Some filters have <i>extent</i>-Properties which are used to calculate, but in themselves distinct from kernel sizes. "
                "The extent describes many voxels a kernel extends from the central voxel in either direction. "
                "In other words: <i>kernel_size</i> = 2 * <i>extent</i> + 1"
        );
    }

    virtual void adjustPropertiesToInput();
    virtual void initialize();
    virtual void deinitialize();
    virtual void dataWillChange(const Port* source);

private:
    void storeInstance(InteractiveListProperty::Instance&);
    void restoreInstance(InteractiveListProperty::Instance&);
    bool hasConfiguredFilters() const;

    void onFilterListChange();
    void onFilterPropertyChange(Property* property);
    void inputOutputChannelCheck();

    template<typename Settings>
    void addFilter();
    void disableTracking(Property* property);

    /// Filter list
    std::vector<std::unique_ptr<FilterProperties>> filterProperties_;

    VolumePort inport_;
    VolumePort outport_;

    // General properties
    BoolProperty enabled_;
    TempPathProperty outputVolumeFilePath_;
    IntProperty outputVolumeDeflateLevel_;

    InteractiveListProperty filterList_;
    boost::optional<InteractiveListProperty::Instance> selectedInstance_;
    size_t numInstances_;

    bool skipPropertySync_;

    // Disable properties when processor is not enabled:
    PropertyDisabler propertyDisabler_;

    static const std::string loggerCat_; ///< category used in logging
};

template<typename Settings>
void VolumeFilterList::addFilter() {
    FilterProperties* filterProperties = new TemplateFilterProperties<Settings>();

    filterList_.addItem(filterProperties->getVolumeFilterName());
    filterProperties_.push_back(std::unique_ptr<FilterProperties>(filterProperties));
    for(Property* property : filterProperties->getProperties()) {
        addProperty(property);
        disableTracking(property);
        property->setGroupID(filterProperties->getVolumeFilterName());
        ON_CHANGE_LAMBDA((*property), ([this,property] () { this->onFilterPropertyChange(property); }));
    }
    filterProperties->storeVisibility();
    setPropertyGroupGuiName(filterProperties->getVolumeFilterName(), filterProperties->getVolumeFilterName());
    setPropertyGroupVisible(filterProperties->getVolumeFilterName(), false);
}

}

#endif // VRN_VOLUMEFILTERLIST_H
