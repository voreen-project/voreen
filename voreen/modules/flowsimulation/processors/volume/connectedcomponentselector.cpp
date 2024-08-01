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

#include "connectedcomponentselector.h"

#include "voreen/core/ports/conditions/portconditionvolumetype.h"
#include "voreen/core/datastructures/volume/volumeminmax.h"

namespace {

template<typename T>
voreen::VolumeAtomic<T>* selectComponent(const voreen::VolumeAtomic<T>* components, const std::vector<int>& selectedComponents) {
    voreen::VolumeAtomic<T>* output = components->clone();

    std::unordered_set<T> selectedIds(selectedComponents.begin(), selectedComponents.end());
    const T emptyId = static_cast<T>(0);
    const T minusOne = static_cast<T>(-1); // Member indices start counting at 0, components at 1.

    for(size_t i=0; i<output->getNumVoxels(); i++) {
        if(selectedIds.find(output->voxel(i) + minusOne) == selectedIds.end()) {
            output->voxel(i) = emptyId;
        }
    }

    return output;
}

}

namespace voreen {

const std::string ConnectedComponentSelector::loggerCat_("voreen.ensembleanalysis.ConnectedComponentSelector");

ConnectedComponentSelector::ConnectedComponentSelector()
    : Processor()
    , inport_(Port::INPORT, "volume.input", "Connected Components")
    , outport_(Port::OUTPORT, "volume.output", "Selected Component")
    , components_("components", "Selected Components")
{
    addPort(inport_);
        inport_.addCondition(new PortConditionVolumeTypeInteger());
        inport_.addCondition(new PortConditionVolumeChannelCount(1));
    addPort(outport_);

    addProperty(components_);
}

ConnectedComponentSelector::~ConnectedComponentSelector() {
}

Processor* ConnectedComponentSelector::create() const {
    return new ConnectedComponentSelector();
}

void ConnectedComponentSelector::process() {

    VolumeRAMRepresentationLock componentsVolume(inport_.getData());

    VolumeRAM* output = nullptr;
    if(const auto* components = dynamic_cast<const VolumeRAM_UInt8*>(*componentsVolume)) {
        output = selectComponent(components, components_.get());
    }
    else if(const auto* components = dynamic_cast<const VolumeRAM_UInt16*>(*componentsVolume)) {
        output = selectComponent(components, components_.get());
    }
    else if(const auto* components = dynamic_cast<const VolumeRAM_UInt32*>(*componentsVolume)) {
        output = selectComponent(components, components_.get());
    }
    else {
        LERROR("Unsupported input volume");
        return;
    }

    Volume* volume = new Volume(output, inport_.getData());
    outport_.setData(volume, true);
}

void ConnectedComponentSelector::adjustPropertiesToInput() {
    if(!inport_.isReady()) {
        return;
    }

    components_.reset();
    VolumeMinMax* vmm = inport_.getData()->getDerivedData<VolumeMinMax>();

    // Remove the "zero" aka "empty" component. We assume only positive ids!
    int firstId = std::max(vmm->getMin(), 1.0f);
    int lastId = vmm->getMax();
    std::vector<int> selected;
    for(int id = firstId; id <= lastId; id++) {
        components_.addRow(std::to_string(id));
        selected.push_back(static_cast<int>(selected.size()));
    }

    // Select all components per default.
    components_.setSelectedRowIndices(selected);
}

} // namespace
