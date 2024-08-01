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

#include "volumelistmodalityfilter.h"
#include "voreen/core/datastructures/callback/memberfunctioncallback.h"

namespace voreen {

const std::string VolumeListModalityFilter::loggerCat_ = "VolumeListModalityFilter";

VolumeListModalityFilter::VolumeListModalityFilter()
    : Processor(),
    inport_(Port::INPORT, "volumecollection", "VolumeList Input"),
    outport_(Port::OUTPORT, "volumecollection.filtered", "VolumeList Output"),
    modalityProp_("modality", "modality: ")
{
    addPort(inport_);
    ON_CHANGE(inport_, VolumeListModalityFilter, updateModalityOptions)
    addPort(outport_);

    addProperty(modalityProp_);
    updateModalityOptions();
}

Processor* VolumeListModalityFilter::create() const {
    return new VolumeListModalityFilter();
}

void VolumeListModalityFilter::process() {

    const VolumeList* collection = inport_.getData();
    if (!collection || collection->empty()) {
        outport_.setData(nullptr);
        return;
    }

    if(modalityProp_.getValue() == Modality::MODALITY_ANY) {
        outport_.setData(collection, false);
        return;
    }

    VolumeList* filteredList = new VolumeList();
    for (size_t i = 0; i < collection->size(); ++i) {
        if (collection->at(i)->getModality() == modalityProp_.getValue()) {
            filteredList->add(collection->at(i));
        }
    }
    outport_.setData(filteredList, true);
}

void VolumeListModalityFilter::updateModalityOptions() {

    // HACK: The way modalities are designed, we have to call their constructor in order to add them to the internal list.
    if(auto* volumes = inport_.getData()) {
        for(size_t i=0; i<volumes->size(); i++) {
            volumes->at(i)->getModality(); // This adds them to the list.
        }
    }

    bool empty = modalityProp_.getOptions().empty();
    Modality selectedModality = empty ? Modality::MODALITY_ANY : modalityProp_.getValue();

    modalityProp_.setOptions({});
    const std::vector<Modality>& modalities = Modality::getModalities();
    for (size_t i = 0; i < modalities.size(); ++i) {
        modalityProp_.addOption(modalities[i].getName(), modalities[i].getName(), modalities[i]);
    }

    modalityProp_.set(selectedModality.getName());
};

}   // namespace voreen
