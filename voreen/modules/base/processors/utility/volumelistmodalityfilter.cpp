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

#include "volumelistmodalityfilter.h"
#include "voreen/core/datastructures/callback/memberfunctioncallback.h"

namespace voreen {

const std::string VolumeListModalityFilter::loggerCat_ = "VolumeListModalityFilter";

VolumeListModalityFilter::VolumeListModalityFilter()
    : Processor(),
    inport_(Port::INPORT, "volumecollection", "VolumeList Input"),
    outport_(Port::OUTPORT, "volumecollection.filtered", "VolumeList Output"),
    modalityProp_("modality", "modality: "),
    currentModality_(Modality::MODALITY_ANY),
    filteredList_()
{
    addPort(inport_);
    addPort(outport_);

    std::vector<Option<Modality> > options;
    const std::vector<Modality>& modalities = Modality::getModalities();
    for (size_t i = 0; i < modalities.size(); ++i) {
        modalityProp_.addOption(modalities[i].getName(), modalities[i].getName(), modalities[i]);
    }
    modalityProp_.set(Modality::MODALITY_ANY.getName());
    modalityProp_.onChange(MemberFunctionCallback<VolumeListModalityFilter>(this,
        &VolumeListModalityFilter::adjustFilteredList));
    addProperty(modalityProp_);
}

Processor* VolumeListModalityFilter::create() const {
    return new VolumeListModalityFilter();
}

void VolumeListModalityFilter::initialize() {
    Processor::initialize();
    adjustFilteredList();
}

void VolumeListModalityFilter::deinitialize() {
    outport_.setData(0);

    Processor::deinitialize();
}

void VolumeListModalityFilter::invalidate(int /*inv*/) {
    adjustFilteredList();
}

void VolumeListModalityFilter::process() {
    // nothing
}

// private methods
//

void VolumeListModalityFilter::adjustFilteredList() {
    const VolumeList* collection = inport_.getData();
    if ((collection == 0) || (collection->empty() == true)) {
        filteredList_.clear();
        outport_.setData(0);
        return;
    }

    if (currentModality_.getName() != modalityProp_.get()) {
        currentModality_ = modalityProp_.getValue();
        filteredList_.clear();
        if (currentModality_ != Modality::MODALITY_ANY) {
            for (size_t i = 0; i < collection->size(); ++i) {
                if (collection->at(i)->getModality() == currentModality_)
                    filteredList_.add(collection->at(i));
            }
            outport_.setData(&filteredList_, false);
        } else
            outport_.setData(const_cast<VolumeList*>(inport_.getData()), false);
    }
}

}   // namespace voreen
