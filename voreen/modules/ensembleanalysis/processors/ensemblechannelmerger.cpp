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

#include "ensemblechannelmerger.h"

#include "voreen/core/datastructures/callback/lambdacallback.h"
#include "voreen/core/datastructures/volume/volumediskmultichanneladapter.h"

#include "../utils/ensemblehash.h"

namespace voreen {

EnsembleChannelMerger::EnsembleChannelMerger()
    : Processor()
    , ensembleInport_(Port::INPORT, "ensembledatastructurein", "Ensemble Datastructure Input", false)
    , ensembleOutport_(Port::OUTPORT, "ensembledatastructureout", "Ensemble Datastructure Output", false)
    , fields_("fields", "Fields", false, Processor::VALID, Property::LOD_DEFAULT, true)
    , apply_("apply", "Apply")
    , reset_("reset", "Reset")
    , autoUpdate_("autoUpdate", "Auto. Update", false)
{
    addPort(ensembleInport_);
    ON_CHANGE(ensembleInport_, EnsembleChannelMerger, adjustToEnsemble);
    addPort(ensembleOutport_);

    addProperty(fields_);
    ON_CHANGE_LAMBDA(fields_, [this] {
        if(!ensembleInport_.hasData()) {
            return;
        }

        const auto* ensemble = ensembleInport_.getData();

        auto& instances = fields_.getInstances();
        if(!instances.empty()) {
            auto referenceDimension = ensemble->getFieldMetaData(instances.front().getName()).dimensions_;
            for(auto& instance : instances) {
                instance.setActive(ensemble->getFieldMetaData(instance.getName()).dimensions_ == referenceDimension);
            }
        }
    });
    addProperty(apply_);
    ON_CHANGE(apply_, EnsembleChannelMerger, merge);
    addProperty(reset_);
    ON_CHANGE_LAMBDA(reset_, [this] { fields_.reset(); });
    addProperty(autoUpdate_);
}

EnsembleChannelMerger::~EnsembleChannelMerger() {
}

Processor* EnsembleChannelMerger::create() const {
    return new EnsembleChannelMerger();
}

void EnsembleChannelMerger::process() {
    if(autoUpdate_.get()) {
        merge();
    }
}

void EnsembleChannelMerger::merge() {

    ensembleOutport_.clear();
    volumes_.clear();

    // Trivial case.
    if(fields_.getInstances().size() <= 1) {
        ensembleOutport_.setData(ensembleInport_.getData(), false);
        return;
    }

    std::unique_ptr<EnsembleDataset> ensemble(new EnsembleDataset());
    const auto& members = ensembleInport_.getData()->getMembers();
    for(size_t i=0; i<members.size(); i++) {

        std::vector<TimeStep> timeSteps;

        for(size_t j=0; j<members[i].getTimeSteps().size(); j++) {
            std::map<std::string, const VolumeBase*> volumes;

            // Add untouched volumes.
            for(const auto& item: fields_.getItems()) {
                if(!fields_.hasInstance(item)) {
                    volumes[item] = members[i].getTimeSteps()[j].getVolume(item);
                }
            }

            std::string fieldName = "";

            std::vector<const VolumeBase*> channels;
            for(const auto& instance: fields_.getInstances()) {
                const VolumeBase* channel = members[i].getTimeSteps()[j].getVolume(instance.getName());
                channels.push_back(channel);

                fieldName += instance.getName();
            }

            // Gather min and max values.
            std::vector<float> minValues;
            std::vector<float> maxValues;
            float minMagnitude = 0.0f;
            float maxMagnitude = 0.0f;
            for(size_t k=0; k<channels.size(); k++) {
                auto vmm = channels[k]->getDerivedData<VolumeMinMax>();
                minValues.push_back(vmm->getMin());
                maxValues.push_back(vmm->getMax());
                minMagnitude += vmm->getMin() * vmm->getMin();
                maxMagnitude += vmm->getMax() * vmm->getMax();
            }

            // For the magnitude, this just represents the theoretically lowest and highest
            // achievable magnitude, this might not actually occur in one of the voxels.
            // But without this approximation, these values would need to be calculated
            // which would be very slow.
            minMagnitude = std::sqrt(minMagnitude);
            maxMagnitude = std::sqrt(maxMagnitude);

            VolumeDisk* vd = new VolumeDiskMultiChannelAdapter(channels);
            VolumeBase* volume = new Volume(vd, channels.front());
            volume->addDerivedData(new VolumeMinMax(minValues, maxValues, minValues, maxValues));
            volume->addDerivedData(new VolumeMinMaxMagnitude(minMagnitude, maxMagnitude, minMagnitude, maxMagnitude));
            volumes_.emplace_back(std::unique_ptr<VolumeBase>(volume));

            volumes[fieldName] = volume;

            TimeStep timeStep(volumes, members[i].getTimeSteps()[j].getTime(), false);
            timeSteps.emplace_back(timeStep);
        }

        ensemble->addMember(EnsembleMember(members[i].getName(), members[i].getColor(), timeSteps));
    }

    ensembleOutport_.setData(ensemble.release(), true);
}

void EnsembleChannelMerger::adjustToEnsemble() {

    const auto* ensemble = ensembleInport_.getData();
    if(ensemble) {

        std::string hash = EnsembleHash(*ensemble).getHash();
        if(hash != hash_) {

            fields_.clear();

            for (auto field: ensemble->getCommonFieldNames()) {
                if(ensemble->getFieldMetaData(field).hasHomogeneousDimensions()) {
                    fields_.addItem(field);
                }
            }

            hash_ = hash;
        }
    }
}

void EnsembleChannelMerger::serialize(Serializer& s) const {
    Processor::serialize(s);
    s.serialize("hash", hash_);
}

void EnsembleChannelMerger::deserialize(Deserializer& s) {
    Processor::deserialize(s);
    s.optionalDeserialize("hash", hash_, std::string(""));
}

} // namespace
