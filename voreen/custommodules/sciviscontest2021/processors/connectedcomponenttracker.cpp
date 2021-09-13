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

#include "connectedcomponenttracker.h"

#include "voreen/core/ports/conditions/portconditionvolumetype.h"
#include "voreen/core/ports/conditions/portconditionvolumelist.h"
#include "voreen/core/datastructures/volume/volumeminmax.h"
#include "voreen/core/io/volumeserializerpopulator.h"
#include "voreen/core/io/volumeserializer.h"
#include "voreen/core/io/volumewriter.h"
#include "voreen/core/io/progressbar.h"

#include <queue>

namespace voreen {

const std::string ConnectedComponentTracker::loggerCat_("voreen.ensembleanalysis.ConnectedComponentTracker");

ConnectedComponentTracker::ConnectedComponentTracker()
    : Processor()
    , inport_(Port::INPORT, "volume.input", "Connected Components")
    , startButton_("startButton", "Start")
    , overlapThreshold_("overlapThreshold", "Overlap Threshold", 0.01f, 0.0f, 1.0f)
    , outputpath_("outputPath", "Output Path", "Output", "", "", FileDialogProperty::DIRECTORY)
{
    addPort(inport_);
    inport_.addCondition(new PortConditionVolumeListEnsemble());
    inport_.addCondition(new PortConditionVolumeListAdapter(new PortConditionVolumeChannelCount(1)));

    addProperty(startButton_);
    ON_CHANGE(startButton_, ConnectedComponentTracker, calculate);

    addProperty(overlapThreshold_);
    addProperty(outputpath_);
}

ConnectedComponentTracker::~ConnectedComponentTracker() {
}

Processor* ConnectedComponentTracker::create() const {
    return new ConnectedComponentTracker();
}

void ConnectedComponentTracker::calculate() {

    using T = uint32_t;
    const T BACKGROUND = 0;

    auto list = inport_.getData();

    if(list->size() < 2) {
        return;
    }


    double overlapThreshold = overlapThreshold_.get();
    std::unique_ptr<ProgressBar> progressBar(VoreenApplication::app()->createProgressDialog());
    progressBar->setTitle("Progress");
    progressBar->show();

    auto curr = list->first();
    auto currLock = VolumeRAMRepresentationLock(curr);

    auto rwm = curr->getRealWorldMapping();

    auto getId = [&rwm] (const VolumeRAM* volume, size_t idx) {
        float value = rwm.normalizedToRealWorld(volume->getVoxelNormalized(idx));
        T id = std::round(value);
        return id;
    };

    auto setId = [&rwm] (VolumeRAM* volume, size_t idx, T id) {
        float value = rwm.realWorldToNormalized(static_cast<float>(id));
        volume->setVoxelNormalized(value, idx);
    };

    T idCounter = T(0);
    std::queue<T> availableIds;

    struct Feature {

        struct Overlap {
            T originalId_;
            size_t numVoxels_;
        };

        T originalId_;
        T remappedId_;
        bool wasRemapped_;
        std::vector<size_t> voxels_;
        std::vector<Overlap> overlaps_;
//        size_t remappedByOriginalId_;

        Feature(T originalId, T id)
            : originalId_(originalId)
            , remappedId_(id)
            , wasRemapped_(false)
//            , remappedByOriginalId_(0)
        {
        }
    };

    auto generateId = [&] {
        T id;
        if(!availableIds.empty()) {
            id = availableIds.front();
            availableIds.pop();
        }
        else {
            id = ++idCounter;
        }
        return id;
    };

    auto initFeatures = [&] (const VolumeRAM* volume) {
        std::map<T, Feature> ids;

        for(size_t i = 0; i<volume->getNumVoxels(); i++) {
            T id = getId(volume, i);
            if(id == BACKGROUND) {
                continue;
            }

            if(ids.find(id) == ids.end()) {
                ids.insert(std::make_pair(id, Feature(id, generateId())));
            }

            ids.at(id).voxels_.push_back(i);
        }

        return ids;
    };

    auto currFeatures = initFeatures(*currLock);

    for (size_t i = 1; i < list->size(); i++) {

        auto next = list->at(i);
        auto nextLock = VolumeRAMRepresentationLock(next);
        auto nextFeatures = initFeatures(*nextLock);

        std::vector<Feature> currFeaturesSorted;
        currFeaturesSorted.reserve(currFeatures.size());
        for(auto& elem : currFeatures) {
            currFeaturesSorted.push_back(elem.second);
        }

        std::sort(currFeaturesSorted.begin(), currFeaturesSorted.end(), [] (const Feature& lhs, const Feature& rhs) {
            return lhs.voxels_.size() > rhs.voxels_.size();
        });

        for(auto& currFeature : currFeaturesSorted) {

            auto& feature = currFeatures.at(currFeature.originalId_);

            const auto& currVoxels = currFeature.voxels_;

            for(auto& nextFeature : nextFeatures) {

                const auto& nextVoxels = nextFeature.second.voxels_;

                std::vector<size_t> intersection;
                std::set_intersection(currVoxels.begin(), currVoxels.end(),
                                      nextVoxels.begin(), nextVoxels.end(),
                                      std::back_inserter(intersection));

                //double overlap = static_cast<double>(intersection.size()) / static_cast<double>(currVoxels.size() + nextVoxels.size() - intersection.size());
                //if(overlap > overlapThreshold) {
                if(!intersection.empty()) {
                    feature.overlaps_.emplace_back(Feature::Overlap{nextFeature.second.originalId_, intersection.size()});
                    //nextFeature.second.overlaps_.push_back(currFeature.second.id_);
                }
            }

            if(!feature.overlaps_.empty()) {

                // Largest keeps id.
                std::sort(feature.overlaps_.begin(), feature.overlaps_.end(),
                          [&nextFeatures] (const Feature::Overlap& a, const Feature::Overlap& b) {

                    //return nextFeatures.at(a.originalId_).voxels_.size() > nextFeatures.at(b.originalId_).voxels_.size();
                    return a.numVoxels_ > b.numVoxels_; // Another possibility...
                });

                //std::cout << "Front: " << currFeature.second.overlaps_.front().numVoxels_ << std::endl;
                //std::cout << "Back: " << currFeature.second.overlaps_.back().numVoxels_ << std::endl;

                // Remap id of next feature only.
                T largestRegionId = feature.overlaps_.front().originalId_;
                if(!nextFeatures.at(largestRegionId).wasRemapped_) {
                    nextFeatures.at(largestRegionId).remappedId_ = currFeature.remappedId_;
                    nextFeatures.at(largestRegionId).wasRemapped_ = true;
                }

            }
        }

        // Fix ids.
        auto correctedIdVolume = nextLock->clone();
        for(size_t j=0; j<correctedIdVolume->getNumVoxels(); j++) {
            T oldId = getId(correctedIdVolume, j);
            if(oldId == BACKGROUND) {
                continue;
            }

            T correctedId = nextFeatures.at(oldId).remappedId_;
            setId(correctedIdVolume, j, correctedId);
        }

        // Save.
        std::unique_ptr<Volume> vol(new Volume(correctedIdVolume, next));
        vol->setOrigin(next->getOrigin());
        auto url = vol->getOrigin().getPath();
        auto path = outputpath_.get() + '/' + tgt::FileSystem::fileName(url);
        LINFO("Writing " << path);
        VolumeSerializerPopulator().getVolumeSerializer()->write(path, vol.get());

        // Finalize.
        std::swap(curr, next);
        std::swap(currLock, nextLock);
        std::swap(currFeatures, nextFeatures);

        progressBar->setProgress(1.0f * i / list->size());
    }
}


} // namespace
