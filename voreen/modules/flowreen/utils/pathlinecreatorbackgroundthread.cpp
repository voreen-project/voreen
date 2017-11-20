/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2017 University of Muenster, Germany.                        *
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

#include "pathlinecreatorbackgroundthread.h"

#include "voreen/core/ports/conditions/portconditionvolumetype.h"
#include "voreen/core/datastructures/volume/volumeminmaxmagnitude.h"

#include "modules/flowreen/datastructures/streamlinebundle.h"
#include "modules/flowreen/processors/streamline/streamlinecreator.h"

namespace voreen {

PathlineCreatorBackgroundThread::PathlineCreatorBackgroundThread(PathlineCreator* processor, int seedTime,
    const VolumeList* flow, StreamlineList* output,
    int maxNumPathlines, tgt::ivec2 streamlineLengthThreshold,
    tgt::vec2 absoluteMagnitudeThreshold, PathlineCreator::FilterMode filterMode)
    : ProcessorBackgroundThread<PathlineCreator>(processor)
    , rnd(std::bind(std::uniform_real_distribution<float>(0.f, 1.f), std::mt19937(seedTime)))
    , flow_(flow), output_(output)
    , maxNumPathlines_(maxNumPathlines), streamlineLengthThreshold_(tgt::svec2(streamlineLengthThreshold))
    , absoluteMagnitudeThreshold_(absoluteMagnitudeThreshold), filterMode_(filterMode)
{
    //progress is used to get Error-Messages
    setProgress("", 0.f);

    // Try to extract first timestep since it is needed for multiple purposes.
    tgtAssert(!flow->empty(), "No timestep passed");
    firstTimeStepFlow_ = dynamic_cast<const VolumeRAM_3xFloat*>(flow->at(0)->getRepresentation<VolumeRAM>());
    tgtAssert(firstTimeStepFlow_, "RAM representation of first time step not available");
}

PathlineCreatorBackgroundThread::~PathlineCreatorBackgroundThread() {
}

void PathlineCreatorBackgroundThread::threadMain() {
    //--------------------- init seeds ---------------------------------------
    std::vector<tgt::vec3> validPositions;
    for (size_t z = 0; z < firstTimeStepFlow_->getDimensions().z; z++) {
        for (size_t y = 0; y < firstTimeStepFlow_->getDimensions().y; y++) {
            for (size_t x = 0; x < firstTimeStepFlow_->getDimensions().x; x++) {
                float len = tgt::length(firstTimeStepFlow_->voxel(x, y, z));
                if (len >= absoluteMagnitudeThreshold_.x && len <= absoluteMagnitudeThreshold_.y)
                    validPositions.push_back(tgt::vec3(x, y, z));
            }
        }
    }
    if (validPositions.empty()) {
        setProgress("No valid seeding points found. Decrease lower magnitude threshold!", 0.f);
        return;
    }
    interruptionPoint();
    processor_->setProgress(0.1f);

    std::unique_ptr<tgt::vec3[]> seedingPositions(new tgt::vec3[maxNumPathlines_]);
    for (size_t i = 0; i < maxNumPathlines_; ++i) {
        tgt::vec3 basePos = validPositions.at(static_cast<size_t>(rnd() * (validPositions.size() - 1)));
        seedingPositions[i] = tgt::clamp(basePos + tgt::vec3(rnd(), rnd(), rnd()) - tgt::vec3(0.5f), tgt::vec3::zero, tgt::vec3(firstTimeStepFlow_->getDimensions()));
    }
    interruptionPoint();
    processor_->setProgress(0.2f);

    //--------------------- calculate pathlines ----------------------------
    const size_t maxNumTries = 5; // HACK: tries per pathline

    // Setup seeding points.
    std::vector<Streamline> pathlines(maxNumPathlines_);
    for (size_t i = 0; i < pathlines.size(); i++) {
        Streamline::StreamlineElement startElement;

        for (size_t j = 0; j < maxNumTries; j++) {

            // Retrieve velocity.
            startElement.velocity_ = getVelocityAt(firstTimeStepFlow_, seedingPositions[i]);

            // in case of flow at start being zero or with its magnitude not fitting
            // into the range defined by thresholds, the random position leads to no
            // useful streamline so that another position has to be taken
            float len = tgt::length(startElement.velocity_);
            if (len >= absoluteMagnitudeThreshold_.x || len <= absoluteMagnitudeThreshold_.y)
                break;
            
            reseedPosition(seedingPositions.get(), i);
        }

        startElement.position_ = seedingPositions[i];
        pathlines[i].addElementAtEnd(startElement);
    }

    const float progressPerTimeStep = 0.8f / flow_->size();
    for (size_t t = 0; t < flow_->size(); t++) {

        // Request volume for current time step.
        const VolumeRAM_3xFloat* volume = dynamic_cast<const VolumeRAM_3xFloat*>(flow_->at(t)->getRepresentation<VolumeRAM>());
        if (!volume) {
            setProgress("RAM representation request failed. You might have run out of memory.", 0.f);
            return;
        }

        for (size_t i = pathlines.size(); i > 0; i--) {

            // Retrieved offset (shifted by 1 because of reversed loop on unsigned type.
            Streamline& pathline = pathlines[i - 1];

            if (calculateRungeKuttaStep(volume, pathline)) {
                if (pathline.getNumElements() >= streamlineLengthThreshold_.x ||
                    pathline.getNumElements() <= streamlineLengthThreshold_.y) {

                    //we have a valid streamline
                    output_->addStreamline(pathline);
                }

                pathlines.pop_back();
            }

            interruptionPoint();
        }

        // Update progress.
        processor_->setProgress(processor_->getProgress() + progressPerTimeStep);
    }

    processor_->setProgress(1.0f);
}

//---------------------------------------------------------------------------
//          Helpers
//---------------------------------------------------------------------------
void PathlineCreatorBackgroundThread::reseedPosition(tgt::vec3* seedingPositions, size_t currentPosition)
{
    tgt::vec3 randVec = tgt::vec3(rnd(), rnd(), rnd());
    tgt::vec3 dimAsVec3 = tgt::vec3(firstTimeStepFlow_->getDimensions() - tgt::svec3::one);

    // Use a "die" to determine wether a completely new random position
    // will be taken or wether an exisiting one will be used.
    // When the probability for a new position is low, the seeding positions
    // seem be prone to cluster at single location.
    if ((rnd() < 0.5) || (currentPosition <= 1)) {
        seedingPositions[currentPosition] = randVec * dimAsVec3;
        return;
    }

    // if there are already random positions which lead to a vector-field value not being
    // zero or which falls within the limits defined by thresholds, take this
    // position to generate another.
    // Use the position and add some random offset to it.
    size_t index = tgt::iround(rnd() * (currentPosition - 1));
    randVec *= rnd() * 9.f + 1.f;
    seedingPositions[currentPosition] = tgt::clamp((seedingPositions[index] + randVec), tgt::vec3::zero, dimAsVec3);
}

bool PathlineCreatorBackgroundThread::calculateRungeKuttaStep(const VolumeRAM_3xFloat* volume, Streamline& pathline) {

    const tgt::vec3 dimAsVec3 = tgt::vec3(firstTimeStepFlow_->getDimensions() - tgt::svec3::one);
    const float h = 0.5f; //stepwidth;

    tgt::vec3 r = pathline.getLastElement().position_;
    tgt::vec3 velR = pathline.getLastElement().velocity_;

    //no magnitude
    if (velR == tgt::vec3::zero)
        return true;

    tgt::vec3 k1 = tgt::normalize(velR) * h; //v != zero
    tgt::vec3 k2 = getVelocityAt(volume, r + (k1 / 2.0f));
    if (k2 != tgt::vec3::zero) k2 = tgt::normalize(k2) * h;
    tgt::vec3 k3 = getVelocityAt(volume, r + (k2 / 2.0f));
    if (k3 != tgt::vec3::zero) k3 = tgt::normalize(k3) * h;
    tgt::vec3 k4 = getVelocityAt(volume, r + k3);
    if (k4 != tgt::vec3::zero) k4 = tgt::normalize(k4) * h;
    r += ((k1 / 6.0f) + (k2 / 3.0f) + (k3 / 3.0f) + (k4 / 6.0f));

    //is new r valid?
    if ((r != tgt::clamp(r, tgt::vec3::zero, dimAsVec3)) || r == pathline.getLastElement().position_) // in case of no progress on streamline in this direction...
        return true;

    velR = getVelocityAt(volume, r);
    float magnitudeR = tgt::length(velR);
    if ((magnitudeR < absoluteMagnitudeThreshold_.x) || (magnitudeR > absoluteMagnitudeThreshold_.y))
        return true;
    
    pathline.addElementAtEnd(Streamline::StreamlineElement(r, velR));
    if (pathline.getNumElements() == streamlineLengthThreshold_.y)
        return true;

    return false;
}

tgt::vec3 PathlineCreatorBackgroundThread::getVelocityAt(const VolumeRAM_3xFloat* volume, const tgt::vec3& pos) const {
    return (filterMode_ == PathlineCreator::NEAREST ? volume->getVoxelNearest(pos, 0, true)
        : volume->getVoxelLinear(pos, 0, true));
}

}   // namespace
