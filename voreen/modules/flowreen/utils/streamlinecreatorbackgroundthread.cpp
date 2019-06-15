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

#include "streamlinecreatorbackgroundthread.h"

#include "voreen/core/ports/conditions/portconditionvolumetype.h"
#include "voreen/core/datastructures/volume/volumeminmaxmagnitude.h"

#include "modules/flowreen/datastructures/streamlinebundle.h"
#include "modules/flowreen/processors/streamline/streamlinecreator.h"

namespace voreen {

StreamlineCreatorBackgroundThread::StreamlineCreatorBackgroundThread(StreamlineCreator* processor, int seedTime,
    const VolumeBase* flow, StreamlineList* output,
    int maxNumStreamlines, tgt::ivec2 streamlineLengthThreshold,
    tgt::vec2 absoluteMagnitudeThreshold, StreamlineCreator::FilterMode filterMode)
    : ProcessorBackgroundThread<StreamlineCreator>(processor)
    , rnd(std::bind(std::uniform_real_distribution<float>(0.f, 1.f), std::mt19937(seedTime)))
    , flow_(flow), representation_(flow), output_(output)
    , maxNumStreamlines_(maxNumStreamlines), streamlineLengthThreshold_(tgt::svec2(streamlineLengthThreshold))
    , absoluteMagnitudeThreshold_(absoluteMagnitudeThreshold), filterMode_(filterMode)
    , seedingPositions_(0)
{
    //progress is used to get Error-Messages
    setProgress("", 0.f);
}

StreamlineCreatorBackgroundThread::~StreamlineCreatorBackgroundThread() {
    delete[] seedingPositions_;
    seedingPositions_ = 0;
}

void StreamlineCreatorBackgroundThread::handleInterruption() {
    delete[] seedingPositions_;
    seedingPositions_ = 0;
}

void StreamlineCreatorBackgroundThread::threadMain() {

    //--------------------- init seeds ---------------------------------------
    std::vector<tgt::vec3> validPositions;
    for (size_t z = 0; z < flow_->getDimensions().z; z++) {
        for (size_t y = 0; y < flow_->getDimensions().y; y++) {
            for (size_t x = 0; x < flow_->getDimensions().x; x++) {
                tgt::vec3 voxel = getVelocityAt(tgt::vec3(x, y, z));
                float len = tgt::length(voxel);
                if (len >= absoluteMagnitudeThreshold_.x && len <= absoluteMagnitudeThreshold_.y)
                    validPositions.push_back(voxel);
            }
        }
    }
    if (validPositions.empty()) {
        setProgress("No valid seeding points found. Decrease lower magnitude threshold!", 0.f);
        return;
    }
    interruptionPoint();
    processor_->setProgress(0.1f);

    seedingPositions_ = new tgt::vec3[maxNumStreamlines_];
    for (size_t i = 0; i < maxNumStreamlines_; ++i) {
        tgt::vec3 basePos = validPositions.at(static_cast<size_t>(rnd() * (validPositions.size() - 1)));
        seedingPositions_[i] = tgt::clamp(basePos + tgt::vec3(rnd(), rnd(), rnd()) - tgt::vec3(0.5f), tgt::vec3::zero, tgt::vec3(flow_->getDimensions()));
    }
    interruptionPoint();
    processor_->setProgress(0.2f);

    //--------------------- calculate streamlines ----------------------------
    size_t i = 0, numTries = 0;
    const size_t maxNumTries = maxNumStreamlines_ * 5; // HACK: tries per streamline
    size_t updateProcess = std::min<size_t>(maxNumStreamlines_ / 7 + 1, 100);

    for (; ((i < maxNumStreamlines_) && (numTries < maxNumTries)); ++numTries) {
        const tgt::vec3& start = seedingPositions_[i];
        float len = tgt::length(getVelocityAt(start));
        // in case of flow at start being zero or with its magnitude not fitting
        // into the range defined by thresholds, the random position leads to no
        // useful streamline so that another position has to be taken
        if (len < absoluteMagnitudeThreshold_.x || len > absoluteMagnitudeThreshold_.y) {
            reseedPosition(i);
            continue;
        }

        //the main streamline construction
        Streamline line = computeStreamlineRungeKutta(start);

        //check streamline length
        if (line.getNumElements() < streamlineLengthThreshold_.x ||
            line.getNumElements() > streamlineLengthThreshold_.y) {
            reseedPosition(i);
            continue;
        }

        //we have a valid streamline
        output_->addStreamline(line);

        i++; //go to next seed point
        interruptionPoint();
        if (i % updateProcess == 0) {
            processor_->setProgress(0.2f + ((float)i / (float)maxNumStreamlines_)*0.75f);
        }
    }

    //------------------------------------------------------------------------
    delete[] seedingPositions_; seedingPositions_ = 0;
}

//---------------------------------------------------------------------------
//          Helpers
//---------------------------------------------------------------------------
void StreamlineCreatorBackgroundThread::reseedPosition(const size_t currentPosition)
{
    tgt::vec3 randVec = tgt::vec3(rnd(), rnd(), rnd());
    tgt::vec3 dimAsVec3 = tgt::vec3(flow_->getDimensions() - tgt::svec3::one);

    // Use a "die" to determine wether a completely new random position
    // will be taken or wether an exisiting one will be used.
    // When the probability for a new position is low, the seeding positions
    // seem be prone to cluster at single location.
    if ((rnd() < 0.5) || (currentPosition <= 1)) {
        seedingPositions_[currentPosition] = randVec * dimAsVec3;
        return;
    }

    // if there are already random positions which lead to a vector-field value not being
    // zero or which falls within the limits defined by thresholds, take this
    // position to generate another.
    // Use the position and add some random offset to it.
    size_t index = tgt::iround(rnd() * (currentPosition - 1));
    randVec *= rnd() * 9.f + 1.f;
    seedingPositions_[currentPosition] = tgt::clamp((seedingPositions_[index] + randVec), tgt::vec3::zero, dimAsVec3);
}

Streamline StreamlineCreatorBackgroundThread::computeStreamlineRungeKutta(const tgt::vec3& start) {

    const tgt::vec3 dimAsVec3 = tgt::vec3(flow_->getDimensions() - tgt::svec3::one);
    const float h = 0.5f; //stepwidth;
    const size_t maxNumElements = streamlineLengthThreshold_.y;

    tgt::vec3 r(start), r_(start);
    tgt::vec3 k1(0.0f), k2(0.0f), k3(0.0f), k4(0.0f);
    tgt::vec3 k1_(0.0f), k2_(0.0f), k3_(0.0f), k4_(0.0f);
    tgt::vec3 velR = getVelocityAt(r), velR_ = getVelocityAt(r_);

    Streamline line;            // the streamline, which will be returned
    line.addElementAtEnd(Streamline::StreamlineElement(start, velR));

    bool lookupPos = true;  // integrate along the streamline in positive direction?
    bool lookupNeg = true;  // integrate along the streamline in negative direction?

    //we look in both directions at the same time
    while (line.getNumElements() < maxNumElements && (lookupPos || lookupNeg)) {

        //look in direction
        if (lookupPos) {
            //no magnitude
            if (velR == tgt::vec3::zero) {
                lookupPos = false;
            }
            else {
                // get next point
                k1 = tgt::normalize(velR) * h; //v != zero
                k2 = getVelocityAt(r + (k1 / 2.0f));
                if (k2 != tgt::vec3::zero) k2 = tgt::normalize(k2) * h;
                k3 = getVelocityAt(r + (k2 / 2.0f));
                if (k3 != tgt::vec3::zero) k3 = tgt::normalize(k3) * h;
                k4 = getVelocityAt(r + k3);
                if (k4 != tgt::vec3::zero) k4 = tgt::normalize(k4) * h;
                r += ((k1 / 6.0f) + (k2 / 3.0f) + (k3 / 3.0f) + (k4 / 6.0f));

                //is new r valid?
                if ((r != tgt::clamp(r, tgt::vec3::zero, dimAsVec3)) || r == line.getLastElement().position_) { // in case of no progress on streamline in this direction...
                    lookupPos = false;
                }
                else {//check length
                    velR = getVelocityAt(r);
                    float magnitudeR = tgt::length(velR);
                    if ((magnitudeR < absoluteMagnitudeThreshold_.x) ||
                        (magnitudeR > absoluteMagnitudeThreshold_.y)) {
                        lookupPos = false;
                    }
                    else {
                        line.addElementAtEnd(Streamline::StreamlineElement(r, velR));
                        if (line.getNumElements() == maxNumElements)
                            break;
                    }
                }
            }
        }

        // look previous direction
        if (lookupNeg) {
            //no magnitude
            if (velR_ == tgt::vec3::zero) {
                lookupNeg = false;
            }
            else {
                // get next point
                k1_ = tgt::normalize(velR_) * h; // velR_ != zero
                k2_ = getVelocityAt(r_ - (k1_ / 2.0f));
                if (k2_ != tgt::vec3::zero) k2_ = tgt::normalize(k2_) * h;
                k3_ = getVelocityAt(r_ - (k2_ / 2.0f));
                if (k3_ != tgt::vec3::zero) k3_ = tgt::normalize(k3_) * h;
                k4_ = getVelocityAt(r_ - k3_);
                if (k4_ != tgt::vec3::zero) k4_ = tgt::normalize(k4_) * h;
                r_ -= ((k1_ / 6.0f) + (k2_ / 3.0f) + (k3_ / 3.0f) + (k4_ / 6.0f));

                //is new r valid?
                if ((r_ != tgt::clamp(r_, tgt::vec3::zero, dimAsVec3)) || (r_ == line.getFirstElement().position_)) { // in case of no progress on streamline in this direction...
                    lookupNeg = false;
                }
                else { //check length
                    velR_ = getVelocityAt(r_);
                    float magnitudeR_ = tgt::length(velR_);
                    if ((magnitudeR_ < absoluteMagnitudeThreshold_.x) ||
                        (magnitudeR_ > absoluteMagnitudeThreshold_.y)) {
                        lookupNeg = false;
                    }
                    else {
                        line.addElementAtFront(Streamline::StreamlineElement(r_, velR_));
                        //if (line.getNumElements() == maxNumElements)
                        //    break;
                    }
                }
            }
        }
    }

    return line;
}

const tgt::vec3 StreamlineCreatorBackgroundThread::getVelocityAt(const tgt::vec3& pos) {

    RealWorldMapping rwm = flow_->getRealWorldMapping();

    tgt::vec3 voxel;
    if(filterMode_ == StreamlineCreator::NEAREST) {
        for (size_t channel = 0; channel < flow_->getNumChannels(); channel++) {
            voxel[channel] = rwm.normalizedToRealWorld(
                    representation_->getVoxelNormalized(pos, channel));
        }
    }
    else {
        for (size_t channel = 0; channel < flow_->getNumChannels(); channel++) {
            voxel[channel] = rwm.normalizedToRealWorld(
                    representation_->getVoxelNormalizedLinear(pos, channel));
        }
    }
    return voxel;
}

}   // namespace
