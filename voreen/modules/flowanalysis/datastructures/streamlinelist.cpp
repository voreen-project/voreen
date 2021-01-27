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

#include "streamlinelist.h"

#include "voreen/core/datastructures/volume/volumebase.h"

namespace voreen {

StreamlineList::StreamlineList(const VolumeBase* volume)
    : StreamlineListBase()
    , dimensions_(1)
    , spacing_(1.0f)
    , magnitudeRange_(tgt::vec2::zero)
    , temporalRange_(tgt::vec2::zero)
    , listTransformMatrix_(tgt::mat4::identity)
    , velocityTransformMatrix_(tgt::mat4::identity)
{
    if(volume) {
        dimensions_             = volume->getDimensions();
        spacing_                = volume->getSpacing();
        worldBounds_            = volume->getBoundingBox().getBoundingBox();
        voxelToWorldMatrix_     = volume->getVoxelToWorldMatrix();
        worldToVoxelMatrix_     = volume->getWorldToVoxelMatrix();
        voxelBounds_.addPoint(worldToVoxelMatrix_ * worldBounds_.getLLF());
        voxelBounds_.addPoint(worldToVoxelMatrix_ * worldBounds_.getURB());
    }
}

StreamlineList::~StreamlineList() {
}

StreamlineListBase* StreamlineList::clone() const{
    StreamlineList* result           = new StreamlineList();

    result->streamlines_             = this->streamlines_;
    result->dimensions_              = this->dimensions_;
    result->spacing_                 = this->spacing_;
    result->worldBounds_             = this->worldBounds_;
    result->voxelBounds_             = this->voxelBounds_;
    result->voxelToWorldMatrix_      = this->voxelToWorldMatrix_;
    result->worldToVoxelMatrix_      = this->worldToVoxelMatrix_;
    result->magnitudeRange_          = this->magnitudeRange_;
    result->temporalRange_           = this->temporalRange_;
    result->listTransformMatrix_     = this->listTransformMatrix_;
    result->velocityTransformMatrix_ = this->velocityTransformMatrix_;

    return result;
}

    //------------------------
    //  Streamline Handling
    //------------------------
void StreamlineList::addStreamline(const Streamline& line) {

    //ignore degenerated lines
    if(line.getNumElements() < 2) {
        return;
    }

    notifyPendingDataInvalidation();

    bool wasEmpty = streamlines_.empty();
    streamlines_.push_back(line);

    if (wasEmpty || magnitudeRange_.x > line.getMinMagnitude()) {
        magnitudeRange_.x = line.getMinMagnitude();
    }
    if (wasEmpty || magnitudeRange_.y < line.getMaxMagnitude()) {
        magnitudeRange_.y = line.getMaxMagnitude();
    }
    if (wasEmpty || temporalRange_.x > line.getTemporalRange().x) {
        temporalRange_.x = line.getTemporalRange().x;
    }
    if (wasEmpty || temporalRange_.y < line.getTemporalRange().y) {
        temporalRange_.y = line.getTemporalRange().y;
    }
}

void StreamlineList::addStreamlineList(const StreamlineListBase& list) {
    //return if list is empty
    if(list.getStreamlines().empty()) return;

    notifyPendingDataInvalidation();

    bool wasEmpty = streamlines_.empty();
    if (wasEmpty || magnitudeRange_.x > list.getMinMagnitude()) {
        magnitudeRange_.x = list.getMinMagnitude();
    }
    if (wasEmpty || magnitudeRange_.y < list.getMaxMagnitude()) {
        magnitudeRange_.y = list.getMaxMagnitude();
    }
    if (wasEmpty || temporalRange_.x > list.getTemporalRange().x) {
        temporalRange_.x = list.getTemporalRange().x;
    }
    if (wasEmpty || temporalRange_.y < list.getTemporalRange().y) {
        temporalRange_.y = list.getTemporalRange().y;
    }

    //copy streamlines
    streamlines_.insert(streamlines_.end(), list.getStreamlines().begin(), list.getStreamlines().end());
}

void StreamlineList::removeStreamline(size_t pos) {
    tgtAssert(pos < streamlines_.size(), "Index out of bounds.");
    notifyPendingDataInvalidation();
    std::swap(streamlines_[pos], streamlines_[streamlines_.size()-1]);
    streamlines_.pop_back();
}

void StreamlineList::clearStreamlines() {
    notifyPendingDataInvalidation();
    streamlines_.clear();
}

const std::vector<Streamline>& StreamlineList::getStreamlines() const {
    return streamlines_;
}

//----------------
//  Meta
//----------------
const tgt::svec3& StreamlineList::getOriginalDimensions() const {
    return dimensions_;
}

const tgt::vec3& StreamlineList::getOriginalSpacing() const {
    return spacing_;
}

const tgt::Bounds& StreamlineList::getOriginalVoxelBounds() const {
    return voxelBounds_;
}

const tgt::Bounds& StreamlineList::getOriginalWorldBounds() const {
    return worldBounds_;
}

const tgt::mat4& StreamlineList::getOriginalVoxelToWorldMatrix() const {
    return voxelToWorldMatrix_;
}

const tgt::mat4& StreamlineList::getOriginalWorldToVoxelMatrix() const {
    return worldToVoxelMatrix_;
}

float StreamlineList::getMinMagnitude() const {
    return magnitudeRange_.x;
}

float StreamlineList::getMaxMagnitude() const {
    return magnitudeRange_.y;
}

const tgt::vec2& StreamlineList::getTemporalRange() const {
    return temporalRange_;
}

const tgt::mat4& StreamlineList::getListTransformMatrix() const {
    return listTransformMatrix_;
}

const tgt::mat4& StreamlineList::getVelocityTransformMatrix() const {
    return velocityTransformMatrix_;
}

const tgt::mat4 StreamlineList::getVoxelToWorldMatrix() const {
    return listTransformMatrix_ * voxelToWorldMatrix_;
}

void StreamlineList::setTransformMatrices(const tgt::mat4& listMatrix, const tgt::mat4& velocityMatrix) {
    notifyPendingDataInvalidation();
    listTransformMatrix_ = listMatrix;
    velocityTransformMatrix_ = velocityMatrix;
}

    //----------------
    //  Storage
    //----------------
std::string StreamlineList::metaToCSVString() const {
    std::stringstream output;
    output << streamlines_.size() << std::endl;
    output << magnitudeRange_.x << ", " << magnitudeRange_.y << std::endl;
    output << temporalRange_.x << ", " << temporalRange_.y << std::endl;
    output << dimensions_.x << ", " << dimensions_.y << ", " << dimensions_.z << std::endl;
    output << spacing_.x << ", " << spacing_.y << ", " << spacing_.z << std::endl;
    output << worldBounds_.getLLF().x << ", " << worldBounds_.getLLF().y << ", " << worldBounds_.getLLF().z << std::endl;
    output << worldBounds_.getURB().x << ", " << worldBounds_.getURB().y << ", " << worldBounds_.getURB().z << std::endl;
    return output.str();
}

void StreamlineList::serialize(Serializer& s) const {
    tgt::ivec3 tmpDim = dimensions_;
    s.serialize("OriginalDimensions", tmpDim);
    s.serialize("OriginalSpacing", spacing_);
    s.serialize("OriginalWorldBounds", worldBounds_);
    s.serialize("OriginalVoxelBounds", voxelBounds_);
    s.serialize("VoxelToWorld", voxelToWorldMatrix_);
    s.serialize("WorldToVoxel", worldToVoxelMatrix_);
    s.serialize("magnitudeRange", magnitudeRange_);
    s.serialize("temporalRange", temporalRange_);
    s.serialize("listTransformMatrix", listTransformMatrix_);
    s.serialize("velocityTransformMatrix", velocityTransformMatrix_);

    s.serialize("NumStreamlines", streamlines_.size());
    for(size_t i = 0; i < streamlines_.size(); i++) {
        s.serialize("Streamline" + std::to_string(i), streamlines_[i]);
    }
}

void StreamlineList::deserialize(Deserializer& s) {
    tgt::ivec3 tmpDim = dimensions_;
    s.deserialize("OriginalDimensions",tmpDim);
    dimensions_ = tmpDim;
    s.deserialize("OriginalSpacing", spacing_);
    s.deserialize("OriginalWorldBounds",worldBounds_);
    s.deserialize("OriginalVoxelBounds",voxelBounds_);
    s.deserialize("VoxelToWorld",voxelToWorldMatrix_);
    s.deserialize("WorldToVoxel",worldToVoxelMatrix_);
    s.deserialize("magnitudeRange", magnitudeRange_);
    s.deserialize("temporalRange", temporalRange_);
    s.optionalDeserialize("listTransformMatrix", listTransformMatrix_, tgt::mat4::identity);
    s.optionalDeserialize("velocityTransformMatrix", velocityTransformMatrix_, tgt::mat4::identity);

    size_t numStreamlines = 0;
    s.deserialize("NumStreamlines", numStreamlines);
    streamlines_.clear();
    streamlines_.resize(numStreamlines);
    for(size_t i = 0; i < numStreamlines; i++) {
        s.deserialize("Streamline" + std::to_string(i), streamlines_[i]);
    }
}

}   // namespace
