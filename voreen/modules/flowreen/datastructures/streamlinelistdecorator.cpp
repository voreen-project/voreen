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

#include "streamlinelistdecorator.h"

#include "streamlinelist.h"

#include "tgt/assert.h"

namespace voreen {

StreamlineListDecoratorIdentity::StreamlineListDecoratorIdentity(StreamlineListBase* base)
    : StreamlineListBase()
    , basePointer_(base)
{
    tgtAssert(base, "Decorator without base pointer!!!");
    base->addObserver(this);
}

StreamlineListDecoratorIdentity::~StreamlineListDecoratorIdentity() {
}

StreamlineListBase* StreamlineListDecoratorIdentity::clone() const {
    if(basePointer_) {
        return basePointer_->clone();
    } else {
        tgtAssert(false,"missing base pointer");
        return 0;
    }
}

    //------------------------
    //  Observer
    //------------------------
void StreamlineListDecoratorIdentity::beforeStreamlineListDelete(const StreamlineListBase* source) {
    tgtAssert(source == basePointer_, "Notified by other StreamlineList than decorated!");

    // not ideal, but the best we can do here: notify observers about deletion of data
    basePointer_->removeObserver(this);
    std::vector<StreamlineListObserver*> observers = getObservers();
    for (size_t i = 0; i < observers.size(); ++i)
        observers[i]->beforeStreamlineListDecoratorPointerReset(this);
    basePointer_ = 0;
}

    //------------------------
    //  Streamline Handling
    //------------------------
void StreamlineListDecoratorIdentity::addStreamline(const Streamline& line) {
    if(basePointer_) {
        basePointer_->addStreamline(line);
    } else {
        tgtAssert(false,"missing base pointer");
    }
}

void StreamlineListDecoratorIdentity::addStreamlineList(const StreamlineListBase& list) {
    if(basePointer_) {
        basePointer_->addStreamlineList(list);
    } else {
        tgtAssert(false,"missing base pointer");
    }
}

const std::vector<Streamline>& StreamlineListDecoratorIdentity::removeStreamline(size_t pos) {
    tgtAssert(basePointer_, "missing base pointer");
    return basePointer_->removeStreamline(pos);
}

const std::vector<Streamline>& StreamlineListDecoratorIdentity::getStreamlines() const {
    tgtAssert(basePointer_, "missing base pointer");
    return basePointer_->getStreamlines();
}


    //------------------------------
    //  Streamline Bundle Handling
    //------------------------------
void StreamlineListDecoratorIdentity::addStreamlineBundle(const StreamlineBundle& bundle) {
    tgtAssert(basePointer_, "missing base pointer");
    basePointer_->addStreamlineBundle(bundle);
}

const std::vector<StreamlineBundle>& StreamlineListDecoratorIdentity::removeStreamlineBundle(size_t pos) {
    tgtAssert(basePointer_, "missing base pointer");
    return basePointer_->removeStreamlineBundle(pos);
}

const std::vector<StreamlineBundle>& StreamlineListDecoratorIdentity::getStreamlineBundles() const {
    tgtAssert(basePointer_, "missing base pointer");
    return basePointer_->getStreamlineBundles();
}

void StreamlineListDecoratorIdentity::setStreamlineNoiseFlag(size_t pos) {
    tgtAssert(basePointer_, "missing base pointer");
    basePointer_->setStreamlineNoiseFlag(pos);
}

const std::vector<size_t>& StreamlineListDecoratorIdentity::getStreamlineNoise() const {
    tgtAssert(basePointer_, "missing base pointer");
    return basePointer_->getStreamlineNoise();
}

    //----------------
    //  Meta
    //----------------
const tgt::svec3& StreamlineListDecoratorIdentity::getOriginalDimensions() const {
    tgtAssert(basePointer_, "missing base pointer");
    return basePointer_->getOriginalDimensions();
}

const tgt::vec3& StreamlineListDecoratorIdentity::getOriginalSpacing() const {
    tgtAssert(basePointer_, "missing base pointer");
    return basePointer_->getOriginalSpacing();
}

const tgt::Bounds& StreamlineListDecoratorIdentity::getOriginalVoxelBounds() const {
    tgtAssert(basePointer_, "missing base pointer");
    return basePointer_->getOriginalVoxelBounds();
}

const tgt::Bounds& StreamlineListDecoratorIdentity::getOriginalWorldBounds() const {
    tgtAssert(basePointer_, "missing base pointer");
    return basePointer_->getOriginalWorldBounds();
}

const tgt::mat4& StreamlineListDecoratorIdentity::getOriginalVoxelToWorldMatrix() const {
    tgtAssert(basePointer_, "missing base pointer");
    return basePointer_->getOriginalVoxelToWorldMatrix();
}

const tgt::mat4& StreamlineListDecoratorIdentity::getOriginalWorldToVoxelMatrix() const {
    tgtAssert(basePointer_, "missing base pointer");
    return basePointer_->getOriginalWorldToVoxelMatrix();
}

float StreamlineListDecoratorIdentity::getMinMagnitude() const {
    tgtAssert(basePointer_, "missing base pointer");
    return basePointer_->getMinMagnitude();
}

float StreamlineListDecoratorIdentity::getMaxMagnitude() const {
    tgtAssert(basePointer_, "missing base pointer");
    return basePointer_->getMaxMagnitude();
}

const tgt::mat4& StreamlineListDecoratorIdentity::getListTransformMatrix() const {
    tgtAssert(basePointer_, "missing base pointer");
    return basePointer_->getListTransformMatrix();
}

const tgt::mat4& StreamlineListDecoratorIdentity::getVelocityTransformMatrix() const {
    tgtAssert(basePointer_, "missing base pointer");
    return basePointer_->getVelocityTransformMatrix();
}

const tgt::mat4 StreamlineListDecoratorIdentity::getVoxelToWorldMatrix() const {
    tgtAssert(basePointer_, "missing base pointer");
    return basePointer_->getVoxelToWorldMatrix();
}

void StreamlineListDecoratorIdentity::setTransformMatrices(const tgt::mat4& listMatrix, const tgt::mat4& velocityMatrix) {
    if(basePointer_) {
        basePointer_->setTransformMatrices(listMatrix, velocityMatrix);
    } else {
        tgtAssert(false,"missing base pointer");
    }
}

    //----------------
    //  Storage
    //----------------
std::string StreamlineListDecoratorIdentity::metaToCSVString() const {
    if(basePointer_) {
        return basePointer_->metaToCSVString();
    } else {
        tgtAssert(false,"missing base pointer");
        return "";
    }
}

void StreamlineListDecoratorIdentity::serialize(Serializer& s) const {
    if(basePointer_) {
        basePointer_->serialize(s);
    } else {
        tgtAssert(false,"missing base pointer");
    }
}

void StreamlineListDecoratorIdentity::deserialize(Deserializer& s) {
    tgtAssert(false,"Decorator can not be deserializedr");
}



//--------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------
//  StreamlineListDecoratorReplaceTransformation
//--------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------
StreamlineListDecoratorReplaceTransformation::StreamlineListDecoratorReplaceTransformation(
        StreamlineListBase* base, const tgt::mat4& listMatrix, const tgt::mat4& velocityMatrix)
    : StreamlineListDecoratorIdentity(base)
    , decoratorListTransformMatrix_(listMatrix)
    , decoratorVelocityTransformMatrix_(velocityMatrix)
{}

StreamlineListBase* StreamlineListDecoratorReplaceTransformation::clone() const {
    if(basePointer_) {
        StreamlineListBase* tmp = basePointer_->clone();
        tmp->setTransformMatrices(decoratorListTransformMatrix_, decoratorVelocityTransformMatrix_);
        return tmp;
    } else {
        tgtAssert(false,"missing base pointer");
        return 0;
    }
}

const tgt::mat4& StreamlineListDecoratorReplaceTransformation::getListTransformMatrix() const {
    return decoratorListTransformMatrix_;
}

const tgt::mat4& StreamlineListDecoratorReplaceTransformation::getVelocityTransformMatrix() const {
    return decoratorVelocityTransformMatrix_;
}

const tgt::mat4 StreamlineListDecoratorReplaceTransformation::getVoxelToWorldMatrix() const {
    return decoratorListTransformMatrix_ * getOriginalVoxelToWorldMatrix();
}

void StreamlineListDecoratorReplaceTransformation::setTransformMatrices(const tgt::mat4& listMatrix, const tgt::mat4& velocityMatrix) {
    decoratorListTransformMatrix_ = listMatrix;
    decoratorVelocityTransformMatrix_ = velocityMatrix;
}

void StreamlineListDecoratorReplaceTransformation::serialize(Serializer& s) const {
    if(basePointer_) {
        tgt::mat4 tmpTransformMatrix, tmpVelocityMatrix;
        tmpTransformMatrix = basePointer_->getListTransformMatrix();
        tmpVelocityMatrix  = basePointer_->getVelocityTransformMatrix();
        basePointer_->setTransformMatrices(decoratorListTransformMatrix_, decoratorVelocityTransformMatrix_);

        basePointer_->serialize(s);

        basePointer_->setTransformMatrices(tmpTransformMatrix, tmpVelocityMatrix);
    } else {
        tgtAssert(false,"missing base pointer");
    }
};



}   // namespace
