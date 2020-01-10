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

#include "streamlinelistdecorator.h"

#include "streamlinelist.h"

#include "tgt/assert.h"

namespace voreen {

StreamlineListDecoratorIdentity::StreamlineListDecoratorIdentity(StreamlineListBase* base)
    : StreamlineListBase()
    , decorated_(base)
{
    tgtAssert(base, "Decorator without base pointer");
    decorated_->Observable<StreamlineListObserver>::addObserver(this);
}

StreamlineListDecoratorIdentity::~StreamlineListDecoratorIdentity() {
}

StreamlineListBase* StreamlineListDecoratorIdentity::clone() const {
    if(decorated_) {
        return decorated_->clone();
    } else {
        tgtAssert(false,"missing base pointer");
        return nullptr;
    }
}

    //------------------------
    //  Observer
    //------------------------
void StreamlineListDecoratorIdentity::beforeStreamlineListDelete(const StreamlineListBase* source) {
    tgtAssert(source == decorated_, "Notified by other StreamlineList than decorated!");

    decorated_->Observable<StreamlineListObserver>::removeObserver(this);
    std::vector<StreamlineListObserver*> observers = Observable<StreamlineListObserver>::getObservers();
    for (size_t i = 0; i < observers.size(); ++i)
        observers[i]->beforeStreamlineListDelete(this);
    decorated_ = nullptr;
}

    //------------------------
    //  Streamline Handling
    //------------------------
void StreamlineListDecoratorIdentity::addStreamline(const Streamline& line) {
    tgtAssert(decorated_, "missing base pointer");
    decorated_->addStreamline(line);
}

void StreamlineListDecoratorIdentity::addStreamlineList(const StreamlineListBase& list) {
    tgtAssert(decorated_, "missing base pointer");
    decorated_->addStreamlineList(list);
}

void StreamlineListDecoratorIdentity::removeStreamline(size_t pos) {
    tgtAssert(decorated_, "missing base pointer");
    decorated_->removeStreamline(pos);
}

void StreamlineListDecoratorIdentity::clearStreamlines() {
    tgtAssert(decorated_, "missing base pointer");
    decorated_->clearStreamlines();
}

const std::vector<Streamline>& StreamlineListDecoratorIdentity::getStreamlines() const {
    tgtAssert(decorated_, "missing base pointer");
    return decorated_->getStreamlines();
}

    //----------------
    //  Meta
    //----------------
const tgt::svec3& StreamlineListDecoratorIdentity::getOriginalDimensions() const {
    tgtAssert(decorated_, "missing base pointer");
    return decorated_->getOriginalDimensions();
}

const tgt::vec3& StreamlineListDecoratorIdentity::getOriginalSpacing() const {
    tgtAssert(decorated_, "missing base pointer");
    return decorated_->getOriginalSpacing();
}

const tgt::Bounds& StreamlineListDecoratorIdentity::getOriginalVoxelBounds() const {
    tgtAssert(decorated_, "missing base pointer");
    return decorated_->getOriginalVoxelBounds();
}

const tgt::Bounds& StreamlineListDecoratorIdentity::getOriginalWorldBounds() const {
    tgtAssert(decorated_, "missing base pointer");
    return decorated_->getOriginalWorldBounds();
}

const tgt::mat4& StreamlineListDecoratorIdentity::getOriginalVoxelToWorldMatrix() const {
    tgtAssert(decorated_, "missing base pointer");
    return decorated_->getOriginalVoxelToWorldMatrix();
}

const tgt::mat4& StreamlineListDecoratorIdentity::getOriginalWorldToVoxelMatrix() const {
    tgtAssert(decorated_, "missing base pointer");
    return decorated_->getOriginalWorldToVoxelMatrix();
}

float StreamlineListDecoratorIdentity::getMinMagnitude() const {
    tgtAssert(decorated_, "missing base pointer");
    return decorated_->getMinMagnitude();
}

float StreamlineListDecoratorIdentity::getMaxMagnitude() const {
    tgtAssert(decorated_, "missing base pointer");
    return decorated_->getMaxMagnitude();
}

const tgt::mat4& StreamlineListDecoratorIdentity::getListTransformMatrix() const {
    tgtAssert(decorated_, "missing base pointer");
    return decorated_->getListTransformMatrix();
}

const tgt::mat4& StreamlineListDecoratorIdentity::getVelocityTransformMatrix() const {
    tgtAssert(decorated_, "missing base pointer");
    return decorated_->getVelocityTransformMatrix();
}

const tgt::mat4 StreamlineListDecoratorIdentity::getVoxelToWorldMatrix() const {
    tgtAssert(decorated_, "missing base pointer");
    return decorated_->getVoxelToWorldMatrix();
}

void StreamlineListDecoratorIdentity::setTransformMatrices(const tgt::mat4& listMatrix, const tgt::mat4& velocityMatrix) {
    tgtAssert(decorated_, "missing base pointer");
    notifyPendingDataInvalidation();
    decorated_->setTransformMatrices(listMatrix, velocityMatrix);
}

    //----------------
    //  Storage
    //----------------
std::string StreamlineListDecoratorIdentity::metaToCSVString() const {
    tgtAssert(decorated_, "missing base pointer");
    return decorated_->metaToCSVString();
}

void StreamlineListDecoratorIdentity::serialize(Serializer& s) const {
    tgtAssert(decorated_, "missing base pointer");
    decorated_->serialize(s);
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
    if(decorated_) {
        StreamlineListBase* tmp = decorated_->clone();
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
    notifyPendingDataInvalidation();
    decoratorListTransformMatrix_ = listMatrix;
    decoratorVelocityTransformMatrix_ = velocityMatrix;
}

void StreamlineListDecoratorReplaceTransformation::serialize(Serializer& s) const {
    if(decorated_) {
        tgt::mat4 tmpTransformMatrix, tmpVelocityMatrix;
        tmpTransformMatrix = decorated_->getListTransformMatrix();
        tmpVelocityMatrix  = decorated_->getVelocityTransformMatrix();
        decorated_->setTransformMatrices(decoratorListTransformMatrix_, decoratorVelocityTransformMatrix_);

        decorated_->serialize(s);

        decorated_->setTransformMatrices(tmpTransformMatrix, tmpVelocityMatrix);
    } else {
        tgtAssert(false,"missing base pointer");
    }
}

}   // namespace
