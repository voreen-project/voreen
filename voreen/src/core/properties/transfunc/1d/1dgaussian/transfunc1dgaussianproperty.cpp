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

#include "voreen/core/properties/transfunc/1d/1dgaussian/transfunc1dgaussianproperty.h"

#include "voreen/core/datastructures/volume/volumeminmax.h"

namespace voreen {

TransFunc1DGaussianProperty::TransFunc1DGaussianProperty(const std::string& ident, const std::string& guiText, int invalidationLevel, Property::LevelOfDetail lod)
    : TransFunc1DProperty(ident, guiText, invalidationLevel, lod)
    , transFunc1DGaussian_(0)
{}

TransFunc1DGaussianProperty::TransFunc1DGaussianProperty()
    : TransFunc1DProperty("", "", Processor::INVALID_RESULT, Property::LOD_DEFAULT)
    , transFunc1DGaussian_(0)
{}

TransFunc1DGaussianProperty::~TransFunc1DGaussianProperty() {
    //delete done in base class
    transFunc1DGaussian_ = 0;
}

void TransFunc1DGaussianProperty::set1DGaussian(TransFunc1DGaussian* tf) {
    // tf object already assigned
    if (tf == transFunc1DGaussian_) {
        return;
    }

    // assign new object, but store previous one for deletion
    transFunc1DGaussian_ = tf;

    //notify changes (and deletes old tf)
    TransFunc1DProperty::set1D(tf);
}

TransFunc1DGaussian* TransFunc1DGaussianProperty::get() const {
    return transFunc1DGaussian_;
}

void TransFunc1DGaussianProperty::serialize(Serializer& s) const {
    TransFunc1DProperty::serialize(s);

    s.serialize("TransferFunction", transFunc1DGaussian_);
}

void TransFunc1DGaussianProperty::deserialize(Deserializer& s) {
    TransFunc1DProperty::deserialize(s);

    TransFunc1DGaussian* tf = 0;
    s.deserialize("TransferFunction", tf);
    set1DGaussian(tf);
}

void TransFunc1DGaussianProperty::initialize() {
    TransFunc1DProperty::initialize();

    // create initial transfer function, if it has not been created during deserialization
    if (!transFunc1DGaussian_) {
        set1DGaussian(new TransFunc1DGaussian());
        LGL_ERROR;
    }
}

void TransFunc1DGaussianProperty::deinitialize() {
    TransFunc1DProperty::deinitialize();
    //delete done in base class
    transFunc1DGaussian_ = 0;
}


void TransFunc1DGaussianProperty::applyTextureResizeFromData(int requestedSize) {
    requestedSize = std::min(requestedSize, 4096); //< limit 1D tfs to 12 bit for reducing memory consumption
    transFunc1DGaussian_->resize(requestedSize,1,1);
}

bool TransFunc1DGaussianProperty::isDomainFittingEnabled() const {
    return (domainFittingStrategy_ == FIT_DOMAIN_ALWAYS ||
           (domainFittingStrategy_ == FIT_DOMAIN_INITIAL && transFunc1DGaussian_->getDomain() == tgt::vec2(0.f, 1.f)));
}

void TransFunc1DGaussianProperty::applyDomainFromData() {
    //set gamma value to default
    if(transFunc1DGaussian_->getGammaValue() != 1.f) transFunc1DGaussian_->setGammaValue(1.f);

    tgt::vec2 oldDom = transFunc1DGaussian_->getDomain();

    //set default, if no volume is present
    if(!volume_) {
        if(oldDom != tgt::vec2(0.f,1.f)) {
            transFunc1DGaussian_->setDomain(tgt::vec2(0.f, 1.f));
            invalidate();
        }
        return;
     }

     RealWorldMapping rwm = volume_->getRealWorldMapping();
     if (volume_->hasDerivedData<VolumeMinMax>() &&
        volume_->getDerivedData<VolumeMinMax>()->getNumChannels() > channel_) { //< if volume min/max values already computed, use them
        float min = rwm.normalizedToRealWorld(volume_->getDerivedData<VolumeMinMax>()->getMinNormalized(channel_));
        float max = rwm.normalizedToRealWorld(volume_->getDerivedData<VolumeMinMax>()->getMaxNormalized(channel_));
        if(oldDom != tgt::vec2(min, max)) {
            transFunc1DGaussian_->setDomain(tgt::vec2(min, max));
            invalidate();
        }
    }
    else {
        // if min/max values not available, compute them asynchronously
        // and fit to full normalized range in the mean time
        volume_->getDerivedDataThreaded<VolumeMinMax>();
        //set tmp value
        float min = rwm.normalizedToRealWorld(0.f);
        float max = rwm.normalizedToRealWorld(1.f);
        if(oldDom != tgt::vec2(min, max)) {
             transFunc1DGaussian_->setDomain(tgt::vec2(min, max));
             invalidate();
         }
    }
}

bool TransFunc1DGaussianProperty::isTFMetaPresent() const {
    bool result = false;
    if (volume_ && volume_->hasMetaData("TransTunc1DGaussian")) {
        const TransFunc1DGaussianMetaData* tfMetaData = dynamic_cast<const TransFunc1DGaussianMetaData*>(volume_->getMetaData("TransFunc1DGaussian"));
        if (tfMetaData) {
            if(channel_ < tfMetaData->getNumChannels()) {
                if (tfMetaData->getTransferFunction(channel_)) {
                    result = true;
                }
            }
        }
    }
    return result;
}

void TransFunc1DGaussianProperty::applyTFMetaFromData() {
    if(isTFMetaPresent()) {
        //volume_ and meta data have been checked in isTFMetapresent()
        set1DGaussian(static_cast<TransFunc1DGaussian*>(
            dynamic_cast<const TransFunc1DGaussianMetaData*>(volume_->getMetaData("TransFunc1DGaussian"))->getTransferFunction(channel_)->clone()));
    }
}


} // namespace voreen
