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

#include "voreen/core/properties/transfunc/2d/2dprimitives/transfunc2dprimitivesproperty.h"

#include "voreen/core/datastructures/volume/histogram.h"

#include "tgt/gpucapabilities.h"

namespace voreen {

TransFunc2DPrimitivesProperty::TransFunc2DPrimitivesProperty(const std::string& ident, const std::string& guiText, int invalidationLevel, Property::LevelOfDetail lod)
    : TransFunc2DProperty(ident, guiText, invalidationLevel, lod)
    , transFunc2DPrimitives_(0)
{}

TransFunc2DPrimitivesProperty::TransFunc2DPrimitivesProperty()
    : TransFunc2DProperty("", "", Processor::INVALID_RESULT, Property::LOD_DEFAULT)
    , transFunc2DPrimitives_(0)
{}

TransFunc2DPrimitivesProperty::~TransFunc2DPrimitivesProperty()  {
    //delete done in base class
    transFunc2DPrimitives_ = 0;
}

void TransFunc2DPrimitivesProperty::set2DPrimitives(TransFunc2DPrimitives* tf) {
    // tf object already assigned
    if (tf == transFunc2DPrimitives_) {
        return;
    }

    // assign new object, but store previous one for deletion
    transFunc2DPrimitives_ = tf;

    //notify changes (and deletes old tf)
    TransFunc2DProperty::set2D(tf);
}

TransFunc2DPrimitives* TransFunc2DPrimitivesProperty::get() const {
    return transFunc2DPrimitives_;
}

void TransFunc2DPrimitivesProperty::serialize(Serializer& s) const {
    TransFunc2DProperty::serialize(s);

    s.serialize("TransferFunction", transFunc2DPrimitives_);
}

void TransFunc2DPrimitivesProperty::deserialize(Deserializer& s) {
    TransFunc2DProperty::deserialize(s);

    TransFunc2DPrimitives* tf = 0;
    s.deserialize("TransferFunction", tf);
    set2DPrimitives(tf);
}

void TransFunc2DPrimitivesProperty::initialize() {
    TransFunc2DProperty::initialize();

    // create initial transfer function, if it has not been created during deserialization
    if (!transFunc2DPrimitives_) {
        set2DPrimitives(new TransFunc2DPrimitives());
        LGL_ERROR;
    }
}

void TransFunc2DPrimitivesProperty::deinitialize() {
    TransFunc2DProperty::deinitialize();
    //delete done in base class
    transFunc2DPrimitives_ = 0;
}


void TransFunc2DPrimitivesProperty::applyTextureResizeFromData(int requestedSize) {
    requestedSize = std::min(requestedSize, 1024); //< limit 2D tfs to 10 bit for reducing memory consumption
    transFunc2DPrimitives_->resize(requestedSize, requestedSize,1);
}

bool TransFunc2DPrimitivesProperty::isDomainFittingEnabled() const {
    return (domainFittingStrategy_ == FIT_DOMAIN_ALWAYS ||
               (domainFittingStrategy_ == FIT_DOMAIN_INITIAL &&
                transFunc2DPrimitives_->getDomain(0) == tgt::vec2(0.f, 1.f)  &&
                transFunc2DPrimitives_->getDomain(1) == tgt::vec2(0.f, 1.f)));
}

void TransFunc2DPrimitivesProperty::applyDomainFromData() {
    //set gamma value to default
    if(transFunc2DPrimitives_->getGammaValue() != tgt::vec2(1.f)) transFunc2DPrimitives_->setGammaValue(tgt::vec2(1.f));

    tgt::vec2 oldDomX = transFunc2DPrimitives_->getDomain(0);
    tgt::vec2 oldDomY = transFunc2DPrimitives_->getDomain(1);

    //set default, if no volume is present
    if(!volume_) {
        if(oldDomX != tgt::vec2(0.f,1.f) || oldDomY != tgt::vec2(0.f,1.f)) {
            transFunc2DPrimitives_->setDomain(tgt::vec2(0.f, 1.f),0);
            transFunc2DPrimitives_->setDomain(tgt::vec2(0.f, 1.f),1);
            invalidate();
        }
        return;
     }

     RealWorldMapping rwm = volume_->getRealWorldMapping();
     if (volume_->hasDerivedData<VolumeHistogramIntensityGradient>() &&
        volume_->getDerivedData<VolumeHistogramIntensityGradient>()->getNumChannels() > channel_) {

        tgt::vec2 newDomX(volume_->getDerivedData<VolumeHistogramIntensityGradient>()->getMinValue(0,channel_),
                          volume_->getDerivedData<VolumeHistogramIntensityGradient>()->getMaxValue(0,channel_));
        tgt::vec2 newDomY(volume_->getDerivedData<VolumeHistogramIntensityGradient>()->getMinValue(1,channel_),
                          volume_->getDerivedData<VolumeHistogramIntensityGradient>()->getMaxValue(1,channel_));

        if(oldDomX != newDomX || oldDomY != newDomY) {
            transFunc2DPrimitives_->setDomain(newDomX,0);
            transFunc2DPrimitives_->setDomain(newDomY,1);
            invalidate();
        }
    }
    else {
        // if VolumeHistogramIntensityGradient values not available, compute them asynchronously
        // and fit to full normalized range in the mean time
        volume_->getDerivedDataThreaded<VolumeHistogramIntensityGradient>();
        //set tmp value
        if(oldDomX != tgt::vec2(0.f,1.f) || oldDomY != tgt::vec2(0.f,1.f)) {
            transFunc2DPrimitives_->setDomain(tgt::vec2(0.f, 1.f),0);
            transFunc2DPrimitives_->setDomain(tgt::vec2(0.f, 1.f),1);
            invalidate();
        }
    }
}

bool TransFunc2DPrimitivesProperty::isTFMetaPresent() const {
    bool result = false;
    if (volume_ && volume_->hasMetaData("TransTunc2DPrimitives")) {
        const TransFunc2DPrimitivesMetaData* tfMetaData = dynamic_cast<const TransFunc2DPrimitivesMetaData*>(volume_->getMetaData("TransFunc2DPrimitives"));
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

void TransFunc2DPrimitivesProperty::applyTFMetaFromData() {
    if(isTFMetaPresent()) {
        //volume_ and meta data have been checked in isTFMetapresent()
        set2DPrimitives(static_cast<TransFunc2DPrimitives*>(
            dynamic_cast<const TransFunc2DPrimitivesMetaData*>(volume_->getMetaData("TransFunc2DPrimitives"))->getTransferFunction(channel_)->clone()));
    }
}


} // namespace voreen
