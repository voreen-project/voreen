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

#include "voreen/core/properties/transfunc/transfuncpropertybase.h"

#include "voreen/core/datastructures/transfunc/transfuncbase.h"

#include "voreen/core/voreenapplication.h"

#include "voreen/core/datastructures/volume/volumeminmax.h"
#include "tgt/gpucapabilities.h"


namespace voreen {

TransFuncPropertyBase::TransFuncPropertyBase(const std::string& ident, const std::string& guiText, int invalidationLevel, Property::LevelOfDetail lod)
    : Property(ident, guiText, invalidationLevel, lod)
    , volume_(0)
    , channel_(0)
    , domainFittingStrategy_(FIT_DOMAIN_INITIAL)
    , baseFunction_(0)
    , computeHistogram_(true)
{
}

TransFuncPropertyBase::TransFuncPropertyBase()
    : Property("", "", Processor::INVALID_RESULT, Property::LOD_DEFAULT)
    , volume_(0)
    , channel_(0)
    , domainFittingStrategy_(FIT_DOMAIN_INITIAL)
    , baseFunction_(0)
{}

TransFuncPropertyBase::~TransFuncPropertyBase() {
    baseFunction_ = 0;
}

//-------------------------------------------------
//  Property Functions
//-------------------------------------------------
void TransFuncPropertyBase::deinitialize() {
    //remove old observer
    if(volume_) volume_->Observable<VolumeObserver>::removeObserver(this);

    if (baseFunction_) {
        delete baseFunction_;
        baseFunction_ = 0;
        if(VoreenApplication::app() && VoreenApplication::app()->isInitializedGL()) {
            LGL_ERROR;
        }
    }
    Property::deinitialize();
}

void TransFuncPropertyBase::serialize(Serializer& s) const {
    Property::serialize(s);

    s.serialize("domainFittingStrategy", domainFittingStrategy_);
    s.serialize("computeHistogram", computeHistogram_);
}

void TransFuncPropertyBase::deserialize(Deserializer& s) {
    Property::deserialize(s);

    try {
        int fittingStrategy = domainFittingStrategy_;
        s.deserialize("domainFittingStrategy", fittingStrategy);
        domainFittingStrategy_ = (DomainAutoFittingStrategy)fittingStrategy;
    }
    catch (SerializationNoSuchDataException&) {
        s.removeLastError();
    }

    s.optionalDeserialize("computeHistogram", computeHistogram_, true);
}

void TransFuncPropertyBase::set(TransFuncBase* tf) {

    // assign new object, but store previous one for deletion
    TransFuncBase* oldValue = baseFunction_;
    baseFunction_ = tf;

    if (oldValue && tf && *tf == *oldValue) {
        // tf is equal -> only inform widgets about new object
        updateWidgets();
    }
    else {
        // tf different -> notify change
        invalidate();
    }
}

void TransFuncPropertyBase::reset() {
    if(baseFunction_)
        baseFunction_->reset();
}

    //-------------------------------------------------
    //  Domain Fitting Handling
    //-------------------------------------------------
void TransFuncPropertyBase::setDomainFittingStrategy(TransFuncPropertyBase::DomainAutoFittingStrategy strategy) {
    domainFittingStrategy_ = strategy;
}

TransFuncPropertyBase::DomainAutoFittingStrategy TransFuncPropertyBase::getDomainFittingStrategy() const {
    return domainFittingStrategy_;
}

    //-------------------------------------------------
    //  Volume Handling
    //-------------------------------------------------
void TransFuncPropertyBase::setVolume(const VolumeBase* volume, size_t channel) {

    // this is necessary to fix problems with onChange callbacks on ports
    if (!isInitialized())
        return;

    //check, if anything has changed
    if (volume_ != volume || channel_ != channel) {
        //remove old observer
        if(volume_) volume_->Observable<VolumeObserver>::removeObserver(this);

        volume_ = volume;

        if (volume_) {
            if(volume_->getNumChannels() <= channel) {
                LWARNINGC("TransFuncPropertyBase","Volume has not requested channel " << itos(channel) << ". Channel 0 will be used instead!");
                channel = 0;
            }
            channel_ = channel;

            volume_->Observable<VolumeObserver>::addObserver(this);

            // Resize texture of tf according to bitdepth of volume
            size_t bits = (volume_->getBytesPerVoxel() * 8) / volume_->getNumChannels();
            if (bits > 16)
                bits = 16; // handle float data as if it was 16 bit to prevent overflow

            int max = static_cast<int>(pow(2.f, (int)bits));
            if (tgt::Singleton<tgt::GpuCapabilities>::isInited())
                max = std::min(max, GpuCaps.getMaxTextureSize());

            applyTextureResizeFromData(max);

            //apply domain fitting if needed
            if(isDomainFittingEnabled())
                applyDomainFromData();
        }

        updateWidgets();
    }
}

const VolumeBase* TransFuncPropertyBase::getVolume() const {
    return volume_;
}

const size_t TransFuncPropertyBase::getVolumeChannel() const {
    return channel_;
}

void TransFuncPropertyBase::volumeDelete(const VolumeBase* source) {
    if (volume_ == source) {
        volume_->Observable<VolumeObserver>::removeObserver(this);
        volume_ = 0;
        channel_ = 0;
        updateWidgets();
    }
}

void TransFuncPropertyBase::volumeChange(const VolumeBase* source) {
    if (volume_ == source) {
        updateWidgets();
    }
}

void TransFuncPropertyBase::derivedDataThreadFinished(const VolumeBase* source) {
    if (source->hasDerivedData<VolumeMinMax>() && volume_ == source &&
        isDomainFittingEnabled()) {
        applyDomainFromData();
    }
    //update widget (histogram callback)
    updateWidgets();
}

bool TransFuncPropertyBase::getComputeHistogram() const {
    return computeHistogram_;
}

void TransFuncPropertyBase::setComputeHistogram(bool autoHist /*= true*/) {
    computeHistogram_ = autoHist;
    updateWidgets();
}


}

