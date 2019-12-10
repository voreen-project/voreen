/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2019 University of Muenster, Germany,                        *
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

#include "voreen/core/properties/transfunc/transfunctypeproperty.h"

#include "voreen/core/properties/transfunc/1d/1dkeys/transfunc1dkeysproperty.h"
#include "voreen/core/properties/transfunc/1d/1dgaussian/transfunc1dgaussianproperty.h"
#include "voreen/core/properties/transfunc/2d/2dprimitives/transfunc2dprimitivesproperty.h"

namespace voreen {

TransFuncTypeProperty::TransFuncTypeProperty(const std::string& id, const std::string& guiText, int invalidationLevel, Property::LevelOfDetail lod)
    : OptionProperty<TransFuncType>(id, guiText, invalidationLevel, false, lod)
{}

TransFuncTypeProperty::TransFuncTypeProperty()
    : OptionProperty<TransFuncType>()
{}

TransFuncTypeProperty::~TransFuncTypeProperty() {
}


void TransFuncTypeProperty::select(const std::string& key) {
    OptionProperty<TransFuncType>::select(key);
    applyTypeSettings();
}

void TransFuncTypeProperty::selectByKey(const std::string& key) {
    OptionProperty<TransFuncType>::selectByKey(key);
    applyTypeSettings();
}

void TransFuncTypeProperty::invalidate() {
    OptionProperty<TransFuncType>::invalidate();
    applyTypeSettings();
}

void TransFuncTypeProperty::applyTypeSettings() {
    if(!options_.empty()) {
        switch(getValue()) {
        case TFT_1DKEYS:
            for(std::set<TransFuncPropertyBase*>::const_iterator it = properties_.begin(); it != properties_.end(); it++)
                (*it)->setVisibleFlag(dynamic_cast<TransFunc1DKeysProperty*>(*it));
            break;
        case TFT_1DGAUSSIAN:
            for(std::set<TransFuncPropertyBase*>::const_iterator it = properties_.begin(); it != properties_.end(); it++)
                (*it)->setVisibleFlag(dynamic_cast<TransFunc1DGaussianProperty*>(*it));
            break;
        case TFT_2DPRIMITIVES:
            for(std::set<TransFuncPropertyBase*>::const_iterator it = properties_.begin(); it != properties_.end(); it++)
                (*it)->setVisibleFlag(dynamic_cast<TransFunc2DPrimitivesProperty*>(*it));
            break;
        default:
            tgtAssert(false,"Unsupported option!");
        }
    }
}

void TransFuncTypeProperty::addTransFuncProperty(TransFuncPropertyBase* transFuncProperty) {
    //check type of property
    TransFuncType currentType = TFT_UNSUPPORTED;
    if (dynamic_cast<TransFunc1DKeysProperty*>(transFuncProperty))
        currentType = TFT_1DKEYS;
    else if (dynamic_cast<TransFunc1DGaussianProperty*>(transFuncProperty))
        currentType = TFT_1DGAUSSIAN;
    else if(dynamic_cast<TransFunc2DPrimitivesProperty*>(transFuncProperty))
        currentType = TFT_2DPRIMITIVES;

    //insert and add option
    if(currentType != TFT_UNSUPPORTED) {
           if((properties_.insert(transFuncProperty)).second) {
               addTransFuncTypeOption(currentType);
           } else {
                tgtAssert(false,"TF Property already added!");
                LERRORC("TransFuncTypeProperty","TF Property already added!");
           }
    } else {
        tgtAssert(false,"Unsupported TF Property!");
        LERRORC("TransFuncTypeProperty","Unsupported TF property!");
    }
    applyTypeSettings();
}

void TransFuncTypeProperty::addTransFuncTypeOption(TransFuncType option) {
    if(option == TFT_UNSUPPORTED) {
        tgtAssert(false,"TFT_UNSUPPORTED is not a valid option!");
        LERRORC("TransFuncTypeProperty","TFT_UNSUPPORTED is not a valid option!");
        return;
    }
    std::string optionString;
    switch(option) {
    case TFT_1DKEYS:
        optionString = "1D Keys Color Map";
        break;
    case TFT_1DGAUSSIAN:
        optionString = "1D Gaussian Color Map";
        break;
    case TFT_2DPRIMITIVES:
        optionString = "2D Primitves Color Map";
        break;
    default:
        tgtAssert(false,"Unsupported option!");
    }
    if(!hasKey(optionString))
        addOption(optionString,optionString,option);
}

}   // namespace
