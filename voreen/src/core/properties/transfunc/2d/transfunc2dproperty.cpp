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

#include "voreen/core/properties/transfunc/2d/transfunc2dproperty.h"

namespace voreen {

TransFunc2DProperty::TransFunc2DProperty(const std::string& ident, const std::string& guiText, int invalidationLevel, Property::LevelOfDetail lod)
    : TransFuncPropertyBase(ident, guiText, invalidationLevel, lod)
    , transFunc2D_(0)
{}

TransFunc2DProperty::TransFunc2DProperty()
    : TransFuncPropertyBase("", "", Processor::INVALID_RESULT, Property::LOD_DEFAULT)
    , transFunc2D_(0)
{}

TransFunc2DProperty::~TransFunc2DProperty() {
    //delete done in base class
    transFunc2D_ = 0;
}

void TransFunc2DProperty::deinitialize() {
    TransFuncPropertyBase::deinitialize();

    //delete done in base class
    transFunc2D_ = 0;
}

void TransFunc2DProperty::set2D(TransFunc2D* tf) {
    // tf object already assigned
    if (tf == transFunc2D_) {
        return;
    }

    // assign new object, but store previous one for deletion
    TransFunc2D* oldValue = transFunc2D_;
    transFunc2D_ = tf;

    //notify changes
    TransFuncPropertyBase::set(tf);

    //delete old function (deletes base function too)
    delete oldValue;
}

TransFunc2D* TransFunc2DProperty::get() const {
    return transFunc2D_;
}

} // namespace
