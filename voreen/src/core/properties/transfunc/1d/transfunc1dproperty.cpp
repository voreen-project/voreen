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

#include "voreen/core/properties/transfunc/1d/transfunc1dproperty.h"

namespace voreen {

TransFunc1DProperty::TransFunc1DProperty(const std::string& ident, const std::string& guiText, int invalidationLevel, Property::LevelOfDetail lod)
    : TransFuncPropertyBase(ident, guiText, invalidationLevel, lod)
    , transFunc1D_(0)
{}

TransFunc1DProperty::TransFunc1DProperty()
    : TransFuncPropertyBase("", "", Processor::INVALID_RESULT, Property::LOD_DEFAULT)
    , transFunc1D_(0)
{}

TransFunc1DProperty::~TransFunc1DProperty() {
    //delete done in base class
    transFunc1D_ = 0;
}

void TransFunc1DProperty::deinitialize() {
    TransFuncPropertyBase::deinitialize();

    //delete done in base class
    transFunc1D_ = 0;
}

void TransFunc1DProperty::set1D(TransFunc1D* tf) {
    // tf object already assigned
    if (tf == transFunc1D_) {
        return;
    }

    // assign new object, but store previous one for deletion
    TransFunc1D* oldValue = transFunc1D_;
    transFunc1D_ = tf;

    //notify changes
    TransFuncPropertyBase::set(tf);

    //delete old function (deletes base function too)
    delete oldValue;
}

TransFunc1D* TransFunc1DProperty::get() const {
    return transFunc1D_;
}

} // namespace

