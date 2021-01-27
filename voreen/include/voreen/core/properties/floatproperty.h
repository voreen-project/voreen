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

#ifndef VRN_FLOATPROPERTY_H
#define VRN_FLOATPROPERTY_H

#include "voreen/core/properties/numericproperty.h"

namespace voreen {

#ifdef DLL_TEMPLATE_INST
template class VRN_CORE_API NumericProperty<float>;
#endif

class VRN_CORE_API FloatProperty : public NumericProperty<float> {
public:
    FloatProperty(const std::string& id, const std::string& guiText,
                float value = 0.f, float minValue = 0.f, float maxValue = 1.f,
                int invalidationLevel=Processor::INVALID_RESULT,
                NumericProperty<float>::BoundaryUpdatePolicy bup = NumericProperty<float>::STATIC,
                Property::LevelOfDetail lod = Property::LOD_DEFAULT);
    FloatProperty();

    virtual Property* create() const;

    virtual std::string getClassName() const       { return "FloatProperty"; }
    virtual std::string getTypeDescription() const { return "Float"; }

    /*
     * Adapt stepping and numDecimals to the current range (max - min),
     * so that numSignificantDecimals decimals will be adjustable.
     */
    void adaptDecimalsToRange(size_t numSignificantDecimals);

};

}   // namespace

#endif
