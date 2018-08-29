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

#include "voreen/core/properties/floatproperty.h"

namespace voreen {

FloatProperty::FloatProperty(const std::string& id, const std::string& guiText,
                     float value, float minValue, float maxValue,
                     int invalidationLevel, NumericProperty<float>::BoundaryUpdatePolicy bup,
                     Property::LevelOfDetail lod)
    : NumericProperty<float>(id, guiText, value, minValue, maxValue, 0.01f,
                               invalidationLevel, bup, lod)
{}

FloatProperty::FloatProperty()
    : NumericProperty<float>("", "", 0.f, 0.f, 100.f, 0.01f, Processor::INVALID_RESULT)
{}

Property* FloatProperty::create() const {
    return new FloatProperty();
}

void FloatProperty::adaptDecimalsToRange(size_t numSignificantDecimals) {
    const float range = getMaxValue() - getMinValue();
    if(range <= 0) {
        return;
    }
    // We can only omit digits behind the period, digits has to be >= 0.
    int digits = std::max(0, tgt::iround(-log10(range)) + static_cast<int>(numSignificantDecimals));
    setNumDecimals(digits);
    setStepping(pow(10.0f, -digits));

    // Hack: If we set min or max to a very specific value before,
    // the GUI representation might be slightly off, because it is
    // affected by the stepping, but the actual min/max values are not.
    // So now that we have set the stepping appropriately, we force
    // an update by changing min/max around and back to their original
    // values.
    float min = getMinValue();
    float max = getMaxValue();
    setMinValue(std::numeric_limits<float>::lowest());
    setMaxValue(std::numeric_limits<float>::max());
    setMinValue(min);
    setMaxValue(max);
}


}   // namespace
