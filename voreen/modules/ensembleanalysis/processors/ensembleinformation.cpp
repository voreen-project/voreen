/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2024 University of Muenster, Germany,                        *
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

#include "ensembleinformation.h"

#include "voreen/core/datastructures/callback/memberfunctioncallback.h"

namespace voreen {

EnsembleInformation::EnsembleInformation()
    : Processor()
    , ensemblePort_(Port::INPORT, "ensemble.inport", "Ensemble Input")
    , fieldNames_("fieldNames", "Field")
    , minValue_("minValue", "Min Value", 0.f, -1e16f, 1e16f)
    , maxValue_("maxValue", "Max Value", 0.f, -1e16f, 1e16f)
    , minMagnitude_("minMagnitude", "Min Magnitude", 0.0f, 0.0f, 1e16f)
    , maxMagnitude_("maxMagnitude", "Max Magnitude", 0.0f, 0.0f, 1e16f)
    , numChannels_("numChannels", "Num Channels", 1, 1, 16)
{
    addPort(ensemblePort_);

    addProperty(fieldNames_);

    addProperty(minValue_);
    minValue_.setNumDecimals(4);
    minValue_.setReadOnlyFlag(true);
    addProperty(maxValue_);
    maxValue_.setNumDecimals(4);
    maxValue_.setReadOnlyFlag(true);

    addProperty(minMagnitude_);
    minMagnitude_.setNumDecimals(4);
    minMagnitude_.setReadOnlyFlag(true);
    addProperty(maxMagnitude_);
    maxMagnitude_.setNumDecimals(4);
    maxMagnitude_.setReadOnlyFlag(true);

    addProperty(numChannels_);
    numChannels_.setReadOnlyFlag(true);
}

Processor* EnsembleInformation::create() const {
    return new EnsembleInformation();
}

void EnsembleInformation::adjustPropertiesToInput() {
    fieldNames_.setOptions(std::deque<Option<std::string>>());

    if(ensemblePort_.hasData()) {
        for(const std::string& fieldName : ensemblePort_.getData()->getUniqueFieldNames()) {
            fieldNames_.addOption(fieldName, fieldName);
        }
    }

    fieldNames_.setReadOnlyFlag(fieldNames_.getOptions().empty());
}

void EnsembleInformation::process() {

    auto setPropertyValue = [] (FloatProperty& property, float value) {
        property.setMinValue(value);
        property.setMaxValue(value);
        property.set(value);
        property.adaptDecimalsToRange(4);
    };

    const EnsembleDataset* ensemble = ensemblePort_.getData();
    tgtAssert(ensemble, "ensemble null");

    if(fieldNames_.getOptions().empty()) {
        return;
    }

    std::string fieldName = fieldNames_.getValue();

    auto minMax = ensemble->getValueRange(fieldName);
    setPropertyValue(minValue_, minMax.x);
    setPropertyValue(maxValue_, minMax.y);

    auto minMaxMag = ensemble->getMagnitudeRange(fieldName);
    setPropertyValue(minMagnitude_, minMaxMag.x);
    setPropertyValue(maxMagnitude_, minMaxMag.y);

    numChannels_.set(ensemble->getNumChannels(fieldName));
}

} // namespace
