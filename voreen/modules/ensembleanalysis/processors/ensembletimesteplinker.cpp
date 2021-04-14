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

#include "ensembletimesteplinker.h"

#include "voreen/core/datastructures/callback/memberfunctioncallback.h"

#include "modules/ensembleanalysis/utils/utils.h"

namespace voreen {

EnsembleTimeStepLinker::EnsembleTimeStepLinker()
    : Processor()
    , ensemblePort_(Port::INPORT, "ensemble.inport", "Ensemble Input")
    , inTimeStep_("inTimeStep", "Time Step (Incoming Link)", 0, 1, 0, Processor::VALID, IntProperty::DYNAMIC)
    , outFromCommonTimeRange_("outFromCommonTimeRange", "To common Time Range (Outgoing Link)", tgt::vec2(0.0f, 1.0f), 0.0f, 1.0f, 0.0f, std::numeric_limits<float>::max(), Processor::VALID)
    , outFromGlobalTimeRange_("outFromGlobalTimeRange", "To global Time Range (Outgoing Link)", tgt::vec2(0.0f, 1.0f), 0.0f, 1.0f, 0.0f, std::numeric_limits<float>::max(), Processor::VALID)
{
    addPort(ensemblePort_);
    ON_CHANGE_LAMBDA(ensemblePort_, [this] {
        const EnsembleDataset* ensemble = ensemblePort_.getData();
        if(!ensemble) {
            return;
        }
        outFromCommonTimeRange_.setMinValue(ensemble->getStartTime());
        outFromCommonTimeRange_.setMaxValue(ensemble->getEndTime());
        outFromGlobalTimeRange_.setMinValue(ensemble->getStartTime());
        outFromGlobalTimeRange_.setMaxValue(ensemble->getEndTime());
    });

    addProperty(inTimeStep_);
    inTimeStep_.setVisibleFlag(false);
    ON_CHANGE(inTimeStep_, EnsembleTimeStepLinker, onTimeStepChange);
    addProperty(outFromCommonTimeRange_);
    outFromCommonTimeRange_.setVisibleFlag(false);
    addProperty(outFromGlobalTimeRange_);
    outFromGlobalTimeRange_.setVisibleFlag(false);
}

Processor* EnsembleTimeStepLinker::create() const {
    return new EnsembleTimeStepLinker();
}

void EnsembleTimeStepLinker::onTimeStepChange() {
    const EnsembleDataset* ensemble = ensemblePort_.getData();
    if(!ensemble) {
        return;
    }

    tgt::vec2 commonTimeRange = ensemble->getCommonTimeInterval();
    float t = mapRange(inTimeStep_.get(), inTimeStep_.getMinValue(), inTimeStep_.getMaxValue(), commonTimeRange.x, commonTimeRange.y);
    outFromCommonTimeRange_.set(tgt::vec2(t, t+ensemble->getMaxTimeStepDuration()));

    tgt::vec2 globalTimeRange = tgt::vec2(ensemble->getStartTime(), ensemble->getEndTime());
    t = mapRange(inTimeStep_.get(), inTimeStep_.getMinValue(), inTimeStep_.getMaxValue(), globalTimeRange.x, globalTimeRange.y);
    outFromGlobalTimeRange_.set(tgt::vec2(t, t+ensemble->getMaxTimeStepDuration()));
}

} // namespace
