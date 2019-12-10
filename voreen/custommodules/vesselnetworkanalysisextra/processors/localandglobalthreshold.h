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

#ifndef VRN_LOCAL_AND_GLOBAL_THRESHOLD_H
#define VRN_LOCAL_AND_GLOBAL_THRESHOLD_H

#include "voreen/core/processors/asynccomputeprocessor.h"

#include "modules/base/processors/geometry/geometryprocessor.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/properties/vectorproperty.h"
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/stringproperty.h"
#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/buttonproperty.h"
#include "voreen/core/properties/numeric/intervalproperty.h"

#include "voreen/core/ports/volumeport.h"

namespace voreen {

struct LAGTInput {
    const VolumeBase& volume;
    float backgroundUpperBoundNormalized;
    float foregroundLowerBoundNormalized;
    uint8_t windowExtent;

    LAGTInput(const VolumeBase& volume, float backgroundUpperBoundNormalized, float foregroundLowerBoundNormalized, uint8_t windowExtent)
        : volume(volume)
        , backgroundUpperBoundNormalized(backgroundUpperBoundNormalized)
        , foregroundLowerBoundNormalized(foregroundLowerBoundNormalized)
        , windowExtent(windowExtent)
    {
    }
    LAGTInput(LAGTInput&& old)
        : volume(old.volume)
        , backgroundUpperBoundNormalized(old.backgroundUpperBoundNormalized)
        , foregroundLowerBoundNormalized(old.foregroundLowerBoundNormalized)
        , windowExtent(old.windowExtent)
    {
    }
};

struct LAGTOutput {
    std::unique_ptr<Volume> volume;
};


class VRN_CORE_API LocalAndGlobalThreshold : public AsyncComputeProcessor<LAGTInput, LAGTOutput>  {
public:
    LocalAndGlobalThreshold();

    virtual ~LocalAndGlobalThreshold();

    virtual std::string getCategory() const { return "Volume Processing"; }
    virtual std::string getClassName() const { return "LocalAndGlobalThreshold"; }
    virtual Processor* create() const { return new LocalAndGlobalThreshold(); }

    virtual void adjustPropertiesToInput();

    virtual LAGTInput prepareComputeInput();
    virtual LAGTOutput compute(LAGTInput input, ProgressReporter& progressReporter) const;
    virtual void processComputeOutput(LAGTOutput output);

protected:
    virtual void setDescriptions() {
        setDescription("Volume segmentation by region growing.");
    }

    FloatIntervalProperty localApplicationRange_;     ///< intensity thresholds which are applied when thresholdFilling_ is true
    IntProperty windowExtent_;
    StringProperty windowSize_;

    VolumePort inport_;
    VolumePort outport_;

    static const std::string loggerCat_; ///< category used in logging
};

} // namespace voreen

#endif // VRN_LOCAL_AND_GLOBAL_THRESHOLD_H
