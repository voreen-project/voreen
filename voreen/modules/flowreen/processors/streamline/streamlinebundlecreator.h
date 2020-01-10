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

#ifndef VRN_STREAMLINEBUNDLECREATOR_H
#define VRN_STREAMLINEBUNDLECREATOR_H

#include "voreen/core/processors/asynccomputeprocessor.h"

#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/buttonproperty.h"
#include "voreen/core/properties/progressproperty.h"
#include "voreen/core/properties/numeric/intervalproperty.h"

#include "../../ports/streamlinelistport.h"

namespace voreen {

struct StreamlineBundleCreatorInput {
    PortDataPointer<StreamlineListBase> streamlines;
    int resampleSize;
    float distanceThreshold;
    float minNumStreamlinesThreshold;
};

struct StreamlineBundleCreatorOutput {
    std::unique_ptr<StreamlineListBase> streamlineBundles;
    std::unique_ptr<StreamlineListBase> streamlineNoise;
};

/**
* This processor is used to create streamline bundles from a list of streamlines.
* It can be used with the StreamlineRenderer3D.
*/
class VRN_CORE_API StreamlineBundleCreator : public AsyncComputeProcessor<StreamlineBundleCreatorInput, StreamlineBundleCreatorOutput> {
public:
    StreamlineBundleCreator();

    virtual Processor* create() const { return new StreamlineBundleCreator(); }

    virtual std::string getCategory() const { return "Streamline Processing"; }
    virtual std::string getClassName() const { return "StreamlineBundleCreator"; }
    virtual Processor::CodeState getCodeState() const { return CODE_STATE_STABLE; }

    virtual ComputeInput prepareComputeInput();
    virtual ComputeOutput compute(ComputeInput input, ProgressReporter& progressReporter) const;
    virtual void processComputeOutput(ComputeOutput output);

protected:

    virtual bool isReady() const;
    virtual void adjustPropertiesToInput();

    virtual void setDescriptions() {
        setDescription("Bundles streamlines using the QuickBundles algorithm.");
        maxAverageDistanceThresholdProp_.setDescription("In principle, this parameter specifies an approximate diameter of the main streams/flows which shall be detected. "
                                                        "Increasing this value leads to less but bigger bundles, decreasing leads to more but smaller bundles. "
                                                        "As a general rule of thumb, for the first results, a value of the same magnitude as the diameter of the vessels "
                                                        "etc. in the actual dataset should be chosen, if existent.");
        minNumStreamlinesPerBundleProp_.setDescription("Each bundle has to contain at least this amount of streamlines "
                                                       "(in percent according to the total amount) in order to be not filtered out as noise.");
        resampleSizeProp_.setDescription("The underlying algorithm needs each streamline to be resampled to a fixed number of elements/segments. "
                                         "A higher value could improve the result on a dataset with very curvy or twisty streamlines, but slows down the process.");
    }

private:

    StreamlineListPort streamlineInport_;
    StreamlineListPort streamlineBundleOutport_;
    StreamlineListPort streamlineNoiseOutport_;

    //streamline bundle settings
    FloatProperty maxAverageDistanceThresholdProp_;         ///< distance threshold for the bundle algorithm
    FloatProperty minNumStreamlinesPerBundleProp_;          ///< bundle must contain more than this percentage of streamlines
    IntProperty resampleSizeProp_;                          ///< streamlines are resampled to this value
};

}   // namespace

#endif
