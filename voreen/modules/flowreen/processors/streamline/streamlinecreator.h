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

#ifndef VRN_STREAMLINECREATOR_H
#define VRN_STREAMLINECREATOR_H

#include "voreen/core/processors/asynccomputeprocessor.h"

#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/buttonproperty.h"
#include "voreen/core/properties/progressproperty.h"
#include "voreen/core/properties/numeric/intervalproperty.h"

#include "voreen/core/ports/volumeport.h"
#include "../../ports/streamlinelistport.h"

namespace voreen {

struct StreamlineCreatorInput {
    tgt::ivec2 streamlineLengthThreshold;
    tgt::vec2 absoluteMagnitudeThreshold;
    float stopIntegrationAngleThreshold;
    VolumeRAM::Filter filterMode;
    PortDataPointer<VolumeBase> flowVolume;
    PortDataPointer<VolumeBase> seedMask; // Might be used later on to restrict integration.
    std::vector<tgt::vec3> seedPoints;
    std::unique_ptr<StreamlineListBase> output;
};

struct StreamlineCreatorOutput {
    std::unique_ptr<StreamlineListBase> streamlines;
};

/**
 * This processor is used to create streamlines from a vec3 volume.
 * It can be used with the StreamlineRenderer3D. At the moment only RAM volumes are supported.
 *
 * @Note: It uses a background thread to handle changed parameters during calculation.
 */
class VRN_CORE_API StreamlineCreator : public AsyncComputeProcessor<StreamlineCreatorInput, StreamlineCreatorOutput> {
public:
    StreamlineCreator();

    virtual Processor* create() const { return new StreamlineCreator(); }

    virtual std::string getCategory() const { return "Streamline Processing"; }
    virtual std::string getClassName() const { return "StreamlineCreator"; }
    virtual Processor::CodeState getCodeState() const { return CODE_STATE_STABLE; }

    virtual ComputeInput prepareComputeInput();
    virtual ComputeOutput compute(ComputeInput input, ProgressReporter& progressReporter) const;
    virtual void processComputeOutput(ComputeOutput output);

protected:

    virtual bool isReady() const;
    virtual void adjustPropertiesToInput();
    virtual std::vector<std::reference_wrapper<Port>> getCriticalPorts();

    virtual void setDescriptions() {
        setDescription("This processor is used to create streamlines from a vec3 float volume. The resulting streamlines can be visualized or modified " \
                        "by other processors of the <i>Flowreen</i> module.");
        numSeedPoints_.setDescription("Can be used to determine the number of streamlines, which should be created. It can be used as a perfromance parameter.");
        seedTime_.setDescription("It is used as debug output to see the current generator. See the next description for more details.");
        streamlineLengthThresholdProp_.setDescription("Streamlines, which are to short will be discarded. Streamlines, which are to long, will be clipped "\
                                                      "to the maximum threshold.");
        absoluteMagnitudeThresholdProp_.setDescription("Flow data points outside the threshold intervall will not be used for streamline construction.");
        relativeMagnitudeThresholdProp_.setDescription("Can be used to adjust the absolut magnitude correctly.");
        /*
        detectStreamlineBundlesProp_.setDescription("If checked, at any time a streamline calculation has finished a bundle detection process is triggered. " \
                                                    "This may take another few minutes. In the meantime observing the calculated streamlines is possible.");
        maxAverageDistanceThresholdProp_.setDescription("In principle, this parameter specifies an approximate diameter of the main streams/flows which shall be detected. " \
                                                        "Increasing this value leads to less but bigger bundles, decreasing leads to more but smaller bundles. " \
                                                        "As a general rule of thumb, for the first results, a value of the same magnitude as the diameter of the vessels " \
                                                        "etc. in the actual dataset should be chosen, if existent.");
        minNumStreamlinesPerBundleProp_.setDescription("Each bundle having been detected has to contain at least this amount of streamlines " \
                                                       "(in percent according to the total amount) in order to be not filtered out as noise.");
        resampleSizeProp_.setDescription("The underlying algorithm needs each streamline to be resampled to a fixed number of elements/segments. " \
                                         "A higher value could improve the result on a dataset with very curvy or twisty streamlines, but slows down the process.");
                                     */
    }

private:

    struct IntegrationInput {
        const VolumeRAM* representation;
        RealWorldMapping rwm;
        VolumeRAM::Filter filterMode;
        tgt::ivec2 streamlineLengthThreshold;
        tgt::vec2 absoluteMagnitudeThreshold;
        float stopIntegrationAngleThreshold;
    };

    /** Adjusts the relative threshold according to the absolute one. */
    void adjustRelativeThreshold();

    Streamline computeStreamlineRungeKutta(const tgt::vec3& start, const IntegrationInput& input) const;
    tgt::vec3 getVelocity(const tgt::vec3& pos, const IntegrationInput& input) const;


    VolumePort volInport_;
    VolumePort seedMask_;
    StreamlineListPort streamlineOutport_;

    // seeding
    IntProperty numSeedPoints_;                             ///< number of seed points
    IntProperty seedTime_;                                  ///< seed

    // streamline settings
    IntIntervalProperty streamlineLengthThresholdProp_;     ///< streamline length must be in this interval
    FloatIntervalProperty absoluteMagnitudeThresholdProp_;  ///< only magnitudes in this interval are used
    BoolProperty fitAbsoluteMagnitudeProp_;                 ///< fit magnitude on input change?
    FloatIntervalProperty relativeMagnitudeThresholdProp_;  ///< debug output
    IntProperty stopIntegrationAngleThresholdProp_;         ///< stop integration when exceeding threshold?
    OptionProperty<VolumeRAM::Filter> filterModeProp_;      ///< filtering inside the dataset
/*
    //streamline bundle settings
    BoolProperty detectStreamlineBundlesProp_;              ///< determines if bundles shall be detected
    ProgressProperty bundleDetectionProgressProp_;          ///< progress bar for the bundle-detection thread
    FloatProperty maxAverageDistanceThresholdProp_;         ///< distance threshold for the bundle algorithm
    FloatProperty minNumStreamlinesPerBundleProp_;          ///< bundle must contain more than this percentage of streamlines
    IntProperty resampleSizeProp_;                          ///< streamlines are resampled to this value
*/
};

}   // namespace

#endif  // VRN_STREAMLINECREATOR_H
