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

#ifndef VRN_PATHLINECREATOR_H
#define VRN_PATHLINECREATOR_H

#include "voreen/core/processors/asynccomputeprocessor.h"

#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/buttonproperty.h"
#include "voreen/core/properties/progressproperty.h"
#include "voreen/core/properties/numeric/intervalproperty.h"

#include "voreen/core/ports/volumeport.h"
#include "../../ports/streamlinelistport.h"

namespace voreen {

class SpatioTemporalSampler;

struct PathlineCreatorInput {
    tgt::svec2 pathlineLengthThreshold;
    tgt::vec2 absoluteMagnitudeThreshold;
    float stopIntegrationAngleThreshold;
    float velocityUnitConversion;
    float temporalResolution;
    int temporalIntegrationSteps;
    VolumeRAM::Filter filterMode;
    PortDataPointer<VolumeList> flowVolumes;
    PortDataPointer<VolumeBase> seedMask; // Might be used later on to restrict integration.
    std::list<Streamline> pathlines;
    std::unique_ptr<StreamlineListBase> output;
};

struct PathlineCreatorOutput {
    std::unique_ptr<StreamlineListBase> pathlines;
};

/**
 * This processor is used to create pathlines from a vec3 volume.
 * It can be used with the StreamlineRenderer3D. At the moment only RAM volumes are supported.
 */
class VRN_CORE_API PathlineCreator : public AsyncComputeProcessor<PathlineCreatorInput, PathlineCreatorOutput> {
public:
    PathlineCreator();

    virtual Processor* create() const { return new PathlineCreator(); }

    virtual std::string getCategory() const { return "Pathline Processing"; }
    virtual std::string getClassName() const { return "PathlineCreator"; }
    virtual Processor::CodeState getCodeState() const { return CODE_STATE_STABLE; }

    virtual ComputeInput prepareComputeInput();
    virtual ComputeOutput compute(ComputeInput input, ProgressReporter& progressReporter) const;
    virtual void processComputeOutput(ComputeOutput output);

protected:

    virtual bool isReady() const;
    virtual void adjustPropertiesToInput();
    virtual std::vector<std::reference_wrapper<Port>> getCriticalPorts();

    virtual void setDescriptions() {
        setDescription("This processor is used to create pathlines from a sequence of vec3 volume. The resulting pathlines can be visualized or modified " \
                    "by other processors of the <i>Flowreen</i> module.");
        numSeedPoints_.setDescription("Can be used to determine the number of pathlines, which should be created. It can be used as a perfromance parameter.");
        seedTime_.setDescription("It is used as debug output to see the current generator. See the next description for more details.");
        pathlineLengthThreshold_.setDescription("Pathlines, which are to short will be discarded. Pathlines, which are to long, will be clipped "\
                                                  "to the maximum threshold.");
        absoluteMagnitudeThreshold_.setDescription("Flow data points outside the threshold intervall will not be used for pathline construction.");
        relativeMagnitudeThreshold_.setDescription("Can be used to adjust the absolut magnitude correctly.");
    }

private:

    struct IntegrationInput {
        float stepSize;
        tgt::Bounds bounds;
        tgt::vec2 absoluteMagnitudeThreshold;
        float stopIntegrationAngleThreshold;
    };

    /// Adjusts the relative threshold according to the absolute one.
    void adjustRelativeThreshold();

    /// Perform a single integration step for the specified pathline.
    bool integrationStep(Streamline& pathline, const SpatioTemporalSampler& sampler, const IntegrationInput& input) const;

    VolumeListPort volumeListInport_;
    VolumePort seedMask_;
    StreamlineListPort pathlineOutport_;

    // seeding
    IntProperty numSeedPoints_;                             ///< number of seed points
    IntProperty seedTime_;                                  ///< seed

    // pathline settings
    IntIntervalProperty pathlineLengthThreshold_;       ///< pathline length must be in this interval
    FloatIntervalProperty absoluteMagnitudeThreshold_;  ///< only magnitudes in this interval are used
    BoolProperty fitAbsoluteMagnitude_;                 ///< fit magnitude on input change?
    FloatIntervalProperty relativeMagnitudeThreshold_;  ///< debug output
    IntProperty stopIntegrationAngleThreshold_;         ///< stop integration when exceeding threshold?
    FloatProperty temporalResolution_;                  ///< (global) temporal resolution between time steps
    OptionProperty<VolumeRAM::Filter> filterMode_;      ///< filtering inside the dataset

    // debug
    FloatOptionProperty velocityUnitConversion_;
    IntProperty temporalIntegrationSteps_;

    static const std::string loggerCat_;
};

}   // namespace

#endif  // VRN_STREAMLINECREATOR_H
