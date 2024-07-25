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

#ifndef VRN_FLOWSIMULATION_H
#define VRN_FLOWSIMULATION_H

#include "voreen/core/processors/asynccomputeprocessor.h"

#include "voreen/core/ports/geometryport.h"
#include "voreen/core/ports/volumeport.h"
#include "../../ports/flowsimulationconfigport.h"

#include "modules/ensembleanalysis/ports/ensembledatasetport.h"

namespace voreen {

struct FlowSimulationInput {
    std::unique_ptr<FlowSimulationConfig> config;
    size_t selectedParametrization;
    std::string simulationResultPath;
    bool deleteOldSimulations;
    int numCuboids;
};

struct FlowSimulationOutput {
};

/**
 * This processor performs simulations using the incoming parameter set on the given geometry.
 */
class VRN_CORE_API FlowSimulation : public AsyncComputeProcessor<FlowSimulationInput, FlowSimulationOutput> {
public:
    FlowSimulation();
    virtual ~FlowSimulation();

    virtual Processor* create() const         { return new FlowSimulation();    }
    virtual std::string getClassName() const  { return "FlowSimulation";        }
    virtual std::string getCategory() const   { return "Simulation";            }
    virtual CodeState getCodeState() const    { return CODE_STATE_EXPERIMENTAL; }

    virtual bool isReady() const;

protected:
    virtual void setDescriptions() {
        setDescription("This processor performs simulations using the incoming parameter set on the given geometry.");
    }

    virtual void adjustPropertiesToInput();
    virtual void clearOutports();

    virtual ComputeInput prepareComputeInput();
    virtual ComputeOutput compute(ComputeInput input, ProgressReporter& progressReporter) const;
    virtual void processComputeOutput(ComputeOutput output);

    virtual void serialize(Serializer& s) const;
    virtual void deserialize(Deserializer& s);

private:

    std::map<float, std::string> checkAndConvertVolumeList(const VolumeList* volumes, tgt::mat4 transformation) const;
    void runSimulation(const FlowSimulationInput& input, ProgressReporter& progressReporter) const;
    void enqueueInsituResult(const std::string& filename, VolumePort& port, std::unique_ptr<Volume> volume = nullptr) const;

    GeometryPort geometryDataPort_;
    VolumeListPort geometryVolumeDataPort_;
    VolumeListPort measuredDataPort_;
    FlowSimulationConfigPort parameterPort_;

    mutable VolumePort debugMaterialsPort_;
    mutable VolumePort debugVelocityPort_;
    mutable VolumePort debugPressurePort_;

    FileDialogProperty simulationResults_;
    BoolProperty deleteOldSimulations_;

    BoolProperty simulateAllParametrizations_;
    IntProperty selectedParametrization_;

    IntProperty numCuboids_;

    static const std::string loggerCat_;
};

}   //namespace

#endif
