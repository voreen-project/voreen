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

#ifndef VRN_FLOWSIMULATIONRESULT_H
#define VRN_FLOWSIMULATIONRESULT_H

#include "voreen/core/processors/asynccomputeprocessor.h"

#include "voreen/core/ports/volumeport.h"
#include "../../ports/flowsimulationconfigport.h"

#include "voreen/core/properties/buttonproperty.h"
#include "voreen/core/properties/stringproperty.h"
#include "voreen/core/properties/optionproperty.h"
#include "voreen/core/properties/progressproperty.h"
#include "voreen/core/properties/temppathproperty.h"

namespace voreen {

struct FlowSimulationResultInput {
    VolumeURL url;
};

struct FlowSimulationResultOutput {
    std::unique_ptr<Volume> volume;
};

/**
 * This processor loads .pvd files as result from an OpenLB simulation.
 * This processor therefore assumes the .pvd structure as in openlb 1.6.
 */
class VRN_CORE_API FlowSimulationResult : public AsyncComputeProcessor<FlowSimulationResultInput, FlowSimulationResultOutput> {
public:
    FlowSimulationResult();
    virtual ~FlowSimulationResult();
    virtual Processor* create() const         { return new FlowSimulationResult();     }

    virtual std::string getClassName() const  { return "FlowSimulationResult";         }
    virtual std::string getCategory() const   { return "Simulation";                   }
    virtual CodeState getCodeState() const    { return CODE_STATE_EXPERIMENTAL;        }

    virtual bool isReady() const override;

    virtual ComputeInput prepareComputeInput();
    virtual ComputeOutput compute(ComputeInput input, ProgressReporter& progressReporter) const;
    virtual void processComputeOutput(ComputeOutput output);

protected:
    virtual void setDescriptions() {
        setDescription("This processor loads .pvd files as result from an OpenLB simulation. "
                       "This processor therefore assumes the .pvd structure as in openlb 1.6.");
    }

    void serialize(Serializer& s) const override;
    void deserialize(Deserializer& d) override;

private:

    void onFileChange();

    std::vector<std::string> loadPvdFile(std::string filename) const;

    VolumePort outport_;

    FileDialogProperty inputFile_;
    IntProperty timeStep_;
    StringOptionProperty fields_;
    BoolProperty selectMostRecentTimeStep_;

    std::vector<std::string> timeStepPaths_;
    std::string selectedField_;

    static const std::string loggerCat_;
};

}   //namespace

#endif
