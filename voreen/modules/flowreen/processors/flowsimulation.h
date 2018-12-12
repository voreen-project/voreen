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

namespace voreen {


struct FlowSimulationInput {

};

struct FlowSimulationOutput {

};

/**
 * This processor TODO.
 */
class VRN_CORE_API FlowSimulation : public AsyncComputeProcessor<FlowSimulationInput, FlowSimulationOutput> {
public:
    FlowSimulation();
    virtual ~FlowSimulation();
    virtual Processor* create() const         { return new FlowSimulation();    }

    virtual std::string getClassName() const  { return "FlowSimulation";        }
    virtual std::string getCategory() const   { return "Simulation";            }
    virtual CodeState getCodeState() const    { return CODE_STATE_EXPERIMENTAL; }

protected:
    virtual void setDescriptions() {
        setDescription("This processor TODO.");
    }

    virtual ComputeInput prepareComputeInput();
    virtual ComputeOutput compute(ComputeInput input, ProgressReporter& progressReporter) const;
    virtual void processComputeOutput(ComputeOutput output);

private:
    GeometryPort inport_;
    VolumeListPort outport_;

    static const std::string loggerCat_;
};

}   //namespace

#endif
