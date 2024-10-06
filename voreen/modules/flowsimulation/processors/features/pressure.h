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

#ifndef VRN_PRESSURE_H
#define VRN_PRESSURE_H

#include "voreen/core/processors/asynccomputeprocessor.h"
#include "voreen/core/ports/volumeport.h"
#include "voreen/core/ports/geometryport.h"
#include "voreen/core/datastructures/geometry/glmeshgeometry.h"

#include "../../processors/simulation/flowsimulation.h"

namespace voreen {

class VRN_CORE_API Pressure : public AsyncComputeProcessor<FlowSimulationInput, FlowSimulationOutput> {
public:

    Pressure();
    virtual ~Pressure();
    virtual Processor* create() const;

    virtual std::string getClassName() const { return "Pressure"; }
    virtual std::string getCategory() const { return "Volume Processing"; }
    virtual CodeState getCodeState() const { return CODE_STATE_EXPERIMENTAL; }

    virtual bool isReady() const;

protected:
    virtual void setDescriptions() {
        setDescription("Calculates the Pressure for a flow Input and a lumen or surface. Output unit is Pa.");
    }

    virtual ComputeInput prepareComputeInput() override;
    virtual ComputeOutput compute(ComputeInput input, ProgressReporter& progressReporter) const override;
    virtual void processComputeOutput(ComputeOutput output) override;

private:
    GeometryPort geometryDataPort_;
    VolumePort geometryVolumeDataPort_;
    VolumePort measuredDataPort_;
    mutable VolumePort pressurePort_;

    FloatProperty kinematicViscosity_;
    FloatProperty density_;

    FlowSimulation flowSimulation_;

    GeometryPort geometryDataPortInternal_;
    VolumeListPort geometryVolumeDataPortInternal_;
    VolumeListPort measuredDataPortInternal_;
    FlowSimulationConfigPort configPortInternal_;

};

}

#endif