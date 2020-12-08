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

#ifndef VRN_VORTEXCRITERION_H
#define VRN_VORTEXCRITERION_H

#include "voreen/core/processors/asynccomputeprocessor.h"
#include "voreen/core/ports/geometryport.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"

namespace voreen {

struct VortexCriterionIO {
    PortDataPointer<VolumeBase> inputJacobian;
    std::unique_ptr<VolumeRAM_Float> outputCriterion;
    std::function<float(const tgt::mat3&)> criterion;
};

class VortexCriterion : public AsyncComputeProcessor<VortexCriterionIO, VortexCriterionIO> {
public:
    VortexCriterion();

    virtual Processor* create() const;
    virtual std::string getClassName() const      { return "VortexCriterion";       }
    virtual std::string getCategory() const       { return "Volume Processing";     }
    virtual CodeState getCodeState() const        { return CODE_STATE_EXPERIMENTAL; }

protected:

    virtual void setDescriptions() {
        setDescription("Calculates a voxel-wise vortex criterion");
        inputJacobian_.setDescription("Jacobi matrix of a 3D vector-valued volume");
        criterion_.setDescription("Choose between on of the currently implemented vortex criteria.");
    }

    virtual ComputeInput prepareComputeInput();
    virtual ComputeOutput compute(ComputeInput input, ProgressReporter& progressReporter) const;
    virtual void processComputeOutput(ComputeOutput output);

private:

    // Ports
    VolumePort inputJacobian_;
    VolumePort outputCriterion_;

    StringOptionProperty criterion_;

    static const std::string loggerCat_;
};

} // namespace voreen

#endif