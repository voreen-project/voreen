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

#ifndef VRN_FLOWSIMULATIONGEOMETRY_H
#define VRN_FLOWSIMULATIONGEOMETRY_H

#include "voreen/core/processors/processor.h"
#include "voreen/core/ports/geometryport.h"
#include "voreen/core/properties/matrixproperty.h"

#include "../../datastructures/flowparameters.h"
#include "../../ports/flowparametrizationport.h"

namespace voreen {

enum FlowSimulationGeometryType {
    FSGT_CYLINDER,
};

/**
 * This processor is being used to generate simple synthetic simulation geometries.
 */
class VRN_CORE_API FlowSimulationGeometry : public Processor {
public:
    FlowSimulationGeometry();
    virtual Processor* create() const         { return new FlowSimulationGeometry();    }

    virtual std::string getClassName() const  { return "FlowSimulationGeometry";        }
    virtual std::string getCategory() const   { return "Simulation";                    }
    virtual CodeState getCodeState() const    { return CODE_STATE_EXPERIMENTAL;         }

    virtual void process();

protected:

    virtual void setDescriptions() {
        setDescription("This processor is being used to generate simple synthetic simulation geometries");
    }

private:

    FlowParametrizationPort flowParametrizationInport_;
    FlowParametrizationPort flowParametrizationOutport_;
    GeometryPort geometryPort_;

    OptionProperty<FlowSimulationGeometryType> geometryType_;
    FloatProperty ratio_;
    FloatMat4Property transformation_;
    OptionProperty<FlowProfile> flowProfile_;

    static const std::string loggerCat_;
};

}   //namespace

#endif
