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

#ifndef VRN_FLOWINDICATORFLUID_H
#define VRN_FLOWINDICATORFLUID_H

#include "voreen/core/processors/processor.h"
#include "../../ports/flowsimulationconfigport.h"

namespace voreen {

/**
 * This processor is being used to select in and out flow.
 */
class VRN_CORE_API FlowIndicatorFluid : public Processor {
public:
    FlowIndicatorFluid();
    virtual Processor* create() const         { return new FlowIndicatorFluid();        }

    virtual std::string getClassName() const  { return "FlowIndicatorFluid";            }
    virtual std::string getCategory() const   { return "Simulation";                    }
    virtual CodeState getCodeState() const    { return CODE_STATE_EXPERIMENTAL;         }

protected:

    virtual void process();

    virtual void setDescriptions() {
        setDescription("This processor is used to create a flow indicator covering the entire fluid domain");
    }

private:

    FlowSimulationConfigPort configInport_;
    FlowSimulationConfigPort configOutport_;

    FloatProperty startDuration_;

    static const std::string loggerCat_;
};

}   //namespace

#endif
