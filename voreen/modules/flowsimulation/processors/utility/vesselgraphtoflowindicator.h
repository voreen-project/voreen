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

#ifndef VRN_VESSELGRAPHTOFLOWINDICATOR_H
#define VRN_VESSELGRAPHTOFLOWINDICATOR_H

#include "voreen/core/processors/processor.h"

#include "../../ports/flowsimulationconfigport.h"
#include "modules/vesselnetworkanalysis/ports/vesselgraphport.h"

namespace voreen {

/**
 * This processor is being used to generate simple synthetic simulation geometries.
 */
class VRN_CORE_API VesselGraphToFlowIndicator : public Processor {
public:
    VesselGraphToFlowIndicator();
    virtual Processor* create() const         { return new VesselGraphToFlowIndicator();    }

    virtual std::string getClassName() const  { return "VesselGraphToFlowIndicator";        }
    virtual std::string getCategory() const   { return "utility";                           }
    virtual CodeState getCodeState() const    { return CODE_STATE_EXPERIMENTAL;             }

protected:

    virtual void setDescriptions() {
        setDescription("This processor is used to populate an existing flow simulation parametrization with a volume "
                       "made from flow indicators covering the entire vessel graph.");
    }

    virtual void process();

private:

    VesselGraphPort vesselGraphInport_;
    FlowSimulationConfigPort flowParametrizationInport_;
    FlowSimulationConfigPort flowParametrizationOutport_;

    static const std::string loggerCat_;
};

}   //namespace

#endif
