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

#ifndef VRN_FLOWPARAMETRIZATIONENSEMBLE_H
#define VRN_FLOWPARAMETRIZATIONENSEMBLE_H

#include "voreen/core/processors/processor.h"

#include "voreen/core/properties/stringproperty.h"
#include "voreen/core/properties/floatproperty.h"
#include "voreen/core/properties/intproperty.h"

#include "../../ports/flowsimulationconfigport.h"

#include "modules/base/properties/interactivelistproperty.h"

namespace voreen {

/**
 * This processor is being used to parametrize a simulation ensemble.
 */
class VRN_CORE_API FlowParametrizationEnsemble : public Processor {
public:
    FlowParametrizationEnsemble();
    virtual Processor* create() const         { return new FlowParametrizationEnsemble(); }

    virtual std::string getClassName() const  { return "FlowParametrizationEnsemble";     }
    virtual std::string getCategory() const   { return "Simulation";                      }
    virtual CodeState getCodeState() const    { return CODE_STATE_EXPERIMENTAL;           }

    virtual void process();

protected:

    virtual void setDescriptions() {
        setDescription("This processor is being used to parameterize a simulation ensemble.");
        outputResolution_.setDescription("This restricts the output resolution. A simulation might require a high "
                                         "resolution for accurately modelling a physical system, however the output "
                                         "of the simulation might not be resolution critical. "
                                         "A value of 0 indicates no restriction.");
    }

private:

    void addFeature(const std::string& name, int id);

    FlowSimulationConfigPort outport_;

    StringProperty ensembleName_;
    FloatProperty simulationTime_;
    IntProperty numTimeSteps_;
    IntProperty outputResolution_;
    StringOptionProperty outputFileFormat_;

    InteractiveListProperty flowFeatures_;
    std::vector<int> flowFeatureIds_;

    static const std::string loggerCat_;
};

}   //namespace

#endif
