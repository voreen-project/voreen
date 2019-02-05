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

#ifndef VRN_FLOWPARAMETRIZATION_H
#define VRN_FLOWPARAMETRIZATION_H

#include "voreen/core/processors/processor.h"

#include "voreen/core/properties/stringproperty.h"
#include "voreen/core/properties/numeric/intervalproperty.h"
#include "voreen/core/properties/string/stringtableproperty.h"

#include "modules/flowreen/ports/flowparametrizationport.h"

namespace voreen {

/**
 * This processor is being used to parametrize a simulation.
 */
class VRN_CORE_API FlowParametrization : public Processor {
public:
    FlowParametrization();
    virtual Processor* create() const         { return new FlowParametrization(); }

    virtual std::string getClassName() const  { return "FlowParametrization";     }
    virtual std::string getCategory() const   { return "Simulation";              }
    virtual CodeState getCodeState() const    { return CODE_STATE_EXPERIMENTAL;  }

    virtual void process();

protected:
    virtual void setDescriptions() {
        setDescription("This processor is being used to parameterize a simulation.");
    }

private:

    void addParametrization();
    void removeParametrization();
    void clearParametrizations();

    std::vector<FlowParameters> flowParameters_;

    FlowParametrizationPort outport_;

    StringProperty parametrizationName_;
    FloatProperty simulationTime_;
    FloatProperty temporalResolution_;
    FloatProperty characteristicLength_;
    FloatProperty viscosity_;
    FloatProperty density_;
    BoolProperty  bouzidi_;

    ButtonProperty addParametrization_;
    ButtonProperty removeParametrization_;
    ButtonProperty clearParametrizations_;
/*
    FloatIntervalProperty simulationTimeInterval_;
    FloatIntervalProperty temporalResolutionInterval_;
    FloatIntervalProperty characteristicLengthInterval_;
    FloatIntervalProperty viscosityInterval_;
    FloatIntervalProperty densityInterval_;
    BoolProperty bouzidi_;
*/
    StringProperty ensembleName_;
    StringTableProperty parametrizations_;

    static const std::string loggerCat_;
};

}   //namespace

#endif
