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

#ifndef VRN_FLOWPARAMETRIZATIONRUN_H
#define VRN_FLOWPARAMETRIZATIONRUN_H

#include "voreen/core/processors/processor.h"

#include "voreen/core/properties/stringproperty.h"
#include "voreen/core/properties/numeric/intervalproperty.h"
#include "voreen/core/properties/string/stringtableproperty.h"

#include "../../ports/flowsimulationconfigport.h"

namespace voreen {

/**
 * This processor is being used to parametrize simulation runs.
 */
class VRN_CORE_API FlowParametrizationRun : public Processor {
public:
    FlowParametrizationRun();
    virtual Processor* create() const         { return new FlowParametrizationRun(); }

    virtual std::string getClassName() const  { return "FlowParametrizationRun";     }
    virtual std::string getCategory() const   { return "Simulation";                 }
    virtual CodeState getCodeState() const    { return CODE_STATE_EXPERIMENTAL;      }

    virtual void process();
    virtual void adjustPropertiesToInput();

    virtual void serialize(Serializer& s) const;
    virtual void deserialize(Deserializer& s);

protected:

    virtual void setDescriptions() {
        setDescription("This processor is being used to parameterize simulation runs.");
    }

private:

    enum Fluid {
        FLUID_ARBITRARY,
        FLUID_WATER,
        FLUID_BLOOD,
    };

    void fluidChanged();
    void addParametrizations();
    void removeParametrization();
    void clearParametrizations();

    std::vector<Parameters> flowParameters_;

    FlowSimulationConfigPort inport_;
    FlowSimulationConfigPort outport_;

    StringProperty parametrizationName_;
    IntIntervalProperty spatialResolution_;
    FloatIntervalProperty relaxationTime_;
    FloatProperty characteristicLength_;
    FloatProperty characteristicVelocity_;
    OptionProperty<Fluid> fluid_;
    FloatIntervalProperty viscosity_;
    FloatIntervalProperty density_;
    OptionProperty<FlowTurbulenceModel> turbulenceModel_;
    FloatIntervalProperty smagorinskyConstant_;
    OptionProperty<FlowBoundaryCondition> wallBoundaryCondition_;
    BoolProperty latticePerturbation_;
    FloatIntervalProperty inletVelocityMultiplier_;
    IntProperty samples_;

    ButtonProperty addParametrization_;
    ButtonProperty removeParametrization_;
    ButtonProperty clearParametrizations_;

    StringTableProperty parametrizations_;
    BoolProperty addInvalidParametrizations_;

    static const std::string loggerCat_;
};

}   //namespace

#endif
