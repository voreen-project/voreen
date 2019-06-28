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

#include "flowparametrization.h"

#define PARAMETER_DISCRETIZATION_BEGIN(property) \
    int discretization ## property = discretization_.get(); \
    if(property ## _.get().x == property ## _.get().y) { \
        discretization ## property = 1; \
    } \
    std::string tmp = name; \
    for(int property ## i = 0; property ## i < discretization ## property; property ## i++) { \
        float property = property ## _.get().x; \
        if(discretization ## property > 1) { \
            property += (property ## _.get().y - property ## _.get().x) \
                            * property ## i / (discretization ## property - 1); \
        } \
        std::string name = tmp + static_cast<char>('A' + property ## i);
//std::string name = tmp + std::string(#property).substr(0, 3) + "=" + std::to_string(property);

#define PARAMETER_DISCRETIZATION_END }


namespace voreen {

const std::string FlowParametrization::loggerCat_("voreen.flowsimulation.FlowParametrization");

FlowParametrization::FlowParametrization()
    : Processor()
    , inport_(Port::INPORT, "inport", "Parameter Inport")
    , outport_(Port::OUTPORT, "outport", "Parameter Inport")
    , parametrizationName_("parametrizationName", "Parametrization Name", "test_parametrization")
    , simulationTime_("simulationTime", "Simulation Time (s)", 2.0f, 0.1f, 20.0f)
    , temporalResolution_("temporalResolution", "Temporal Resolution", 0.1, 0.001, 1.0f) //TODO: define proper semantic
    , spatialResolution_("spatialResolution", "Spatial Resolution", 32, 16, 512)
    , numTimeSteps_("numTimeSteps", "Num. Output Time Steps", 50, 1, 1000)
    , outputResolution_("outputResolution", "Spatial Resolution of Output", 128, 32, 1024)
    , startPhaseFunction_("startPhaseFunction", "Start Phase Function")
    , startPhaseDuration_("startPhaseDuration", "Start Phase Duration")
    , characteristicLength_("characteristicLength", "Characteristic Length (mm)", 10.0f, 0.1f, 1000.0f)
    , characteristicVelocity_("characteristicVelocity", "Characteristic Velocity (mm/s)", 10.0f, 0.0f, 1000.0f)
    , fluid_("fluid", "Fluid")
    , viscosity_("viscosity", "Kinematic Viscosity (x10^-3 m^2/s)", 3.5, 3, 4)
    , density_("density", "Density (kg/m^3)", 1000.0f, 1000.0f, 1100.0f)
    , smagorinskyConstant_("smagorinskyConstant", "Smagorinsky Contant", 0.1f, 0.01f, 2.0f)
    , bouzidi_("bouzidi", "Bouzidi", true)
    , discretization_("discretization", "Discretization", 3, 1, 26)
    , addParametrization_("addParametrizations", "Add Parametrizations")
    , removeParametrization_("removeParametrization", "Remove Parametrization")
    , clearParametrizations_("clearParametrizations", "Clear Parametrizations")
    , ensembleName_("ensembleName", "Ensemble Name", "test_ensemble")
    , parametrizations_("parametrizations", "Parametrizations", 7, Processor::VALID)

{
    addPort(inport_);
    addPort(outport_);

    addProperty(ensembleName_);
        ensembleName_.setGroupID("ensemble");
    addProperty(simulationTime_);
        simulationTime_.setGroupID("ensemble");
    addProperty(temporalResolution_);
        temporalResolution_.adaptDecimalsToRange(3);
        temporalResolution_.setGroupID("ensemble");
    addProperty(spatialResolution_);
        spatialResolution_.setGroupID("ensemble");
    addProperty(numTimeSteps_);
        numTimeSteps_.setGroupID("ensemble");
    addProperty(outputResolution_);
        outputResolution_.setGroupID("ensemble");
    addProperty(startPhaseFunction_);
        startPhaseFunction_.addOption("none", "NONE", FlowFunction::FF_NONE); // get's selected automatically
        startPhaseFunction_.addOption("constant", "CONSTANT", FlowFunction ::FF_CONSTANT);
        startPhaseFunction_.addOption("sinus", "SINUS", FlowFunction::FF_SINUS);
        startPhaseFunction_.setGroupID("ensemble");
    addProperty(startPhaseDuration_);
        startPhaseDuration_.setGroupID("ensemble");
    setPropertyGroupGuiName("ensemble", "Ensemble");

    addProperty(parametrizationName_);
        parametrizationName_.setGroupID("parameters");
    addProperty(characteristicLength_);
        characteristicLength_.setGroupID("parameters");
    addProperty(characteristicVelocity_);
        characteristicVelocity_.setGroupID("parameters");
    addProperty(fluid_);
        ON_CHANGE(fluid_, FlowParametrization, fluidChanged);
        fluid_.addOption("water", "Water", FLUID_WATER);
        fluid_.addOption("blood", "Blood", FLUID_BLOOD);
        fluid_.setGroupID("parameters");
        fluidChanged(); // Init proper values.
    addProperty(viscosity_);
        viscosity_.setGroupID("parameters");
    addProperty(density_);
        density_.setGroupID("parameters");
    addProperty(smagorinskyConstant_);
        smagorinskyConstant_.setGroupID("parameters");
    addProperty(bouzidi_);
        bouzidi_.setGroupID("parameters");
    addProperty(discretization_);
        discretization_.setGroupID("parameters");
    setPropertyGroupGuiName("parameters", "Parameters");

    addProperty(addParametrization_);
    ON_CHANGE(addParametrization_, FlowParametrization, addParametrizations);
    addProperty(removeParametrization_);
    ON_CHANGE(removeParametrization_, FlowParametrization, removeParametrization);
    addProperty(clearParametrizations_);
    ON_CHANGE(clearParametrizations_, FlowParametrization, clearParametrizations);

    addProperty(parametrizations_);
    parametrizations_.setColumnLabel(0, "Name");
    parametrizations_.setColumnLabel(1, "Char. Len.");
    parametrizations_.setColumnLabel(2, "Char. Vel.");
    parametrizations_.setColumnLabel(3, "Viscosity");
    parametrizations_.setColumnLabel(4, "Density");
    parametrizations_.setColumnLabel(5, "Smagorinsky Const.");
    parametrizations_.setColumnLabel(6, "Bouzidi");
}

void FlowParametrization::fluidChanged() {
    switch(fluid_.getValue()) {
    case FLUID_WATER:
        viscosity_.setMinValue(0.79722f);
        viscosity_.setMaxValue(1.35f);
        viscosity_.set(1.0016f); // at room temperature
        density_.setMinValue(988.1f);
        density_.setMaxValue(1000.0f);
        density_.set(998.21f); // at room temperature
        break;
    case FLUID_BLOOD:
        viscosity_.setMinValue(3.0f);
        viscosity_.setMaxValue(4.0f);
        viscosity_.set(4.0f); // literature value
        density_.setMinValue(1043.0f);
        density_.setMaxValue(1057.0f);
        density_.set(1055.0f); // literature value
        break;
    default:
        tgtAssert(false, "Unhandled fluid");
        break;
    }
}

void FlowParametrization::addParametrizations() {

    std::string name = parametrizationName_.get();
    for(const FlowParameters& params : flowParameters_) {
        if(params.getName().find(name) != std::string::npos) {
            LERROR("Already parametrization with prefix " << name);
            return;;
        }
    }

    PARAMETER_DISCRETIZATION_BEGIN(characteristicLength)
    PARAMETER_DISCRETIZATION_BEGIN(characteristicVelocity)
    PARAMETER_DISCRETIZATION_BEGIN(viscosity)
    PARAMETER_DISCRETIZATION_BEGIN(density)
    PARAMETER_DISCRETIZATION_BEGIN(smagorinskyConstant)
    bool bouzidi = bouzidi_.get();
//   for (bool bouzidi : {true, false})
    {
        FlowParameters parameters(name);
        parameters.setCharacteristicLength(characteristicLength);
        parameters.setCharacteristicVelocity(characteristicVelocity);
        parameters.setViscosity(viscosity);
        parameters.setDensity(density);
        parameters.setSmagorinskyConstant(smagorinskyConstant);
        parameters.setBouzidi(bouzidi);
        flowParameters_.push_back(parameters);

        std::vector<std::string> row(parametrizations_.getNumColumns());
        row[0] = parameters.getName();
        row[1] = std::to_string(parameters.getCharacteristicLength());
        row[2] = std::to_string(parameters.getCharacteristicVelocity());
        row[3] = std::to_string(parameters.getViscosity());
        row[4] = std::to_string(parameters.getDensity());
        row[5] = std::to_string(parameters.getSmagorinskyConstant());
        row[6] = std::to_string(parameters.getBouzidi());
        parametrizations_.addRow(row);
    }
    PARAMETER_DISCRETIZATION_END
    PARAMETER_DISCRETIZATION_END
    PARAMETER_DISCRETIZATION_END
    PARAMETER_DISCRETIZATION_END
    PARAMETER_DISCRETIZATION_END
}

void FlowParametrization::removeParametrization() {
    if(parametrizations_.getNumRows() > 0 && parametrizations_.getSelectedRowIndex() >= 0) {
        parametrizations_.removeRow(parametrizations_.getSelectedRowIndex());
        flowParameters_.erase(flowParameters_.begin() + parametrizations_.getSelectedRowIndex());
    }
    else {
        LWARNING("No parametrization selected");
    }
}

void FlowParametrization::clearParametrizations() {
    parametrizations_.reset();
    flowParameters_.clear();
}

void FlowParametrization::adjustPropertiesToInput() {
    setPropertyGroupVisible("ensemble", !inport_.hasData());
}

bool FlowParametrization::isReady() const {
    // Ignore inport!
    return outport_.isReady();
}

void FlowParametrization::process() {

    FlowParametrizationList* flowParametrizationList = nullptr;

    if(inport_.isReady()) {
        flowParametrizationList = new FlowParametrizationList(*inport_.getData());
    }
    else {
        flowParametrizationList = new FlowParametrizationList(ensembleName_.get());
        flowParametrizationList->setSimulationTime(simulationTime_.get());
        flowParametrizationList->setTemporalResolution(temporalResolution_.get());
        flowParametrizationList->setSpatialResolution(spatialResolution_.get());
        flowParametrizationList->setNumTimeSteps(numTimeSteps_.get());
        flowParametrizationList->setOutputResolution(outputResolution_.get());
        flowParametrizationList->setStartPhaseFunction(startPhaseFunction_.getValue());
    }

    for (const FlowParameters& flowParameters : flowParameters_) {
        flowParametrizationList->addFlowParameters(flowParameters);
    }

    outport_.setData(flowParametrizationList);
}

void FlowParametrization::serialize(Serializer& s) const {
    Processor::serialize(s);
    s.serialize("flowParameters", flowParameters_);
}

void FlowParametrization::deserialize(Deserializer& s) {
    Processor::deserialize(s);
    s.deserialize("flowParameters", flowParameters_);
}

}   // namespace
