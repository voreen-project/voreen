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

#include "flowparametrizationrun.h"

#define PARAMETER_DISCRETIZATION_BEGIN(PROPERTY, TYPE) \
    int discretization ## PROPERTY = discretization_.get(); \
    if(PROPERTY ## _.get().x == PROPERTY ## _.get().y) { \
        discretization ## PROPERTY = 1; \
    } \
    std::string tmp = name; \
    for(int PROPERTY ## i = 0; PROPERTY ## i < discretization ## PROPERTY; PROPERTY ## i++) { \
        TYPE PROPERTY = PROPERTY ## _.get().x; \
        std::string name = tmp; \
        if(discretization ## PROPERTY > 1) { \
            PROPERTY += (PROPERTY ## _.get().y - PROPERTY ## _.get().x) \
                            * PROPERTY ## i / (discretization ## PROPERTY - 1); \
            name += static_cast<char>('A' + PROPERTY ## i); \
        } \
//std::string name = tmp + std::string(#PROPERTY).substr(0, 3) + "=" + std::to_string(PROPERTY);

#define PARAMETER_DISCRETIZATION_END }


namespace voreen {

const std::string FlowParametrizationRun::loggerCat_("voreen.flowsimulation.FlowParametrizationRun");

FlowParametrizationRun::FlowParametrizationRun()
    : Processor()
    , inport_(Port::INPORT, "inport", "Parameter Inport")
    , outport_(Port::OUTPORT, "outport", "Parameter Inport")
    , parametrizationName_("parametrizationName", "Parametrization Name", "test_parametrization")
    , spatialResolution_("spatialResolution", "Spatial Resolution", 32, 8, 512)
    , relaxationTime_("relaxationTime", "Relaxation Time", 1.0f, 0.5f, 4.0f)
    , characteristicLength_("characteristicLength", "Characteristic Length (mm)", 10.0f, 0.1f, 10000.0f)
    , characteristicVelocity_("characteristicVelocity", "Characteristic Velocity (m/s)", 1.0f, 0.0f, 10.0f)
    , fluid_("fluid", "Fluid")
    , viscosity_("viscosity", "Kinematic Viscosity (x10^(-6) m^2/s)", 3.5, 3, 4)
    , density_("density", "Density (kg/m^3)", 1000.0f, 1000.0f, 1100.0f)
    , smagorinskyConstant_("smagorinskyConstant", "Smagorinsky Contant", 0.1f, 0.1f, 5.0f)
    , bouzidi_("bouzidi", "Bouzidi", true)
    , discretization_("discretization", "Discretization", 3, 1, 26)
    , addParametrization_("addParametrizations", "Add Parametrizations")
    , removeParametrization_("removeParametrization", "Remove Parametrization")
    , clearParametrizations_("clearParametrizations", "Clear Parametrizations")
    , parametrizations_("parametrizations", "Parametrizations", 11, Processor::VALID)
    , addInvalidParametrizations_("addInvalidParametrizations", "Add invalid Parametrizations", false)
{
    addPort(inport_);
    addPort(outport_);

    addProperty(parametrizationName_);
        parametrizationName_.setGroupID("parameters");
    addProperty(spatialResolution_);
        spatialResolution_.setGroupID("parameters");
    addProperty(relaxationTime_);
        relaxationTime_.setNumDecimals(3);
        relaxationTime_.setGroupID("parameters");
    addProperty(characteristicLength_);
        characteristicLength_.setGroupID("parameters");
    addProperty(characteristicVelocity_);
        characteristicVelocity_.setNumDecimals(3);
        characteristicVelocity_.setGroupID("parameters");
    addProperty(fluid_);
        ON_CHANGE(fluid_, FlowParametrizationRun, fluidChanged);
        fluid_.addOption("arbitrary", "Arbitrary", FLUID_ARBITRARY);
        fluid_.addOption("water", "Water", FLUID_WATER);
        fluid_.addOption("blood", "Blood", FLUID_BLOOD);
        fluid_.setGroupID("parameters");
        fluidChanged(); // Init proper values.
    addProperty(viscosity_);
        viscosity_.setNumDecimals(4);
        viscosity_.setGroupID("parameters");
    addProperty(density_);
        density_.setNumDecimals(0);
        density_.setGroupID("parameters");
    addProperty(smagorinskyConstant_);
        smagorinskyConstant_.setNumDecimals(3);
        smagorinskyConstant_.setGroupID("parameters");
    addProperty(bouzidi_);
        bouzidi_.setGroupID("parameters");
    addProperty(discretization_);
        discretization_.setGroupID("parameters");
    setPropertyGroupGuiName("parameters", "Parameters");

    addProperty(addParametrization_);
    ON_CHANGE(addParametrization_, FlowParametrizationRun, addParametrizations);
    addProperty(removeParametrization_);
    ON_CHANGE(removeParametrization_, FlowParametrizationRun, removeParametrization);
    addProperty(clearParametrizations_);
    ON_CHANGE(clearParametrizations_, FlowParametrizationRun, clearParametrizations);

    addProperty(parametrizations_);
    parametrizations_.setColumnLabel(0, "Valid");
    parametrizations_.setColumnLabel(1, "Name");
    parametrizations_.setColumnLabel(2, "Reynolds");
    parametrizations_.setColumnLabel(3, "Res.");
    parametrizations_.setColumnLabel(4, "Relax.");
    parametrizations_.setColumnLabel(5, "Char. Len.");
    parametrizations_.setColumnLabel(6, "Char. Vel.");
    parametrizations_.setColumnLabel(7, "viscosity");
    parametrizations_.setColumnLabel(8, "density");
    parametrizations_.setColumnLabel(9, "Smagorinsky");
    parametrizations_.setColumnLabel(10, "Bouzidi");

    addProperty(addInvalidParametrizations_);
}

void FlowParametrizationRun::fluidChanged() {
    switch(fluid_.getValue()) {
    case FLUID_ARBITRARY:
        viscosity_.setMinValue(0.1f);
        viscosity_.setMaxValue(100.0f);
        viscosity_.set(1.0f);
        density_.setMinValue(0.1f);
        density_.setMaxValue(10000.0f);
        density_.set(1000.0f);
        break;
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

void FlowParametrizationRun::addParametrizations() {

    std::string name = parametrizationName_.get();
    for(const FlowParameterSet& params : flowParameters_) {
        if(params.getName().find(name) != std::string::npos) {
            VoreenApplication::app()->showMessageBox("Warning", "Already parametrization added with the same prefix");
            LWARNING("Already parametrization with prefix " << name);
            return;
        }
    }

    PARAMETER_DISCRETIZATION_BEGIN(spatialResolution, int)
    PARAMETER_DISCRETIZATION_BEGIN(relaxationTime, float)
    PARAMETER_DISCRETIZATION_BEGIN(viscosity, float)
    PARAMETER_DISCRETIZATION_BEGIN(density, float)
    PARAMETER_DISCRETIZATION_BEGIN(smagorinskyConstant, float)
    float characteristicLength = characteristicLength_.get();
    float characteristicVelocity = characteristicVelocity_.get();
    bool bouzidi = bouzidi_.get();
//   for (bool bouzidi : {true, false})
    {
        FlowParameterSet parameters(name);
        parameters.setSpatialResolution(spatialResolution);
        parameters.setRelaxationTime(relaxationTime);
        parameters.setCharacteristicLength(characteristicLength * 0.001f); // [mm] to [m]
        parameters.setCharacteristicVelocity(characteristicVelocity);
        parameters.setViscosity(viscosity * 10e-6f); // Due to interface value range.
        parameters.setDensity(density);
        parameters.setSmagorinskyConstant(smagorinskyConstant);
        parameters.setBouzidi(bouzidi);
        flowParameters_.emplace_back(parameters);

        std::vector<std::string> row;
        row.push_back(parameters.isValid() ? "✓" : "✗");
        row.push_back(parameters.getName());
        row.push_back(std::to_string(parameters.getReynoldsNumber()));
        row.push_back(std::to_string(parameters.getSpatialResolution()));
        row.push_back(std::to_string(parameters.getRelaxationTime()));
        row.push_back(std::to_string(parameters.getCharacteristicLength()));
        row.push_back(std::to_string(parameters.getCharacteristicVelocity()));
        row.push_back(std::to_string(parameters.getViscosity()));
        row.push_back(std::to_string(parameters.getDensity()));
        row.push_back(std::to_string(parameters.getSmagorinskyConstant()));
        row.push_back(std::to_string(parameters.getBouzidi()));
        parametrizations_.addRow(row);
    }
    //PARAMETER_DISCRETIZATION_END
    PARAMETER_DISCRETIZATION_END
    PARAMETER_DISCRETIZATION_END
    PARAMETER_DISCRETIZATION_END
    PARAMETER_DISCRETIZATION_END
    PARAMETER_DISCRETIZATION_END
}

void FlowParametrizationRun::removeParametrization() {
    if(parametrizations_.getNumRows() > 0 && parametrizations_.getSelectedRowIndex() >= 0) {
        flowParameters_.erase(flowParameters_.begin() + parametrizations_.getSelectedRowIndex());
        parametrizations_.removeRow(parametrizations_.getSelectedRowIndex());
    }
    else {
        VoreenApplication::app()->showMessageBox("Warning", "No parametrization selected");
        LWARNING("No parametrization selected");
    }
}

void FlowParametrizationRun::clearParametrizations() {
    parametrizations_.reset();
    flowParameters_.clear();
}

void FlowParametrizationRun::process() {

    FlowParameterSetEnsemble* flowParametrizationList = new FlowParameterSetEnsemble(*inport_.getData());

    for (const FlowParameterSet& flowParameters : flowParameters_) {
        if(addInvalidParametrizations_.get() || flowParameters.isValid()) {
            flowParametrizationList->addFlowParameterSet(flowParameters);
        }
    }

    outport_.setData(flowParametrizationList);
}

void FlowParametrizationRun::adjustPropertiesToInput() {
    if(!inport_.hasData()) {
        return;
    }

    float minCharacteristicVelocity = 0.0f;
    for(const FlowIndicator& indicator : inport_.getData()->getFlowIndicators()) {
        minCharacteristicVelocity = std::max(minCharacteristicVelocity, indicator.velocityCurve_.getMaxVelocity());
    }
    characteristicVelocity_.setMinValue(minCharacteristicVelocity);
}

void FlowParametrizationRun::serialize(Serializer& s) const {
    Processor::serialize(s);
    s.serialize("flowParameters", flowParameters_);
}

void FlowParametrizationRun::deserialize(Deserializer& s) {
    Processor::deserialize(s);
    s.deserialize("flowParameters", flowParameters_);
}

}   // namespace
