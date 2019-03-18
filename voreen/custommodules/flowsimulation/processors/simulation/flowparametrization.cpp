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

namespace voreen {

const std::string FlowParametrization::loggerCat_("voreen.flowreen.FlowParametrization");

FlowParametrization::FlowParametrization()
    : Processor()
    , inport_(Port::INPORT, "inport", "Parameter Inport")
    , outport_(Port::OUTPORT, "outport", "Parameter Inport")
    , parametrizationName_("parametrizationName", "Parametrization Name", "test_parametrization")
    , simulationTime_("simulationTime", "Simulation Time (s)", 2.0f, 0.1f, 20.0f)
    , temporalResolution_("temporalResolution", "Temporal Resolution (ms)", 3.1f, 1.0f, 30.0f)
    , spatialResolution_("spatialResolution", "Spatial Resolution", 128, 32, 1024)
    , flowFunction_("flowFunction", "Flow Function")
    , characteristicLength_("characteristicLength", "Characteristic Length (mm)", 10.0f, 0.1f, 1000.0f)
    , characteristicVelocity_("characteristicVelocity", "Characteristic Velocity (mm/s)", 10.0f, 0.0f, 1000.0f)
    , viscosity_("viscosity", "Viscosity (e-6 m^2/s)", 3.5, 3, 4)
    , density_("density", "Density (kg/m^3)", 1000.0f, 1000.0f, 1100.0f)
    , bouzidi_("bouzidi", "Bouzidi", true)
    , addParametrization_("addParametrization", "Add Parametrization")
    , removeParametrization_("removeParametrization", "Remove Parametrization")
    , clearParametrizations_("clearParametrizations", "Clear Parametrizations")
    , autoGenerateParametrizations_("autoGenerateParametrizations", "auto-generate parametrizations")
    , ensembleName_("ensembleName", "Ensemble Name", "test_ensemble")
    , parametrizations_("parametrizations", "Parametrizations", 6, Processor::VALID)

{
    addPort(inport_);
    addPort(outport_);

    addProperty(ensembleName_);
        ensembleName_.setGroupID("ensemble");
    addProperty(simulationTime_);
        simulationTime_.setGroupID("ensemble");
    addProperty(temporalResolution_);
        temporalResolution_.setGroupID("ensemble");
    addProperty(spatialResolution_);
        spatialResolution_.setGroupID("ensemble");
    addProperty(flowFunction_);
        flowFunction_.addOption("none", "NONE", FlowFunction::FF_NONE); // get's selected automatically
        flowFunction_.addOption("constant", "CONSTANT", FlowFunction ::FF_CONSTANT);
        flowFunction_.addOption("sinus", "SINUS", FlowFunction::FF_SINUS);
        flowFunction_.setGroupID("ensemble");
    setPropertyGroupGuiName("ensemble", "Ensemble");

    addProperty(parametrizationName_);
        parametrizationName_.setGroupID("parameters");
    addProperty(characteristicLength_);
        characteristicLength_.setGroupID("parameters");
    addProperty(characteristicVelocity_);
        characteristicVelocity_.setGroupID("parameters");
    addProperty(viscosity_);
        viscosity_.setGroupID("parameters");
    addProperty(density_);
        density_.setGroupID("parameters");
    addProperty(bouzidi_);
        bouzidi_.setGroupID("parameters");
    setPropertyGroupGuiName("parameters", "Parameters");

    addProperty(addParametrization_);
    ON_CHANGE(addParametrization_, FlowParametrization, addParametrization);
    addProperty(removeParametrization_);
    ON_CHANGE(removeParametrization_, FlowParametrization, removeParametrization);
    addProperty(clearParametrizations_);
    ON_CHANGE(clearParametrizations_, FlowParametrization, clearParametrizations);

    addProperty(autoGenerateParametrizations_);
    ON_CHANGE(autoGenerateParametrizations_, FlowParametrization, autoGenerateEnsemble);

    addProperty(parametrizations_);
    parametrizations_.setColumnLabel(0, "Name");
    parametrizations_.setColumnLabel(1, "Char. Len.");
    parametrizations_.setColumnLabel(2, "Char. Vel.");
    parametrizations_.setColumnLabel(3, "Viscosity");
    parametrizations_.setColumnLabel(4, "Density");
    parametrizations_.setColumnLabel(5, "Bouzidi");
}

void FlowParametrization::addParametrization() {

    // Check for duplicates (name is ID).
    for(const FlowParameters& parameters : flowParameters_) {
        if(parameters.getName() == parametrizationName_.get()) {
            LERROR("The name has already been used within the ensemble.");
            return;
        }
    }

    FlowParameters parameters(parametrizationName_.get());
    parameters.setCharacteristicLength(characteristicLength_.get());
    parameters.setCharacteristicVelocity(characteristicVelocity_.get());
    parameters.setViscosity(viscosity_.get());
    parameters.setDensity(density_.get());
    parameters.setBouzidi(bouzidi_.get());
    flowParameters_.push_back(parameters);

    std::vector <std::string> row(6);
    row[0] = parameters.getName();
    row[1] = std::to_string(parameters.getCharacteristicLength());
    row[2] = std::to_string(parameters.getCharacteristicVelocity());
    row[3] = std::to_string(parameters.getViscosity());
    row[4] = std::to_string(parameters.getDensity());
    row[5] = std::to_string(parameters.getBouzidi());
    parametrizations_.addRow(row);
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

void FlowParametrization::autoGenerateEnsemble() {
    // TODO: implement intelligent heuristic!
//    const size_t split = 3;
//    for(size_t i=0; i<split; i++) {
//        float characteristicLength = characteristicLength_.getMinValue() + (characteristicLength_.getMaxValue() - characteristicLength_.getMinValue()) * i / split;
//        for(size_t l=0; l<split; l++) {
//            float characteristicVelocity = characteristicVelocity_.getMinValue() + (characteristicVelocity_.getMaxValue() - characteristicVelocity_.getMinValue()) * l / split;
            float characteristicLength = characteristicLength_.get();
            float characteristicVelocity = characteristicVelocity_.get();
            for(size_t j=0; j<6; j++) {
                float viscosity = viscosity_.getMinValue() + (viscosity_.getMaxValue() - viscosity_.getMinValue()) * j / 5;
                for (size_t k = 0; k < 5; k++) {
                    float density = density_.getMinValue() + (density_.getMaxValue() - density_.getMinValue()) * k / 4;

                    for (bool bouzidi : {true, false}) {
                        std::string name = "config"
//                                "_len=" + std::to_string(characteristicLength) + "_vel=" + std::to_string(characteristicVelocity) +
                                "_v=" + std::to_string(viscosity) + "_d=" + std::to_string(density) + "_b=" + std::to_string(bouzidi);

                        FlowParameters parameters(name);
                        parameters.setCharacteristicLength(characteristicLength);
                        parameters.setCharacteristicVelocity(characteristicVelocity);
                        parameters.setViscosity(viscosity);
                        parameters.setDensity(density);
                        parameters.setBouzidi(bouzidi);
                        flowParameters_.push_back(parameters);

                        std::vector<std::string> row(6);
                        row[0] = parameters.getName();
                        row[1] = std::to_string(parameters.getCharacteristicLength());
                        row[2] = std::to_string(parameters.getCharacteristicVelocity());
                        row[3] = std::to_string(parameters.getViscosity());
                        row[4] = std::to_string(parameters.getDensity());
                        row[5] = std::to_string(parameters.getBouzidi());
                        parametrizations_.addRow(row);
                    }
                }
            }
//        }
//    }
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
        flowParametrizationList->setFlowFunction(flowFunction_.getValue());
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
    s.deserialize("flowParameters", flowParameters_,
                  XmlSerializationConstants::ITEMNODE,
                  std::function<FlowParameters()>([]{ return FlowParameters(""); }));
}

}   // namespace