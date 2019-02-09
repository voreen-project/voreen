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
    , outport_(Port::OUTPORT, "outport", "Parameter Port")
    , parametrizationName_("parametrizationName", "Parametrization Name", "test_parametrization")
    , simulationTime_("simulationTime", "Simulation Time (s)", 2.0f, 0.1f, 10.0f)
    , temporalResolution_("temporalResolution", "Temporal Resolution (ms)", 3.1f, 1.0f, 30.0f)
    , characteristicLength_("characteristicLength", "Characteristic Length (mm)", 22.46f, 1.0f, 100.0f)
    , viscosity_("viscosity", "Viscosity (e-6 m^2/s)", 3.5, 3, 4)
    , density_("density", "Density (kg/m^3)", 1000.0f, 1000.0f, 1100.0f)
    , bouzidi_("bounzidi", "Bounzidi", true)
    , addParametrization_("addParametrization", "Add Parametrization")
    , removeParametrization_("removeParametrization", "Remove Parametrization")
    , clearParametrizations_("clearParametrizations", "Clear Parametrizations")
    , autoGenerateParametrizations_("autoGenerateParametrizations", "auto-generate parametrizations")
    , ensembleName_("ensembleName", "Ensemble Name", "test_ensemble")
    , parametrizations_("parametrizations", "Parametrizations", 5, Processor::VALID)
{
    addPort(outport_);

    addProperty(parametrizationName_);
        parametrizationName_.setGroupID("parameters");
    addProperty(simulationTime_);
        simulationTime_.setGroupID("parameters");
    addProperty(temporalResolution_);
        temporalResolution_.setGroupID("parameters");
    addProperty(characteristicLength_);
        characteristicLength_.setGroupID("parameters");
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

    addProperty(ensembleName_);
    addProperty(parametrizations_);
    parametrizations_.setColumnLabel(0, "Name");
    parametrizations_.setColumnLabel(1, "Char. Len.");
    parametrizations_.setColumnLabel(2, "Viscosity");
    parametrizations_.setColumnLabel(3, "Density");
    parametrizations_.setColumnLabel(4, "Bouzidi");
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
    parameters.setSimulationTime(simulationTime_.get());
    parameters.setTemporalResolution(temporalResolution_.get());
    parameters.setCharacteristicLength(characteristicLength_.get());
    parameters.setViscosity(viscosity_.get());
    parameters.setDensity(density_.get());
    parameters.setBouzidi(bouzidi_.get());
    flowParameters_.push_back(parameters);

    std::vector <std::string> row(5);
    row[0] = parameters.getName();
    row[1] = std::to_string(parameters.getCharacteristicLength());
    row[2] = std::to_string(parameters.getViscosity());
    row[3] = std::to_string(parameters.getDensity());
    row[4] = std::to_string(parameters.getBouzidi());
    parametrizations_.addRow(row);
}
void FlowParametrization::removeParametrization() {
    if(parametrizations_.getNumRows() > 0 && parametrizations_.getSelectedRowIndex() != -1) {
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
    const size_t split = 3;

    std::vector<float> characteristicLengths;
    for(size_t i=0; i<split; i++) {
        float characteristicLength = characteristicLength_.getMinValue() + (characteristicLength_.getMaxValue() - characteristicLength_.getMinValue()) * i / split;
        for(size_t j=0; j<split; j++) {
            float viscosity = viscosity_.getMinValue() + (viscosity_.getMaxValue() - viscosity_.getMinValue()) * j / split;
            for(size_t k=0; k<split; k++) {
                float density = density_.getMinValue() + (density_.getMaxValue() - density_.getMinValue()) * k / split;

                for(bool bouzidi : {true, false}) {
                    std::string name = "run_l=" + std::to_string(characteristicLength) + "_v=" + std::to_string(viscosity) + "_d=" + std::to_string(density) + "_b=" + std::to_string(bouzidi);
                    FlowParameters parameters(name);
                    parameters.setSimulationTime(simulationTime_.get());
                    parameters.setTemporalResolution(temporalResolution_.get());
                    parameters.setCharacteristicLength(characteristicLength);
                    parameters.setViscosity(viscosity);
                    parameters.setDensity(density);
                    parameters.setBouzidi(bouzidi);
                    flowParameters_.push_back(parameters);

                    std::vector <std::string> row(5);
                    row[0] = parameters.getName();
                    row[1] = std::to_string(parameters.getCharacteristicLength());
                    row[2] = std::to_string(parameters.getViscosity());
                    row[3] = std::to_string(parameters.getDensity());
                    row[4] = std::to_string(parameters.getBouzidi());
                    parametrizations_.addRow(row);
                }
            }
        }
    }
}

void FlowParametrization::process() {

    FlowParametrizationList* flowParametrizationList = new FlowParametrizationList(ensembleName_.get());

    for(const FlowParameters& flowParameters : flowParameters_) {
        flowParametrizationList->addFlowParameters(flowParameters);
    }

    outport_.setData(flowParametrizationList);

}

}   // namespace
