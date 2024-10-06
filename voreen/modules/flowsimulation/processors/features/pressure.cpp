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

#include "pressure.h"

#include "voreen/core/datastructures/volume/volumeminmaxmagnitude.h"
#include "voreen/core/ports/conditions/portconditionvolumetype.h"


namespace voreen {

Pressure::Pressure()
    : AsyncComputeProcessor()
    , geometryDataPort_(Port::INPORT, "geometryDataPort", "Geometry Input", false)
    , geometryVolumeDataPort_(Port::INPORT, "geometryVolumeDataPort", "Segmentation Input", false)
    , measuredDataPort_(Port::INPORT, "measuredDataPort", "Measured Data Input", false)
    , pressurePort_(Port::OUTPORT, "pressurePort", "Pressure Output", false, Processor::VALID)
    , kinematicViscosity_("kinematicViscosityProp", "Kinematic viscosity (x10^(-6) m^2/s)", 1.35f, 0.1f, 100.0f)
    , density_("densityProp", "Density (kg/m^3)", 988.21f, 0.1f, 10000.0f)
    , geometryDataPortInternal_(Port::OUTPORT, "geometry", "", false)
    , geometryVolumeDataPortInternal_(Port::OUTPORT, "geometryVolumeData", "", false, Processor::VALID)
    , measuredDataPortInternal_(Port::OUTPORT, "volumelist", "", false, Processor::VALID)
    , configPortInternal_(Port::OUTPORT, "config", "", false)
{
    addPort(geometryDataPort_);
    addPort(geometryVolumeDataPort_);
    addPort(measuredDataPort_);
    measuredDataPort_.addCondition(new PortConditionVolumeChannelCount(3));
    addPort(pressurePort_);

    addProperty(kinematicViscosity_);
    addProperty(density_);

    auto simulationOutputPath = VoreenApplication::app() ? VoreenApplication::app()->getUniqueTmpFilePath("") : "simulation";
    flowSimulation_.setSimulationResultPath(simulationOutputPath);
    flowSimulation_.setNumCuboids(1);

    geometryDataPortInternal_.setProcessor(const_cast<Pressure*>(this));
    geometryDataPortInternal_.connect(flowSimulation_.getPort("geometryDataPort"));

    geometryVolumeDataPortInternal_.setProcessor(const_cast<Pressure*>(this));
    geometryVolumeDataPortInternal_.connect(flowSimulation_.getPort("geometryVolumeDataPort"));

    measuredDataPortInternal_.setProcessor(const_cast<Pressure*>(this));
    measuredDataPortInternal_.connect(flowSimulation_.getPort("measuredDataPort"));

    configPortInternal_.setProcessor(const_cast<Pressure*>(this));
    configPortInternal_.connect(flowSimulation_.getPort("parameterPort"));
}

Pressure::~Pressure() {
    if(VoreenApplication* app = VoreenApplication::app()) {
        app->getCommandQueue()->removeAll(this);
    }
}

Processor* Pressure::create() const {
    return new Pressure();
}

FlowSimulationInput Pressure::prepareComputeInput() {

    auto* measuredData = measuredDataPort_.getData();
    if(!measuredData) {
        throw InvalidInputException("No velocity data available", InvalidInputException::S_ERROR);
    }

    tgt::vec3 spacing = measuredData->getSpacing();
    tgt::svec3 dim = measuredData->getDimensions();

    int resolution = dim.x;
    float characteristicLength = spacing.x * resolution;
    
    auto* vmmm = measuredData->getDerivedData<VolumeMinMaxMagnitude>();
    float maxVelocityMagnitude = vmmm->getMaxMagnitude();

    // Set up an indicator that covers the entire fluid domain.
    FlowIndicator indicator;
    indicator.type_ = FIT_VELOCITY;
    indicator.flowProfile_ = FP_VOLUME;
    indicator.id_ = MAT_FLUID;
    indicator.velocityCurve_ = VelocityCurve::createConstantCurve(1.0f);

    Parameters parameters;
    parameters.name_ = "pressure";
    parameters.density_ = density_.get();
    parameters.viscosity_ = kinematicViscosity_.get() * 1e-6f;
    parameters.turbulenceModel_ = FTM_BGK;
    // parameters.smagorinskyConstant_ = 0.0f; // unused.
    parameters.characteristicLength_ = characteristicLength * 0.001f; // Convert to meters.
    parameters.characteristicVelocity_ = maxVelocityMagnitude * 1.1f; // Grant some wiggle room.
    parameters.inletVelocityMultiplier_ = 1.0f;
    parameters.relaxationTime_ = 1.0f;
    parameters.latticePerturbation_ = false;
    parameters.wallBoundaryCondition_ = FBC_BOUNCE_BACK;
    parameters.spatialResolution_ = resolution;

    auto* config = new FlowSimulationConfig("pressure");
    config->addFlowIndicator(indicator, true);
    config->addFlowParameterSet(parameters);
    config->setFlowFeatures(FF_PRESSURE); // To prevent simulation to emit warning.
    config->setNumTimeSteps(1);
    config->setSimulationTime(0.0f);
    config->setOutputResolution(resolution);
    config->setOutputFileFormat("vti");

    geometryDataPortInternal_.setData(geometryDataPort_.getData(), false);

    auto geometryVolumes = new VolumeList();
    if (geometryVolumeDataPort_.hasData()) {
        geometryVolumes->add(const_cast<VolumeBase*>(geometryVolumeDataPort_.getData()));
    }
    geometryVolumeDataPortInternal_.setData(geometryVolumes, true);

    auto measuredDataVolumes = new VolumeList();
    measuredDataVolumes->add(const_cast<VolumeBase*>(measuredData));
    measuredDataPortInternal_.setData(measuredDataVolumes, true);

    configPortInternal_.setData(config, true);

    auto input = flowSimulation_.prepareComputeInput();
    input.deleteOldSimulations = true;
    return input;
}

FlowSimulationOutput Pressure::compute(FlowSimulationInput input, ProgressReporter& progressReporter) const {

    try {
        flowSimulation_.compute(std::move(input), progressReporter);

        if(auto* app = VoreenApplication::app()) {
            app->getCommandQueue()->enqueue(this, LambdaFunctionCallback([this] {
                if(auto* port = dynamic_cast<VolumePort*>(flowSimulation_.getPort("debugPressurePort"))) {
                    pressurePort_.setData(port->getData(), false);
                }
            }));
        }
    }
    catch (InvalidInputException& e) {
        LERRORC("Pressure", e.what());
        return {};
    }

    return {};
}

void Pressure::processComputeOutput(FlowSimulationOutput output) {
}

bool Pressure::isReady() const {
    if(!isInitialized()) {
        setNotReadyErrorMessage("Not initialized");
        return false;
    }

    if(!measuredDataPort_.isReady()) {
        setNotReadyErrorMessage("No velocity data available");
        return false;
    }

    return true;
}

}