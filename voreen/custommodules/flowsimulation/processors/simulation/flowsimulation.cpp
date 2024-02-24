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

#include "flowsimulation.h"

#include "voreen/core/datastructures/geometry/geometrysequence.h"
#include "voreen/core/datastructures/geometry/glmeshgeometry.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/datastructures/volume/volumeelement.h"
#include "voreen/core/datastructures/volume/volumeminmaxmagnitude.h"
#include "voreen/core/ports/conditions/portconditionvolumelist.h"

#include "modules/flowanalysis/utils/flowutils.h"
#include "modules/hdf5/io/hdf5volumereader.h"
#include "modules/hdf5/io/hdf5volumewriter.h"

#include "../../ext/openlb/voreen/simulation_core.h"

#include <thread>

namespace utils {

using namespace voreen;

// Caution: The returned volume is bound to the life-time of the input volume.
template <typename VolumeType, typename BaseType=typename VolumeElement<typename VolumeType::VoxelType>::BaseType>
std::unique_ptr<Volume> wrapSimpleIntoVoreenVolume(SimpleVolume<BaseType>& volume) {

    auto dimensions = tgt::ivec3::fromPointer(volume.dimensions.data());
    auto spacing = tgt::Vector3<T>::fromPointer(volume.spacing.data());
    auto offset = tgt::Vector3<T>::fromPointer(volume.offset.data());

    auto* data = reinterpret_cast<typename VolumeType::VoxelType*>(volume.data.data());
    auto* representation = new VolumeType(data, dimensions, false);

    std::unique_ptr<Volume> output(new Volume(representation, spacing, offset));
    output->setRealWorldMapping(RealWorldMapping::createDenormalizingMapping<BaseType>());
    //output->addDerivedData(new VolumeMinMax(volume.minValues, volume.maxValues, volume.minValues, volume.maxValues));
    //output->addDerivedData(new VolumeMinMaxMagnitude(volume.minMagnitude, volume.maxMagnitude, volume.minMagnitude, volume.maxMagnitude));
    return output;
}

void indicateErroneousVoxels(SimpleVolume<uint8_t>& volume) {
    const auto errorValue = volume.maxValues[0] + 1;
    auto dim = volume.dimensions;
#ifdef VRN_MODULE_OPENMP
#pragma omp parallel for
#endif
    for (int z = 1; z < dim[2] - 1; z++) {
        for (int y = 1; y < dim[1] - 1; y++) {
            for (int x = 1; x < dim[0] - 1; x++) {
                if (volume.getValue(x, y, z) == 0) {
                    for(int dz : {-1, 1}) {
                        for(int dy : {-1, 1}) {
                            for(int dx : {-1, 1}) {
                                if(volume.getValue(x+dx, y+dy, z+dz) == 1) {
                                    volume.setValue(errorValue, x+dx, y+dy, z+dz);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

}


namespace voreen {

const std::string FlowSimulation::loggerCat_("voreen.flowsimulation.FlowSimulation");

FlowSimulation::FlowSimulation()
    : AsyncComputeProcessor<ComputeInput, ComputeOutput>()
    , geometryDataPort_(Port::INPORT, "geometryDataPort", "Geometry Input", false)
    , segmentationDataPort_(Port::INPORT, "segmentationDataPort", "Segmentation Input", false)
    , measuredDataPort_(Port::INPORT, "measuredDataPort", "Measured Data Input", false)
    , parameterPort_(Port::INPORT, "parameterPort", "Parameterization", false)
    , debugMaterialsPort_(Port::OUTPORT, "debugMaterialsPort", "Debug Materials Port", false, Processor::VALID)
    , debugVelocityPort_(Port::OUTPORT, "debugVelocityPort", "Debug Velocity Port", false, Processor::VALID)
    , debugPressurePort_(Port::OUTPORT, "debugPressurePort", "Debug Pressure Port", false, Processor::VALID)
    , simulationResults_("simulationResults", "Simulation Results", "Simulation Results", "", "", FileDialogProperty::DIRECTORY, Processor::VALID, Property::LOD_DEFAULT, VoreenFileWatchListener::ALWAYS_OFF)
    , deleteOldSimulations_("deleteOldSimulations", "Delete old Simulations", false)
    , simulateAllParametrizations_("simulateAllParametrizations", "Simulate all Parametrizations", false)
    , selectedParametrization_("selectedSimulation", "Selected Parametrization", 0, 0, 0)
    , numCuboids_("numCuboids", "Number of Cuboids", 1, 1, std::thread::hardware_concurrency(), Processor::INVALID_RESULT, IntProperty::STATIC, Property::LOD_DEBUG)
{
    addPort(geometryDataPort_);
    addPort(segmentationDataPort_);
    segmentationDataPort_.addCondition(new PortConditionVolumeListEnsemble());
    segmentationDataPort_.addCondition(new PortConditionVolumeListAdapter(new PortConditionVolumeChannelCount(1)));
    addPort(measuredDataPort_);
    measuredDataPort_.addCondition(new PortConditionVolumeListEnsemble());
    measuredDataPort_.addCondition(new PortConditionVolumeListAdapter(new PortConditionVolumeChannelCount(3)));
    addPort(parameterPort_);
    addPort(debugMaterialsPort_);
    addPort(debugVelocityPort_);
    addPort(debugPressurePort_);

    addProperty(simulationResults_);
    simulationResults_.setGroupID("results");
    addProperty(deleteOldSimulations_);
    deleteOldSimulations_.setGroupID("results");

    addProperty(simulateAllParametrizations_);
    simulateAllParametrizations_.setGroupID("results");
    ON_CHANGE_LAMBDA(simulateAllParametrizations_, [this]{
        selectedParametrization_.setReadOnlyFlag(simulateAllParametrizations_.get());
    });
    addProperty(selectedParametrization_);
    selectedParametrization_.setGroupID("results");

    setPropertyGroupGuiName("results", "Results");

    addProperty(numCuboids_);
    numCuboids_.setGroupID("debug");
    setPropertyGroupGuiName("debug", "Debug");
}

FlowSimulation::~FlowSimulation() {
    if(VoreenApplication* app = VoreenApplication::app()) {
        app->getCommandQueue()->removeAll(this);
    }
}

bool FlowSimulation::isReady() const {
    if(!isInitialized()) {
        setNotReadyErrorMessage("Not initialized.");
        return false;
    }

    if(!parameterPort_.isReady()) {
        setNotReadyErrorMessage("Parameter Port not ready.");
        return false;
    }

    // Node: debug ports are optional.

    return true;
}

void FlowSimulation::adjustPropertiesToInput() {
    const FlowSimulationConfig* config = parameterPort_.getData();
    if(!config || config->empty()) {
        selectedParametrization_.setMinValue(-1);
        selectedParametrization_.setMaxValue(-1);
    }
    else {
        selectedParametrization_.setMinValue(0);
        selectedParametrization_.setMaxValue(static_cast<int>(config->size()) - 1);
        selectedParametrization_.set(0);
    }
}

void FlowSimulation::clearOutports() {
    // We do not clear debug data.
    // It may be available even if the processor is not ready.
}

FlowSimulationInput FlowSimulation::prepareComputeInput() {

    std::unique_ptr<GlMeshGeometryBase> geometry;
    if(const auto* geometryData = dynamic_cast<const GlMeshGeometryBase*>(geometryDataPort_.getData())) {
        geometry.reset(dynamic_cast<GlMeshGeometryBase*>(geometryData->clone().release()));
    }

    // Currently, we prefer the geometry files over the segmentation data.
    std::map<float, std::string> segmentationDataUrls;
    if (!geometry) {
        if (const auto* segmentationData = segmentationDataPort_.getData()) {
            for (size_t i = 0; i < segmentationData->size(); i++) {
                auto* volume = segmentationData->at(i);
                segmentationDataUrls[volume->getTimestep()] = volume->getOrigin().getURL();
            }
        }

        if(segmentationDataUrls.empty()) {
            throw InvalidInputException("No boundary geometry input", InvalidInputException::S_ERROR);
        }
    }
    else if (segmentationDataPort_.hasData()) {
        LWARNING("Both geometry and segmentation input input present. Using geometry input.");
    }

    tgtAssert(measuredDataPort_.isDataInvalidationObservable(), "VolumeListPort must be DataInvalidationObservable!");
    const VolumeList* measuredData = measuredDataPort_.getThreadSafeData();

    tgtAssert(parameterPort_.isDataInvalidationObservable(), "FlowParametrizationPort must be DataInvalidationObservable!");
    auto config = parameterPort_.getThreadSafeData();
    if(!config || config->empty()) {
        throw InvalidInputException("No parameterization", InvalidInputException::S_ERROR);
    }

    if(config->getFlowFeatures() == FF_NONE) {
        throw InvalidInputException("No flow feature selected", InvalidInputException::S_WARNING);
    }

    std::map<float, std::string> measuredDataUrls;
    if(measuredData && !measuredData->empty()) {
        LINFO("Configuring a steered simulation");
        measuredDataUrls[measuredData->first()->getTimestep()] = measuredData->first()->getOrigin().getURL();
        // Check for volume compatibility
        for (size_t i = 1; i < measuredData->size(); i++) {
            VolumeBase* volumeTi = measuredData->at(i);
            measuredDataUrls[volumeTi->getTimestep()] = volumeTi->getOrigin().getURL();

            if (measuredData->at(i-1)->getTimestep() >= volumeTi->getTimestep()) {
                throw InvalidInputException("Time Steps of are not ordered", InvalidInputException::S_ERROR);
            }
        }
    }
    else {
        for(const auto& indicator : config->getFlowIndicators()) {
            if(indicator.type_ == FIT_VELOCITY && indicator.flowProfile_ == FP_VOLUME) {
                throw InvalidInputException("Volume input required", InvalidInputException::S_ERROR);
            }
        }

        measuredData = nullptr;
        LINFO("Configuring an unsteered simulation");
    }

    std::string geometryPath = VoreenApplication::app()->getUniqueTmpFilePath(".stl");
    try {
        std::ofstream file(geometryPath);
        geometry->setTransformationMatrix(geometry->getTransformationMatrix() * config->getTransformationMatrix());
        geometry->exportAsStl(file);
        file.close();
    }
    catch (std::exception&) {
        throw InvalidInputException("Geometry could not be exported", InvalidInputException::S_ERROR);
    }


    if(simulationResults_.get().empty()) {
        throw InvalidInputException("No output directory selected", InvalidInputException::S_WARNING);
    }

    std::string simulationPath = simulationResults_.get() + "/" + config->getName() + "/";
    if (!tgt::FileSystem::createDirectoryRecursive(simulationPath)) {
        throw InvalidInputException("Output directory could not be created", InvalidInputException::S_ERROR);
    }

    size_t selectedParametrization = FlowSimulationConfig::ALL_PARAMETER_SETS;
    if(!simulateAllParametrizations_.get()) {
        selectedParametrization = static_cast<size_t>(selectedParametrization_.get());
    }

    // Clear debug data, so their file handles of old volumes get closed and can be overwritten.
    debugMaterialsPort_.clear();
    debugVelocityPort_.clear();
    debugPressurePort_.clear();

    return FlowSimulationInput{
            geometryPath,
            segmentationDataUrls,
            measuredDataUrls,
            config,
            selectedParametrization,
            simulationPath,
            deleteOldSimulations_.get(),
            numCuboids_.get()
    };
}

FlowSimulationOutput FlowSimulation::compute(FlowSimulationInput input, ProgressReporter& progressReporter) const {

    // Needs to be initialized in each new thread to be used.
    olb::olbInit(nullptr, nullptr);

    // Run either all or just a single simulation.
    if(input.selectedParametrization == FlowSimulationConfig::ALL_PARAMETER_SETS) {
        progressReporter.setProgress(0.0f);
        size_t numRuns = input.config->size();
        for(size_t i=0; i<numRuns; i++) {
            // Define run input.
            FlowSimulationInput runInput = input;
            runInput.selectedParametrization = i;

            // Run the ith simulation.
            SubtaskProgressReporter runProgressReporter(progressReporter, tgt::vec2(i, i+1)/tgt::vec2(numRuns));
            runSimulation(runInput, runProgressReporter);
        }
        progressReporter.setProgress(1.0f);
    }
    else {
        // Pass through.
        runSimulation(input, progressReporter);
    }

    // Done.
    return FlowSimulationOutput{};
}

void FlowSimulation::processComputeOutput(FlowSimulationOutput output) {
    // Nothing to do (yet).
}

void FlowSimulation::serialize(Serializer& s) const {
    AsyncComputeProcessor::serialize(s);

    if(debugMaterialsPort_.hasData()) {
        s.serialize("debugMaterialsPortDataPath", debugMaterialsPort_.getData()->getOrigin());
    }
    if(debugVelocityPort_.hasData()) {
        s.serialize("debugVelocityPortDataPath", debugVelocityPort_.getData()->getOrigin());
    }
    if(debugPressurePort_.hasData()) {
        s.serialize("debugPressurePortDataPath", debugPressurePort_.getData()->getOrigin());
    }
}

void FlowSimulation::deserialize(Deserializer& s) {
    AsyncComputeProcessor::deserialize(s);

    auto deserializeDebugData = [&] (const std::string& key, VolumePort& port) {
        VolumeURL debugPortDataPath;

        try {
            s.deserialize(key, debugPortDataPath);
        } catch (SerializationException& e) {
            s.removeLastError();
            return;
        }

        enqueueInsituResult(debugPortDataPath.getURL(), port);
    };

    deserializeDebugData("debugMaterialsPortDataPath", debugMaterialsPort_);
    deserializeDebugData("debugVelocityPortDataPath", debugVelocityPort_);
    deserializeDebugData("debugPressurePortDataPath", debugPressurePort_);
}

void FlowSimulation::enqueueInsituResult(const std::string& filename, VolumePort& port, std::unique_ptr<Volume> volume) const {

    auto* app = VoreenApplication::app();
    if(!app) {
        return;
    }

    // Write volume to disk, if provided.
    if(volume) {
        try {
            HDF5VolumeWriter().write(filename, volume.get());
        } catch (tgt::FileException& e) {
            LERROR(e.what());
            return;
        }
    }

    // Enqueue setting of port data for synchronous execution.
    app->getCommandQueue()->enqueue(this, LambdaFunctionCallback([filename, &port] {
        port.clear();
        try {
            std::unique_ptr<VolumeList> volumes(HDF5VolumeReaderOriginal().read(filename));
            port.setData(volumes->at(0), true);
        } catch(tgt::FileException& e) {
            LERROR(e.what());
        }
    }));
}

void FlowSimulation::runSimulation(const FlowSimulationInput& input,
                                   ProgressReporter& progressReporter) const {

    const auto segmentationDataUrls = input.segmentationDataUrls_;
    const auto measuredDataUrls = input.measuredDataUrls_;
    const FlowSimulationConfig& config = *input.config;
    const Parameters& parameters = config.at(input.selectedParametrization);

    LINFO("Starting simulation run: " << parameters.name_);
    progressReporter.setProgress(0.0f);

    std::string simulationResultPath = input.simulationResultPath + parameters.name_ + "/";
    if (input.deleteOldSimulations && tgt::FileSystem::dirExists(simulationResultPath)) {
        tgt::FileSystem::deleteDirectoryRecursive(simulationResultPath);
    }
    if (!tgt::FileSystem::createDirectory(simulationResultPath)) {
        LERROR("Output directory could not be created. It may already exist.");
        return;
    }

    singleton::directories().setOutputDir(simulationResultPath);

    UnitConverterFromResolutionAndRelaxationTime<T, DESCRIPTOR> converter(
            parameters.spatialResolution_,              // Resolution that charPhysLength is resolved by.
            (T) parameters.relaxationTime_,             // Relaxation time
            (T) parameters.characteristicLength_,       // charPhysLength: reference length of simulation geometry
            (T) parameters.characteristicVelocity_,     // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
            (T) parameters.viscosity_,                  // physViscosity: physical kinematic viscosity in __m^2 / s__
            (T) parameters.density_                     // physDensity: physical density in __kg / m^3__
    );

    // Prints the converter log as console output
    converter.print();

    // === 2nd Step: Prepare Geometry ===

    // Voxelization. For now, we only support a single geometry.
    STLreader<T> stlReader(input.geometryPath, converter.getConversionFactorLength(), VOREEN_LENGTH_TO_SI, 1);

    // Instantiation.
    IndicatorLayer3D<T> extendedDomain(stlReader, converter.getConversionFactorLength());
    const int noOfCuboids = input.numCuboids;
    CuboidGeometry3D<T> cuboidGeometry(extendedDomain, converter.getConversionFactorLength(), noOfCuboids);
    HeuristicLoadBalancer<T> loadBalancer(cuboidGeometry);
    SuperGeometry<T,3> superGeometry(cuboidGeometry, loadBalancer, 2);
    SuperLattice<T, DESCRIPTOR> lattice(superGeometry);

    // Geometry preparation.
    LINFO("Preparing Geometry ...");
    bool success = prepareGeometry(converter, extendedDomain, stlReader, superGeometry, config.getFlowIndicators());

    // Now already write geometry debug data so that we can see what went wrong..
    SuperLatticeGeometry3D<T, DESCRIPTOR> geometry( lattice, superGeometry );

    // And pass it to the outport.
    {
        const Vector<T, 3> diag = extendedDomain.getMax() - extendedDomain.getMin();
        const int len = std::round(std::max({diag[0], diag[1], diag[2]}) / converter.getConversionFactorLength());
        SimpleVolume<uint8_t> geometryVolume = sampleVolume<uint8_t>(extendedDomain, converter, len, geometry);
        if(!success) utils::indicateErroneousVoxels(geometryVolume);
        auto volume = utils::wrapSimpleIntoVoreenVolume<VolumeRAM_UInt8>(geometryVolume);
        enqueueInsituResult(simulationResultPath + "geometry.h5", debugMaterialsPort_, std::move(volume));
    }

    if(!success) {
        LERROR("The model contains errors! Check resolution and geometry.");
        return;
    }

    interruptionPoint();

    // === 3rd Step: Prepare Lattice ===
    LINFO("Preparing Lattice ...");
    prepareLattice(lattice, converter,
                   stlReader, superGeometry,
                   config.getFlowIndicators(),
                   parameters);

    interruptionPoint();

    // === 4th Step: Main Loop  ===
    util::ValueTracer<T> converge( converter.getLatticeTime(0.5), 1e-5);
    bool swapVelocityFile = false; // See below.

    const int maxIteration = converter.getLatticeTime(config.getSimulationTime());
    auto checkpoint = [&] (int iteration, bool enforce=false) {
        return getResults(
                lattice,
                superGeometry,
                stlReader,
                converter,
                iteration,
                maxIteration,
                simulationResultPath,
                config.getNumTimeSteps(),
                config.getOutputResolution(),
                config.getOutputFileFormat(),
                config.getFlowFeatures(),
                parameters,
                enforce
        );
    };

    for (int iteration = 0; iteration <= maxIteration; iteration++) {

        // === 5th Step: Definition of Initial and Boundary Conditions ===
        setBoundaryValues(lattice, converter, iteration, superGeometry, config.getFlowIndicators(), parameters);

        // === 6th Step: Collide and Stream Execution ===
        lattice.collideAndStream();

        // === 7th Step: Computation and Output of the Results ===
        bool abort = !checkpoint(iteration);
        if(abort) {
            LWARNING("Simulation diverged!");
        }

        // Store last velocity.
        const int outputIter = maxIteration / config.getNumTimeSteps();
        if (abort || iteration % outputIter == 0) {
            {
                // Write velocity file.
                SuperLatticePhysVelocity3D<T, DESCRIPTOR> velocity(lattice, converter);
                auto velocityVolume = sampleVolume<float>(stlReader, converter, config.getOutputResolution(), velocity);
                auto volume = utils::wrapSimpleIntoVoreenVolume<VolumeRAM_3xFloat>(velocityVolume);
                std::string filename = "velocity" + std::to_string(swapVelocityFile) + ".h5";
                enqueueInsituResult(simulationResultPath + filename, debugVelocityPort_, std::move(volume));
            }
            {
                // Write pressure file.
                SuperLatticePhysPressure3D<T, DESCRIPTOR> pressure(lattice, converter);
                auto pressureVolume = sampleVolume<float>(stlReader, converter, config.getOutputResolution(), pressure);
                auto volume = utils::wrapSimpleIntoVoreenVolume<VolumeRAM_Float>(pressureVolume);
                std::string filename = "pressure" + std::to_string(swapVelocityFile) + ".h5";
                enqueueInsituResult(simulationResultPath + filename, debugPressurePort_, std::move(volume));
            }

            // As the old file might still be in use, we need two files in order
            // to ensure that we can write the current time step.
            swapVelocityFile = !swapVelocityFile;
        }

        // === 8th Step: Check for convergence.
        if(!abort) {
            converge.takeValue(lattice.getStatistics().getAverageEnergy(), true);
            if (converge.hasConverged()) {
                LINFO("Simulation converged!");
                abort = true;
            }
        }

        if(abort) {
            checkpoint(iteration, true);
            break;
        }

        float progress = iteration / (maxIteration + 1.0f);
        progressReporter.setProgress(progress);
    }
    progressReporter.setProgress(1.0f);
    LINFO("Finished simulation run: " << parameters.name_);
}

}   // namespace
