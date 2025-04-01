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

#include "flowsimulation.h"

#include "voreen/core/datastructures/geometry/geometrysequence.h"
#include "voreen/core/datastructures/geometry/glmeshgeometry.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/datastructures/volume/volumedecorator.h"
#include "voreen/core/datastructures/volume/volumeelement.h"
#include "voreen/core/datastructures/volume/volumeminmaxmagnitude.h"
#include "voreen/core/ports/conditions/portconditionvolumelist.h"

#include "modules/flowanalysis/utils/flowutils.h"
#include "modules/hdf5/io/hdf5volumereader.h"
#include "modules/hdf5/io/hdf5volumewriter.h"
#ifdef VRN_MODULE_VTK
#include "modules/vtk/io/vtivolumewriter.h"
#endif

#ifndef WIN32

#include "../../ext/openlb/voreen/shared/simulation_core.h"

#include <thread>

namespace utils {

using namespace voreen;

// Caution: The returned volume is bound to the life-time of the input volume.
template <typename VolumeType, typename BaseType=typename VolumeElement<typename VolumeType::VoxelType>::BaseType>
std::unique_ptr<Volume> wrapSimpleIntoVoreenVolume(SimpleVolume<BaseType>& volume, tgt::mat4 physicalToWorldMatrix = tgt::mat4::identity) {

    auto dimensions = tgt::ivec3::fromPointer(volume.dimensions.data());
    auto spacing = tgt::Vector3<T>::fromPointer(volume.spacing.data());
    auto offset = tgt::Vector3<T>::fromPointer(volume.offset.data());

    auto* data = reinterpret_cast<typename VolumeType::VoxelType*>(volume.data.data());
    auto* representation = new VolumeType(data, dimensions, false);

    std::unique_ptr<Volume> output(new Volume(representation, spacing, offset));
    output->setPhysicalToWorldMatrix(physicalToWorldMatrix);
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

#endif

namespace voreen {

const std::string FlowSimulation::loggerCat_("voreen.flowsimulation.FlowSimulation");

FlowSimulation::FlowSimulation()
    : AsyncComputeProcessor<ComputeInput, ComputeOutput>()
    , geometryDataPort_(Port::INPORT, "geometryDataPort", "Geometry Input", false)
    , geometryVolumeDataPort_(Port::INPORT, "geometryVolumeDataPort", "Segmentation Input", false)
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
    addPort(geometryVolumeDataPort_);
    geometryVolumeDataPort_.addCondition(new PortConditionVolumeListEnsemble());
    geometryVolumeDataPort_.addCondition(new PortConditionVolumeListAdapter(new PortConditionVolumeChannelCount(1)));
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

#ifndef WIN32

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

#else

bool FlowSimulation::isReady() const {
    setNotReadyErrorMessage("Processor not available on Windows");
    return false;
}

#endif

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

void FlowSimulation::setNumCuboids(int num) {
    numCuboids_.set(num);
}

int FlowSimulation::getNumCuboids() const {
    return numCuboids_.get();
}

void FlowSimulation::setSimulationResultPath(const std::string &path) {
    simulationResults_.set(path);
}

std::string FlowSimulation::getSimulationResultPath() const {
    return simulationResults_.get();
}

#ifndef WIN32

FlowSimulationInput FlowSimulation::prepareComputeInput() {

    tgtAssert(parameterPort_.isDataInvalidationObservable(), "FlowParametrizationPort must be DataInvalidationObservable!");
    auto originalConfig = parameterPort_.getThreadSafeData();
    if(!originalConfig || originalConfig->empty()) {
        throw InvalidInputException("No parameterization", InvalidInputException::S_ERROR);
    }

    if(originalConfig->getFlowFeatures() == FF_NONE) {
        throw InvalidInputException("No flow feature selected", InvalidInputException::S_WARNING);
    }

    std::unique_ptr<FlowSimulationConfig> config(new FlowSimulationConfig(*originalConfig));

    auto configTransformationMatrix = config->getTransformationMatrix();

    if(const auto* geometryData = dynamic_cast<const GlMeshGeometryBase*>(geometryDataPort_.getData())) {
        std::unique_ptr<GlMeshGeometryBase> geometry(dynamic_cast<GlMeshGeometryBase*>(geometryData->clone().release()));
        std::string geometryPath = VoreenApplication::app()->getUniqueTmpFilePath(".stl");

        try {
            std::ofstream file(geometryPath);
            geometry->setTransformationMatrix(configTransformationMatrix * geometry->getTransformationMatrix());
            geometry->exportAsStl(file);
            file.close();

            std::map<float, std::string> geometryFiles;
            geometryFiles[0.0f] = geometryPath;
            config->setGeometryFiles(geometryFiles, true);
        }
        catch (std::exception&) {
            throw InvalidInputException("Geometry could not be exported", InvalidInputException::S_ERROR);
        }

        if(geometryVolumeDataPort_.hasData()) {
            LWARNING("Both geometry and segmentation input input present. Using geometry input.");
        }
    }
    else {

        if (const auto* segmentationData = geometryVolumeDataPort_.getData()) {
            auto geometryFiles = checkAndConvertVolumeList(segmentationData, configTransformationMatrix);
            config->setGeometryFiles(geometryFiles, false);
        }
    }

    if(config->getGeometryFiles().empty()) {
        throw InvalidInputException("No boundary geometry input", InvalidInputException::S_ERROR);
    }

    tgtAssert(measuredDataPort_.isDataInvalidationObservable(), "VolumeListPort must be DataInvalidationObservable!");
    const VolumeList* measuredData = measuredDataPort_.getThreadSafeData();

    std::map<float, std::string> measuredDataUrls;
    if(measuredData && !measuredData->empty()) {
        LINFO("Configuring a steered simulation");
        auto measuredDataFiles = checkAndConvertVolumeList(measuredData, configTransformationMatrix);
        config->setMeasuredDataFiles(measuredDataFiles);
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
            std::move(config),
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
            // Set run id.
            input.selectedParametrization = i;

            // Run the ith simulation.
            SubtaskProgressReporter runProgressReporter(progressReporter, tgt::vec2(i, i+1)/tgt::vec2(numRuns));
            runSimulation(input, runProgressReporter);
        }
        progressReporter.setProgress(1.0f);
    }
    else {
        // Pass through.
        runSimulation(input, progressReporter);
    }

    // Done.
    return {};
}

#else

FlowSimulationInput FlowSimulation::prepareComputeInput() {
    throw InvalidInputException("This processor is not available on windows", InvalidInputException::S_ERROR);
}

FlowSimulationOutput FlowSimulation::compute(FlowSimulationInput, ProgressReporter&) const {
    return {};
}

#endif

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

#ifndef WIN32

void FlowSimulation::runSimulation(const FlowSimulationInput& input,
                                   ProgressReporter& progressReporter) const {

    auto& config = *input.config;

    const Parameters& parameters = config.at(input.selectedParametrization);

    auto physicalToWorldMatrix = config.getInvertedTransformationMatrix();

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

    VolumeTimeSeries measuredDataTimeSeries(config.getMeasuredDataFiles());

    // === 2nd Step: Prepare Geometry ===

    // Voxelization. For now, we only support a single geometry.
    std::unique_ptr<IndicatorF3D<T>> boundaryGeometry;
    std::unique_ptr<VolumeTimeSeries> geometryVolumeTimeSeries;
    if(config.isGeometryMesh()) {
        std::string geometryFileName = config.getGeometryFiles().begin()->second;
        boundaryGeometry.reset(new STLreader<T>(geometryFileName, converter.getConversionFactorLength(), VOREEN_LENGTH_TO_SI, 1));
    }
    else {
        geometryVolumeTimeSeries.reset(new VolumeTimeSeries(config.getGeometryFiles()));
        boundaryGeometry.reset(new VolumeDataMapperIndicator(geometryVolumeTimeSeries->createSampler(0.0f)));
    }
    IndicatorLayer3D<T> extendedDomain(*boundaryGeometry, converter.getConversionFactorLength());
    const int noOfCuboids = input.numCuboids;
    CuboidGeometry3D<T> cuboidGeometry(extendedDomain, converter.getConversionFactorLength(), noOfCuboids);
    HeuristicLoadBalancer<T> loadBalancer(cuboidGeometry);
    SuperGeometry<T,3> superGeometry(cuboidGeometry, loadBalancer, 2);
    SuperLattice<T, DESCRIPTOR> lattice(superGeometry);

    // Geometry preparation.
    LINFO("Preparing Geometry ...");
    bool success = prepareGeometry(converter, extendedDomain, *boundaryGeometry, superGeometry, config.getFlowIndicators(true));

    // Now already write geometry debug data so that we can see what went wrong..
    SuperLatticeGeometry3D<T, DESCRIPTOR> geometry( lattice, superGeometry );

    // And pass it to the outport.
    {
        const Vector<T, 3> diag = extendedDomain.getMax() - extendedDomain.getMin();
        const int len = std::round(std::max({diag[0], diag[1], diag[2]}) / converter.getConversionFactorLength());
        SimpleVolume<uint8_t> geometryVolume = sampleVolume<uint8_t>(extendedDomain, converter, len, geometry);
        if(!success) utils::indicateErroneousVoxels(geometryVolume);
        auto volume = utils::wrapSimpleIntoVoreenVolume<VolumeRAM_UInt8>(geometryVolume, physicalToWorldMatrix);
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
                   *boundaryGeometry, superGeometry,
                   config.getFlowIndicators(true),
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
                *boundaryGeometry,
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

    // Writes velocity and pressure to their respective debug ports.
    auto writeDebugData = [&] () {
        {
            // Write velocity file.
            SuperLatticePhysVelocity3D<T, DESCRIPTOR> velocity(lattice, converter);
            auto velocityVolume = sampleVolume<float>(*boundaryGeometry, converter, config.getOutputResolution(), velocity);
            auto volume = utils::wrapSimpleIntoVoreenVolume<VolumeRAM_3xFloat>(velocityVolume, physicalToWorldMatrix);
            std::string filename = "velocity" + std::to_string(swapVelocityFile) + ".h5";
            enqueueInsituResult(simulationResultPath + filename, debugVelocityPort_, std::move(volume));
        }
        {
            // Write pressure file.
            SuperLatticePhysPressure3D<T, DESCRIPTOR> pressure(lattice, converter);
            auto pressureVolume = sampleVolume<float>(*boundaryGeometry, converter, config.getOutputResolution(), pressure);
            auto volume = utils::wrapSimpleIntoVoreenVolume<VolumeRAM_Float>(pressureVolume, physicalToWorldMatrix);
            std::string filename = "pressure" + std::to_string(swapVelocityFile) + ".h5";
            enqueueInsituResult(simulationResultPath + filename, debugPressurePort_, std::move(volume));
        }
    };

    // The simulation loop.
    for (int iteration = 0; iteration <= maxIteration; iteration++) {

        // === 5th Step: Definition of Initial and Boundary Conditions ===
        setBoundaryValues(lattice, converter, iteration, superGeometry, config.getFlowIndicators(true), parameters, &measuredDataTimeSeries);

        // Store debug data.
        const int outputIter = maxIteration / config.getNumTimeSteps();
        if (outputIter == 0 || iteration % outputIter == 0) {
            writeDebugData();

            // As the old file might still be in use, we need two files in order
            // to ensure that we can write the current time step.
            swapVelocityFile = !swapVelocityFile;
        }

        // === 6th Step: Collide and Stream Execution ===
        lattice.collideAndStream();

        // === 7th Step: Computation and Output of the Results ===
        bool abort = !checkpoint(iteration);
        if(abort) {
            LWARNING("Simulation diverged!");
        }
        // === 8th Step: Check for convergence.
        else {
            converge.takeValue(lattice.getStatistics().getAverageEnergy(), true);
            if (converge.hasConverged()) {
                LINFO("Simulation converged!");
                abort = true;
            }
        }

        // Write last valid time step.
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

#else

void FlowSimulation::runSimulation(const FlowSimulationInput&, ProgressReporter&) const {}

#endif

std::map<float, std::string> FlowSimulation::checkAndConvertVolumeList(const VolumeList* volumes, tgt::mat4 transformation) const {
    std::map<float, std::string> result;
    for (size_t i = 0; i < volumes->size(); i++) {
        VolumeBase* volume = volumes->at(i);
        auto path = volume->getOrigin().getPath();
        auto extension = tgt::FileSystem::fileExtension(path, true);
        if (extension != "vti" || transformation != tgt::mat4::identity) {
            path = VoreenApplication::app()->getUniqueTmpFilePath(".vti");
#ifdef VRN_MODULE_VTK
            LWARNING("Need to convert to vti..");
            std::unique_ptr<VolumeDecoratorReplace> transformedVolume(
                    new VolumeDecoratorReplaceTransformation(volume, transformation * volume->getPhysicalToWorldMatrix())
            );

            VTIVolumeWriter().writeVectorField(path, transformedVolume.get());
#else
            throw InvalidInputException("Need to convert to vti, but VTI module not enabled", InvalidInputException::S_ERROR);
#endif
        }
        std::string url = path + "?fieldName=" + volume->getModality().getName();
        result[volume->getTimestep()] = url;
    }

    return result;
}

}   // namespace
