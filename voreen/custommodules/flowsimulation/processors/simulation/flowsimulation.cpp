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
#include "voreen/core/datastructures/volume/volumefactory.h"
#include "voreen/core/datastructures/volume/volumeminmaxmagnitude.h"
#include "voreen/core/ports/conditions/portconditionvolumelist.h"

#include "voreen/core/utils/stringutils.h"

#include "modules/flowanalysis/utils/flowutils.h"

#include "../../ext/openlb/voreen/simulation_core.h"

#include <thread>

namespace {

using namespace voreen;

// Allow for simulation lattice initialization by volume data.
class MeasuredDataMapper : public AnalyticalF3D<T, T> {
public:
    MeasuredDataMapper(UnitConverter<T, DESCRIPTOR> const& converter, const VelocitySampler& sampler, float multiplier = 1.0f)
        : AnalyticalF3D<T, T>(3)
        , converter_(converter)
        , sampler_(sampler)
        , multiplier_(multiplier)
    {
    }
    virtual bool operator() (T output[], const T input[]) {
        // Store simulation positions in world coordinates.
        tgt::vec3 pos = tgt::vec3(tgt::Vector3<T>::fromPointer(input)) / VOREEN_LENGTH_TO_SI;
        tgt::vec3 vel = sampler_.sample(pos) * multiplier_;

        for(size_t i=0; i<3; i++) {
            output[i] = converter_.getLatticeVelocity(vel[i]);
        }

        return true;
    }

private:
    const UnitConverter<T, DESCRIPTOR>& converter_;
    const VelocitySampler& sampler_;
    const float multiplier_;
};

template<typename T>
struct TypedGeometryPair {

    TypedGeometryPair(const Geometry* lhs, const Geometry* rhs)
        : lhs_(dynamic_cast<const T*>(lhs)), rhs_(dynamic_cast<const T*>(rhs)) {}

    operator bool () {
        return lhs_ != nullptr && rhs_ != nullptr;
    }

    const T* lhs_;
    const T* rhs_;
};

template<typename T>
std::unique_ptr<GlMeshGeometryBase> mergeGeometriesTyped(const T* lhs, const T* rhs) {
    tgtAssert(lhs, "lhs null");
    tgtAssert(rhs, "rhs null");

    auto merged = lhs->clone().release();
    auto mergedTyped = static_cast<T*>(merged);

    // Transform all vertex position to world space.
    auto vertices = mergedTyped->getVertices();
    for(auto& vertex : vertices) {
        vertex.pos_ = mergedTyped->getTransformationMatrix() * vertex.pos_;
    }
    mergedTyped->setTransformationMatrix(tgt::mat4::identity);

    // Copy over and transform vertices.
    for(auto vertex : rhs->getVertices()) {
        vertex.pos_ = rhs->getTransformationMatrix() * vertex.pos_;
        vertices.emplace_back(vertex);
    }
    mergedTyped->setVertices(vertices);

    // Add and adjust indices, if required.
    if (lhs->usesIndexedDrawing()) {
        auto indexOffset = lhs->getNumVertices();
        for (auto index : rhs->getIndices()) {
            mergedTyped->addIndex(index + indexOffset);
        }
    }

    return std::unique_ptr<GlMeshGeometryBase>(mergedTyped);
}


std::unique_ptr<GlMeshGeometryBase> mergeGeometries(const GlMeshGeometryBase* lhs, const GlMeshGeometryBase* rhs) {
    if(lhs->getPrimitiveType() != rhs->getPrimitiveType()) {
        return nullptr;
    }

    if(lhs->getVertexLayout() != rhs->getVertexLayout()) {
        return nullptr;
    }

    if(lhs->getIndexType() != rhs->getIndexType()) {
        return nullptr;
    }

    if(lhs->usesIndexedDrawing() != rhs->usesIndexedDrawing()) {
        return nullptr;
    }

    if(auto typedGeometry = TypedGeometryPair<GlMeshGeometryUInt32Normal>(lhs, rhs)) {
        return mergeGeometriesTyped(typedGeometry.lhs_, typedGeometry.rhs_);
    }

    return nullptr;
}
/*
template <typename I>
Volume* convertSimpleVolumeToVoreenVolume(const SimpleVolume<I>& volume) {
    auto dimensions = tgt::ivec3::fromPointer(volume.dimensions.data());
    auto spacing = tgt::Vector3<T>::fromPointer(volume.spacing.data());
    auto offset = tgt::Vector3<T>::fromPointer(volume.offset.data());

    std::string baseType = getBaseTypeFromType<I>();
    std::string format = getFormatFromBaseTypeAndChannels(baseType, volume.numChannels);
    auto* representation = VolumeFactory().create(format, dimensions);
    for(size_t i=0, numVoxels = tgt::hmul(dimensions); i < numVoxels; i++) {
        for (int channel = 0; channel < volume.numChannels; channel++) {
            representation->voxel(i)[channel] = volume.getValue(i, channel), i, channel);
        }
    }

    auto* output = new Volume(representation, spacing, offset);
    output->setRealWorldMapping(RealWorldMapping::createDenormalizingMapping(baseType));
    return output;
}
*/
}


namespace voreen {

const std::string FlowSimulation::loggerCat_("voreen.flowsimulation.FlowSimulation");

FlowSimulation::FlowSimulation()
    : AsyncComputeProcessor<ComputeInput, ComputeOutput>()
    , geometryDataPort_(Port::INPORT, "geometryDataPort", "Geometry Input", false)
    , measuredDataPort_(Port::INPORT, "measuredDataPort", "Measured Data Input", false)
    , parameterPort_(Port::INPORT, "parameterPort", "Parameterization", false)
    , debugMaterialsPort_(Port::OUTPORT, "debugMaterialsPort", "Debug Materials Port", false, Processor::VALID)
    , debugVelocityPort_(Port::OUTPORT, "debugVelocityPort", "Debug Velocity Port", false, Processor::VALID)
    , simulationResults_("simulationResults", "Simulation Results", "Simulation Results", VoreenApplication::app()->getTemporaryPath("simulation"), "", FileDialogProperty::DIRECTORY, Processor::VALID, Property::LOD_DEFAULT, VoreenFileWatchListener::ALWAYS_OFF)
    , deleteOldSimulations_("deleteOldSimulations", "Delete old Simulations", false)
    , simulateAllParametrizations_("simulateAllParametrizations", "Simulate all Parametrizations", false)
    , selectedParametrization_("selectedSimulation", "Selected Parametrization", 0, 0, 0)
{
    addPort(geometryDataPort_);
    addPort(measuredDataPort_);
    measuredDataPort_.addCondition(new PortConditionVolumeListEnsemble());
    measuredDataPort_.addCondition(new PortConditionVolumeListAdapter(new PortConditionVolumeChannelCount(3)));
    addPort(parameterPort_);
    addPort(debugMaterialsPort_);
    addPort(debugVelocityPort_);

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
    if(!geometryDataPort_.isReady()) {
        setNotReadyErrorMessage("Geometry Port not ready.");
        return false;
    }

    // Note: measuredDataPort is optional!
    if(measuredDataPort_.hasData() && !measuredDataPort_.isReady()) {
        setNotReadyErrorMessage("Measured Data Port not ready.");
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

FlowSimulationInput FlowSimulation::prepareComputeInput() {

    PortDataPointer<GlMeshGeometryBase> geometry(nullptr, false);

    if(auto geometrySequence = dynamic_cast<const GeometrySequence*>(geometryDataPort_.getData())) {

        std::vector<const GlMeshGeometryBase*> geometries;

        for(size_t i=0; i<geometrySequence->getNumGeometries(); i++) {
            auto geometry = dynamic_cast<const GlMeshGeometryBase*>(geometrySequence->getGeometry(i));
            if(geometry) {
                geometries.push_back(geometry);
            }
            else {
                throw InvalidInputException("GeometrySequence contains non-GlMeshGeometry", InvalidInputException::S_WARNING);
            }
        }

        if(!geometries.empty()) {
            geometry = PortDataPointer<GlMeshGeometryBase>(geometries.front(), false);
            for (size_t i=1; i<geometries.size(); i++) {
                auto combined = mergeGeometries(geometry, geometries.at(i));

                if(!combined) {
                    throw InvalidInputException("Encountered Mesh of unexpected Type", InvalidInputException::S_ERROR);
                }

                geometry = PortDataPointer<GlMeshGeometryBase>(combined.release(), true);
            }
        }
    }
    else if(auto geometryData = dynamic_cast<const GlMeshGeometryBase*>(geometryDataPort_.getData())) {
        geometry = PortDataPointer<GlMeshGeometryBase>(geometryData, false);
    }

    if (!geometry) {
        throw InvalidInputException("No GlMeshGeometry input", InvalidInputException::S_ERROR);
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

    if(measuredData && !measuredData->empty()) {
        LINFO("Configuring a steered simulation");
        // Check for volume compatibility
        for (size_t i = 1; i < measuredData->size(); i++) {
            VolumeBase* volumeTi = measuredData->at(i);

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

    return FlowSimulationInput{
            geometryPath,
            measuredData,
            config,
            selectedParametrization,
            simulationPath,
            deleteOldSimulations_.get()
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

void FlowSimulation::runSimulation(const FlowSimulationInput& input,
                                   ProgressReporter& progressReporter) const {

    const VolumeList* measuredData = input.measuredData;
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

    VolumeSampler volumeSampler = [&] (UnitConverter<T, DESCRIPTOR> const& converter, float time, std::function<void(AnalyticalF3D<T,T>&)>& target, float multiplier=1.0f) {
        // If a single time step is attached, we use it throughout the entire simulation.
        if(measuredData->size() == 1) {
            const VolumeBase* volume = measuredData->first();
            VolumeRAMRepresentationLock data(measuredData->first());
            SpatialSampler sampler(*data, volume->getRealWorldMapping(), VolumeRAM::LINEAR, volume->getWorldToVoxelMatrix());
            MeasuredDataMapper mapper(converter, sampler, multiplier);
            target(mapper);
        }
        else {
            // If we have multiple time steps, we need to update the volume we sample from.
            tgtAssert(measuredData->size() > 1, "expected more than 1 volume");

            // Query time.
            // Note that we periodically sample the volumes if the simulation time exceeds measurement time.
            // TODO: Make adjustable.
            float start = measuredData->first()->getTimestep();
            float end   = measuredData->at(measuredData->size() - 1)->getTimestep();
            time        = std::fmod(time - start, end - start);

            // Find the volume whose time step is right before the current time.
            size_t idx = 0;
            while (idx < measuredData->size() - 1 && measuredData->at(idx + 1)->getTimestep() < time) idx++;

            const VolumeBase* volume0 = measuredData->at(idx + 0);
            VolumeRAMRepresentationLock data0(volume0);

            const VolumeBase* volume1 = measuredData->at(idx + 1);
            VolumeRAMRepresentationLock data1(volume1);

            float alpha = (time - volume0->getTimestep()) / (volume1->getTimestep() - volume0->getTimestep());
            SpatioTemporalSampler sampler(*data0, *data1, alpha, volume0->getRealWorldMapping(),
                                          VolumeRAM::LINEAR, volume0->getWorldToVoxelMatrix());

            MeasuredDataMapper mapper(converter, sampler, multiplier);
            target(mapper);
        }
    };

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

    // Instantiation of the STLreader class
    // file name, voxel size in meter, stl unit in meter, outer voxel no., inner voxel no.
    STLreader<T> stlReader(input.geometryPath, converter.getConversionFactorLength(), VOREEN_LENGTH_TO_SI, 1);
    IndicatorLayer3D<T> extendedDomain(stlReader, converter.getConversionFactorLength());

    // Instantiation of a cuboidGeometry with weights
    const int noOfCuboids = 1;//std::thread::hardware_concurrency();
    CuboidGeometry3D<T> cuboidGeometry(extendedDomain, converter.getConversionFactorLength(), noOfCuboids);

    // Instantiation of a loadBalancer
    HeuristicLoadBalancer<T> loadBalancer(cuboidGeometry);

    interruptionPoint();

    // Instantiation of a superGeometry
    LINFO("Preparing Geometry ...");
    SuperGeometry<T,3> superGeometry(cuboidGeometry, loadBalancer, 2);

    bool success = prepareGeometry(converter, extendedDomain, stlReader, superGeometry, config.getFlowIndicators());
    if(!success) {
        LERROR("The model contains errors! Check resolution and geometry.");
        return;
    }

    interruptionPoint();

    // === 3rd Step: Prepare Lattice ===
    LINFO("Preparing Lattice ...");
    SuperLattice<T, DESCRIPTOR> lattice(superGeometry);
    prepareLattice(lattice, converter,
                   stlReader, superGeometry,
                   config.getFlowIndicators(),
                   parameters);

    interruptionPoint();

    // Always write geometry debug data..
    SuperVTMwriter3D<T> vtmWriter( "results" );
    SuperLatticeGeometry3D<T, DESCRIPTOR> geometry( lattice, superGeometry );
    vtmWriter.write( geometry );
    vtmWriter.createMasterFile();
/*
    // And pass it to the outport.
    if(VoreenApplication* app = VoreenApplication::app()) {
        app->getCommandQueue()->removeAll(this);

        auto* volume = convertSimpleVolumeToVoreenVolume(sampleVolume<int>(stlReader, converter, config.getOutputResolution(), geometry));

        app->getCommandQueue()->enqueue(this, LambdaFunctionCallback([&] {
            debugMaterialsPort_.setData(volume, true);
        }));
    }
*/
    // And pass it to the outport.
    if(VoreenApplication* app = VoreenApplication::app()) {
        SimpleVolume<int> geometryVolume = sampleVolume<int>(extendedDomain, converter, config.getOutputResolution(), geometry);

        auto dimensions = tgt::ivec3::fromPointer(geometryVolume.dimensions.data());
        auto spacing = tgt::Vector3<T>::fromPointer(geometryVolume.spacing.data());
        auto offset = tgt::Vector3<T>::fromPointer(geometryVolume.offset.data());

        auto* representation = new VolumeRAM_UInt8(dimensions);
        for(size_t i=0, numVoxels = tgt::hmul(dimensions); i < numVoxels; i++) {
            representation->voxel(i) = geometryVolume.data[i];
        }

        auto* volume = new Volume(representation, spacing, offset);
        volume->setRealWorldMapping(RealWorldMapping::createDenormalizingMapping<uint8_t>());

        app->getCommandQueue()->enqueue(this, LambdaFunctionCallback([&] {
            debugMaterialsPort_.setData(volume, true);
        }));
    }

    interruptionPoint();

    // === 4th Step: Main Loop  ===
    const int maxIteration = converter.getLatticeTime(config.getSimulationTime());
    util::ValueTracer<T> converge( converter.getLatticeTime(0.5), 1e-5);
    for (int iteration = 0; iteration <= maxIteration; iteration++) {

        // === 5th Step: Definition of Initial and Boundary Conditions ===
        setBoundaryValues(lattice, converter, iteration, superGeometry, config.getFlowIndicators(), parameters, volumeSampler);

        // === 6th Step: Collide and Stream Execution ===
        lattice.collideAndStream();

        // === 7th Step: Computation and Output of the Results ===
        bool success = getResults(lattice, converter, iteration, maxIteration, superGeometry, stlReader,
                                  simulationResultPath,
                                  config.getNumTimeSteps(),
                                  config.getOutputResolution(),
                                  config.getOutputFileFormat(),
                                  config.getFlowFeatures(),
                                  parameters);

        /*
        SuperLatticePhysVelocity3D<T, DESCRIPTOR> velocity(lattice, converter);
        auto volume = sampleVolume<float>(stlReader, converter, config.getOutputResolution(), velocity);
        someFunction(velocity, volume, debugVelocityPort_);

        // Store last velocity.
        if(VoreenApplication* app = VoreenApplication::app()) {
            app->getCommandQueue()->removeAll(this);

            SuperLatticePhysVelocity3D<T, DESCRIPTOR> velocity(lattice, converter);

            auto* volume = convertSimpleVolumeToVoreenVolume();

            app->getCommandQueue()->enqueue(this, LambdaFunctionCallback([&] {
                debugVelocityPort_.setData(volume, true);
            }));
        }
         */

        // Store last velocity.
        const int outputIter = maxIteration / config.getNumTimeSteps();
        if (iteration % outputIter == 0) {
            if (VoreenApplication * app = VoreenApplication::app()) {
                app->getCommandQueue()->removeAll(this);

                SuperLatticePhysVelocity3D<T, DESCRIPTOR> velocity(lattice, converter);

                SimpleVolume<float> velocityVolume = sampleVolume<float>(extendedDomain, converter,
                                                                         config.getOutputResolution(), velocity);

                auto dimensions = tgt::ivec3::fromPointer(velocityVolume.dimensions.data());
                auto spacing = tgt::Vector3<T>::fromPointer(velocityVolume.spacing.data());
                auto offset = tgt::Vector3<T>::fromPointer(velocityVolume.offset.data());

                auto* representation = new VolumeRAM_3xFloat(dimensions);
                for (size_t i = 0, numVoxels = tgt::hmul(dimensions); i < numVoxels; i++) {
                    for (int channel = 0; channel < velocityVolume.numChannels; channel++) {
                        representation->setVoxelNormalized(velocityVolume.getValue(i, channel), i, channel);
                    }
                }

                auto* volume = new Volume(representation, spacing, offset);
                volume->setRealWorldMapping(RealWorldMapping());

                app->getCommandQueue()->enqueue(this, LambdaFunctionCallback([&] {
                    debugVelocityPort_.setData(volume, true);
                }));
            }
        }

        if(!success) {
            break;
        }

        // === 8th Step: Check for convergence.
        converge.takeValue(lattice.getStatistics().getAverageEnergy(), true);
        if(converge.hasConverged()) {
            LINFO("Simulation converged!");
            break;
        }

        float progress = iteration / (maxIteration + 1.0f);
        progressReporter.setProgress(progress);
    }
    progressReporter.setProgress(1.0f);
    LINFO("Finished simulation run: " << parameters.name_);
}

}   // namespace
