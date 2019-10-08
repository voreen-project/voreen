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

#include "voreen/core/datastructures/volume/volumeatomic.h"

#include "flowsimulation.h"

#include "voreen/core/datastructures/geometry/glmeshgeometry.h"

#include "voreen/core/datastructures/volume/volumeminmaxmagnitude.h"
#include "voreen/core/ports/conditions/portconditionvolumelist.h"

#ifndef OLB_PRECOMPILED
#include "olb3D.hh"
#endif

namespace voreen {

const T FlowSimulation::VOREEN_LENGTH_TO_SI = 0.001;
const T FlowSimulation::VOREEN_TIME_TO_SI = 0.001;
const std::string FlowSimulation::loggerCat_("voreen.flowsimulation.FlowSimulation");


FlowSimulation::MeasuredDataMapper::MeasuredDataMapper(const VolumeBase* volume)
    : AnalyticalF3D<T, T>(3)
    , volume_(volume)
{
    tgtAssert(volume_, "No volume");
    tgtAssert(volume_->getNumChannels() == 3, "Num channels != 3");
    bounds_ = volume_->getBoundingBox(false).getBoundingBox(false);
    representation_.reset(new VolumeRAMRepresentationLock(volume_));
    physicalToVoxelMatrix_ = volume_->getPhysicalToVoxelMatrix();
}

bool FlowSimulation::MeasuredDataMapper::operator() (T output[], const T input[]) {
    tgt::vec3 rwPos = tgt::Vector3<T>::fromPointer(input) / VOREEN_LENGTH_TO_SI;
    if (!bounds_.containsPoint(rwPos)) {
        return false;
    }

    auto* initialState = dynamic_cast<const VolumeRAM_3xFloat*>(**representation_);
    if(!initialState)
        return false;

    tgt::vec3 voxel = initialState->getVoxelLinear(physicalToVoxelMatrix_ * rwPos, 0, true);
    for(size_t i=0; i < initialState->getNumChannels(); i++) {
        output[i] = voxel[i] * VOREEN_LENGTH_TO_SI;
    }

    return true;
}


FlowSimulation::FlowSimulation()
    : AsyncComputeProcessor<ComputeInput, ComputeOutput>()
    , geometryDataPort_(Port::INPORT, "geometryDataPort", "Geometry Input", false)
    , measuredDataPort_(Port::INPORT, "measuredDataPort", "Measured Data Input", false)
    , parameterPort_(Port::INPORT, "parameterPort", "Parameterization", false)
    , simulationResults_("simulationResults", "Simulation Results", "Simulation Results", VoreenApplication::app()->getTemporaryPath("simulation"), "", FileDialogProperty::DIRECTORY, Processor::VALID, Property::LOD_DEFAULT, VoreenFileWatchListener::ALWAYS_OFF)
    , deleteOldSimulations_("deleteOldSimulations", "Delete old Simulations", false)
    , simulateAllParametrizations_("simulateAllParametrizations", "Simulate all Parametrizations", false)
    , selectedParametrization_("selectedSimulation", "Selected Parametrization", 0, 0, 0)
{
    addPort(geometryDataPort_);
    addPort(measuredDataPort_);
    //measuredDataPort_.addCondition(new PortConditionVolumeListEnsemble());
    measuredDataPort_.addCondition(new PortConditionVolumeListAdapter(new PortConditionVolumeType3xFloat()));
    addPort(parameterPort_);

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

    // Note: measuredDataPort ist optional!
    if(measuredDataPort_.hasData() && !measuredDataPort_.isReady()) {
        setNotReadyErrorMessage("Measured Data Port not ready.");
        return false;
    }

    if(!parameterPort_.isReady()) {
        setNotReadyErrorMessage("Parameter Port not ready.");
        return false;
    }

    return true;
}

void FlowSimulation::adjustPropertiesToInput() {
    const FlowParametrizationList* flowParameterList = parameterPort_.getData();
    if(!flowParameterList || flowParameterList->empty()) {
        selectedParametrization_.setMinValue(-1);
        selectedParametrization_.setMaxValue(-1);
    }
    else {
        selectedParametrization_.setMinValue(0);
        selectedParametrization_.setMaxValue(static_cast<int>(flowParameterList->size())-1);
        selectedParametrization_.set(0);
    }
}

FlowSimulationInput FlowSimulation::prepareComputeInput() {
    const GlMeshGeometryBase* geometryData = dynamic_cast<const GlMeshGeometryBase*>(geometryDataPort_.getData());
    if (!geometryData) {
        throw InvalidInputException("Invalid simulation geometry", InvalidInputException::S_WARNING);
    }

    tgtAssert(measuredDataPort_.isDataInvalidationObservable(), "VolumeListPort must be DataInvalidationObservable!");
    const VolumeList* measuredData = measuredDataPort_.getThreadSafeData();

    tgtAssert(parameterPort_.isDataInvalidationObservable(), "FlowParametrizationPort must be DataInvalidationObservable!");
    const FlowParametrizationList* flowParameterList = parameterPort_.getThreadSafeData();
    if(!flowParameterList || flowParameterList->empty()) {
        throw InvalidInputException("No parameterization", InvalidInputException::S_ERROR);
    }

    if(flowParameterList->getFlowFeatures() == FF_NONE) {
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
        measuredData = nullptr;
        LINFO("Configuring an unsteered simulation");
    }

    std::string geometryPath = VoreenApplication::app()->getUniqueTmpFilePath(".stl");
    try {
        std::ofstream file(geometryPath);
        geometryData->exportAsStl(file);
        file.close();
    }
    catch (std::exception&) {
        throw InvalidInputException("Geometry could not be exported", InvalidInputException::S_ERROR);
    }

    if(simulationResults_.get().empty()) {
        throw InvalidInputException("No output directory selected", InvalidInputException::S_WARNING);
    }

    std::string simulationPath = simulationResults_.get() + "/" + flowParameterList->getName() + "/";
    if (!tgt::FileSystem::createDirectoryRecursive(simulationPath)) {
        throw InvalidInputException("Output directory could not be created", InvalidInputException::S_ERROR);
    }

    size_t selectedParametrization = FlowParametrizationList::ALL_PARAMETRIZATIONS;
    if(!simulateAllParametrizations_.get()) {
        selectedParametrization = static_cast<size_t>(selectedParametrization_.get());
    }

    return FlowSimulationInput{
            geometryPath,
            measuredData,
            flowParameterList,
            selectedParametrization,
            simulationPath,
            deleteOldSimulations_.get()
    };
}

FlowSimulationOutput FlowSimulation::compute(FlowSimulationInput input, ProgressReporter& progressReporter) const {

    // Needs to be initialized in each new thread to be used.
    olb::olbInit(nullptr, nullptr);

    // Run either all or just a single simulation.
    if(input.selectedParametrization == FlowParametrizationList::ALL_PARAMETRIZATIONS) {
        progressReporter.setProgress(0.0f);
        size_t numRuns = input.parametrizationList->size();
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
    // Nothing to do.
}

void FlowSimulation::runSimulation(const FlowSimulationInput& input,
                                   ProgressReporter& progressReporter) const {

    const VolumeList* measuredData = input.measuredData;
    const FlowParametrizationList& parametrizationList = *input.parametrizationList;
    const FlowParameters& parameters = parametrizationList.at(input.selectedParametrization);

    LINFO("Starting simulation run: " << parameters.getName());
    progressReporter.setProgress(0.0f);

    std::string simulationResultPath = input.simulationResultPath + parameters.getName() + "/";
    if (input.deleteOldSimulations && tgt::FileSystem::dirExists(simulationResultPath)) {
        tgt::FileSystem::deleteDirectoryRecursive(simulationResultPath);
    }
    if (!tgt::FileSystem::createDirectory(simulationResultPath)) {
        LERROR("Output directory could not be created. It may already exist.");
        return;
    }

    const int N = parameters.getSpatialResolution();
    UnitConverter<T, DESCRIPTOR> converter(
            (T) parameters.getCharacteristicLength() * VOREEN_LENGTH_TO_SI / N,     // resolution for charPhysLength
            (T) parameters.getTemporalResolution() * VOREEN_TIME_TO_SI,             // relaxation time
            (T) parameters.getCharacteristicLength() * VOREEN_LENGTH_TO_SI,         // charPhysLength: reference length of simulation geometry
            (T) parameters.getCharacteristicVelocity() * VOREEN_LENGTH_TO_SI,       // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
            (T) parameters.getViscosity() * 0.001 / parameters.getDensity(),        // physViscosity: physical kinematic viscosity in __m^2 / s__
            (T) parameters.getDensity()                                             // physDensity: physical density in __kg / m^3__
    );
    // Prints the converter log as console output
    converter.print();

    // === 2nd Step: Prepare Geometry ===

    // Instantiation of the STLreader class
    // file name, voxel size in meter, stl unit in meter, outer voxel no., inner voxel no.
    STLreader<T> stlReader(input.geometryPath, converter.getConversionFactorLength(), VOREEN_LENGTH_TO_SI, 1);
    IndicatorLayer3D<T> extendedDomain(stlReader, converter.getConversionFactorLength());

    interruptionPoint();

    std::vector<FlowIndicatorMaterial> flowIndicators;
    for(const FlowIndicator& indicator : parametrizationList.getFlowIndicators()) {
        FlowIndicatorMaterial indicatorMaterial;
        indicatorMaterial.direction_            = indicator.direction_;
        indicatorMaterial.flowProfile_          = indicator.flowProfile_;
        indicatorMaterial.startPhaseFunction_   = indicator.startPhaseFunction_;
        indicatorMaterial.startPhaseDuration_   = indicator.startPhaseDuration_;
        indicatorMaterial.center_               = indicator.center_;
        indicatorMaterial.normal_               = indicator.normal_;
        indicatorMaterial.radius_               = indicator.radius_;
        indicatorMaterial.materialId_           = MAT_EMPTY;
        flowIndicators.push_back(indicatorMaterial);
    }

    // Instantiation of a cuboidGeometry with weights
    const int noOfCuboids = 1;
    CuboidGeometry3D<T> cuboidGeometry(extendedDomain, converter.getConversionFactorLength(), noOfCuboids);

    // Instantiation of a loadBalancer
    HeuristicLoadBalancer<T> loadBalancer(cuboidGeometry);

    // Instantiation of a superGeometry
    SuperGeometry3D<T> superGeometry(cuboidGeometry, loadBalancer, 2);

    interruptionPoint();

    prepareGeometry(converter, extendedDomain, stlReader, superGeometry,
                    parametrizationList, input.selectedParametrization,
                    flowIndicators);

    interruptionPoint();

    // === 3rd Step: Prepare Lattice ===
    SuperLattice3D<T, DESCRIPTOR> sLattice(superGeometry);

    SmagorinskyBGKdynamics<T, DESCRIPTOR> bulkDynamics(converter.getLatticeRelaxationFrequency(),
                                                       instances::getBulkMomenta<T, DESCRIPTOR>(),
                                                       parameters.getSmagorinskyConstant());

    // choose between local and non-local boundary condition
    sOnLatticeBoundaryCondition3D<T, DESCRIPTOR> sBoundaryCondition(sLattice);
    createInterpBoundaryCondition3D<T, DESCRIPTOR>(sBoundaryCondition);
    // createLocalBoundaryCondition3D<T,DESCRIPTOR>(sBoundaryCondition);

    sOffLatticeBoundaryCondition3D<T, DESCRIPTOR> sOffBoundaryCondition(sLattice);
    createBouzidiBoundaryCondition3D<T, DESCRIPTOR>(sOffBoundaryCondition);

    interruptionPoint();

    prepareLattice(sLattice, converter, bulkDynamics,
                   sBoundaryCondition, sOffBoundaryCondition,
                   stlReader, superGeometry,
                   measuredData,
                   parametrizationList, input.selectedParametrization,
                   flowIndicators);

    interruptionPoint();

    // === 4th Step: Main Loop  ===
    const int tmax = converter.getLatticeTime(parametrizationList.getSimulationTime());
    util::ValueTracer<T> converge( converter.getLatticeTime(0.5), 1e-5);
    for (int ti = 0; ti <= tmax; ti++) {

        // === 5th Step: Definition of Initial and Boundary Conditions ===
        setBoundaryValues(sLattice, sOffBoundaryCondition, converter, ti, superGeometry,
                          parametrizationList, input.selectedParametrization, flowIndicators);

        // === 6th Step: Collide and Stream Execution ===
        sLattice.collideAndStream();

        // === 7th Step: Computation and Output of the Results ===
        bool success = getResults(sLattice, converter, ti, tmax, bulkDynamics, superGeometry, stlReader,
                                  parametrizationList, input.selectedParametrization, flowIndicators,
                                  simulationResultPath);
        if(!success) {
            break;
        }

        // Check for convergence.
        converge.takeValue(sLattice.getStatistics().getAverageEnergy(), true);
        if(converge.hasConverged()) {
            LINFO("Simulation converged!");
            break;
        }

        float progress = ti / (tmax + 1.0f);
        progressReporter.setProgress(progress);
    }
    progressReporter.setProgress(1.0f);
    LINFO("Finished simulation run: " << parameters.getName());
}

// Stores data from stl file in geometry in form of material numbers
void FlowSimulation::prepareGeometry( UnitConverter<T,DESCRIPTOR> const& converter, IndicatorF3D<T>& indicator,
                                      STLreader<T>& stlReader, SuperGeometry3D<T>& superGeometry,
                                      const FlowParametrizationList& parametrizationList,
                                      size_t selectedParametrization,
                                      std::vector<FlowIndicatorMaterial>& flowIndicators) const {

    LINFO("Prepare Geometry ...");

    superGeometry.rename( MAT_EMPTY, MAT_WALL,  indicator );
    superGeometry.rename( MAT_WALL,  MAT_FLUID, stlReader );

    superGeometry.clean();

    int materialId = MAT_COUNT;

    for (size_t i = 0; i < flowIndicators.size(); i++) {

        const tgt::vec3& center = flowIndicators[i].center_;
        const tgt::vec3& normal = flowIndicators[i].normal_;
        T radius = flowIndicators[i].radius_;

        // Set material number for inflow
        IndicatorCircle3D<T> flow(center[0]*VOREEN_LENGTH_TO_SI, center[1]*VOREEN_LENGTH_TO_SI, center[2]*VOREEN_LENGTH_TO_SI,
                                  normal[0], normal[1], normal[2], radius*VOREEN_LENGTH_TO_SI);
        IndicatorCylinder3D<T> layerFlow(flow, 2. * converter.getConversionFactorLength());
        superGeometry.rename(MAT_WALL, materialId, MAT_FLUID, layerFlow);
        flowIndicators[i].materialId_ = materialId;
        materialId++;
    }

    // Removes all not needed boundary voxels outside the surface
    superGeometry.clean();
    // Removes all not needed boundary voxels inside the surface
    superGeometry.innerClean(MAT_COUNT);
    superGeometry.checkForErrors();

    LINFO("Prepare Geometry ... OK");
}

// Set up the geometry of the simulation
void FlowSimulation::prepareLattice( SuperLattice3D<T, DESCRIPTOR>& lattice,
                                     UnitConverter<T,DESCRIPTOR> const& converter,
                                     Dynamics<T, DESCRIPTOR>& bulkDynamics,
                                     sOnLatticeBoundaryCondition3D<T, DESCRIPTOR>& bc,
                                     sOffLatticeBoundaryCondition3D<T,DESCRIPTOR>& offBc,
                                     STLreader<T>& stlReader, SuperGeometry3D<T>& superGeometry,
                                     const VolumeList* measuredData,
                                     const FlowParametrizationList& parametrizationList,
                                     size_t selectedParametrization,
                                     std::vector<FlowIndicatorMaterial>& flowIndicators) const {

    LINFO("Prepare Lattice ...");

    const bool bouzidiOn = parametrizationList.at(selectedParametrization).getBouzidi();
    const T omega = converter.getLatticeRelaxationFrequency();

    // material=0 --> do nothing
    lattice.defineDynamics(superGeometry, MAT_EMPTY, &instances::getNoDynamics<T, DESCRIPTOR>());

    // material=1 --> bulk dynamics
    lattice.defineDynamics(superGeometry, MAT_FLUID, &bulkDynamics);

    if (bouzidiOn) {
        // material=2 --> no dynamics + bouzidi zero velocity
        lattice.defineDynamics(superGeometry, MAT_WALL, &instances::getNoDynamics<T, DESCRIPTOR>());
        offBc.addZeroVelocityBoundary(superGeometry, MAT_WALL, stlReader);
    } else {
        // material=2 --> bounceBack dynamics
        lattice.defineDynamics(superGeometry, MAT_WALL, &instances::getBounceBack<T, DESCRIPTOR>());
    }

    for(const FlowIndicatorMaterial& indicator : flowIndicators) {
        if(indicator.direction_ == FD_IN) {
            if(bouzidiOn) {
                // no dynamics + bouzidi velocity (inflow)
                lattice.defineDynamics(superGeometry, indicator.materialId_, &instances::getNoDynamics<T, DESCRIPTOR>());
                offBc.addVelocityBoundary(superGeometry, indicator.materialId_, stlReader);
            }
            else {
                // bulk dynamics + velocity (inflow)
                lattice.defineDynamics(superGeometry, indicator.materialId_, &bulkDynamics);
                bc.addVelocityBoundary(superGeometry, indicator.materialId_, omega);
            }
        }
        else if(indicator.direction_ == FD_OUT) {
            lattice.defineDynamics(superGeometry, indicator.materialId_, &bulkDynamics);
            bc.addPressureBoundary(superGeometry, indicator.materialId_, omega);
        }
    }

    // Unsteered simulation.
    if(!measuredData) {

        // Initial conditions
        AnalyticalConst3D<T, T> rhoF(1);
        std::vector<T> velocity(3, T());
        AnalyticalConst3D<T, T> uF(velocity);

        lattice.defineRhoU(superGeometry, MAT_FLUID, rhoF, uF);
        lattice.iniEquilibrium(superGeometry, MAT_FLUID, rhoF, uF);

        // Initialize all values of distribution functions to their local equilibrium
        for (const FlowIndicatorMaterial& indicator : flowIndicators) {
            lattice.defineRhoU(superGeometry, indicator.materialId_, rhoF, uF);
            lattice.iniEquilibrium(superGeometry, indicator.materialId_, rhoF, uF);
        }
    }
    // Steered simulation - currently only initializes the first time step!
    else {
        MeasuredDataMapper mapper(measuredData->first());
        lattice.defineU(superGeometry, MAT_FLUID, mapper);
    }

    // Lattice initialize
    lattice.initialize();

    LINFO("Prepare Lattice ... OK");
}

// Generates a slowly increasing sinuidal inflow
void FlowSimulation::setBoundaryValues( SuperLattice3D<T, DESCRIPTOR>& sLattice,
                                        sOffLatticeBoundaryCondition3D<T,DESCRIPTOR>& offBc,
                                        UnitConverter<T,DESCRIPTOR> const& converter, int iT,
                                        SuperGeometry3D<T>& superGeometry,
                                        const FlowParametrizationList& parametrizationList,
                                        size_t selectedParametrization,
                                        std::vector<FlowIndicatorMaterial>& flowIndicators) const {
    // No of time steps for smooth start-up
    int iTupdate = 50;
    bool bouzidiOn = parametrizationList.at(selectedParametrization).getBouzidi();

    if (iT % iTupdate == 0) {
        for(const FlowIndicatorMaterial& indicator : flowIndicators) {
            if (indicator.direction_ == FD_IN) {

                int iTvec[1] = {iT};
                T maxVelocity[1] = {T()};

                switch(indicator.startPhaseFunction_) {
                case FSP_SINUS:
                {
                    int iTperiod = converter.getLatticeTime(indicator.startPhaseDuration_);
                    if(iT < iTperiod) {
                        SinusStartScale<T, int> nSinusStartScale(iTperiod, converter.getCharLatticeVelocity());
                        nSinusStartScale(maxVelocity, iTvec);
                        break;
                    }
                    // Else: fallthrough
                }
                case FSP_CONSTANT:
                {
                    AnalyticalConst1D<T, int> nConstantStartScale(converter.getCharLatticeVelocity());
                    nConstantStartScale(maxVelocity, iTvec);
                    break;
                }
                case FSP_NONE:
                default:
                    // Skip!
                    continue;
                }

                // This function applies the velocity profile to the boundary condition and the lattice.
                auto applyFlowProfile = [&] (AnalyticalF3D<T,T>& profile) {
                    if (bouzidiOn) {
                        offBc.defineU(superGeometry, indicator.materialId_, profile);
                    } else {
                        sLattice.defineU(superGeometry, indicator.materialId_, profile);
                    }
                };

                switch(indicator.flowProfile_) {
                case FP_POISEUILLE:
                {
                    //*
                    // TODO: The following line tends to crash, if there is no overlap between the region and the wall.
                    CirclePoiseuille3D<T> profile(superGeometry, indicator.materialId_, maxVelocity[0]);
                    /*/
                    // TODO: the following lines, however, lead to a strange behavior at in- and outlets.
                    const tgt::vec3& center = indicator.center_;
                    const tgt::vec3& normal = indicator.normal_;
                    T radius = indicator.radius_;
                    CirclePoiseuille3D<T> profile(center[0]*VOREEN_LENGTH_TO_SI, center[1]*VOREEN_LENGTH_TO_SI, center[2]*VOREEN_LENGTH_TO_SI,
                                                  normal[0], normal[1], normal[2], radius * VOREEN_LENGTH_TO_SI, maxVelocity[0]);
                    //*/
                    applyFlowProfile(profile);
                    break;
                }
                case FP_POWERLAW:
                {
                    T n = 1.03 * std::log(converter.getReynoldsNumber()) - 3.6; // Taken from OLB documentation.
                    CirclePowerLawTurbulent3D<T> profile(superGeometry, indicator.materialId_, maxVelocity[0], n);
                    applyFlowProfile(profile);
                    break;
                }
                case FP_CONSTANT:
                {
                    AnalyticalConst3D<T, T> profile(maxVelocity[0]);
                    applyFlowProfile(profile);
                    break;
                }
                case FP_NONE:
                default:
                    // Skip!
                    continue;
                }
            }
        }
    }
}

// Computes flux at inflow and outflow
bool FlowSimulation::getResults( SuperLattice3D<T, DESCRIPTOR>& sLattice,
                                 UnitConverter<T,DESCRIPTOR>& converter, int ti, int tmax,
                                 Dynamics<T, DESCRIPTOR>& bulkDynamics,
                                 SuperGeometry3D<T>& superGeometry,
                                 STLreader<T>& stlReader,
                                 const FlowParametrizationList& parametrizationList,
                                 size_t selectedParametrization,
                                 std::vector<FlowIndicatorMaterial>& flowIndicators,
                                 const std::string& simulationOutputPath) const {

    const int outputIter = tmax / parametrizationList.getNumTimeSteps();

    if ( ti == 0 ) {
        SuperVTMwriter3D<T> vtmWriter( "debug" );
        SuperLatticeGeometry3D<T, DESCRIPTOR> geometry( sLattice, superGeometry );
        SuperLatticeCuboid3D<T, DESCRIPTOR> cuboid( sLattice );
        SuperLatticeRank3D<T, DESCRIPTOR> rank( sLattice );
        vtmWriter.write( geometry );
        vtmWriter.write( cuboid );
        vtmWriter.write( rank );
        vtmWriter.createMasterFile();
    }

    if (ti % outputIter == 0) {

        if(parametrizationList.getFlowFeatures() & FF_VELOCITY) {
            SuperLatticePhysVelocity3D<T, DESCRIPTOR> velocity(sLattice, converter);
            writeResult(stlReader, converter, ti, tmax, velocity, parametrizationList, selectedParametrization,
                        simulationOutputPath, "velocity");
        }

        if(parametrizationList.getFlowFeatures() & FF_MAGNITUDE) {
            SuperLatticePhysVelocity3D<T, DESCRIPTOR> velocity(sLattice, converter);
            SuperEuklidNorm3D<T, DESCRIPTOR> magnitude(velocity);
            writeResult(stlReader, converter, ti, tmax, magnitude, parametrizationList, selectedParametrization,
                        simulationOutputPath, "magnitude");
        }

        if(parametrizationList.getFlowFeatures() & FF_PRESSURE) {
            SuperLatticePhysPressure3D<T, DESCRIPTOR> pressure(sLattice, converter);
            writeResult(stlReader, converter, ti, tmax, pressure, parametrizationList, selectedParametrization,
                        simulationOutputPath, "pressure");
        }

#ifndef OLB_PRECOMPILED
        if(parametrizationList.getFlowFeatures() & FF_WALLSHEARSTRESS) {
            SuperLatticePhysWallShearStress3D<T, DESCRIPTOR> wallShearStress(sLattice, superGeometry, MAT_WALL,
                                                                             converter, stlReader);
            writeResult(stlReader, converter, ti, tmax, wallShearStress, parametrizationList,
                        selectedParametrization,
                        simulationOutputPath, "wallShearStress");
        }
#endif

        // Lattice statistics console output
        LINFO("step="     << ti << "; " <<
              "t="        << converter.getPhysTime(ti) << "; " <<
              "uMax="     << sLattice.getStatistics().getMaxU() << "; " <<
              "avEnergy=" << sLattice.getStatistics().getAverageEnergy() << "; " <<
              "avRho="    << sLattice.getStatistics().getAverageRho()
              );
        sLattice.getStatistics().print(ti, converter.getPhysTime(ti));
    }

    T tau = converter.getLatticeRelaxationFrequency();
    T threshold = tau < 0.55 ? 0.125*(tau - 0.5) : 0.4;
    if (sLattice.getStatistics().getMaxU() >= threshold) {
        LERROR("uMax=" << sLattice.getStatistics().getMaxU() << " above threshold=" << threshold);
        return false;
    }

    return true;
}

void FlowSimulation::writeResult(STLreader<T>& stlReader,
                                 UnitConverter<T,DESCRIPTOR>& converter, int ti, int tmax,
                                 SuperLatticeF3D<T, DESCRIPTOR>& feature,
                                 const FlowParametrizationList& parametrizationList,
                                 size_t selectedParametrization,
                                 const std::string& simulationOutputPath,
                                 const std::string& name) const {

    const Vector<T, 3>& min = stlReader.getMin();
    const Vector<T, 3>& max = stlReader.getMax();

    const Vector<T, 3> len = (max - min);
    const T maxLen = std::max({len[0], len[1], len[2]});
    const int gridResolution = static_cast<int>(std::round(maxLen / converter.getConversionFactorLength()));
    const int resolution = std::min(parametrizationList.getOutputResolution(), gridResolution);

    Vector<T, 3> offset = min + (len - maxLen) * 0.5;
    Vector<T, 3> spacing(maxLen / (resolution-1));

    // Determine format.
    // This could be done in a more dynamic way, but the code should be easily portable to the cluster.
    std::string format;
    switch (feature.getTargetDim()) {
        case 1:
            format = "float";
            break;
        case 3:
            format = "Vector3(float)";
            break;
        default:
            LERROR("Unhandled target dimensions");
            return;
    }

    std::vector<float> rawFeatureData;
    rawFeatureData.reserve(static_cast<size_t>(resolution * resolution * resolution * feature.getTargetDim()));
    AnalyticalFfromSuperF3D<T> interpolateFeature(feature, true);

    std::vector<T> minValue(feature.getTargetDim(), std::numeric_limits<T>::max());
    std::vector<T> maxValue(feature.getTargetDim(), std::numeric_limits<T>::lowest());

    T minMagnitude = std::numeric_limits<T>::max();
    T maxMagnitude = 0;

    for (int z = 0; z < resolution; z++) {
        for (int y = 0; y < resolution; y++) {
            for (int x = 0; x < resolution; x++) {

                T pos[3] = {offset[0] + x * spacing[0], offset[1] + y * spacing[1], offset[2] + z * spacing[2]};
                std::vector<T> val(feature.getTargetDim(), 0.0f);

                if (pos[0] >= min[0] && pos[1] >= min[1] && pos[2] >= min[2] &&
                    pos[0] <= max[0] && pos[1] <= max[1] && pos[2] <= max[2]) {
                    interpolateFeature(&val[0], pos);

                    // Update min/max.
                    T magnitude = 0;
                    for (int i = 0; i < feature.getTargetDim(); i++) {
                        minValue[i] = std::min(minValue[i], val[i]);
                        maxValue[i] = std::max(maxValue[i], val[i]);
                        magnitude += val[i] * val[i];
                    }
                    minMagnitude = std::min(minMagnitude, magnitude);
                    maxMagnitude = std::max(maxMagnitude, magnitude);
                }

                // Downgrade to float.
                for (int i = 0; i < feature.getTargetDim(); i++) {
                    rawFeatureData.push_back(static_cast<float>(val[i] / VOREEN_LENGTH_TO_SI));
                }
            }
        }
    }

    // Adapt to Voreen units.
    offset  = offset  * (1/VOREEN_LENGTH_TO_SI);
    spacing = spacing * (1/VOREEN_LENGTH_TO_SI);

    for (int i = 0; i < feature.getTargetDim(); i++) {
        minValue[i] /= VOREEN_LENGTH_TO_SI;
        maxValue[i] /= VOREEN_LENGTH_TO_SI;
    }

    minMagnitude = std::sqrt(minMagnitude) / VOREEN_LENGTH_TO_SI;
    maxMagnitude = std::sqrt(maxMagnitude) / VOREEN_LENGTH_TO_SI;

    // Set output names.
    int tmaxLen = static_cast<int>(std::to_string(tmax).length());
    std::ostringstream suffix;
    suffix << std::setw(tmaxLen) << std::setfill('0') << ti;
    std::string featureFilename = name + "_" + suffix.str();
    std::string rawFilename = simulationOutputPath + featureFilename + ".raw";
    std::string vvdFilename = simulationOutputPath + featureFilename + ".vvd";

    const FlowParameters& parameters = parametrizationList.at(selectedParametrization);
    const LatticeStatistics<T>& statistics = feature.getSuperLattice().getStatistics();
    std::fstream vvdFeatureFile(vvdFilename.c_str(), std::ios::out);
    vvdFeatureFile
            // Header.
            << "<?xml version=\"1.0\" ?>"
            << "<VoreenData version=\"1\">"
            << "<Volumes>"
            << "<Volume>"
            // Data.
            << "<RawData filename=\"" << featureFilename << ".raw\" format=\"" << format << "\" x=\"" << resolution
            << "\" y=\"" << resolution << "\" z=\"" << resolution << "\" />"
            // Mandatory Meta data.
            << "<MetaData>"
            << "<MetaItem name=\"" << VolumeBase::META_DATA_NAME_OFFSET << "\" type=\"Vec3MetaData\">"
            << "<value x=\"" << offset[0] << "\" y=\"" << offset[1] << "\" z=\"" << offset[2] << "\" />"
            << "</MetaItem>"
            << "<MetaItem name=\"" << VolumeBase::META_DATA_NAME_SPACING << "\" type=\"Vec3MetaData\">"
            << "<value x=\"" << spacing[0] << "\" y=\"" << spacing[1] << "\" z=\"" << spacing[2] << "\" />"
            << "</MetaItem>"
            << "<MetaItem name=\"" << VolumeBase::META_DATA_NAME_TIMESTEP << "\" type=\"FloatMetaData\" value=\"" << converter.getPhysTime(ti) << "\" />"
            << "<MetaItem name=\"" << VolumeBase::META_DATA_NAME_REAL_WORLD_MAPPING << "\" type=\"RealWorldMappingMetaData\"><value scale=\"1\" offset=\"0\" unit=\"\" /></MetaItem>"
            << "<MetaItem name=\"" << "name" << "\" type=\"StringMetaData\" value=\"" << name << "\" />"
            // Parameters.
            << "<MetaItem name=\"" << "ParameterSpatialResolution" << "\" type=\"IntMetaData\" value=\"" << parameters.getSpatialResolution() << "\" />"
            << "<MetaItem name=\"" << "ParameterTemporalResolution" << "\" type=\"FloatMetaData\" value=\"" << parameters.getTemporalResolution() << "\" />"
            << "<MetaItem name=\"" << "ParameterCharacteristicLength" << "\" type=\"FloatMetaData\" value=\"" << parameters.getCharacteristicLength() << "\" />"
            << "<MetaItem name=\"" << "ParameterCharacteristicVelocity" << "\" type=\"FloatMetaData\" value=\"" << parameters.getCharacteristicVelocity() << "\" />"
            << "<MetaItem name=\"" << "ParameterViscosity" << "\" type=\"FloatMetaData\" value=\"" << parameters.getViscosity() << "\" />"
            << "<MetaItem name=\"" << "ParameterDensity" << "\" type=\"FloatMetaData\" value=\"" << parameters.getDensity() << "\" />"
            << "<MetaItem name=\"" << "ParameterSmagorinskyConstant" << "\" type=\"FloatMetaData\" value=\"" << parameters.getSmagorinskyConstant() << "\" />"
            << "<MetaItem name=\"" << "ParameterBouzidi" << "\" type=\"BoolMetaData\" value=\"" << (parameters.getBouzidi() ? "true" : "false") << "\" />"
            // Additional meta data.
            << "<MetaItem name=\"" << "StatisticsMaxVelocity" << "\" type=\"FloatMetaData\" value=\"" << statistics.getMaxU() << "\" />"
            << "<MetaItem name=\"" << "StatisticsAvgEnergy" << "\" type=\"FloatMetaData\" value=\"" << statistics.getAverageEnergy() << "\" />"
            << "<MetaItem name=\"" << "StatisticsMaxRho" << "\" type=\"FloatMetaData\" value=\"" << statistics.getAverageRho() << "\" />"
            << "</MetaData>"
            // Derived data.
            << "<DerivedData>";
            // * VolumeMinMaxMagnitude
    if (feature.getTargetDim() > 1) {
        vvdFeatureFile
            << "<DerivedItem type=\"VolumeMinMaxMagnitude\" minMagnitude=\"" << minMagnitude << "\" maxMagnitude=\"" << maxMagnitude << "\" minNormalizedMagnitude=\"" << minMagnitude << "\" maxNormalizedMagnitude=\"" << maxMagnitude << "\" />";
    }
            // * VolumeMinMax
    vvdFeatureFile << "<DerivedItem type=\"VolumeMinMax\"><minValues>";
    for(int i=0; i<feature.getTargetDim(); i++) {
        vvdFeatureFile << "<channel value=\"" << minValue[i] << "\" />";
    }
    vvdFeatureFile << "</minValues><maxValues>";
    for(int i=0; i<feature.getTargetDim(); i++) {
        vvdFeatureFile << "<channel value=\"" << maxValue[i] << "\" />";
    }
    vvdFeatureFile << "</maxValues><minNormValues>";
    for(int i=0; i<feature.getTargetDim(); i++) {
        vvdFeatureFile << "<channel value=\"" << minValue[i] << "\" />";
    }
    vvdFeatureFile << "</minNormValues><maxNormValues>";
    for(int i=0; i<feature.getTargetDim(); i++) {
        vvdFeatureFile << "<channel value=\"" << maxValue[i] << "\" />";
    }
    vvdFeatureFile
            << "</maxNormValues></DerivedItem>"
            << "</DerivedData>"
            // Footer.
            << "</Volume>"
            << "</Volumes>"
            << "</VoreenData>";

    vvdFeatureFile.close();

    std::fstream rawFeatureFile(rawFilename.c_str(), std::ios::out | std::ios::binary);
    size_t numBytes = rawFeatureData.size() * sizeof(float) / sizeof(char);
    rawFeatureFile.write(reinterpret_cast<const char*>(rawFeatureData.data()), numBytes);
    if (!rawFeatureFile.good()) {
        LERROR("Could not write " << name << " file");
    }
}

}   // namespace
