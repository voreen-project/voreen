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

#include "voreen/core/datastructures/geometry/glmeshgeometry.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/datastructures/volume/volumefactory.h"
#include "voreen/core/datastructures/volume/volumeminmaxmagnitude.h"
#include "voreen/core/ports/conditions/portconditionvolumetype.h"

#include "../../utils/geometryconverter.h"

#include "modules/hdf5/io/hdf5volumereader.h"
#include "modules/hdf5/io/hdf5volumewriter.h"

namespace voreen {

const T FlowSimulation::VOREEN_LENGTH_TO_SI = 0.001;
const T FlowSimulation::VOREEN_TIME_TO_SI = 0.001;
const std::string FlowSimulation::loggerCat_("voreen.flowreen.FlowSimulation");


FlowSimulation::MeasuredDataMapper::MeasuredDataMapper(const VolumeBase* volume)
    : AnalyticalF3D<T, T>(3)
    , volume_(volume)
{
    tgtAssert(volume_, "No volume");
    tgtAssert(volume_->getNumChannels() == 3, "Num channels != 3");
    bounds_ = volume_->getBoundingBox(true).getBoundingBox(true);
}

bool FlowSimulation::MeasuredDataMapper::operator() (T output[], const T input[]) {
    tgt::vec3 rwPos = tgt::Vector3<T>::fromPointer(input) / VOREEN_LENGTH_TO_SI;
    if (!bounds_.containsPoint(rwPos)) {
        return false;
    }

    const VolumeRAM_3xFloat* initialState = dynamic_cast<const VolumeRAM_3xFloat*>(volume_->getRepresentation<VolumeRAM>());
    if(!initialState)
        return false;

    for(size_t i=0; i < volume_->getNumChannels(); i++) {
        tgt::vec3 voxel = initialState->getVoxelLinear(rwPos, 0, false);
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
    , simulateAllParametrizations_("simulateAllParametrizations", "Simulate all Parametrizations?", false)
    , selectedParametrization_("selectedSimulation", "Selected Parametrization", 0, 0, 0)
{
    addPort(geometryDataPort_);
    addPort(measuredDataPort_);
    measuredDataPort_.addCondition(new PortConditionVolumeList(new PortConditionVolumeType3xFloat()));
    addPort(parameterPort_);

    addProperty(simulationResults_);

    //addProperty(simulateAllParametrizations_);
    ON_CHANGE_LAMBDA(simulateAllParametrizations_, [this]{
        selectedParametrization_.setReadOnlyFlag(simulateAllParametrizations_.get());
    });
    addProperty(selectedParametrization_);
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
    auto geometryData = geometryDataPort_.getThreadSafeData();
    if (!geometryData) {
        throw InvalidInputException("No simulation geometry", InvalidInputException::S_WARNING);
    }

    tgtAssert(measuredDataPort_.isDataInvalidationObservable(), "VolumeListPort must be DataInvalidationObservable!");
    const VolumeList* measuredData = measuredDataPort_.getThreadSafeData();

    tgtAssert(parameterPort_.isDataInvalidationObservable(), "FlowParametrizationPort must be DataInvalidationObservable!");
    const FlowParametrizationList* flowParameterList = parameterPort_.getThreadSafeData();
    if(!flowParameterList || flowParameterList->empty()) {
        throw InvalidInputException("No parameterization", InvalidInputException::S_ERROR);
    }

    if(measuredData && !measuredData->empty()) {
        LINFO("Configuring a steered simulation");
        // Check for volume compatibility
        VolumeBase* volumeT0 = measuredData->first();
        if (volumeT0->getDimensions() != tgt::svec3(volumeT0->getDimensions().x)) {
            throw InvalidInputException("Measured data must have dimensions: n x n x n",
                                        InvalidInputException::S_ERROR);
        }
        if (volumeT0->getSpacing() != tgt::vec3(volumeT0->getSpacing().x)) {
            throw InvalidInputException("Measured data must have spacing: n x n x n", InvalidInputException::S_ERROR);
        }

        for (size_t i = 1; i < measuredData->size(); i++) {
            VolumeBase* volumeTi = measuredData->at(i);

            if (measuredData->at(i-1)->getTimestep() >= volumeTi->getTimestep()) {
                throw InvalidInputException("Time Steps of are not ordered", InvalidInputException::S_ERROR);
            }

            if (volumeT0->getDimensions() != volumeTi->getDimensions() ||
                volumeT0->getSpacing() != volumeTi->getSpacing()) {
                throw InvalidInputException("Measured data contains different kinds of volumes.",
                                            InvalidInputException::S_ERROR);
            }
        }
    }
    else {
        measuredData = nullptr;
        LINFO("Configuring an unsteered simulation");
    }

    std::string geometryPath = VoreenApplication::app()->getUniqueTmpFilePath(".stl");
    if(!exportGeometryToSTL(geometryData, geometryPath)) {
        throw InvalidInputException("Geometry could not be initialized", InvalidInputException::S_ERROR);
    }

    if(simulationResults_.get().empty()) {
        throw InvalidInputException("No output directory selected", InvalidInputException::S_WARNING);
    }

    size_t selectedParametrization = static_cast<size_t>(selectedParametrization_.get());
    std::string simulationPath = simulationResults_.get() + "/" + flowParameterList->getName() + "/";
    simulationPath += flowParameterList->at(selectedParametrization).getName() + "/";
    if(!tgt::FileSystem::createDirectoryRecursive(simulationPath)) {
        throw InvalidInputException("Output directory could not be created", InvalidInputException::S_ERROR);
    }

    return FlowSimulationInput{
            geometryPath,
            measuredData,
            flowParameterList,
            selectedParametrization,
            simulationPath
    };
}

FlowSimulationOutput FlowSimulation::compute(FlowSimulationInput input, ProgressReporter& progressReporter) const {

    const VolumeList* measuredData = input.measuredData;
    const FlowParametrizationList& parametrizationList = *input.parametrizationList;
    const FlowParameters& parameters = parametrizationList.at(input.selectedParametrization);

    progressReporter.setProgress(0.0f);

    // Needs to be initialized in each new thread to be used.
    olb::olbInit(nullptr, nullptr);

    const int N = parametrizationList.getSpatialResolution();
    UnitConverter<T, DESCRIPTOR> converter(
            (T) parameters.getCharacteristicLength() * VOREEN_LENGTH_TO_SI / N, // physDeltaX: spacing between two lattice cells in __m__
            (T) parametrizationList.getTemporalResolution() * VOREEN_TIME_TO_SI,// physDeltaT: time step in __s__
            (T) parameters.getCharacteristicLength() * VOREEN_LENGTH_TO_SI,     // charPhysLength: reference length of simulation geometry
            (T) parameters.getCharacteristicVelocity() * VOREEN_LENGTH_TO_SI,   // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
            (T) parameters.getViscosity() * 1e-6,                               // physViscosity: physical kinematic viscosity in __m^2 / s__
            (T) parameters.getDensity()                                         // physDensity: physical density in __kg / m^3__
    );
    // Prints the converter log as console output
    converter.print();

    // === 2nd Step: Prepare Geometry ===

    // Instantiation of the STLreader class
    // file name, voxel size in meter, stl unit in meter, outer voxel no., inner voxel no.
    STLreader<T> stlReader(input.geometryPath, converter.getConversionFactorLength(), VOREEN_LENGTH_TO_SI, 1);
    IndicatorLayer3D<T> extendedDomain(stlReader, converter.getConversionFactorLength());

    std::vector<FlowIndicatorMaterial> flowIndicators;
    for(const FlowIndicator& indicator : parametrizationList.getFlowIndicators()) {
        FlowIndicatorMaterial indicatorMaterial;
        indicatorMaterial.direction_ = indicator.direction_;
        indicatorMaterial.function_  = indicator.function_;
        indicatorMaterial.center_ = indicator.center_;
        indicatorMaterial.normal_ = indicator.normal_;
        indicatorMaterial.radius_ = indicator.radius_;
        indicatorMaterial.materialId_ = 0;
        flowIndicators.push_back(indicatorMaterial);
    }

    // Instantiation of a cuboidGeometry with weights
    const int noOfCuboids = 1;
    CuboidGeometry3D<T> cuboidGeometry(extendedDomain, converter.getConversionFactorLength(), noOfCuboids);

    // Instantiation of a loadBalancer
    HeuristicLoadBalancer<T> loadBalancer(cuboidGeometry);

    // Instantiation of a superGeometry
    SuperGeometry3D<T> superGeometry(cuboidGeometry, loadBalancer, 2);

    prepareGeometry(converter, extendedDomain, stlReader, superGeometry,
                    parametrizationList, input.selectedParametrization,
                    flowIndicators);

    // === 3rd Step: Prepare Lattice ===
    SuperLattice3D<T, DESCRIPTOR> sLattice(superGeometry);

    SmagorinskyBGKdynamics<T, DESCRIPTOR> bulkDynamics(converter.getLatticeRelaxationFrequency(),
                                                       instances::getBulkMomenta<T, DESCRIPTOR>(), 0.1);

    // choose between local and non-local boundary condition
    sOnLatticeBoundaryCondition3D<T, DESCRIPTOR> sBoundaryCondition(sLattice);
    createInterpBoundaryCondition3D<T, DESCRIPTOR>(sBoundaryCondition);
    // createLocalBoundaryCondition3D<T,DESCRIPTOR>(sBoundaryCondition);

    sOffLatticeBoundaryCondition3D<T, DESCRIPTOR> sOffBoundaryCondition(sLattice);
    createBouzidiBoundaryCondition3D<T, DESCRIPTOR>(sOffBoundaryCondition);

    prepareLattice(sLattice, converter, bulkDynamics,
                   sBoundaryCondition, sOffBoundaryCondition,
                   stlReader, superGeometry,
                   measuredData,
                   parametrizationList, input.selectedParametrization,
                   flowIndicators);

    // === 4th Step: Main Loop  ===
    util::ValueTracer<T> converge( converter.getLatticeTime(1.0), 1e-5 );
    LINFO("Starting simulation...");
    for (int iT = 0; iT <= converter.getLatticeTime(parametrizationList.getSimulationTime()); iT++) {

        // === 5th Step: Definition of Initial and Boundary Conditions ===
        setBoundaryValues(sLattice, sOffBoundaryCondition, converter, iT, superGeometry,
                          parametrizationList, input.selectedParametrization, flowIndicators);

        // === 6th Step: Collide and Stream Execution ===
        sLattice.collideAndStream();

        // === 7th Step: Computation and Output of the Results ===
        bool success = getResults(sLattice, converter, iT, bulkDynamics, superGeometry, stlReader,
                                  parametrizationList, input.selectedParametrization, flowIndicators,
                                  input.simulationResultPath);
        if(!success) {
            break;
        }

        // Check for convergence.
        converge.takeValue(sLattice.getStatistics().getAverageEnergy(), true);
        if(converge.hasConverged()) {
            //LINFO("Simulation converged!");
            //break;
        }

        float progress = iT / (converter.getLatticeTime( parametrizationList.getSimulationTime() ) + 1.0f);
        progressReporter.setProgress(progress);
    }
    progressReporter.setProgress(1.0f);
    LINFO("Finished simulation...");

    // Done.
    return FlowSimulationOutput{};
}

void FlowSimulation::processComputeOutput(FlowSimulationOutput output) {
}

// Stores data from stl file in geometry in form of material numbers
void FlowSimulation::prepareGeometry( UnitConverter<T,DESCRIPTOR> const& converter, IndicatorF3D<T>& indicator,
                                      STLreader<T>& stlReader, SuperGeometry3D<T>& superGeometry,
                                      const FlowParametrizationList& parametrizationList,
                                      size_t selectedParametrization,
                                      std::vector<FlowIndicatorMaterial>& flowIndicators) const {

    LINFO("Prepare Geometry ...");

    superGeometry.rename( MAT_EMPTY, MAT_WALL,  indicator );
    superGeometry.rename( MAT_WALL,  MAT_LIQUID, stlReader );

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
        superGeometry.rename(MAT_WALL, materialId, MAT_LIQUID, layerFlow);
        flowIndicators[i].materialId_ = materialId;
        materialId++;
    }

    // Removes all not needed boundary voxels outside the surface
    superGeometry.clean();
    // Removes all not needed boundary voxels inside the surface
    superGeometry.innerClean(3);
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
    lattice.defineDynamics(superGeometry, MAT_LIQUID, &bulkDynamics);

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
                // material=3 --> no dynamics + bouzidi velocity (inflow)
                lattice.defineDynamics(superGeometry, indicator.materialId_, &instances::getNoDynamics<T, DESCRIPTOR>());
                offBc.addVelocityBoundary(superGeometry, indicator.materialId_, stlReader);
            }
            else {
                // material=3 --> bulk dynamics + velocity (inflow)
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

        lattice.defineRhoU(superGeometry, MAT_LIQUID, rhoF, uF);
        lattice.iniEquilibrium(superGeometry, MAT_LIQUID, rhoF, uF);

        // Initialize all values of distribution functions to their local equilibrium
        for (const FlowIndicatorMaterial &indicator : flowIndicators) {
            lattice.defineRhoU(superGeometry, indicator.materialId_, rhoF, uF);
            lattice.iniEquilibrium(superGeometry, indicator.materialId_, rhoF, uF);
        }
    }
    // Steered simulation.
    else {
        MeasuredDataMapper mapper(measuredData->first());
        lattice.defineU(superGeometry, MAT_LIQUID, mapper);
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
    int iTperiod = converter.getLatticeTime(0.5);
    int iTupdate = 50;
    bool bouzidiOn = parametrizationList.at(selectedParametrization).getBouzidi();

    if (iT % iTupdate == 0) {
        for(const FlowIndicatorMaterial& indicator : flowIndicators) {
            if (indicator.direction_ == FD_IN) {

                int iTvec[1] = {iT};
                T maxVelocity[1] = {T()};

                switch(indicator.function_) {
                case FF_CONSTANT:
                {
                    AnalyticalConst1D<T, int> nConstantStartScale(converter.getCharLatticeVelocity());
                    nConstantStartScale(maxVelocity, iTvec);
                    break;
                }
                case FF_SINUS:
                {
                    SinusStartScale<T, int> nSinusStartScale(iTperiod, converter.getCharLatticeVelocity());
                    nSinusStartScale(maxVelocity, iTvec);
                    break;
                }
                case FF_NONE:
                default:
                    // Skip!
                    continue;
                }

                CirclePoiseuille3D<T> velocity(superGeometry, indicator.materialId_, maxVelocity[0]);
                if (bouzidiOn) {
                    offBc.defineU(superGeometry, indicator.materialId_, velocity);
                } else {
                    sLattice.defineU(superGeometry, indicator.materialId_, velocity);
                }
            }
        }
    }
}

// Computes flux at inflow and outflow
bool FlowSimulation::getResults( SuperLattice3D<T, DESCRIPTOR>& sLattice,
                                 UnitConverter<T,DESCRIPTOR>& converter, int iT,
                                 Dynamics<T, DESCRIPTOR>& bulkDynamics,
                                 SuperGeometry3D<T>& superGeometry,
                                 STLreader<T>& stlReader,
                                 const FlowParametrizationList& parametrizationList,
                                 size_t selectedParametrization,
                                 std::vector<FlowIndicatorMaterial>& flowIndicators,
                                 const std::string& simulationOutputPath) const {

    const int vtkIter = converter.getLatticeTime(.1);
    const int statIter = converter.getLatticeTime(.1);

    if (iT % vtkIter == 0) {
        // Write velocity.
        SuperLatticePhysVelocity3D<T, DESCRIPTOR> velocity(sLattice, converter);
        writeResult(stlReader, converter, iT, velocity, simulationOutputPath, "velocity");

        // Write pressure.
        SuperLatticePhysPressure3D<T, DESCRIPTOR> pressure(sLattice, converter);
        writeResult(stlReader, converter, iT, pressure, simulationOutputPath, "pressure");

        // Write WallShearStress
        // TODO: figure out how and estimate somehow from measured data!
        //SuperLatticePhysWallShearStress3D<T, DESCRIPTOR> wallShearStress(sLattice, superGeometry, MAT_WALL, converter, stlReader);
        //writeResult(stlReader, converter, iT, wallShearStress, simulationOutputPath, "wallShearStress");

        // Write Temperature.
        //SuperLatticePhysTemperature3D<T, DESCRIPTOR, TODO> temperature(sLattice, converter);
        //writeResult(stlReader, converter, iT, temperature, simulationOutputPath, "temperature");
    }

    // Writes output on the console
    if (iT % statIter == 0) {
        // Lattice statistics console output
        LINFO("step="     << iT << "; " <<
              "t="        << converter.getPhysTime(iT) << "; " <<
              "uMax="     << sLattice.getStatistics().getMaxU() << "; " <<
              "avEnergy=" << sLattice.getStatistics().getAverageEnergy() << "; " <<
              "avRho="    << sLattice.getStatistics().getAverageRho()
              );
        sLattice.getStatistics().print(iT, converter.getPhysTime(iT));
    }

    if (sLattice.getStatistics().getMaxU() > 0.3) {
        LERROR("PROBLEM uMax=" << sLattice.getStatistics().getMaxU());
        return false;
    }

    return true;
}

void FlowSimulation::writeResult(STLreader<T>& stlReader,
                                 UnitConverter<T,DESCRIPTOR>& converter, int iT,
                                 SuperLatticePhysF3D<T, DESCRIPTOR>& property,
                                 const std::string& simulationOutputPath,
                                 const std::string& name) const {

    const Vector<T, 3>& min = stlReader.getMin();
    const Vector<T, 3>& max = stlReader.getMax();

    const int resolution = converter.getResolution();
    const Vector<T, 3> len = (max - min);
    const T maxLen = std::max({len[0], len[1], len[2]});
    const Vector<T, 3> offset = min + (len - maxLen) * 0.5;

    std::vector<float> rawPropertyData;
    rawPropertyData.reserve(static_cast<size_t>(resolution * resolution * resolution * property.getTargetDim()));
    AnalyticalFfromSuperF3D<T> interpolateProperty(property, true);

    for(int z=0; z<resolution; z++) {
        for(int y=0; y<resolution; y++) {
            for(int x=0; x<resolution; x++) {

                T pos[3] = {offset[0]+x*maxLen/resolution, offset[1]+y*maxLen/resolution, offset[2]+z*maxLen/resolution};
                std::vector<T> dat(property.getTargetDim(), 0.0f);

                if(pos[0] >= min[0] && pos[1] >= min[1] && pos[2] >= min[2] &&
                   pos[0] <= max[0] && pos[1] <= max[1] && pos[2] <= max[2]) {
                    interpolateProperty(&dat[0], pos);
                }

                // Downgrade to float.
                for(int i = 0; i < property.getTargetDim(); i++) {
                    rawPropertyData.push_back(static_cast<float>(dat[i]/VOREEN_LENGTH_TO_SI));
                }
            }
        }
    }

    std::string propertyFilename = simulationOutputPath + name + "_" + std::to_string(iT) + ".raw"; // TODO: encode simulated time.
    std::fstream propertyFile(propertyFilename.c_str(), std::ios::out | std::ios::binary);
    size_t numBytes = rawPropertyData.size() * sizeof(float) / sizeof(char);
    propertyFile.write(reinterpret_cast<const char*>(rawPropertyData.data()), numBytes);
    if (!propertyFile.good()) {
        LERROR("Could not write " << name << " file");
    }
}

}   // namespace
