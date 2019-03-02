/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2011-2014 Mathias J. Krause
 *  E-mail contact: info@openlb.net
 *  The most recent release of OpenLB can be downloaded at
 *  <http://www.openlb.net/>
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2
 *  of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public
 *  License along with this program; if not, write to the Free
 *  Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
 *  Boston, MA  02110-1301, USA.
 */


#include "olb3D.h"
#ifndef OLB_PRECOMPILED // Unless precompiled version is used,
#include "olb3D.hh"   // include full template code;
#endif

#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <dirent.h>

using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;
using namespace olb::util;

typedef double T;
#define DESCRIPTOR D3Q19Descriptor

enum FlowDirection {
    FD_NONE = -1,
    FD_IN   =  0,
    FD_OUT  =  1,
};

enum FlowFunction {
    FF_NONE     = -1,
    FF_CONSTANT =  0,
    FF_SINUS    =  1,
};

// Indicates flux through an arbitrary, circle-shaped area.
// This code is adapted from the voreen host code.
struct FlowIndicator {
    FlowDirection   direction_{FD_NONE};
    FlowFunction    function_{FF_NONE};
    T               center_[3]{0.0};
    T               normal_[3]{0.0};
    T               radius_{0.0};
    int             materialId_{0};
};

////////// Globals //////////////////
// Meta
static const std::string simulation = "default";
static const std::string base = "/scratch/tmp/s_leis06/simulations/";
static const T VOREEN_LENGTH_TO_SI = 0.001;
static const T VOREEN_TIME_TO_SI = 0.001;

// Config
T simulationTime = 0.0;
T temporalResolution = 0.0;
int spatialResolution = 1;
std::vector<FlowIndicator> flowIndicators;

// Parameters
T characteristicLength = 0.0;
T characteristicVelocity = 0.0;
T viscosity = 0.0;
T density = 0.0;
bool bouzidiOn = false;
//////////////////////////////////////


// Stores data from stl file in geometry in form of material numbers
void prepareGeometry(UnitConverter<T, DESCRIPTOR> const& converter, IndicatorF3D<T>& indicator,
                     STLreader<T>& stlReader, SuperGeometry3D<T>& superGeometry) {

    OstreamManager clout(std::cout, "prepareGeometry");
    clout << "Prepare Geometry ..." << std::endl;

    superGeometry.rename(0, 2, indicator);
    superGeometry.rename(2, 1, stlReader);

    superGeometry.clean();

    int materialId = 3; // 0=empty, 1=liquid, 2=walls

    for (size_t i = 0; i < flowIndicators.size(); i++) {

        const T* center = flowIndicators[i].center_;
        const T* normal = flowIndicators[i].normal_;
        T radius = flowIndicators[i].radius_;

        // Set material number for inflow
        IndicatorCircle3D<T> flow(center[0]*VOREEN_LENGTH_TO_SI, center[1]*VOREEN_LENGTH_TO_SI, center[2]*VOREEN_LENGTH_TO_SI,
                                  normal[0], normal[1], normal[2],
                                  radius*VOREEN_LENGTH_TO_SI);
        IndicatorCylinder3D<T> layerFlow(flow, 2. * converter.getConversionFactorLength());
        superGeometry.rename(2, materialId, 1, layerFlow);
        flowIndicators[i].materialId_ = materialId;
        materialId++;
    }

    // Removes all not needed boundary voxels outside the surface
    superGeometry.clean();
    // Removes all not needed boundary voxels inside the surface
    superGeometry.innerClean(3);
    superGeometry.checkForErrors();

    superGeometry.print();
    clout << "Prepare Geometry ... OK" << std::endl;
}

// Set up the geometry of the simulation
void prepareLattice(SuperLattice3D<T, DESCRIPTOR>& lattice,
                    UnitConverter<T, DESCRIPTOR> const& converter, Dynamics<T, DESCRIPTOR>& bulkDynamics,
                    sOnLatticeBoundaryCondition3D<T, DESCRIPTOR>& bc,
                    sOffLatticeBoundaryCondition3D<T, DESCRIPTOR>& offBc,
                    STLreader<T>& stlReader, SuperGeometry3D<T>& superGeometry) {

    OstreamManager clout(std::cout, "prepareLattice");
    clout << "Prepare Lattice ..." << std::endl;

    const T omega = converter.getLatticeRelaxationFrequency();

    // material=0 --> do nothing
    lattice.defineDynamics(superGeometry, 0, &instances::getNoDynamics<T, DESCRIPTOR>());

    // material=1 --> bulk dynamics
    lattice.defineDynamics(superGeometry, 1, &bulkDynamics);

    if (bouzidiOn) {
        // material=2 --> no dynamics + bouzidi zero velocity
        lattice.defineDynamics(superGeometry, 2, &instances::getNoDynamics<T, DESCRIPTOR>());
        offBc.addZeroVelocityBoundary(superGeometry, 2, stlReader);
    } else {
        // material=2 --> bounceBack dynamics
        lattice.defineDynamics(superGeometry, 2, &instances::getBounceBack<T, DESCRIPTOR>());
    }

    for(const FlowIndicator& indicator : flowIndicators) {
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

    // Initial conditions
    AnalyticalConst3D<T, T> rhoF(1);
    std::vector<T> velocity(3, T());
    AnalyticalConst3D<T, T> uF(velocity);

    lattice.defineRhoU( superGeometry,1,rhoF,uF );
    lattice.iniEquilibrium( superGeometry,1,rhoF,uF );

    // Initialize all values of distribution functions to their local equilibrium
    for(const FlowIndicator& indicator : flowIndicators) {
        lattice.defineRhoU(superGeometry, indicator.materialId_, rhoF, uF);
        lattice.iniEquilibrium(superGeometry, indicator.materialId_, rhoF, uF);
    }

    // Lattice initialize
    lattice.initialize();

    clout << "Prepare Lattice ... OK" << std::endl;
}

// Generates a slowly increasing sinuidal inflow
void setBoundaryValues(SuperLattice3D<T, DESCRIPTOR>& sLattice,
                       sOffLatticeBoundaryCondition3D<T, DESCRIPTOR>& offBc,
                       UnitConverter<T, DESCRIPTOR> const& converter, int iT,
                       SuperGeometry3D<T>& superGeometry) {

    // No of time steps for smooth start-up
    int iTperiod = converter.getLatticeTime(0.5);
    int iTupdate = 50;

    if (iT % iTupdate == 0) {
        for(const FlowIndicator& indicator : flowIndicators) {
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
void getResults(SuperLattice3D<T, DESCRIPTOR>& sLattice,
                UnitConverter<T, DESCRIPTOR>& converter, int iT,
                Dynamics<T, DESCRIPTOR>& bulkDynamics,
                SuperGeometry3D<T>& superGeometry, Timer<T>& timer, STLreader<T>& stlReader) {

    OstreamManager clout(std::cout, "getResults");

    SuperVTMwriter3D<T> vtmWriter(simulation);
    SuperLatticePhysVelocity3D<T, DESCRIPTOR> velocity(sLattice, converter);
    SuperLatticePhysPressure3D<T, DESCRIPTOR> pressure(sLattice, converter);
    vtmWriter.addFunctor(velocity);
    vtmWriter.addFunctor(pressure);

    const int vtkIter = converter.getLatticeTime(.1);
    const int statIter = converter.getLatticeTime(.1);

    int rank = 0;
#ifdef PARALLEL_MODE_MPI
    //rank = singleton::mpi().getRank();
#endif
    if (rank == 0 && iT % vtkIter == 0) {
        const Vector<T, 3>& min = stlReader.getMin();
        const Vector<T, 3>& max = stlReader.getMax();

        const int resolution = converter.getResolution();
        const Vector<T, 3> len = max - min;//converter.getCharPhysLength();
        std::vector<float> rawVelocityData;
        rawVelocityData.reserve(static_cast<size_t>(resolution * resolution * resolution * 3));

        AnalyticalFfromSuperF3D<T> interpolateVelocity(velocity, true);
        for(int z=0; z<resolution; z++) {
            for(int y=0; y<resolution; y++) {
                for(int x=0; x<resolution; x++) {

                    T pos[3] = {min[0]+x*len[0]/resolution, min[1]+y*len[1]/resolution, min[2]+z*len[2]/resolution};
                    T u[3] = {0.0, 0.0, 0.0};

                    if(pos[0] >= min[0] && pos[1] >= min[1] && pos[2] >= min[2] &&
                       pos[0] <= max[0] && pos[1] <= max[1] && pos[2] <= max[2])
                        interpolateVelocity(u, pos);

                    // Downgrade to float.
                    rawVelocityData.push_back(static_cast<float>(u[0]/VOREEN_LENGTH_TO_SI));
                    rawVelocityData.push_back(static_cast<float>(u[1]/VOREEN_LENGTH_TO_SI));
                    rawVelocityData.push_back(static_cast<float>(u[2]/VOREEN_LENGTH_TO_SI));
                }
            }
        }

        std::string velocityFilename = singleton::directories().getVtkOutDir() + "/velocity_" + std::to_string(iT) + ".raw";
        std::fstream velocityFile(velocityFilename.c_str(), std::ios::out | std::ios::binary);
        size_t numBytes = rawVelocityData.size() * sizeof(float) / sizeof(char);
        velocityFile.write(reinterpret_cast<const char*>(rawVelocityData.data()), numBytes);
        if (!velocityFile.good()) {
            clout << "Could not write velocity file" << std::endl;
        }
    }

    // Writes output on the console
    if (iT % statIter == 0) {
        // Lattice statistics console output
        sLattice.getStatistics().print(iT, converter.getPhysTime(iT));
    }

    if (sLattice.getStatistics().getMaxU() > 0.3) {
        clout << "PROBLEM uMax=" << sLattice.getStatistics().getMaxU() << std::endl;
        std::exit(0);
    }
}

int main(int argc, char* argv[]) {

    // === 1st Step: Initialization ===
    olbInit(&argc, &argv);
    OstreamManager clout(std::cout, "main");
    // don't display messages from every single mpi process
    clout.setMultiOutput(false);

    if(argc != 3) {
        clout << "Invalid number of arguments!" << std::endl;
        return EXIT_FAILURE;
    }

    //std::string simulation = argv[0];
    std::string ensemble = argv[1];
    std::string run = argv[2];

    clout << "Running: " << simulation << std::endl;
    clout << "Ensemble:" << ensemble << std::endl;
    clout << "Run: " << run << std::endl;

    __mode_t mode = ACCESSPERMS;
    std::string output = base;
    if (DIR* dir = opendir(output.c_str())) {
        closedir(dir);
    } else {
        output += simulation + "/";
        mkdir(output.c_str(), mode); // ignore result
        output += ensemble + "/";
        mkdir(output.c_str(), mode); // ignore result
        output += run + "/";
        if (mkdir(output.c_str(), mode) != 0) {
            clout << "Could not create output directory: '" << output << "'" << std::endl;
            return EXIT_FAILURE;
        }
    }
    singleton::directories().setOutputDir(output.c_str());

    XMLreader config("config.xml");
    simulationTime = std::atof(config["simulationTime"].getAttribute("value").c_str());
    temporalResolution = std::atof(config["temporalResolution"].getAttribute("value").c_str());
    spatialResolution = std::atoi(config["spatialResolution"].getAttribute("value").c_str());

    XMLreader parameters = config["flowParameters"];
    characteristicLength = std::atof(parameters["characteristicLength"].getAttribute("value").c_str());
    characteristicVelocity = std::atof(parameters["characteristicVelocity"].getAttribute("value").c_str());
    viscosity = std::atof(parameters["viscosity"].getAttribute("value").c_str());
    density = std::atof(parameters["density"].getAttribute("value").c_str());
    bouzidiOn = parameters["bouzidi"].getAttribute("value") == "true";

    XMLreader indicators = config["flowIndicators"];
    for(auto iter : indicators) {
        FlowIndicator indicator;
        indicator.direction_ = static_cast<FlowDirection>(std::atoi((*iter)["direction"].getAttribute("value").c_str()));
        indicator.function_  = static_cast<FlowFunction>(std::atoi((*iter)["function"].getAttribute("value").c_str()));
        indicator.center_[0] = std::atof((*iter)["center"].getAttribute("x").c_str());
        indicator.center_[1] = std::atof((*iter)["center"].getAttribute("y").c_str());
        indicator.center_[2] = std::atof((*iter)["center"].getAttribute("z").c_str());
        indicator.normal_[0] = std::atof((*iter)["normal"].getAttribute("x").c_str());
        indicator.normal_[1] = std::atof((*iter)["normal"].getAttribute("y").c_str());
        indicator.normal_[2] = std::atof((*iter)["normal"].getAttribute("z").c_str());
        indicator.radius_ = std::atof((*iter)["radius"].getAttribute("value").c_str());
        flowIndicators.push_back(indicator);
    }
    std::cout << "Found " << flowIndicators.size() << " Flow Indicators";

    const int N = spatialResolution;
    UnitConverter<T, DESCRIPTOR> converter(
            (T) characteristicLength * VOREEN_LENGTH_TO_SI / N,  // physDeltaX: spacing between two lattice cells in __m__
            (T) temporalResolution * VOREEN_TIME_TO_SI,          // physDeltaT: time step in __s__
            (T) characteristicLength * VOREEN_LENGTH_TO_SI,      // charPhysLength: reference length of simulation geometry
            (T) characteristicVelocity * VOREEN_LENGTH_TO_SI,    // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
            (T) viscosity * 1e-6,                                // physViscosity: physical kinematic viscosity in __m^2 / s__
            (T) density                                          // physDensity: physical density in __kg / m^3__
    );
    // Prints the converter log as console output
    converter.print();
    // Writes the converter log in a file
    converter.write(simulation.c_str());

    // === 2nd Step: Prepare Geometry ===

    // Instantiation of the STLreader class
    // file name, voxel size in meter, stl unit in meter, outer voxel no., inner voxel no.
    std::string geometryFileName = "../geometry/geometry.stl";
    STLreader<T> stlReader(geometryFileName.c_str(), converter.getConversionFactorLength(), VOREEN_LENGTH_TO_SI, 0, true);
    IndicatorLayer3D<T> extendedDomain(stlReader, converter.getConversionFactorLength());

    // Instantiation of a cuboidGeometry with weights
#ifdef PARALLEL_MODE_MPI
    const int noOfCuboids = std::min( 16*N,2*singleton::mpi().getSize() );
#else
    const int noOfCuboids = 2;
#endif
    CuboidGeometry3D<T> cuboidGeometry(extendedDomain, converter.getConversionFactorLength(), noOfCuboids);

    // Instantiation of a loadBalancer
    HeuristicLoadBalancer<T> loadBalancer(cuboidGeometry);

    // Instantiation of a superGeometry
    SuperGeometry3D<T> superGeometry(cuboidGeometry, loadBalancer, 2);

    prepareGeometry(converter, extendedDomain, stlReader, superGeometry);

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

    Timer<T> timer1(converter.getLatticeTime(simulationTime), superGeometry.getStatistics().getNvoxel());
    timer1.start();

    prepareLattice(sLattice, converter, bulkDynamics,
                   sBoundaryCondition, sOffBoundaryCondition,
                   stlReader, superGeometry);

    timer1.stop();
    timer1.printSummary();

    // === 4th Step: Main Loop with Timer ===
    clout << "starting simulation..." << std::endl;
    util::ValueTracer<T> converge(converter.getLatticeTime(1.0), 1e-5);
    Timer<T> timer(converter.getLatticeTime(simulationTime), superGeometry.getStatistics().getNvoxel());
    timer.start();

    for (int iT = 0; iT <= converter.getLatticeTime(simulationTime); iT++) {

        // === 5th Step: Definition of Initial and Boundary Conditions ===
        setBoundaryValues(sLattice, sOffBoundaryCondition, converter, iT, superGeometry);

        // === 6th Step: Collide and Stream Execution ===
        sLattice.collideAndStream();

        // === 7th Step: Computation and Output of the Results ===
        getResults(sLattice, converter, iT, bulkDynamics, superGeometry, timer, stlReader);

        // === 8th Step: Check for convergence.
        converge.takeValue(sLattice.getStatistics().getAverageEnergy(), true);
        if(converge.hasConverged()) {
            clout << "Simulation converged!";
            break;
        }
    }

    timer.stop();
    timer.printSummary();
}
