#include "shared/simulation_core.h"
#include "shared/openlb_parameters.cpp" // We treat the definition file as header, we anyway only include it once.

////////// Globals //////////////////
// Meta
const std::string simulation = "simulation_cluster";
const std::string base = "/scratch/tmp/s_leis06/simulations/";

enum ExitCodes {
    EXIT_CODE_SUCCESS = EXIT_SUCCESS,
    EXIT_CODE_FAILURE = EXIT_FAILURE,
    EXIT_CODE_FINISHED = 6,
    EXIT_CODE_CONVERGED = 7,
    EXIT_CODE_DIVERGED = 8,
};


template<typename T>
T convertTo(const std::string& string);

template<>
std::string convertTo(const std::string& string) {
    return string;
}

template<>
float convertTo(const std::string& string) {
    return std::atof(string.c_str());
}

template<typename R, typename S>
std::map<R, S> deserializeMap(const XMLreader& items) {
    std::map<R, S> map;

    auto iter = items.begin();
    while(iter != items.end()) {

        // Read key.
        XMLreader* keyItem = *iter;
        if(keyItem->getName() != "key") {
            std::cout << "Expected key, aborting..." << std::endl;
            return {};
        }

        auto key = convertTo<R>(keyItem->getAttribute("value"));

        // Go to next entry (which is expected to be the value for the key).
        if(++iter == items.end()) {
            std::cout << "No matching value for key" << std::endl;
            return {};
        }

        // Read value.
        XMLreader* valueItem = *iter;
        if(valueItem->getName() != "value") {
            std::cout << "Expected value, aborting..." << std::endl;
            return {};
        }

        auto value = convertTo<S>(valueItem->getAttribute("value"));

        // Add key-value pair.
        map[key] = value;

        // Next key-value pair.
        iter++;
    }

    return map;
}

VelocityCurve deserializeVelocityCurve(const XMLreader& reader) {

    auto values = deserializeMap<float, float>(reader["peakVelocities"]);
    auto curve = VelocityCurve::createFromMap(values);

    bool periodic = reader["periodic"].getAttribute("value") == "true";
    curve.setPeriodic(periodic);

    float scale = std::atof(reader["scale"].getAttribute("value").c_str());
    curve.setScale(scale);

    return curve;
}

Parameters deserializeParameters(const XMLreader& reader) {
    Parameters parameters;
    parameters.name_ = reader["name"].getAttribute("value");
    parameters.spatialResolution_ = std::atoi(reader["spatialResolution"].getAttribute("value").c_str());
    parameters.relaxationTime_ = std::atof(reader["relaxationTime"].getAttribute("value").c_str());
    parameters.characteristicLength_ = std::atof(reader["characteristicLength"].getAttribute("value").c_str());
    parameters.characteristicVelocity_ = std::atof(reader["characteristicVelocity"].getAttribute("value").c_str());
    parameters.viscosity_ = std::atof(reader["viscosity"].getAttribute("value").c_str());
    parameters.density_ = std::atof(reader["density"].getAttribute("value").c_str());
    parameters.turbulenceModel_ = static_cast<FlowTurbulenceModel>(std::atof(reader["turbulenceModel"].getAttribute("value").c_str()));
    parameters.smagorinskyConstant_ = std::atof(reader["smagorinskyConstant"].getAttribute("value").c_str());
    parameters.wallBoundaryCondition_ = static_cast<FlowBoundaryCondition>(std::atof(reader["wallBoundaryCondition"].getAttribute("value").c_str())); // Replaces: bouzidiOn = parameters["bouzidi"].getAttribute("value") == "true";
    parameters.inletVelocityMultiplier_ = std::atof(reader["inletVelocityMultiplier"].getAttribute("value").c_str());
    parameters.latticePerturbation_ = reader["latticePerturbation"].getAttribute("value") == "true";
    return parameters;
}

std::vector<FlowIndicator> deserializeFlowIndicators(const XMLreader& reader) {
    std::vector<FlowIndicator> indicators;
    for(auto iter : reader) {
        FlowIndicator indicator;
        indicator.type_                 = static_cast<FlowIndicatorType>(std::atoi((*iter)["type_"].getAttribute("value").c_str()));
        indicator.id_                   = std::atoi((*iter)["id_"].getAttribute("value").c_str());
        indicator.center_.x             = std::atof((*iter)["center"].getAttribute("x").c_str());
        indicator.center_.y             = std::atof((*iter)["center"].getAttribute("y").c_str());
        indicator.center_.z             = std::atof((*iter)["center"].getAttribute("z").c_str());
        indicator.normal_.x             = std::atof((*iter)["normal"].getAttribute("x").c_str());
        indicator.normal_.y             = std::atof((*iter)["normal"].getAttribute("y").c_str());
        indicator.normal_.z             = std::atof((*iter)["normal"].getAttribute("z").c_str());
        indicator.radius_               = std::atof((*iter)["radius"].getAttribute("value").c_str());
        indicator.flowProfile_          = static_cast<FlowProfile>(std::atoi((*iter)["flowProfile"].getAttribute("value").c_str()));
        indicator.velocityCurve_        = deserializeVelocityCurve((*iter)["velocityCurve"]);
        indicators.push_back(indicator);
    }
    return indicators;
}

int main(int argc, char* argv[]) {

    // === 1st Step: Initialization ===
    olbInit(&argc, &argv);

    if(argc != 4) {
        std::cout << "Invalid number of arguments! Usage:" << std::endl;
        std::cout << "./" << simulation << " <ensemble name> <run name> <output directory>" << std::endl;
        return EXIT_FAILURE;
    }

    //std::string simulation = argv[0];
    std::string ensemble = argv[1];
    std::string run = argv[2];

    //std::string output = base; // hardcoded path
    std::string output = argv[3];
    int rank = 0;
#ifdef PARALLEL_MODE_MPI
    rank = singleton::mpi().getRank();
#endif
    if (rank == 0) {
        __mode_t mode = ACCESSPERMS;
        output += simulation + "/";
        mkdir(output.c_str(), mode); // ignore result
        output += ensemble + "/";
        mkdir(output.c_str(), mode); // ignore result
        output += run + "/";
        struct stat statbuf;
        if (stat(output.c_str(), &statbuf) != 0 && mkdir(output.c_str(), mode) != 0) {
            std::cout << "Could not create output directory: '" << output << "'" << std::endl;
            return EXIT_FAILURE;
        }
    }
    else {
        output += simulation + "/";
        output += ensemble + "/";
        output += run + "/";
    }
    singleton::directories().setOutputDir(output);

    OstreamManager clout(std::cout, "main");
    clout.setMultiOutput(false); // don't display messages from every single mpi process.
    clout << "Set output directory: " << output << std::endl;

    clout << "Running: " << simulation << std::endl;
    clout << "Ensemble: " << ensemble << std::endl;
    clout << "Run: " << run << std::endl;

    // Parse XML simulation config.
    XMLreader config("config.xml");
    T simulationTime             = std::atof(config["simulationTime"].getAttribute("value").c_str());
    int numTimeSteps             = std::atoi(config["numTimeSteps"].getAttribute("value").c_str());
    int outputResolution         = std::atoi(config["outputResolution"].getAttribute("value").c_str());
    std::string outputFileFormat =           config["outputFileFormat"].getAttribute("value");
    int flowFeatures             = std::atoi(config["flowFeatures"].getAttribute("value").c_str());

    auto geometryFiles           = deserializeMap<float, std::string>(config["geometryFiles"]);
    bool geometryIsMesh          = config["geometryIsMesh"].getAttribute("value") == "true";
    auto measuredDataFiles       = deserializeMap<float, std::string>(config["measuredDataFiles"]);

    auto parameters = deserializeParameters(config["flowParameters"]);
    auto indicators = deserializeFlowIndicators(config["transformedFlowIndicators"]);

    clout << "Found " << indicators.size() << " Flow Indicators" << std::endl;

    UnitConverterFromResolutionAndRelaxationTime<T, DESCRIPTOR> converter(
                parameters.spatialResolution_,      // Resolution that charPhysLength is resolved by.
            (T) parameters.relaxationTime_,         // Relaxation time
            (T) parameters.characteristicLength_,   // charPhysLength: reference length of simulation geometry
            (T) parameters.characteristicVelocity_, // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
            (T) parameters.viscosity_,              // physViscosity: physical kinematic viscosity in __m^2 / s__
            (T) parameters.density_                 // physDensity: physical density in __kg / m^3__
    );
    // Prints the converter log as console output
    converter.print();
    // Writes the converter log in a file
    converter.write(simulation.c_str());

    clout << "Loading measured data..." << std::endl;
    VolumeTimeSeries measuredDataTimeSeries(measuredDataFiles);

    // === 2nd Step: Prepare Geometry ===
    clout << "Meshing..." << std::endl;
    olb::util::Timer<T> voxelizationTime(1);
    voxelizationTime.start();

    // Instantiation of the STLreader class
    // file name, voxel size in meter, stl unit in meter, outer voxel no., inner voxel no.

    // TODO: for now, we only support a single geometry.
    std::unique_ptr<IndicatorF3D<T>> boundaryGeometry;
    std::unique_ptr<VolumeTimeSeries> geometryVolumeTimeSeries;
    if(geometryIsMesh) {
        std::string geometryFileName = geometryFiles.begin()->second;
        boundaryGeometry.reset(new STLreader<T>(geometryFileName, converter.getConversionFactorLength(), VOREEN_LENGTH_TO_SI, 1));
    }
    else {
        geometryVolumeTimeSeries.reset(new VolumeTimeSeries(geometryFiles));
        boundaryGeometry.reset(new VolumeDataMapperIndicator(geometryVolumeTimeSeries->createSampler(0.0f)));
    }
    IndicatorLayer3D<T> extendedDomain(*boundaryGeometry, converter.getConversionFactorLength());

    // Instantiation of a cuboidGeometry with weights
#ifdef PARALLEL_MODE_MPI
    const int noOfCuboids = std::min( 16*parameters.spatialResolution_, singleton::mpi().getSize() );
#else
    const int noOfCuboids = 1;
#endif
    CuboidGeometry3D<T> cuboidGeometry(extendedDomain, converter.getConversionFactorLength(), noOfCuboids);

    // Instantiation of a loadBalancer
    HeuristicLoadBalancer<T> loadBalancer(cuboidGeometry);

    // Instantiation of a superGeometry
    SuperGeometry<T,3> superGeometry(cuboidGeometry, loadBalancer, 2);
    bool success = prepareGeometry(converter, extendedDomain, *boundaryGeometry, superGeometry, indicators);

    voxelizationTime.stop();
    voxelizationTime.printSummary();

    if(!success) {
        clout << "The model contains errors! Check resolution and geometry." << std::endl;
        return EXIT_FAILURE;
    }

    // === 3rd Step: Prepare Lattice ===
    clout << "Lattice Initialization..." << std::endl;
    SuperLattice<T, DESCRIPTOR> lattice(superGeometry);

    olb::util::Timer<T> timer(converter.getLatticeTime(simulationTime), superGeometry.getStatistics().getNvoxel());
    timer.start();

    prepareLattice(lattice, converter, *boundaryGeometry, superGeometry, indicators, parameters);

    timer.stop();
    timer.printSummary();

    // === 4th Step: Main Loop with Timer ===
    clout << "Starting simulation..." << std::endl;
    util::ValueTracer<T> converge(converter.getLatticeTime(0.5), 1e-5);
    timer.start();

    const int maxIteration = converter.getLatticeTime(simulationTime);
    auto checkpoint = [&] (int iteration, bool enforce=false) {
        return getResults(
                lattice,
                superGeometry,
                *boundaryGeometry,
                converter,
                iteration,
                maxIteration,
                singleton::directories().getLogOutDir(),
                numTimeSteps,
                outputResolution,
                outputFileFormat,
                flowFeatures,
                parameters,
                enforce
        );
    };

    ExitCodes exitCode = EXIT_CODE_FINISHED;
    for (int iteration = 0; iteration <= maxIteration; iteration++) {

        // === 5th Step: Definition of Initial and Boundary Conditions ===
        setBoundaryValues(lattice, converter, iteration, superGeometry, indicators, parameters);

        // === 6th Step: Collide and Stream Execution ===
        lattice.collideAndStream();

        // Print some update
        if(iteration % 10 == 0) {
            timer.print(iteration);
        }

        // === 7th Step: Computation and Output of the Results ===
        bool abort = !checkpoint(iteration);
        if(abort) {
            clout << "Simulation diverged!" << std::endl;
            exitCode = ExitCodes::EXIT_CODE_DIVERGED;
        }

        // === 8th Step: Check for convergence.
        if(!abort) {
            converge.takeValue(lattice.getStatistics().getAverageEnergy(), true);
            if (converge.hasConverged()) {
                clout << "Simulation converged!" << std::endl;
                exitCode = ExitCodes::EXIT_CODE_CONVERGED;
            }
        }

        if(abort) {
            timer.update(iteration);
            checkpoint(iteration, true);
            break;
        }

        // === 9th Step: Write checkpoint
        // TODO: implement!
        //sLattice.save("simulation.checkpoint");
    }

    timer.stop();
    timer.printSummary();

    return exitCode;
}
