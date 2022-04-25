#include "simulation_core.h"


class VelocityCurve {
public:

    VelocityCurve() {
        peakVelocities_[0.0f] = 0.0f;
    }

    float operator()(float t) const {
        if(periodic_) {
            float begin = peakVelocities_.begin()->first;
            float end = peakVelocities_.rbegin()->first;
            t = std::fmod(t - begin, end - begin);
        }
        else {
            if(t < peakVelocities_.begin()->first) {
                return peakVelocities_.begin()->second;
            }

            if(t > peakVelocities_.rbegin()->first) {
                return peakVelocities_.rbegin()->second;
            }
        }

        struct Comparator {
            bool operator()(const std::pair<float, float>& p, float value) {
                return p.first < value;
            }
        };

        auto upper = std::lower_bound(peakVelocities_. begin(), peakVelocities_.end(), t, Comparator());
        auto lower = upper++;

        float a = (t - lower->first) / (upper->first - lower->first);

        return (1.0f - a) * lower->second + a * upper->second;
    }

    void deserialize(const XMLreader& reader) {
        periodic_ = reader["periodic"].getAttribute("value") == "true";

        peakVelocities_.clear();
        XMLreader items = reader["peakVelocities"];
        auto iter = items.begin();
        while(iter != items.end()) {

            // Read key.
            XMLreader* keyItem = *iter;
            if(keyItem->getName() != "key") {
                std::cout << "VelocityCurve: Expected key, aborting..." << std::endl;
                return;
            }

            float key = std::atof(keyItem->getAttribute("value").c_str());

            // Go to next entry (which is expected to be the value for the key).
            if(++iter == items.end()) {
                std::cout << "VelocityCurve: No matching value for key" << std::endl;
                return;
            }

            // Read value.
            XMLreader* valueItem = *iter;
            if(valueItem->getName() != "value") {
                std::cout << "VelocityCurve: Expected value, aborting..." << std::endl;
                return;
            }

            float value = std::atof(valueItem->getAttribute("value").c_str());

            // Add key-value pair.
            peakVelocities_[key] = value;

            // Next key-value pair.
            iter++;
        }
    }

private:
    std::map<float, float> peakVelocities_;
    bool periodic_;
};

////////// Globals //////////////////
// Meta
const std::string simulation = "default";
const std::string base = "/scratch/tmp/s_leis06/simulations/";

// Config
T simulationTime = 0.0;
int numTimeSteps = 1;
int outputResolution = 1;
std::string outputFileFormat;
int flowFeatures = FF_NONE;
std::vector<FlowIndicator> flowIndicators;
std::vector<MeasuredData> measuredData;

// Parameters
int spatialResolution = 1;
T relaxationTime = 0.0;
T characteristicLength = 0.0;
T characteristicVelocity = 0.0;
T viscosity = 0.0;
T density = 0.0;
FlowTurbulenceModel turbulenceModel = FTM_NONE;
T smagorinskyConstant = 0.0;
FlowBoundaryCondition wallBoundaryCondition = FBC_NONE;
T inletVelocityMultiplier = 1.0;
//////////////////////////////////////


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
    simulationTime          = std::atof(config["simulationTime"].getAttribute("value").c_str());
    numTimeSteps            = std::atoi(config["numTimeSteps"].getAttribute("value").c_str());
    outputResolution        = std::atoi(config["outputResolution"].getAttribute("value").c_str());
    outputFileFormat        =           config["outputFileFormat"].getAttribute("value");
    flowFeatures            = std::atoi(config["flowFeatures"].getAttribute("value").c_str());

    XMLreader parameters    = config["flowParameters"];
    spatialResolution       = std::atoi(parameters["spatialResolution"].getAttribute("value").c_str());
    relaxationTime          = std::atof(parameters["relaxationTime"].getAttribute("value").c_str());
    characteristicLength    = std::atof(parameters["characteristicLength"].getAttribute("value").c_str());
    characteristicVelocity  = std::atof(parameters["characteristicVelocity"].getAttribute("value").c_str());
    viscosity               = std::atof(parameters["viscosity"].getAttribute("value").c_str());
    density                 = std::atof(parameters["density"].getAttribute("value").c_str());
    turbulenceModel         = static_cast<FlowTurbulenceModel>(std::atof(parameters["turbulenceModel"].getAttribute("value").c_str()));
    smagorinskyConstant     = std::atof(parameters["smagorinskyConstant"].getAttribute("value").c_str());
    wallBoundaryCondition   = static_cast<FlowBoundaryCondition>(std::atof(parameters["wallBoundaryCondition"].getAttribute("value").c_str())); // Replaces: bouzidiOn = parameters["bouzidi"].getAttribute("value") == "true";
    inletVelocityMultiplier = std::atof(parameters["inletVelocityMultiplier"].getAttribute("value").c_str());

    XMLreader indicators = config["flowIndicators"];
    for(auto iter : indicators) {
        FlowIndicator indicator;
        indicator.type_                 = static_cast<FlowIndicatorType>(std::atoi((*iter)["type_"].getAttribute("value").c_str()));
        indicator.id_                   = std::atoi((*iter)["id_"].getAttribute("value").c_str());
        indicator.center_[0]            = std::atof((*iter)["center"].getAttribute("x").c_str());
        indicator.center_[1]            = std::atof((*iter)["center"].getAttribute("y").c_str());
        indicator.center_[2]            = std::atof((*iter)["center"].getAttribute("z").c_str());
        indicator.normal_[0]            = std::atof((*iter)["normal"].getAttribute("x").c_str());
        indicator.normal_[1]            = std::atof((*iter)["normal"].getAttribute("y").c_str());
        indicator.normal_[2]            = std::atof((*iter)["normal"].getAttribute("z").c_str());
        indicator.radius_               = std::atof((*iter)["radius"].getAttribute("value").c_str());
        indicator.flowProfile_          = static_cast<FlowProfile>(std::atoi((*iter)["flowProfile"].getAttribute("value").c_str()));
        indicator.velocityCurve_.deserialize((*iter)["velocityCurve"]);
        flowIndicators.push_back(indicator);
    }
    clout << "Found " << flowIndicators.size() << " Flow Indicators" << std::endl;

    // TODO: implement measured data support!

    const int N = spatialResolution;
    UnitConverterFromResolutionAndRelaxationTime<T, DESCRIPTOR> converter(
            (T) N,                      // Resolution that charPhysLength is resolved by.
            (T) relaxationTime,         // Relaxation time
            (T) characteristicLength,   // charPhysLength: reference length of simulation geometry
            (T) characteristicVelocity, // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
            (T) viscosity,              // physViscosity: physical kinematic viscosity in __m^2 / s__
            (T) density                 // physDensity: physical density in __kg / m^3__
    );
    // Prints the converter log as console output
    converter.print();
    // Writes the converter log in a file
    converter.write(simulation.c_str());

    // === 2nd Step: Prepare Geometry ===

    // Instantiation of the STLreader class
    // file name, voxel size in meter, stl unit in meter, outer voxel no., inner voxel no.
    std::string geometryFileName = "../geometry/geometry.stl";
    STLreader<T> stlReader(geometryFileName.c_str(), converter.getConversionFactorLength(), VOREEN_LENGTH_TO_SI, 1);
    IndicatorLayer3D<T> extendedDomain(stlReader, converter.getConversionFactorLength());

    // Instantiation of a cuboidGeometry with weights
#ifdef PARALLEL_MODE_MPI
    const int noOfCuboids = std::min( 16*spatialResolution, 2*singleton::mpi().getSize() );
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

    Timer<T> timer1(converter.getLatticeTime(simulationTime), superGeometry.getStatistics().getNvoxel());
    timer1.start();

    prepareLattice(sLattice, converter, *bulkDynamics, stlReader, superGeometry);

    timer1.stop();
    timer1.printSummary();

    // === 4th Step: Main Loop with Timer ===
    clout << "starting simulation..." << std::endl;
    util::ValueTracer<T> converge(converter.getLatticeTime(0.5), 1e-5);
    Timer<T> timer(converter.getLatticeTime(simulationTime), superGeometry.getStatistics().getNvoxel());
    timer.start();

    const int maxIteration = converter.getLatticeTime(simulationTime);
    for (int iteration = 0; iteration <= maxIteration; iteration++) {

        // === 5th Step: Definition of Initial and Boundary Conditions ===
        setBoundaryValues(sLattice, converter, iteration, superGeometry);

        // === 6th Step: Collide and Stream Execution ===
        sLattice.collideAndStream();

        // === 7th Step: Computation and Output of the Results ===
        getResults(sLattice, converter, iteration, maxIteration, *bulkDynamics, superGeometry, timer, stlReader);

        // === 8th Step: Check for convergence.
        converge.takeValue(sLattice.getStatistics().getAverageEnergy(), true);
        if(converge.hasConverged()) {
            clout << "Simulation converged!" << std::endl;
            break;
        }

        // === 9th Step: Write checkpoint
        // TODO: implement!
        //sLattice.save("simulation.checkpoint");
    }

    timer.stop();
    timer.printSummary();

    return EXIT_SUCCESS;
}