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

enum Material {
    MAT_EMPTY  = 0,
    MAT_FLUID  = 1,
    MAT_WALL   = 2,
    MAT_COUNT,
};

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
    FlowFunction    startPhaseFunction_{FF_NONE};
    T               startPhaseDuration_{0};
    T               center_[3]{0};
    T               normal_[3]{0};
    T               radius_{0};
    int             materialId_{0};
};

// Measured data.
struct MeasuredData {
    std::vector<float> data;
    Vector<T, 3> min;
    Vector<T, 3> max;
};

class MeasuredDataMapper : public AnalyticalF3D<T, T> {
public:
    MeasuredDataMapper(const MeasuredData& data)
        : AnalyticalF3D<T, T>(3)
        , data(data)
    {
    }
    virtual bool operator() (T output[], const T input[]) {

        if(input[0] >= data.min[0] && input[1] >= data.min[1] && input[2] >= data.min[2] &&
           input[0] <= data.max[0] && input[1] <= data.max[1] && input[2] <= data.max[2]) {
            return false;
        }
/*
        // TODO: implement linear interpolation!
        for(size_t i=0; i < volume_->getNumChannels(); i++) {
            tgt::vec3 voxel = initialState->getVoxelLinear(rwPos, 0, false);
            output[i] = voxel[i] * VOREEN_LENGTH_TO_SI;
        }
*/
        return true;
    }

private:

    const MeasuredData& data;
};

////////// Globals //////////////////
// Meta
const std::string simulation = "default";
const std::string base = "/scratch/tmp/s_leis06/simulations/";
const T VOREEN_LENGTH_TO_SI = 0.001;
const T VOREEN_TIME_TO_SI = 0.001;
const std::string META_DATA_NAME_OFFSET = "Offset";
const std::string META_DATA_NAME_SPACING = "Spacing";
const std::string META_DATA_NAME_TIMESTEP = "Timestep";
const std::string META_DATA_NAME_REAL_WORLD_MAPPING = "RealWorldMapping";

// Config
T simulationTime = 0.0;
T temporalResolution = 0.0;
int spatialResolution = 1;
int numTimeSteps = 1;
int outputResolution = 1;
std::vector<FlowIndicator> flowIndicators;
std::vector<MeasuredData> measuredData;

// Parameters
T characteristicLength = 0.0;
T characteristicVelocity = 0.0;
T viscosity = 0.0;
T density = 0.0;
T smagorinskyConstant = 0.0;
bool bouzidiOn = false;
//////////////////////////////////////


// Stores data from stl file in geometry in form of material numbers
void prepareGeometry(UnitConverter<T, DESCRIPTOR> const& converter, IndicatorF3D<T>& indicator,
                     STLreader<T>& stlReader, SuperGeometry3D<T>& superGeometry) {

    OstreamManager clout(std::cout, "prepareGeometry");
    clout << "Prepare Geometry ..." << std::endl;

    superGeometry.rename(MAT_EMPTY, MAT_WALL,   indicator);
    superGeometry.rename(MAT_WALL,  MAT_FLUID, stlReader);

    superGeometry.clean();

    int materialId = MAT_COUNT;

    for (size_t i = 0; i < flowIndicators.size(); i++) {

        const T* center = flowIndicators[i].center_;
        const T* normal = flowIndicators[i].normal_;
        T radius = flowIndicators[i].radius_;

        // Set material number for inflow
        IndicatorCircle3D<T> flow(center[0]*VOREEN_LENGTH_TO_SI, center[1]*VOREEN_LENGTH_TO_SI, center[2]*VOREEN_LENGTH_TO_SI,
                                  normal[0], normal[1], normal[2],
                                  radius*VOREEN_LENGTH_TO_SI);
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

    // Unsteered simulation.
    if(measuredData.empty()) {

        // Initial conditions
        AnalyticalConst3D<T, T> rhoF(1);
        std::vector<T> velocity(3, T());
        AnalyticalConst3D<T, T> uF(velocity);

        lattice.defineRhoU(superGeometry, MAT_FLUID, rhoF, uF);
        lattice.iniEquilibrium(superGeometry, MAT_FLUID, rhoF, uF);

        // Initialize all values of distribution functions to their local equilibrium
        for (const FlowIndicator& indicator : flowIndicators) {
            lattice.defineRhoU(superGeometry, indicator.materialId_, rhoF, uF);
            lattice.iniEquilibrium(superGeometry, indicator.materialId_, rhoF, uF);
        }
    }
    // Steered simulation.
    else {
        MeasuredDataMapper mapper(measuredData.front());
        lattice.defineU(superGeometry, MAT_FLUID, mapper);
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
    int iTupdate = 50;

    if (iT % iTupdate == 0) {
        for(const FlowIndicator& indicator : flowIndicators) {
            if (indicator.direction_ == FD_IN) {

                int iTvec[1] = {iT};
                T maxVelocity[1] = {T()};

                switch(indicator.startPhaseFunction_) {
                case FF_SINUS:
                {
                    int iTperiod = converter.getLatticeTime(indicator.startPhaseDuration_);
                    if(iT < iTperiod) {
                        SinusStartScale<T, int> nSinusStartScale(iTperiod, converter.getCharLatticeVelocity());
                        nSinusStartScale(maxVelocity, iTvec);
                        break;
                    }
                    // Else: fallthrough
                }
                case FF_CONSTANT:
                {
                    AnalyticalConst1D<T, int> nConstantStartScale(converter.getCharLatticeVelocity());
                    nConstantStartScale(maxVelocity, iTvec);
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

// Writes result
void writeResult(STLreader<T>& stlReader,
                 UnitConverter<T,DESCRIPTOR>& converter, int ti, int tmax,
                 SuperLatticeF3D<T, DESCRIPTOR>& feature,
                 const std::string& name) {

    OstreamManager clout(std::cout, "writeResult");

    const Vector<T, 3>& min = stlReader.getMin();
    const Vector<T, 3>& max = stlReader.getMax();

    const Vector<T, 3> len = (max - min);
    const T maxLen = std::max({len[0], len[1], len[2]});
    const int gridResolution = static_cast<int>(std::round(maxLen / converter.getConversionFactorLength()));
    const int resolution = std::min(outputResolution, gridResolution);

    Vector<T, 3> offset = min + (len - maxLen) * 0.5;
    Vector<T, 3> spacing(maxLen / (resolution-1));

    // Determine format.
    // This could be done in a more dynamic way, but the code should be easily portable to the cluster.
    std::string format;
    switch(feature.getTargetDim()) {
        case 1:
            format = "float";
            break;
        case 3:
            format = "Vector3(float)";
            break;
        default:
            clout << "Unhandled target dimensions" << std::endl;
            return;
    }

    std::vector<float> rawFeatureData;
    rawFeatureData.reserve(static_cast<size_t>(resolution * resolution * resolution * feature.getTargetDim()));
    AnalyticalFfromSuperF3D<T> interpolateProperty(feature, true);

    std::vector<T> minValue(feature.getTargetDim(), std::numeric_limits<T>::max());
    std::vector<T> maxValue(feature.getTargetDim(), std::numeric_limits<T>::lowest());

    T minMagnitude = std::numeric_limits<T>::max();
    T maxMagnitude = 0.0f;

    for(int z=0; z<resolution; z++) {
        for(int y=0; y<resolution; y++) {
            for(int x=0; x<resolution; x++) {

                T pos[3] = {offset[0]+x*maxLen/resolution, offset[1]+y*maxLen/resolution, offset[2]+z*maxLen/resolution};
                std::vector<T> val(feature.getTargetDim(), 0.0f);

                if(pos[0] >= min[0] && pos[1] >= min[1] && pos[2] >= min[2] &&
                   pos[0] <= max[0] && pos[1] <= max[1] && pos[2] <= max[2]) {
                    interpolateProperty(&val[0], pos);

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
                for(int i = 0; i < feature.getTargetDim(); i++) {
                    rawFeatureData.push_back(static_cast<float>(val[i]/VOREEN_LENGTH_TO_SI));
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
    std::string rawFilename = singleton::directories().getLogOutDir() + "/" + featureFilename + ".raw";
    std::string vvdFilename = singleton::directories().getLogOutDir() + "/" + featureFilename + ".vvd";

    const LatticeStatistics<T>& statistics = feature.getSuperLattice().getStatistics();
    std::fstream vvdFeatureFile(vvdFilename.c_str(), std::ios::out);
    vvdFeatureFile
            // Header.
            << "<?xml version=\"1.0\" ?>"
            << "<VoreenData version=\"1\">"
            << "<Volumes>"
            << "<Volume>"
            // Data.
            << "<RawData filename=\"" << featureFilename << ".raw\" format=\"" << format << "\" x=\"" << resolution << "\" y=\""<< resolution << "\" z=\"" << resolution << "\" />"
            // Mandatory Meta data.
            << "<MetaData>"
            << "<MetaItem name=\""<< META_DATA_NAME_OFFSET << "\" type=\"Vec3MetaData\">"
            << "<value x=\"" << offset[0] << "\" y=\"" << offset[1] << "\" z=\"" << offset[2] << "\" />"
            << "</MetaItem>"
            << "<MetaItem name=\"" << META_DATA_NAME_SPACING << "\" type=\"Vec3MetaData\">"
            << "<value x=\"" << spacing[0] << "\" y=\"" << spacing[1] << "\" z=\"" << spacing[2] << "\" />"
            << "</MetaItem>"
            << "<MetaItem name=\"" << META_DATA_NAME_TIMESTEP << "\" type=\"FloatMetaData\" value=\"" << converter.getPhysTime(ti) << "\" />"
            << "<MetaItem name=\"" << META_DATA_NAME_REAL_WORLD_MAPPING << "\" type=\"RealWorldMappingMetaData\"><value scale=\"1\" offset=\"0\" unit=\"\" /></MetaItem>"
            << "<MetaItem name=\"" << "name" << "\" type=\"StringMetaData\" value=\"" << name << "\" />"
            // Parameters.
            << "<MetaItem name=\"" << "ParameterCharacteristicLength" << "\" type=\"FloatMetaData\" value=\"" << characteristicLength << "\" />"
            << "<MetaItem name=\"" << "ParameterCharacteristicVelocity" << "\" type=\"FloatMetaData\" value=\"" << characteristicVelocity << "\" />"
            << "<MetaItem name=\"" << "ParameterViscosity" << "\" type=\"FloatMetaData\" value=\"" << viscosity << "\" />"
            << "<MetaItem name=\"" << "ParameterDensity" << "\" type=\"FloatMetaData\" value=\"" << density << "\" />"
            << "<MetaItem name=\"" << "ParameterSmagorinskyConstant" << "\" type=\"FloatMetaData\" value=\"" << smagorinskyConstant << "\" />"
            << "<MetaItem name=\"" << "ParameterBouzidi" << "\" type=\"BoolMetaData\" value=\"" << (bouzidiOn ? "true" : "false") << "\" />"
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
        clout << "Could not write " << name << " file" << std::endl;
    }
}

// Computes flux at inflow and outflow
void getResults(SuperLattice3D<T, DESCRIPTOR>& sLattice,
                UnitConverter<T, DESCRIPTOR>& converter, int ti, int tmax,
                Dynamics<T, DESCRIPTOR>& bulkDynamics,
                SuperGeometry3D<T>& superGeometry, Timer<T>& timer, STLreader<T>& stlReader) {

    OstreamManager clout(std::cout, "getResults");

    const int outputIter = tmax / numTimeSteps;

    int rank = 0;
#ifdef PARALLEL_MODE_MPI
    rank = singleton::mpi().getRank();
#endif
    if (rank == 0 && ti % outputIter == 0) {
        // Write velocity.
        SuperLatticePhysVelocity3D<T, DESCRIPTOR> velocity(sLattice, converter);
        writeResult(stlReader, converter, ti, tmax, velocity, "velocity");

        // Write magnitude.
        SuperEuklidNorm3D <T, DESCRIPTOR> magnitude(velocity);
        writeResult(stlReader, converter, ti, tmax, magnitude, "magnitude");
/*
        // TODO: Pressure and WSS currently do not have an equivalent in measured data!

        // Write pressure.
        SuperLatticePhysPressure3D<T, DESCRIPTOR> pressure(sLattice, converter);
        writeResult(stlReader, converter, ti, tmax, pressure, "pressure");

#ifndef OLB_PRECOMPILED
        // Write wallShearStress.
        SuperLatticePhysWallShearStress3D <T, DESCRIPTOR> wallShearStress(sLattice, superGeometry, MAT_WALL, converter, stlReader);
        writeResult(stlReader, converter, ti, tmax, wallShearStress, "wallShearStress");
#endif
*/
        // Lattice statistics console output
        sLattice.getStatistics().print(ti, converter.getPhysTime(ti));
    }

    if (sLattice.getStatistics().getMaxU() > 0.3) {
        clout << "PROBLEM uMax=" << sLattice.getStatistics().getMaxU() << std::endl;
        std::exit(EXIT_FAILURE);
    }
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
    simulationTime          = std::atof(config["simulationTime"].getAttribute("value").c_str());
    temporalResolution      = std::atof(config["temporalResolution"].getAttribute("value").c_str());
    spatialResolution       = std::atoi(config["spatialResolution"].getAttribute("value").c_str());
    numTimeSteps            = std::atoi(config["numTimeSteps"].getAttribute("value").c_str());
    outputResolution        = std::atoi(config["outputResolution"].getAttribute("value").c_str());

    XMLreader parameters = config["flowParameters"];
    characteristicLength    = std::atof(parameters["characteristicLength"].getAttribute("value").c_str());
    characteristicVelocity  = std::atof(parameters["characteristicVelocity"].getAttribute("value").c_str());
    viscosity               = std::atof(parameters["viscosity"].getAttribute("value").c_str());
    density                 = std::atof(parameters["density"].getAttribute("value").c_str());
    smagorinskyConstant     = std::atof(parameters["smagorinskyConstant"].getAttribute("value").c_str());
    bouzidiOn               = parameters["bouzidi"].getAttribute("value") == "true";

    XMLreader indicators = config["flowIndicators"];
    for(auto iter : indicators) {
        FlowIndicator indicator;
        indicator.direction_            = static_cast<FlowDirection>(std::atoi((*iter)["direction"].getAttribute("value").c_str()));
        indicator.startPhaseFunction_   = static_cast<FlowFunction>(std::atoi((*iter)["startPhaseFunction"].getAttribute("value").c_str()));
        indicator.startPhaseDuration_   = std::atof((*iter)["startPhaseDuration"].getAttribute("value").c_str());
        indicator.center_[0]            = std::atof((*iter)["center"].getAttribute("x").c_str());
        indicator.center_[1]            = std::atof((*iter)["center"].getAttribute("y").c_str());
        indicator.center_[2]            = std::atof((*iter)["center"].getAttribute("z").c_str());
        indicator.normal_[0]            = std::atof((*iter)["normal"].getAttribute("x").c_str());
        indicator.normal_[1]            = std::atof((*iter)["normal"].getAttribute("y").c_str());
        indicator.normal_[2]            = std::atof((*iter)["normal"].getAttribute("z").c_str());
        indicator.radius_               = std::atof((*iter)["radius"].getAttribute("value").c_str());
        flowIndicators.push_back(indicator);
    }
    clout << "Found " << flowIndicators.size() << " Flow Indicators" << std::endl;

    // TODO: implement measured data support!

    const int N = spatialResolution;
    UnitConverter<T, DESCRIPTOR> converter(
            (T) characteristicLength * VOREEN_LENGTH_TO_SI / N,  // physDeltaX: spacing between two lattice cells in __m__
            (T) temporalResolution * VOREEN_TIME_TO_SI,          // TODO: define proper semantic
            (T) characteristicLength * VOREEN_LENGTH_TO_SI,      // charPhysLength: reference length of simulation geometry
            (T) characteristicVelocity * VOREEN_LENGTH_TO_SI,    // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
            (T) viscosity * 0.001 / density,                     // physViscosity: physical kinematic viscosity in __m^2 / s__
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
    STLreader<T> stlReader(geometryFileName.c_str(), converter.getConversionFactorLength(), VOREEN_LENGTH_TO_SI, 1, true);
    IndicatorLayer3D<T> extendedDomain(stlReader, converter.getConversionFactorLength());

    // Instantiation of a cuboidGeometry with weights
#ifdef PARALLEL_MODE_MPI
    const int noOfCuboids = std::min( 16*spatialResolution,2*singleton::mpi().getSize() );
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
                                                       instances::getBulkMomenta<T, DESCRIPTOR>(),
                                                       smagorinskyConstant);

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

    const int tmax = converter.getLatticeTime(simulationTime);
    for (int ti = 0; ti <= tmax; ti++) {

        // === 5th Step: Definition of Initial and Boundary Conditions ===
        setBoundaryValues(sLattice, sOffBoundaryCondition, converter, ti, superGeometry);

        // === 6th Step: Collide and Stream Execution ===
        sLattice.collideAndStream();

        // === 7th Step: Computation and Output of the Results ===
        getResults(sLattice, converter, ti, tmax, bulkDynamics, superGeometry, timer, stlReader);

        // === 8th Step: Check for convergence.
        converge.takeValue(sLattice.getStatistics().getAverageEnergy(), true);
        if(converge.hasConverged()) {
            clout << "Simulation converged!" << std::endl;
            break;
        }
    }

    timer.stop();
    timer.printSummary();

    return EXIT_SUCCESS;
}
