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
typedef D3Q19<> DESCRIPTOR;

enum Material {
    MAT_EMPTY  = 0,
    MAT_FLUID  = 1,
    MAT_WALL   = 2,
    MAT_COUNT,
};

enum FlowFeatures {
    FF_NONE             = 0,
    FF_VELOCITY         = 1,
    FF_MAGNITUDE        = 2,
    FF_PRESSURE         = 4,
    FF_WALLSHEARSTRESS  = 8,
};

enum FlowIndicatorType {
    FIT_INVALID         = -1,
    FIT_CANDIDATE       =  0,
    FIT_VELOCITY        =  1,
    FIT_PRESSURE        =  2,
    FIT_MEASURE         =  3,
};

enum FlowProfile {
    FP_NONE       = 0,
    FP_POISEUILLE = 1,
    FP_POWERLAW   = 2,
    FP_CONSTANT   = 3,
};

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

// Indicates flux through an arbitrary, circle-shaped area.
// This code is adapted from the voreen host code.
struct FlowIndicator {
    FlowIndicatorType   type_{FIT_INVALID};
    int                 id_{0};
    T                   center_[3]{0};
    T                   normal_[3]{0};
    T                   radius_{0};
    FlowProfile         flowProfile_{FP_NONE};
    VelocityCurve       velocityCurve_;
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
const std::string META_DATA_NAME_MODALITY = "Modality";

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
T smagorinskyConstant = 0.0;
bool bouzidiOn = false;
//////////////////////////////////////


// Stores data from stl file in geometry in form of material numbers
void prepareGeometry(UnitConverter<T, DESCRIPTOR> const& converter,
                     IndicatorF3D<T>& indicator,
                     STLreader<T>& stlReader,
                     SuperGeometry3D<T>& superGeometry) {

    OstreamManager clout(std::cout, "prepareGeometry");
    clout << "Prepare Geometry ..." << std::endl;

    superGeometry.rename(MAT_EMPTY, MAT_WALL,  indicator);
    superGeometry.rename(MAT_WALL,  MAT_FLUID, stlReader);

    superGeometry.clean();

    for (size_t i = 0; i < flowIndicators.size(); i++) {

        Vector<T, 3> normal(flowIndicators[i].normal_);

        // Convert to SI units.
        Vector<T, 3> center(flowIndicators[i].center_);
        center *= VOREEN_LENGTH_TO_SI;

        // Add one voxel to account for rounding errors.
        T radius = flowIndicators[i].radius_ * VOREEN_LENGTH_TO_SI + converter.getConversionFactorLength();

        // Define a local disk volume.
        IndicatorCircle3D<T> flow(center[0]*VOREEN_LENGTH_TO_SI, center[1]*VOREEN_LENGTH_TO_SI, center[2]*VOREEN_LENGTH_TO_SI,
                                  normal[0], normal[1], normal[2],
                                  radius);
        IndicatorCylinder3D<T> layerFlow(flow, 2 * converter.getConversionFactorLength());

        // Rename both, wall and fluid, since the indicator might also be inside the fluid domain.
        superGeometry.rename(MAT_WALL, flowIndicators[i].id_, MAT_FLUID, layerFlow);
        superGeometry.rename(MAT_FLUID, flowIndicators[i].id_,  layerFlow);

        // Exclude area behind inlet and in front of outlet - it will otherwise cause unstable simulations.
        bool isInlet = flowIndicators[i].type_ == FIT_VELOCITY;
        bool isOutlet = flowIndicators[i].type_ == FIT_PRESSURE;
        if(isInlet || isOutlet) {
            T sign = isInlet ? T(-1) : T(1);
            center += sign * normal * (converter.getConversionFactorLength() * 2);
            IndicatorCircle3D<T> capFlow(center[0], center[1], center[2],
                                         normal[0], normal[1], normal[2],
                                         radius);
            IndicatorCylinder3D<T> layerCapFlow(capFlow, 2 * converter.getConversionFactorLength());
            superGeometry.rename(MAT_FLUID, MAT_WALL, layerCapFlow);
        }
    }

    // Removes all not needed boundary voxels outside the surface
    superGeometry.clean();
    // Removes all not needed boundary voxels inside the surface
    superGeometry.innerClean(MAT_COUNT);
    superGeometry.checkForErrors();

    superGeometry.print();
    clout << "Prepare Geometry ... OK" << std::endl;
}

// If any dynamics should be missing, we fix this by using "no dynamics" to prevent a crash.
size_t fixupLattice(SuperLattice3D<T, DESCRIPTOR>& lattice) {
    size_t errors = 0;
    for (int iC = 0; iC < lattice.getLoadBalancer().size(); ++iC) {
        auto& blockLattice = lattice.getExtendedBlockLattice(iC);
        for (int iX = 0; iX < blockLattice.getNx(); ++iX) {
            for (int iY = 0; iY < blockLattice.getNy(); ++iY) {
                for (int iZ = 0; iZ < blockLattice.getNz(); ++iZ) {
                    auto cell = blockLattice.get(iX, iY, iZ);
                    auto dynamics = cell.getDynamics();
                    if (!dynamics) {
                        errors++;
                        //std::cerr << "no dynamics at: " << iX << ", " << iY << ", " << iZ << std::endl;
                        cell.defineDynamics(&instances::getNoDynamics<T, DESCRIPTOR>());
                    }
                }
            }
        }
    }
    return errors;
}

// Set up the geometry of the simulation
void prepareLattice(SuperLattice3D<T, DESCRIPTOR>& lattice,
                    UnitConverter<T, DESCRIPTOR> const& converter,
                    Dynamics<T, DESCRIPTOR>& bulkDynamics,
                    STLreader<T>& stlReader,
                    SuperGeometry3D<T>& superGeometry) {

    OstreamManager clout(std::cout, "prepareLattice");
    clout << "Prepare Lattice ..." << std::endl;

    const T omega = converter.getLatticeRelaxationFrequency();

    // material=0 --> do nothing
    lattice.defineDynamics(superGeometry.getMaterialIndicator(MAT_EMPTY), &instances::getNoDynamics<T, DESCRIPTOR>());

    // material=1 --> bulk dynamics
    lattice.defineDynamics(superGeometry.getMaterialIndicator(MAT_FLUID), &bulkDynamics);

    if (bouzidiOn) {
        // material=2 --> no dynamics + bouzidi zero velocity
        lattice.defineDynamics(superGeometry.getMaterialIndicator(MAT_WALL), &instances::getNoDynamics<T, DESCRIPTOR>());
        setBouzidiZeroVelocityBoundary<T, DESCRIPTOR>(lattice, superGeometry, MAT_WALL, stlReader);
    } else {
        // material=2 --> bounceBack dynamics
        lattice.defineDynamics(superGeometry.getMaterialIndicator(MAT_WALL), &instances::getBounceBack<T, DESCRIPTOR>());
    }

    for(const auto& indicator : flowIndicators) {
        if(indicator.type_ == FIT_VELOCITY) {
            if(bouzidiOn) {
                // no dynamics + bouzidi velocity (inflow)
                lattice.defineDynamics(superGeometry.getMaterialIndicator(indicator.id_), &instances::getNoDynamics<T, DESCRIPTOR>());
                setBouzidiVelocityBoundary<T,DESCRIPTOR>(lattice, superGeometry, indicator.id_, stlReader);
            }
            else {
                // bulk dynamics + velocity (inflow)
                lattice.defineDynamics(superGeometry.getMaterialIndicator(indicator.id_), &bulkDynamics);
                setInterpolatedVelocityBoundary<T,DESCRIPTOR>(lattice, omega, superGeometry.getMaterialIndicator(indicator.id_));
            }
        }
        else if(indicator.type_ == FIT_PRESSURE) {
            lattice.defineDynamics(superGeometry.getMaterialIndicator(indicator.id_), &bulkDynamics);
            setInterpolatedPressureBoundary<T,DESCRIPTOR>(lattice, omega, superGeometry.getMaterialIndicator(indicator.id_));
        }
    }

    // If any dynamics should be missing up to this point, we fix this by using "no dynamics".
    size_t errors = fixupLattice(lattice);
    if(errors > 0) {
        clout << errors << " errors have been fixed" << std::endl;
    }

    // Unsteered simulation.
    if(measuredData.empty()) {

        // Initial conditions
        AnalyticalConst3D<T, T> rhoF(1);
        std::vector<T> velocity(3, T());
        AnalyticalConst3D<T, T> uF(velocity);

        lattice.defineRhoU(superGeometry.getMaterialIndicator(MAT_FLUID), rhoF, uF);
        lattice.iniEquilibrium(superGeometry.getMaterialIndicator(MAT_FLUID), rhoF, uF);

        // Initialize all values of distribution functions to their local equilibrium.
        for (const auto& indicator : flowIndicators) {
            lattice.defineRhoU(superGeometry.getMaterialIndicator(indicator.id_), rhoF, uF);
            lattice.iniEquilibrium(superGeometry.getMaterialIndicator(indicator.id_), rhoF, uF);
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
                       UnitConverter<T, DESCRIPTOR> const& converter,
                       int iteration,
                       SuperGeometry3D<T>& superGeometry) {

    for(const auto& indicator : flowIndicators) {
        if (indicator.type_ == FIT_VELOCITY) {

            T targetPhysVelocity = indicator.velocityCurve_(converter.getPhysTime(iteration));
            T targetLatticeVelocity = converter.getLatticeVelocity(targetPhysVelocity);

            // This function applies the velocity profile to the boundary condition and the lattice.
            auto applyFlowProfile = [&] (AnalyticalF3D<T,T>& profile) {
                if (bouzidiOn) {
                    defineUBouzidi<T, DESCRIPTOR>(sLattice, superGeometry, indicator.id_, profile);
                } else {
                    sLattice.defineU(superGeometry.getMaterialIndicator(indicator.id_), profile);
                }
            };

            // Create shortcuts.
            const Vector<T, 3> center(indicator.center_);
            const Vector<T, 3> normal(indicator.normal_);
            T radius = indicator.radius_ * VOREEN_LENGTH_TO_SI;

            // Apply the indicator's profile.
            switch(indicator.flowProfile_) {
            case FP_POISEUILLE:
            {
//                CirclePoiseuille3D<T> profile(superGeometry, indicator.id_, targetLatticeVelocity); // This is the alternative way, but how does it work?
                    CirclePoiseuille3D<T> profile(center[0]*VOREEN_LENGTH_TO_SI, center[1]*VOREEN_LENGTH_TO_SI, center[2]*VOREEN_LENGTH_TO_SI,
                                                  normal[0], normal[1], normal[2], radius, targetLatticeVelocity);
                applyFlowProfile(profile);
                break;
            }
            case FP_POWERLAW:
            {
                T n = 1.03 * std::log(converter.getReynoldsNumber()) - 3.6; // Taken from OLB documentation.
                CirclePowerLawTurbulent3D<T> profile(center[0]*VOREEN_LENGTH_TO_SI, center[1]*VOREEN_LENGTH_TO_SI, center[2]*VOREEN_LENGTH_TO_SI,
                                                     normal[0], normal[1], normal[2], radius, targetLatticeVelocity, n);
                applyFlowProfile(profile);
                break;
            }
            case FP_CONSTANT:
            {
                AnalyticalConst3D<T, T> profile(normal[0] * targetLatticeVelocity, normal[1] * targetLatticeVelocity, normal[2] * targetLatticeVelocity);
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

// Writes result
void writeVVDFile(STLreader<T>& stlReader,
                  UnitConverter<T,DESCRIPTOR>& converter,
                  int iteration, int maxIteration,
                  SuperLatticeF3D<T, DESCRIPTOR>& feature,
                  const std::string& name) {

    OstreamManager clout(std::cout, "writeVVDFile");

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
                    rawFeatureData.push_back(static_cast<float>(val[i]));
                }
            }
        }
    }

    // Adapt to Voreen units.
    offset  = offset  * (1/VOREEN_LENGTH_TO_SI);
    spacing = spacing * (1/VOREEN_LENGTH_TO_SI);

    minMagnitude = std::sqrt(minMagnitude);
    maxMagnitude = std::sqrt(maxMagnitude);

    // Set output names.
    int maxIterationLen = static_cast<int>(std::to_string(maxIteration).length());
    std::ostringstream suffix;
    suffix << std::setw(maxIterationLen) << std::setfill('0') << iteration;
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
            << "<MetaItem name=\"" << META_DATA_NAME_TIMESTEP << "\" type=\"FloatMetaData\" value=\"" << converter.getPhysTime(iteration) << "\" />"
            << "<MetaItem name=\"" << META_DATA_NAME_REAL_WORLD_MAPPING << "\" type=\"RealWorldMappingMetaData\"><value scale=\"1\" offset=\"0\" unit=\"\" /></MetaItem>"
            << "<MetaItem name=\"" << META_DATA_NAME_MODALITY << "\" type=\"StringMetaData\" value=\"" << name << "\" />"
            // Parameters.
            << "<MetaItem name=\"" << "ParameterSpatialResolution" << "\" type=\"IntMetaData\" value=\"" << spatialResolution << "\" />"
            << "<MetaItem name=\"" << "ParameterRelaxationTime" << "\" type=\"FloatMetaData\" value=\"" << relaxationTime << "\" />"
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
                UnitConverter<T, DESCRIPTOR>& converter,
                int iteration, int maxIteration,
                Dynamics<T, DESCRIPTOR>& bulkDynamics,
                SuperGeometry3D<T>& superGeometry, Timer<T>& timer, STLreader<T>& stlReader) {

    OstreamManager clout(std::cout, "getResults");

    const int outputIter = maxIteration / numTimeSteps;

    int rank = 0;
#ifdef PARALLEL_MODE_MPI
    rank = singleton::mpi().getRank();
#endif
    if (rank == 0 && iteration % outputIter == 0) {

        bool writeVVD = outputFileFormat == ".vvd";
        bool writeVTI = outputFileFormat == ".vti";

        SuperVTMwriter3D<T> vtmWriter( "results" );

        // Always write debug data.
        if(iteration == 0) {
            SuperLatticeGeometry3D<T, DESCRIPTOR> geometry( sLattice, superGeometry );
            vtmWriter.write( geometry );

            SuperLatticeCuboid3D<T, DESCRIPTOR> cuboid( sLattice );
            vtmWriter.write( cuboid );

            SuperLatticeRank3D<T, DESCRIPTOR> rank( sLattice );
            vtmWriter.write( rank );

            vtmWriter.createMasterFile();
        }

        if(flowFeatures & FF_VELOCITY) {
            SuperLatticePhysVelocity3D<T, DESCRIPTOR> velocity(sLattice, converter);
            if(writeVVD) {
                writeVVDFile(stlReader, converter, iteration, maxIteration, velocity, "velocity");
            }
            if(writeVTI) {
                vtmWriter.write(velocity, iteration);
            }
        }

        if(flowFeatures & FF_MAGNITUDE) {
            SuperLatticePhysVelocity3D<T, DESCRIPTOR> velocity(sLattice, converter);
            SuperEuklidNorm3D<T, DESCRIPTOR> magnitude(velocity);
            if(writeVVD) {
                writeVVDFile(stlReader, converter, iteration, maxIteration, magnitude, "magnitude");
            }
            if(writeVTI) {
                vtmWriter.write(magnitude, iteration);
            }
        }

        if(flowFeatures & FF_PRESSURE) {
            SuperLatticePhysPressure3D<T, DESCRIPTOR> pressure(sLattice, converter);
            if(writeVVD) {
                writeVVDFile(stlReader, converter, iteration, maxIteration, pressure, "pressure");
            }
            if(writeVTI) {
                vtmWriter.write(pressure, iteration);
            }
        }

#ifndef OLB_PRECOMPILED
        if(flowFeatures & FF_WALLSHEARSTRESS) {
            SuperLatticePhysWallShearStress3D<T, DESCRIPTOR> wallShearStress(sLattice, superGeometry, MAT_WALL,
                                                                             converter, stlReader);
            if(writeVVD) {
                writeVVDFile(stlReader, converter, iteration, maxIteration, wallShearStress, "wallShearStress");
            }
            if(writeVTI) {
                vtmWriter.write(wallShearStress, iteration);
            }
        }
#endif

        // Lattice statistics console output
        sLattice.getStatistics().print(iteration, converter.getPhysTime(iteration));
    }

    T tau = converter.getLatticeRelaxationFrequency();
    T threshold = tau < 0.55 ? 0.125*(tau - 0.5) : 0.4;
    if (sLattice.getStatistics().getMaxU() >= threshold) {
        clout << "uMax=" << sLattice.getStatistics().getMaxU() << " above threshold=" << threshold;
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
    smagorinskyConstant     = std::atof(parameters["smagorinskyConstant"].getAttribute("value").c_str());
    bouzidiOn               = parameters["bouzidi"].getAttribute("value") == "true";

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
    const int noOfCuboids = std::min( 16*spatialResolution, singleton::mpi().getSize() );
#else
    const int noOfCuboids = 1;
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

    Timer<T> timer1(converter.getLatticeTime(simulationTime), superGeometry.getStatistics().getNvoxel());
    timer1.start();

    prepareLattice(sLattice, converter, bulkDynamics, stlReader, superGeometry);

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
        getResults(sLattice, converter, iteration, maxIteration, bulkDynamics, superGeometry, timer, stlReader);

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