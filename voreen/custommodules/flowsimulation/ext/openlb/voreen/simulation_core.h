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

#include "olb3D.h"
#include "olb3D.hh"

#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <dirent.h>
#include <functional>

#include "openlb_parameters.hh"

using namespace olb;
using namespace olb::descriptors;

#ifdef PLATFORM_GPU_CUDA
using T = float;
#else
using T = double;
#endif
using DESCRIPTOR = D3Q19<>;

using namespace voreen;

namespace {

using VolumeSampler = std::function<void(UnitConverter<T, DESCRIPTOR> const&, float, std::function<void(AnalyticalF3D<T,T>&)>&, float)>;

class LatticePerturber : public AnalyticalF3D<T, T> {
public:
    LatticePerturber(T maxNoise=1e-5)
        : AnalyticalF3D<T, T>(3)
        , rnd_(std::bind(std::uniform_real_distribution<T>(-maxNoise, maxNoise), std::mt19937(time(nullptr))))
    {
    }
    virtual bool operator() (T output[], const T input[]) {
        for(size_t i=0; i<3; i++) {
            output[i] = rnd_();
        }
        return true;
    }

private:
    std::function<T()> rnd_;
};

// Stores data from stl file in geometry in form of material numbers
void prepareGeometry(UnitConverter<T,DESCRIPTOR> const& converter,
                     IndicatorF3D<T>& indicator,
                     STLreader<T>& stlReader,
                     SuperGeometry<T,3>& superGeometry,
                     const std::vector<FlowIndicator>& indicators)
{
    superGeometry.rename( MAT_EMPTY, MAT_WALL,  indicator );
    superGeometry.rename( MAT_WALL,  MAT_FLUID, stlReader );

    superGeometry.clean();

    for (size_t i = 0; i < indicators.size(); i++) {

        Vector<T, 3> normal(indicators[i].normal_.x, indicators[i].normal_.y, indicators[i].normal_.z);

        // Convert to SI units.
        Vector<T, 3> center(indicators[i].center_.x, indicators[i].center_.y, indicators[i].center_.z);
        center *= VOREEN_LENGTH_TO_SI;

        // Add one voxel to account for rounding errors.
        T radius = indicators[i].radius_ * VOREEN_LENGTH_TO_SI + converter.getConversionFactorLength();
        T length = indicators[i].length_ * VOREEN_LENGTH_TO_SI + converter.getConversionFactorLength() * 2;

        // Define a local disk volume.
        IndicatorCircle3D<T> flow(center[0], center[1], center[2],
                                  normal[0], normal[1], normal[2],
                                  radius);
        IndicatorCylinder3D<T> layerFlow(flow, length);

        // Rename both, wall and fluid, since the indicator might also be inside the fluid domain.
        superGeometry.rename(MAT_WALL, indicators[i].id_, layerFlow);
        superGeometry.rename(MAT_FLUID, indicators[i].id_,  layerFlow);

        // Exclude area behind inlet and in front of outlet - it will otherwise cause unstable simulations.
        bool isInlet = indicators[i].type_ == FIT_VELOCITY;
        bool isOutlet = indicators[i].type_ == FIT_PRESSURE;
        if(isInlet || isOutlet) {
            T sign = isInlet ? T(-1) : T(1);
            center += sign * normal * T(length * 0.5 + converter.getConversionFactorLength());
            IndicatorCircle3D<T> capFlowWall(center[0], center[1], center[2],
                                             normal[0], normal[1], normal[2],
                                             radius);

            IndicatorCylinder3D<T> layerCapFlowWall(capFlowWall, 4 * converter.getConversionFactorLength());
            superGeometry.rename(MAT_FLUID, MAT_WALL, layerCapFlowWall);

            IndicatorCircle3D<T> capFlowEmpty(center[0], center[1], center[2],
                                              normal[0], normal[1], normal[2],
                                              radius);
            IndicatorCylinder3D<T> layerCapFlowEmpty(capFlowEmpty, 2 * converter.getConversionFactorLength());
            superGeometry.rename(MAT_WALL, MAT_EMPTY, layerCapFlowEmpty);
        }
    }

    // TODO: clean regions that are isolated from simulation domain.
    // Removes all not needed boundary voxels outside the surface
    //superGeometry.clean();
    // Removes all not needed boundary voxels inside the surface
    superGeometry.innerClean(MAT_COUNT);
    superGeometry.checkForErrors();
}

/*
// If any dynamics should be missing, we fix this by using "no dynamics" to prevent a crash.
size_t fixupLattice(SuperLattice<T, DESCRIPTOR>& lattice) {
    size_t errors = 0;
    for (int iC = 0; iC < lattice.getLoadBalancer().size(); ++iC) {
        auto& blockLattice = lattice.getBlock(iC);
        for (int iX = 0; iX < blockLattice.getNx(); ++iX) {
            for (int iY = 0; iY < blockLattice.getNy(); ++iY) {
                for (int iZ = 0; iZ < blockLattice.getNz(); ++iZ) {
                    auto cell = blockLattice.get(iX, iY, iZ);
                    auto dynamics = cell.getDynamics();
                    if (!dynamics) {
                        errors++;
                        //std::cerr << "no dynamics at: " << iX << ", " << iY << ", " << iZ << std::endl;
                        //cell.defineDynamics(&instances::getNoDynamics<T, DESCRIPTOR>());
                        // TODO: no longer required??
                    }
                }
            }
        }
    }
    return errors;
}
*/

void defineBulkDynamics(FlowTurbulenceModel model,
                        SuperLattice<T, DESCRIPTOR>& lattice,
                        SuperGeometry<T,3>& geometry,
                        int material) {
    auto indicator = geometry.getMaterialIndicator(material);
    switch(model) {
        case FTM_SMAGORINSKY:
            lattice.defineDynamics<SmagorinskyBGKdynamics<T, DESCRIPTOR>>(indicator);
            lattice.setParameter<collision::LES::Smagorinsky>(0.1);
            return;
        case FTM_SMAGORINSKY_SHEAR_IMPROVED: // Does not compile in OpenLB 1.4.
            lattice.defineDynamics<ShearSmagorinskyBGKdynamics<T, DESCRIPTOR>>(indicator);
            lattice.setParameter<collision::LES::Smagorinsky>(0.1);
            return;
        case FTM_SMAGORINSKY_CONSISTENT:
            lattice.defineDynamics<ConSmagorinskyBGKdynamics<T, DESCRIPTOR>>(indicator);
            lattice.setParameter<collision::LES::Smagorinsky>(0.1);
            return;
        case FTM_SMAGORINSKY_CONSISTENT_STRAIN:
            lattice.defineDynamics<ConStrainSmagorinskyBGKdynamics<T, DESCRIPTOR>>(indicator);
            lattice.setParameter<collision::LES::Smagorinsky>(0.1);
            return;
        case FTM_SMAGORINSKY_DYNAMIC:
            lattice.defineDynamics<ConSmagorinskyBGKdynamics<T, DESCRIPTOR>>(indicator);
            lattice.setParameter<collision::LES::Smagorinsky>(0.1);
            return;
        case FTM_BGK:
            lattice.defineDynamics<BGKdynamics<T, DESCRIPTOR>>(indicator);
            return;
        case FTM_NONE:
        default:
            return;
    };
};

// Set up the geometry of the simulation.
void prepareLattice( SuperLattice<T, DESCRIPTOR>& lattice,
                     UnitConverter<T,DESCRIPTOR> const& converter,
                     STLreader<T>& stlReader,
                     SuperGeometry<T,3>& superGeometry,
                     const std::vector<FlowIndicator>& indicators,
                     const Parameters& parameters)
{
    // material=0 --> do nothing
    lattice.defineDynamics<NoDynamics>(superGeometry, MAT_EMPTY);

    // material=1 --> bulk dynamics
    defineBulkDynamics(parameters.turbulenceModel_, lattice, superGeometry, MAT_FLUID);

    switch(parameters.wallBoundaryCondition_) {
        case FBC_BOUZIDI:
            // material=2 --> no dynamics + bouzidi zero velocity
            lattice.defineDynamics<NoDynamics>(superGeometry, MAT_WALL);
            setBouzidiZeroVelocityBoundary<T, DESCRIPTOR>(lattice, superGeometry, MAT_WALL, stlReader);
            break;
        case FBC_BOUNCE_BACK:
            // material=2 --> bounceBack dynamics
            lattice.defineDynamics<BounceBack>(superGeometry, MAT_WALL);
            break;
        case FBC_NONE:
        default:
            lattice.defineDynamics<NoDynamics>(superGeometry, MAT_WALL);
            break;
    }

    const T omega = converter.getLatticeRelaxationFrequency();

    for(const auto& indicator : indicators) {
        if(indicator.type_ == FIT_VELOCITY) {
            switch(parameters.wallBoundaryCondition_) {
                case FBC_BOUZIDI:
                    // no dynamics + bouzidi velocity (inflow)
                    lattice.defineDynamics<NoDynamics>(superGeometry, indicator.id_);
                    setBouzidiZeroVelocityBoundary<T, DESCRIPTOR>(lattice, superGeometry, indicator.id_, stlReader);
                    break;
                case FBC_BOUNCE_BACK:
                    // bulk dynamics + velocity (inflow)
                    defineBulkDynamics(parameters.turbulenceModel_, lattice, superGeometry, indicator.id_);
                    setInterpolatedVelocityBoundary<T,DESCRIPTOR>(lattice, omega, superGeometry, indicator.id_);
                    break;
                case FBC_NONE:
                default:
                    lattice.defineDynamics<NoDynamics>(superGeometry, MAT_WALL);
                    break;
            }
        }
        else if(indicator.type_ == FIT_PRESSURE) {
            defineBulkDynamics(parameters.turbulenceModel_, lattice, superGeometry, indicator.id_);
            setInterpolatedPressureBoundary<T,DESCRIPTOR>(lattice, omega, superGeometry.getMaterialIndicator(indicator.id_));
        }
    }

    /*
    // If any dynamics should be missing up to this point, we fix this by using "no dynamics".
    size_t errors = fixupLattice(lattice);
    if(errors > 0) {
        std::cout << "[fixupLattice] " << errors << " errors have been fixed" << std::endl;
    }
    */

    // Initial conditions
    AnalyticalConst3D<T, T> rhoF(1);
    std::vector<T> velocity(3, T());
    AnalyticalConst3D<T, T> uF(velocity);

    lattice.defineRhoU(superGeometry.getMaterialIndicator(MAT_FLUID), rhoF, uF);
    lattice.iniEquilibrium(superGeometry.getMaterialIndicator(MAT_FLUID), rhoF, uF);

    // Initialize all values of distribution functions to their local equilibrium.
    for (const auto& indicator : indicators) {
        lattice.defineRhoU(superGeometry.getMaterialIndicator(indicator.id_), rhoF, uF);
        lattice.iniEquilibrium(superGeometry.getMaterialIndicator(indicator.id_), rhoF, uF);
    }

    // Add noise to fluid cells to enforce turbulent regime.
    if(parameters.latticePerturbation_) {
        LatticePerturber perturbation;
        lattice.defineU(superGeometry, MAT_FLUID, perturbation);
    }

    // Lattice initialize
    lattice.setParameter<descriptors::OMEGA>(omega);
    lattice.initialize();
}

void setBoundaryValues( SuperLattice<T, DESCRIPTOR>& lattice,
                        UnitConverter<T,DESCRIPTOR> const& converter,
                        int iteration,
                        SuperGeometry<T,3>& superGeometry,
                        const std::vector<FlowIndicator>& indicators,
                        const Parameters& parameters,
                        VolumeSampler volumeSampler = {}
                        )
{
    float time = converter.getPhysTime(iteration);

    for(const auto& indicator : indicators) {
        if (indicator.type_ == FIT_VELOCITY) {

            T targetPhysVelocity = indicator.velocityCurve_(time) * parameters.inletVelocityMultiplier_;
            T targetLatticeVelocity = converter.getLatticeVelocity(targetPhysVelocity);

            // This function applies the velocity profile to the boundary condition and the lattice.
            std::function<void(AnalyticalF3D<T,T>&)> applyFlowProfile = [&] (AnalyticalF3D<T,T>& profile) {
                if(parameters.wallBoundaryCondition_ == FBC_BOUZIDI) {
                    defineUBouzidi<T, DESCRIPTOR>(lattice, superGeometry, indicator.id_, profile);
                } else if(parameters.wallBoundaryCondition_ == FBC_BOUNCE_BACK) {
                    lattice.defineU(superGeometry.getMaterialIndicator(indicator.id_), profile);
                }
            };

            // Create shortcuts.
            const Vector<T, 3> center(indicator.center_.x, indicator.center_.y, indicator.center_.z);
            const Vector<T, 3> normal(indicator.normal_.x, indicator.normal_.y, indicator.normal_.z);
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
                case FP_VOLUME:
                {
                    // TODO: Find better solution, e.g., a sampler shall be returned and applied in this scope.
                    if(volumeSampler) {
                        // For volume indicators, we normalize the velocity curve to [0, 1]
                        // because multiplying the measurement by another velocity is confusing.
                        // We keep, however, the velocity multiplier, so that we can basically
                        // amplify the measurement.
                        auto multiplier = std::min<float>(targetPhysVelocity / indicator.velocityCurve_.getMaxVelocity(), 1);
                        volumeSampler(converter, time, applyFlowProfile, multiplier);
                    }
                    else {
                        std::cerr << "No Volume Sampler defined!" << std::endl;
                    }
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

void writeVVDFile(STLreader<T>& stlReader,
                  UnitConverter<T,DESCRIPTOR>& converter,
                  int iteration, int maxIteration,
                  int outputResolution,
                  const Parameters& parameters,
                  const std::string& simulationOutputPath,
                  const std::string& name,
                  SuperLatticeF3D<T, DESCRIPTOR>& feature) {

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
    switch (feature.getTargetDim()) {
        case 1:
            format = "float";
            break;
        case 3:
            format = "Vector3(float)";
            break;
        default:
            std::cerr << "Unhandled target dimensions" << std::endl;
            return;
    }

    std::vector<float> rawFeatureData(resolution * resolution * resolution * feature.getTargetDim());
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
                size_t index = z*resolution*resolution*feature.getTargetDim() + y*resolution*feature.getTargetDim() + x*feature.getTargetDim();
                for (int i = 0; i < feature.getTargetDim(); i++) {
                    rawFeatureData[index + i] = static_cast<float>(val[i]);
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
    std::string rawFilename = simulationOutputPath + featureFilename + ".raw";
    std::string vvdFilename = simulationOutputPath + featureFilename + ".vvd";

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
            << "<MetaItem name=\"" << META_DATA_NAME_OFFSET << "\" type=\"Vec3MetaData\">"
            << "<value x=\"" << offset[0] << "\" y=\"" << offset[1] << "\" z=\"" << offset[2] << "\" />"
            << "</MetaItem>"
            << "<MetaItem name=\"" << META_DATA_NAME_SPACING << "\" type=\"Vec3MetaData\">"
            << "<value x=\"" << spacing[0] << "\" y=\"" << spacing[1] << "\" z=\"" << spacing[2] << "\" />"
            << "</MetaItem>"
            << "<MetaItem name=\"" << META_DATA_NAME_TIMESTEP << "\" type=\"FloatMetaData\" value=\"" << converter.getPhysTime(iteration) << "\" />"
            << "<MetaItem name=\"" << META_DATA_NAME_REAL_WORLD_MAPPING << "\" type=\"RealWorldMappingMetaData\"><value scale=\"1\" offset=\"0\" unit=\"\" /></MetaItem>"
            << "<MetaItem name=\"" << META_DATA_NAME_MODALITY << "\" type=\"StringMetaData\" value=\"" << name << "\" />"
            // Parameters.
            << "<MetaItem name=\"" << "ParameterSpatialResolution" << "\" type=\"IntMetaData\" value=\"" << parameters.spatialResolution_ << "\" />"
            << "<MetaItem name=\"" << "ParameterRelaxationTime" << "\" type=\"FloatMetaData\" value=\"" << parameters.relaxationTime_ << "\" />"
            << "<MetaItem name=\"" << "ParameterCharacteristicLength" << "\" type=\"FloatMetaData\" value=\"" << parameters.characteristicLength_ << "\" />"
            << "<MetaItem name=\"" << "ParameterCharacteristicVelocity" << "\" type=\"FloatMetaData\" value=\"" << parameters.characteristicVelocity_ << "\" />"
            << "<MetaItem name=\"" << "ParameterViscosity" << "\" type=\"FloatMetaData\" value=\"" << parameters.viscosity_ << "\" />"
            << "<MetaItem name=\"" << "ParameterDensity" << "\" type=\"FloatMetaData\" value=\"" << parameters.density_ << "\" />"
            << "<MetaItem name=\"" << "ParameterTurbulenceModel" << "\" type=\"IntMetaData\" value=\"" << parameters.turbulenceModel_ << "\" />"
            << "<MetaItem name=\"" << "ParameterSmagorinskyConstant" << "\" type=\"FloatMetaData\" value=\"" << parameters.smagorinskyConstant_ << "\" />"
            << "<MetaItem name=\"" << "ParameterWallBoundaryCondition" << "\" type=\"IntMetaData\" value=\"" << parameters.wallBoundaryCondition_ << "\" />"
            << "<MetaItem name=\"" << "ParameterInletVelocityMultiplier" << "\" type=\"FloatMetaData\" value=\"" << parameters.inletVelocityMultiplier_ << "\" />"
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
        std::cerr << "Could not write " << name << " file" << std::endl;
    }
}

// Computes flux at inflow and outflow
bool getResults( SuperLattice<T, DESCRIPTOR>& lattice,
                 UnitConverter<T,DESCRIPTOR>& converter,
                 int iteration, int maxIteration,
                 SuperGeometry<T,3>& superGeometry,
                 STLreader<T>& stlReader,
                 const std::string& simulationOutputPath,
                 int numTimeSteps,
                 int outputResolution,
                 const std::string& outputFormat,
                 int flowFeatures,
                 const Parameters& parameters)
{
    const int outputIter = maxIteration / numTimeSteps;

    if (iteration % outputIter == 0) {

        bool writeVVD = outputFormat == ".vvd";
        bool writeVTI = outputFormat == ".vti";

        SuperVTMwriter3D<T> vtmWriter( "results" );

        // Always write debug data.
        if(iteration == 0) {
            SuperLatticeGeometry3D<T, DESCRIPTOR> geometry( lattice, superGeometry );
            vtmWriter.write( geometry );

            SuperLatticeCuboid3D<T, DESCRIPTOR> cuboid( lattice );
            vtmWriter.write( cuboid );

            SuperLatticeRank3D<T, DESCRIPTOR> rank( lattice );
            vtmWriter.write( rank );

            vtmWriter.createMasterFile();
        }

        if(flowFeatures & FF_VELOCITY) {
            SuperLatticePhysVelocity3D<T, DESCRIPTOR> velocity(lattice, converter);
            if(writeVVD) {
                writeVVDFile(stlReader, converter, iteration, maxIteration, outputResolution, parameters,
                             simulationOutputPath, "velocity", velocity);
            }
            if(writeVTI) {
                vtmWriter.write(velocity, iteration);
            }
        }

        if(flowFeatures & FF_MAGNITUDE) {
            SuperLatticePhysVelocity3D<T, DESCRIPTOR> velocity(lattice, converter);
            SuperEuklidNorm3D<T, DESCRIPTOR> magnitude(velocity);
            if(writeVVD) {
                writeVVDFile(stlReader, converter, iteration, maxIteration, outputResolution, parameters,
                             simulationOutputPath, "magnitude", magnitude);
            }
            if(writeVTI) {
                vtmWriter.write(magnitude, iteration);
            }
        }

        if(flowFeatures & FF_PRESSURE) {
            SuperLatticePhysPressure3D<T, DESCRIPTOR> pressure(lattice, converter);
            if(writeVVD) {
                writeVVDFile(stlReader, converter, iteration, maxIteration, outputResolution, parameters,
                             simulationOutputPath, "pressure", pressure);
            }
            if(writeVTI) {
                vtmWriter.write(pressure, iteration);
            }
        }

        if(flowFeatures & FF_WALLSHEARSTRESS) {
            SuperLatticePhysWallShearStress3D<T, DESCRIPTOR> wallShearStress(lattice, superGeometry, MAT_WALL,
                                                                             converter, stlReader);
            if(writeVVD) {
                writeVVDFile(stlReader, converter, iteration, maxIteration, outputResolution, parameters,
                             simulationOutputPath, "wallShearStress", wallShearStress);
            }
            if(writeVTI) {
                vtmWriter.write(wallShearStress, iteration);
            }
        }

        // Lattice statistics console output
        std::cout << "step="     << iteration << "; " <<
                     "t="        << converter.getPhysTime(iteration) << "; " <<
                     "uMax="     << lattice.getStatistics().getMaxU() << "; " <<
                     "avEnergy=" << lattice.getStatistics().getAverageEnergy() << "; " <<
                     "avRho="    << lattice.getStatistics().getAverageRho() << std::endl;
        lattice.getStatistics().print(iteration, converter.getPhysTime(iteration));
    }

    T tau = converter.getLatticeRelaxationFrequency();
    T threshold = tau < 0.55 ? 0.125*(tau - 0.5) : 0.4;
    if (lattice.getStatistics().getMaxU() >= threshold) {
        std::cerr << "uMax=" << lattice.getStatistics().getMaxU() << " above threshold=" << threshold << std::endl;
        return false;
    }

    return true;
}

}   // namespace
