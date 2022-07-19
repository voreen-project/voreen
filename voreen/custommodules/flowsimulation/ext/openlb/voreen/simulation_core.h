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

#include "openlb_parameters.h"

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

template<typename S>
struct SimpleVolume {
    std::vector<S> data;

    Vector<int, 3> dimensions;
    int numChannels;

    Vector<T, 3> offset;
    Vector<T, 3> spacing;

    std::vector<S> minValues, maxValues;
    S minMagnitude, maxMagnitude;

    SimpleVolume(Vector<int, 3> dimensions, int numChannels)
        : data(dimensions[0]*dimensions[1]*dimensions[2]*numChannels, S(0))
        , dimensions(dimensions)
        , numChannels(numChannels)
        , minValues(numChannels, std::numeric_limits<S>::max())
        , maxValues(numChannels, std::numeric_limits<S>::lowest())
        , minMagnitude(std::numeric_limits<S>::max())
        , maxMagnitude(0) {}

    void setValue(S value, int x, int y, int z, int channel = 0) {
        long index = z*dimensions[1]*dimensions[0]*numChannels + y*dimensions[0]*numChannels + x*numChannels + channel;
        data[index] = value;
    }
    void setValue(S value, int index, int channel = 0) {
        index = index*numChannels + channel;
        data[index] = value;
    }
    S getValue(int x, int y, int z, int channel = 0) const {
        long index = z*dimensions[1]*dimensions[0]*numChannels + y*dimensions[0]*numChannels + x*numChannels + channel;
        return data[index];
    }
    S getValue(int index, int channel = 0) const {
        index = index*numChannels + channel;
        return data[index];
    }

    S getValueLinear(float x, float y, float z, int channel = 0) const {

        x = std::clamp<float>(x, 0, dimensions[0]);
        y = std::clamp<float>(y, 0, dimensions[1]);
        z = std::clamp<float>(z, 0, dimensions[2]);

        vec3 llb{ std::floor(x), std::floor(y), std::floor(z) };
        vec3 urf{ std::ceil(x), std::ceil(y), std::ceil(z) };
        vec3 p{ x - llb.x, y - llb.y, z - llb.z };

        return    getValue(llb.x, llb.y, llb.z, channel) * (1.f-p.x)*(1.f-p.y)*(1.f-p.z) // llB
                + getValue(urf.x, llb.y, llb.z, channel) * (    p.x)*(1.f-p.y)*(1.f-p.z) // lrB
                + getValue(urf.x, urf.y, llb.z, channel) * (    p.x)*(    p.y)*(1.f-p.z) // urB
                + getValue(llb.x, urf.y, llb.z, channel) * (1.f-p.x)*(    p.y)*(1.f-p.z) // ulB
                + getValue(llb.x, llb.y, urf.z, channel) * (1.f-p.x)*(1.f-p.y)*(    p.z) // llF
                + getValue(urf.x, llb.y, urf.z, channel) * (    p.x)*(1.f-p.y)*(    p.z) // lrF
                + getValue(urf.x, urf.y, urf.z, channel) * (    p.x)*(    p.y)*(    p.z) // urF
                + getValue(llb.x, urf.y, urf.z, channel) * (1.f-p.x)*(    p.y)*(    p.z);// ulF
    }
};

using VolumeSampler = std::function<void(UnitConverter<T, DESCRIPTOR> const&, float, std::function<void(AnalyticalF3D<T,T>&)>&, float)>;

class ZeroAnalyticalF3D : public AnalyticalF3D<T, T> {
public:
    ZeroAnalyticalF3D() : AnalyticalF3D<T, T>(3) {}
    virtual bool operator() (T output[], const T input[]) {
        for(size_t i=0; i<3; i++) {
            output[i] = 0;
        }
        return true;
    }
} ZeroFunctor;

class LatticePerturber : public AnalyticalF3D<T, T> {
public:
    LatticePerturber(AnalyticalF3D<T, T>& original = ZeroFunctor, T maxNoise=1e-5)
        : AnalyticalF3D<T, T>(3)
        , original_(original)
        , rnd_(std::bind(std::uniform_real_distribution<T>(-maxNoise, maxNoise), std::mt19937(time(nullptr))))
    {
    }
    virtual bool operator() (T output[], const T input[]) {
        original_(output, input);
        for(size_t i=0; i<3; i++) {
            output[i] += rnd_();
        }
        return true;
    }

private:
    AnalyticalF3D<T, T>& original_;
    std::function<T()> rnd_;
};

// Stores data from stl file in geometry in form of material numbers
bool prepareGeometry(UnitConverter<T,DESCRIPTOR> const& converter,
                     IndicatorF3D<T>& indicator,
                     STLreader<T>& stlReader,
                     SuperGeometry<T,3>& superGeometry,
                     const std::vector<FlowIndicator>& indicators)
{
    superGeometry.rename(MAT_EMPTY, MAT_WALL, indicator);
    superGeometry.rename(MAT_WALL, MAT_FLUID, stlReader);

    superGeometry.clean();

    for (size_t i = 0; i < indicators.size(); i++) {

        Vector<T, 3> normal(indicators[i].normal_.x, indicators[i].normal_.y, indicators[i].normal_.z);

        // Convert to SI units.
        Vector<T, 3> center(indicators[i].center_.x, indicators[i].center_.y, indicators[i].center_.z);
        center *= VOREEN_LENGTH_TO_SI;

        // Add one voxel to account for precision/rounding errors.
        T radius = indicators[i].radius_ * VOREEN_LENGTH_TO_SI + converter.getConversionFactorLength() * 2;
        T length = indicators[i].length_ * VOREEN_LENGTH_TO_SI + converter.getConversionFactorLength() * 8;

        // Define a local disk volume.
        IndicatorCircle3D<T> flow(center, normal, radius);
        IndicatorCylinder3D<T> layerFlow(flow, length);

        // Rename both, wall and fluid, since the indicator might also be inside the fluid domain.
        //superGeometry.rename(MAT_WALL, indicators[i].id_, layerFlow);
        superGeometry.rename(MAT_FLUID, indicators[i].id_, layerFlow);

        // Exclude area behind inlet and in front of outlet - it will otherwise cause unstable simulations.
        bool isInlet = indicators[i].type_ == FIT_VELOCITY;
        bool isOutlet = indicators[i].type_ == FIT_PRESSURE;
        if(isInlet || isOutlet) {
            T sign = isInlet ? T(-1) : T(1);
            center += sign * normal * T(length * 0.5 + converter.getConversionFactorLength());

            IndicatorCircle3D<T> capFlowWall(center, normal, radius);
            IndicatorCylinder3D<T> layerCapFlowWall(capFlowWall, 4 * converter.getConversionFactorLength());
            superGeometry.rename(MAT_FLUID, MAT_WALL, layerCapFlowWall);

            IndicatorCircle3D<T> capFlowEmpty(center, normal, radius);
            IndicatorCylinder3D<T> layerCapFlowEmpty(capFlowEmpty, 4 * converter.getConversionFactorLength());
            superGeometry.rename(MAT_WALL, MAT_EMPTY, layerCapFlowEmpty);

//            // Fill cut-off region.
//            center += sign * normal * converter.getConversionFactorLength();
//            floodRegion(superGeometry, center[0], center[1], center[2], MAT_FLUID);
        }
    }

    // Removes all not needed boundary voxels outside the surface
    //superGeometry.clean();
    // Removes all not needed boundary voxels inside the surface
    superGeometry.innerClean(MAT_COUNT);

    // return !superGeometry.checkForErrors();
    return true; // TODO: do not ignore return value, actually clean geometry.
}

void defineBulkDynamics(FlowTurbulenceModel model,
                        SuperLattice<T, DESCRIPTOR>& lattice,
                        SuperGeometry<T,3>& geometry,
                        int material) {
    auto indicator = geometry.getMaterialIndicator(material);
    switch(model) {
        case FTM_SMAGORINSKY:
            lattice.defineDynamics<SmagorinskyBGKdynamics<T, DESCRIPTOR>>(indicator);
            return;
        case FTM_SMAGORINSKY_SHEAR_IMPROVED:
            lattice.defineDynamics<ShearSmagorinskyBGKdynamics<T, DESCRIPTOR>>(indicator);
            return;
        case FTM_SMAGORINSKY_CONSISTENT:
            lattice.defineDynamics<ConSmagorinskyBGKdynamics<T, DESCRIPTOR>>(indicator);
            return;
        case FTM_SMAGORINSKY_CONSISTENT_STRAIN:
            lattice.defineDynamics<ConStrainSmagorinskyBGKdynamics<T, DESCRIPTOR>>(indicator);
            return;
        case FTM_BGK:
            lattice.defineDynamics<BGKdynamics<T, DESCRIPTOR>>(indicator);
            return;
        case FTM_NONE:
        default:
            lattice.defineDynamics<NoDynamics>(indicator);
            return;
    }
}

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
            //setInterpolatedPressureBoundary<T,DESCRIPTOR>(lattice, omega, superGeometry.getMaterialIndicator(indicator.id_));
            setLocalPressureBoundary<T,DESCRIPTOR>(lattice, omega, superGeometry.getMaterialIndicator(indicator.id_));
        }
    }

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
    lattice.setParameter<collision::LES::Smagorinsky>(0.1);
    lattice.initialize();
}

void setBoundaryValues( SuperLattice<T, DESCRIPTOR>& lattice,
                        UnitConverter<T,DESCRIPTOR> const& converter,
                        int iteration,
                        SuperGeometry<T,3>& superGeometry,
                        const std::vector<FlowIndicator>& indicators,
                        const Parameters& parameters,
                        const VolumeSampler& volumeSampler = {}
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

template<typename S>
SimpleVolume<S> sampleVolume(IndicatorF3D<T>& indicator,
                             UnitConverter<T,DESCRIPTOR>& converter,
                             int maxOutputResolution,
                             SuperF3D<T, T>& feature) {

    AnalyticalFfromSuperF3D<T> interpolateFeature(feature, true);

    const Vector<T, 3>& min = indicator.getMin();
    const Vector<T, 3>& max = indicator.getMax();

    const Vector<T, 3> len = max - min;
    const T maxLen = std::max({len[0], len[1], len[2]});
    const int maxLenN = std::round(maxLen / converter.getConversionFactorLength());
    const T scaling = std::min<T>(maxOutputResolution, maxLenN) / maxLenN;

    const Vector<int, 3> gridResolution(std::round(len[0] / converter.getConversionFactorLength() * scaling)+1,
                                        std::round(len[1] / converter.getConversionFactorLength() * scaling)+1,
                                        std::round(len[2] / converter.getConversionFactorLength() * scaling)+1);

    // Init volume.
    SimpleVolume<S> volume(gridResolution, feature.getTargetDim());

    volume.offset = min;
    volume.spacing = converter.getConversionFactorLength() / scaling;

    for (int z = 0; z < volume.dimensions[2]; z++) {
        for (int y = 0; y < volume.dimensions[1]; y++) {
            for (int x = 0; x < volume.dimensions[0]; x++) {

                T pos[3] = {volume.offset[0] + x * volume.spacing[0],
                            volume.offset[1] + y * volume.spacing[1],
                            volume.offset[2] + z * volume.spacing[2]};
                std::vector<T> val(volume.numChannels, 0.0f);

                // Retrieve data.
                interpolateFeature(val.data(), pos);

                // Update min/max.
                T magnitude = 0;
                for (int i = 0; i < feature.getTargetDim(); i++) {
                    volume.minValues[i] = std::min<S>(volume.minValues[i], val[i]);
                    volume.maxValues[i] = std::max<S>(volume.maxValues[i], val[i]);
                    magnitude += val[i] * val[i];
                }

                // Update min/max magnitude.
                volume.minMagnitude = std::min<S>(volume.minMagnitude, magnitude);
                volume.maxMagnitude = std::max<S>(volume.maxMagnitude, magnitude);

                // Downgrade to required type.
                for (int i = 0; i < volume.numChannels; i++) {
                    T value = val[i];
                    if(std::numeric_limits<S>::is_integer) {
                        value = std::round(value);
                    }
                    volume.setValue(static_cast<S>(value), x, y, z, i);
                }
            }
        }
    }

    // Adapt to Voreen units.
    volume.offset  *= (1/VOREEN_LENGTH_TO_SI);
    volume.spacing *= (1/VOREEN_LENGTH_TO_SI);

    volume.minMagnitude = std::sqrt(volume.minMagnitude);
    volume.maxMagnitude = std::sqrt(volume.maxMagnitude);

    return volume;
}

void writeVVDFile(IndicatorF3D<T>& indicator,
                  UnitConverter<T,DESCRIPTOR>& converter,
                  int iteration, int maxIteration,
                  int outputResolution,
                  const Parameters& parameters,
                  const std::string& simulationOutputPath,
                  const std::string& name,
                  SuperLatticeF3D<T, DESCRIPTOR>& feature) {

    SimpleVolume<float> volume = sampleVolume<float>(indicator, converter, outputResolution, feature);

    // Determine format.
    // This could be done in a more dynamic way, but the code should be easily portable to the cluster.
    std::string format;
    switch (volume.numChannels) {
        case 1:
            format = "float";
            break;
        case 3:
            format = "Vector3(float)";
            break;
        default:
            std::cerr << "Unhandled number of channels" << std::endl;
            return;
    }

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
            << "<RawData filename=\"" << featureFilename << ".raw\" format=\"" << format << "\" x=\"" << volume.dimensions[0]
            << "\" y=\"" << volume.dimensions[1] << "\" z=\"" << volume.dimensions[2] << "\" />"
            // Mandatory Meta data.
            << "<MetaData>"
            << "<MetaItem name=\"" << META_DATA_NAME_OFFSET << "\" type=\"Vec3MetaData\">"
            << "<value x=\"" << volume.offset[0] << "\" y=\"" << volume.offset[1] << "\" z=\"" << volume.offset[2] << "\" />"
            << "</MetaItem>"
            << "<MetaItem name=\"" << META_DATA_NAME_SPACING << "\" type=\"Vec3MetaData\">"
            << "<value x=\"" << volume.spacing[0] << "\" y=\"" << volume.spacing[1] << "\" z=\"" << volume.spacing[2] << "\" />"
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
    if (volume.numChannels > 1) {
        vvdFeatureFile << "<DerivedItem type=\"VolumeMinMaxMagnitude\" minMagnitude=\"" << volume.minMagnitude << "\" maxMagnitude=\"" << volume.maxMagnitude << "\" minNormalizedMagnitude=\"" << volume.minMagnitude << "\" maxNormalizedMagnitude=\"" << volume.maxMagnitude << "\" />";
    }
    // * VolumeMinMax
    vvdFeatureFile << "<DerivedItem type=\"VolumeMinMax\"><minValues>";
    for(int i=0; i<volume.numChannels; i++) {
        vvdFeatureFile << "<channel value=\"" << volume.minValues[i] << "\" />";
    }
    vvdFeatureFile << "</minValues><maxValues>";
    for(int i=0; i<volume.numChannels; i++) {
        vvdFeatureFile << "<channel value=\"" << volume.maxValues[i] << "\" />";
    }
    vvdFeatureFile << "</maxValues><minNormValues>";
    for(int i=0; i<volume.numChannels; i++) {
        vvdFeatureFile << "<channel value=\"" << volume.minValues[i] << "\" />";
    }
    vvdFeatureFile << "</minNormValues><maxNormValues>";
    for(int i=0; i<volume.numChannels; i++) {
        vvdFeatureFile << "<channel value=\"" << volume.maxValues[i] << "\" />";
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
    size_t numBytes = volume.data.size() * sizeof(float) / sizeof(char);
    rawFeatureFile.write(reinterpret_cast<const char*>(volume.data.data()), numBytes);
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
#ifndef VRN_MODULE_FLOWSIMULATION
        if(iteration == 0) {
            SuperLatticeGeometry3D<T, DESCRIPTOR> geometry( lattice, superGeometry );
            vtmWriter.write( geometry );

            SuperLatticeCuboid3D<T, DESCRIPTOR> cuboid( lattice );
            vtmWriter.write( cuboid );

            SuperLatticeRank3D<T, DESCRIPTOR> rank( lattice );
            vtmWriter.write( rank );

            vtmWriter.createMasterFile();
        }
#endif

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
