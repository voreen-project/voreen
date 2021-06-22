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

#include <olb3D.h>
#ifndef OLB_PRECOMPILED
#include "olb3D.hh"
#endif

using namespace olb;
using namespace olb::descriptors;
typedef double T;
typedef D3Q19<> DESCRIPTOR;

namespace {

static const T VOREEN_LENGTH_TO_SI = 0.001;
static const T VOREEN_TIME_TO_SI = 0.001;

// Static material ids (as defined by OpenLB).
enum Material {
    MAT_EMPTY = 0,
    MAT_FLUID = 1,
    MAT_WALL = 2,
    MAT_COUNT,
};

// Allow for simulation lattice initialization by volume data.
class MeasuredDataMapper : public AnalyticalF3D<T, T> {
public:
    MeasuredDataMapper(const voreen::VolumeBase* volume)
            : AnalyticalF3D<T, T>(3)
            , volume_(volume)
    {
        tgtAssert(volume_, "No volume");
        tgtAssert(volume_->getNumChannels() == 3, "Num channels != 3");
        bounds_ = volume_->getBoundingBox().getBoundingBox();
        representation_.reset(new voreen::VolumeRAMRepresentationLock(volume_));
        worldToVoxelMatrix_ = volume_->getWorldToVoxelMatrix();
    }
    virtual bool operator() (T output[], const T input[]) {
        tgt::vec3 pos = tgt::Vector3<T>::fromPointer(input);
        if (!bounds_.containsPoint(pos)) {
            return false;
        }

        pos = worldToVoxelMatrix_ * pos;

        for(size_t i=0; i < (**representation_)->getNumChannels(); i++) {
            output[i] = (**representation_)->getVoxelNormalized(pos, i);
        }

        return true;
    }

private:
    const voreen::VolumeBase* volume_;
    std::unique_ptr<voreen::VolumeRAMRepresentationLock> representation_;
    tgt::Bounds bounds_;
    tgt::mat4 worldToVoxelMatrix_;
};


// Stores data from stl file in geometry in form of material numbers
void prepareGeometry(UnitConverter<T,DESCRIPTOR> const& converter,
                     IndicatorF3D<T>& indicator,
                     STLreader<T>& stlReader,
                     SuperGeometry3D<T>& superGeometry,
                     const voreen::FlowParameterSetEnsemble& parametrizationList,
                     size_t selectedParametrization)
{
    superGeometry.rename( MAT_EMPTY, MAT_WALL,  indicator );
    superGeometry.rename( MAT_WALL,  MAT_FLUID, stlReader );

    superGeometry.clean();

    auto flowIndicators = parametrizationList.getFlowIndicators();
    for (size_t i = 0; i < flowIndicators.size(); i++) {

        tgt::Vector3<T> normal = flowIndicators[i].normal_;

        // Convert to SI units.
        tgt::Vector3<T> center = flowIndicators[i].center_;
        center *= VOREEN_LENGTH_TO_SI;

        // Add one voxel to account for rounding errors.
        T radius = flowIndicators[i].radius_ * VOREEN_LENGTH_TO_SI + converter.getConversionFactorLength();

        // Define a local disk volume.
        IndicatorCircle3D<T> flow(center[0], center[1], center[2],
                                  normal[0], normal[1], normal[2],
                                  radius);
        IndicatorCylinder3D<T> layerFlow(flow, 2 * converter.getConversionFactorLength());

        // Rename both, wall and fluid, since the indicator might also be inside the fluid domain.
        superGeometry.rename(MAT_WALL, flowIndicators[i].id_, MAT_FLUID, layerFlow);
        superGeometry.rename(MAT_FLUID, flowIndicators[i].id_,  layerFlow);

        // Exclude area behind inlet and in front of outlet - it will otherwise cause unstable simulations.
        bool isInlet = flowIndicators[i].type_ == voreen::FIT_VELOCITY;
        bool isOutlet = flowIndicators[i].type_ == voreen::FIT_PRESSURE;
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
}

// Set up the geometry of the simulation.
void prepareLattice( SuperLattice3D<T, DESCRIPTOR>& lattice,
                     UnitConverter<T,DESCRIPTOR> const& converter,
                     Dynamics<T, DESCRIPTOR>& bulkDynamics,
                     STLreader<T>& stlReader,
                     SuperGeometry3D<T>& superGeometry,
                     const voreen::VolumeList* measuredData,
                     const voreen::FlowParameterSetEnsemble& parameterSetEnsemble,
                     size_t selectedParametrization)
{
    const bool bouzidiOn = parameterSetEnsemble.at(selectedParametrization).getBouzidi();
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

    for(const auto& indicator : parameterSetEnsemble.getFlowIndicators()) {
        if(indicator.type_ == voreen::FIT_VELOCITY) {
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
        else if(indicator.type_ == voreen::FIT_PRESSURE) {
            lattice.defineDynamics(superGeometry.getMaterialIndicator(indicator.id_), &bulkDynamics);
            setInterpolatedPressureBoundary<T,DESCRIPTOR>(lattice, omega, superGeometry.getMaterialIndicator(indicator.id_));
        }
    }

    // Unsteered simulation.
    if(!measuredData) {

        // Initial conditions
        AnalyticalConst3D<T, T> rhoF(1);
        std::vector<T> velocity(3, T(0));
        AnalyticalConst3D<T, T> uF(velocity);

        lattice.defineRhoU(superGeometry.getMaterialIndicator(MAT_FLUID), rhoF, uF);
        lattice.iniEquilibrium(superGeometry.getMaterialIndicator(MAT_FLUID), rhoF, uF);

        // Initialize all values of distribution functions to their local equilibrium.
        for (const auto& indicator : parameterSetEnsemble.getFlowIndicators()) {
            lattice.defineRhoU(superGeometry.getMaterialIndicator(indicator.id_), rhoF, uF);
            lattice.iniEquilibrium(superGeometry.getMaterialIndicator(indicator.id_), rhoF, uF);
        }
    }
    // Steered simulation - currently only initializes the first time step!
    else {
        MeasuredDataMapper mapper(measuredData->first());
        lattice.defineU(superGeometry, MAT_FLUID, mapper);
    }

    // Lattice initialize
    lattice.initialize();
}

// Generates a slowly increasing sinuidal inflow.
void setBoundaryValues( SuperLattice3D<T, DESCRIPTOR>& sLattice,
                        UnitConverter<T,DESCRIPTOR> const& converter,
                        int iteration,
                        SuperGeometry3D<T>& superGeometry,
                        const voreen::FlowParameterSetEnsemble& parameterSetEnsemble,
                        size_t selectedParametrization)
{

    bool bouzidiOn = parameterSetEnsemble.at(selectedParametrization).getBouzidi();

    for(const auto& indicator : parameterSetEnsemble.getFlowIndicators()) {
        if (indicator.type_ == voreen::FIT_VELOCITY) {

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
            const tgt::vec3& center = indicator.center_;
            const tgt::vec3& normal = indicator.normal_;
            T radius = indicator.radius_ * VOREEN_LENGTH_TO_SI;

            // Apply the indicator's profile.
            switch(indicator.flowProfile_) {
            case voreen::FP_POISEUILLE:
            {
//                CirclePoiseuille3D<T> profile(superGeometry, indicator.id_, targetLatticeVelocity); // This is the alternative way, but how does it work?
                CirclePoiseuille3D<T> profile(center[0]*VOREEN_LENGTH_TO_SI, center[1]*VOREEN_LENGTH_TO_SI, center[2]*VOREEN_LENGTH_TO_SI,
                                              normal[0], normal[1], normal[2], radius, targetLatticeVelocity);
                applyFlowProfile(profile);
                break;
            }
            case voreen::FP_POWERLAW:
            {
                T n = 1.03 * std::log(converter.getReynoldsNumber()) - 3.6; // Taken from OLB documentation.
                CirclePowerLawTurbulent3D<T> profile(center[0]*VOREEN_LENGTH_TO_SI, center[1]*VOREEN_LENGTH_TO_SI, center[2]*VOREEN_LENGTH_TO_SI,
                                                     normal[0], normal[1], normal[2], radius, targetLatticeVelocity, n);
                applyFlowProfile(profile);
                break;
            }
            case voreen::FP_CONSTANT:
            {
                AnalyticalConst3D<T, T> profile(normal[0] * targetLatticeVelocity, normal[1] * targetLatticeVelocity, normal[2] * targetLatticeVelocity);
                applyFlowProfile(profile);
                break;
            }
            case voreen::FP_NONE:
            default:
                // Skip!
                continue;
            }
        }
    }
}

}


namespace voreen {

void writeVVDFile(STLreader<T>& stlReader,
                  UnitConverter<T,DESCRIPTOR>& converter,
                  int iteration, int maxIteration,
                  const FlowParameterSetEnsemble& parameterSetEnsemble,
                  size_t selectedParametrization,
                  std::string simulationOutputPath,
                  const std::string& name,
                  SuperLatticeF3D<T, DESCRIPTOR>& feature) {

    // Create output directory.
    simulationOutputPath += name + '/';
    if(!tgt::FileSystem::dirExists(simulationOutputPath)) {
        tgt::FileSystem::createDirectory(simulationOutputPath);
    }

    const Vector<T, 3>& min = stlReader.getMin();
    const Vector<T, 3>& max = stlReader.getMax();

    const Vector<T, 3> len = (max - min);
    const T maxLen = std::max({len[0], len[1], len[2]});
    const int gridResolution = static_cast<int>(std::round(maxLen / converter.getConversionFactorLength()));
    const int resolution = std::min(parameterSetEnsemble.getOutputResolution(), gridResolution);

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
            tgtAssert(false,"Unhandled target dimensions");
            return;
    }

    std::vector<float> rawFeatureData(resolution * resolution * resolution * feature.getTargetDim());
    AnalyticalFfromSuperF3D<T> interpolateFeature(feature, true);

    std::vector<T> minValue(feature.getTargetDim(), std::numeric_limits<T>::max());
    std::vector<T> maxValue(feature.getTargetDim(), std::numeric_limits<T>::lowest());

    T minMagnitude = std::numeric_limits<T>::max();
    T maxMagnitude = 0;

#ifdef VRN_MODULE_OPENMP
#pragma omp parallel for
#endif
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

    const FlowParameterSet& parameters = parameterSetEnsemble.at(selectedParametrization);
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
            << "<MetaItem name=\"" << VolumeBase::META_DATA_NAME_TIMESTEP << "\" type=\"FloatMetaData\" value=\"" << converter.getPhysTime(iteration) << "\" />"
            << "<MetaItem name=\"" << VolumeBase::META_DATA_NAME_REAL_WORLD_MAPPING << "\" type=\"RealWorldMappingMetaData\"><value scale=\"1\" offset=\"0\" unit=\"\" /></MetaItem>"
            << "<MetaItem name=\"" << VolumeBase::META_DATA_NAME_MODALITY << "\" type=\"StringMetaData\" value=\"" << name << "\" />"
            // Parameters.
            << "<MetaItem name=\"" << "ParameterSpatialResolution" << "\" type=\"IntMetaData\" value=\"" << parameters.getSpatialResolution() << "\" />"
            << "<MetaItem name=\"" << "ParameterRelaxationTime" << "\" type=\"FloatMetaData\" value=\"" << parameters.getRelaxationTime() << "\" />"
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
        LERRORC("voreen.flowsimulation.FlowSimulation","Could not write " << name << " file");
    }
}

// Computes flux at inflow and outflow
bool getResults( SuperLattice3D<T, DESCRIPTOR>& sLattice,
                                 UnitConverter<T,DESCRIPTOR>& converter,
                                 int iteration, int maxIteration,
                                 Dynamics<T, DESCRIPTOR>& bulkDynamics,
                                 SuperGeometry3D<T>& superGeometry,
                                 STLreader<T>& stlReader,
                                 const FlowParameterSetEnsemble& parameterSetEnsemble,
                                 size_t selectedParametrization,
                                 const std::string& simulationOutputPath)
{
    const int outputIter = maxIteration / parameterSetEnsemble.getNumTimeSteps();

    if (iteration % outputIter == 0) {

        bool writeVVD = parameterSetEnsemble.getOutputFileFormat() == ".vvd";
        bool writeVTI = parameterSetEnsemble.getOutputFileFormat() == ".vti";

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

        if(parameterSetEnsemble.getFlowFeatures() & FF_VELOCITY) {
            SuperLatticePhysVelocity3D<T, DESCRIPTOR> velocity(sLattice, converter);
            if(writeVVD) {
                writeVVDFile(stlReader, converter, iteration, maxIteration, parameterSetEnsemble, selectedParametrization,
                             simulationOutputPath, "velocity", velocity);
            }
            if(writeVTI) {
                vtmWriter.write(velocity, iteration);
            }
        }

        if(parameterSetEnsemble.getFlowFeatures() & FF_MAGNITUDE) {
            SuperLatticePhysVelocity3D<T, DESCRIPTOR> velocity(sLattice, converter);
            SuperEuklidNorm3D<T, DESCRIPTOR> magnitude(velocity);
            if(writeVVD) {
                writeVVDFile(stlReader, converter, iteration, maxIteration, parameterSetEnsemble, selectedParametrization,
                             simulationOutputPath, "magnitude", magnitude);
            }
            if(writeVTI) {
                vtmWriter.write(magnitude, iteration);
            }
        }

        if(parameterSetEnsemble.getFlowFeatures() & FF_PRESSURE) {
            SuperLatticePhysPressure3D<T, DESCRIPTOR> pressure(sLattice, converter);
            if(writeVVD) {
                writeVVDFile(stlReader, converter, iteration, maxIteration, parameterSetEnsemble, selectedParametrization,
                             simulationOutputPath, "pressure", pressure);
            }
            if(writeVTI) {
                vtmWriter.write(pressure, iteration);
            }
        }

#ifndef OLB_PRECOMPILED
        if(parameterSetEnsemble.getFlowFeatures() & FF_WALLSHEARSTRESS) {
            SuperLatticePhysWallShearStress3D<T, DESCRIPTOR> wallShearStress(sLattice, superGeometry, MAT_WALL,
                                                                             converter, stlReader);
            if(writeVVD) {
                writeVVDFile(stlReader, converter, iteration, maxIteration, parameterSetEnsemble, selectedParametrization,
                             simulationOutputPath, "wallShearStress", wallShearStress);
            }
            if(writeVTI) {
                vtmWriter.write(wallShearStress, iteration);
            }
        }
#endif

        // Lattice statistics console output
        LINFOC("voreen.flowsimulation.FlowSimulation","step="     << iteration << "; " <<
                          "t="        << converter.getPhysTime(iteration) << "; " <<
                          "uMax="     << sLattice.getStatistics().getMaxU() << "; " <<
                          "avEnergy=" << sLattice.getStatistics().getAverageEnergy() << "; " <<
                          "avRho="    << sLattice.getStatistics().getAverageRho()
        );
        sLattice.getStatistics().print(iteration, converter.getPhysTime(iteration));
    }

    T tau = converter.getLatticeRelaxationFrequency();
    T threshold = tau < 0.55 ? 0.125*(tau - 0.5) : 0.4;
    if (sLattice.getStatistics().getMaxU() >= threshold) {
        LERRORC("voreen.flowsimulation.FlowSimulation", "uMax=" << sLattice.getStatistics().getMaxU() << " above threshold=" << threshold);
        return false;
    }

    return true;
}


const std::string FlowSimulation::loggerCat_("voreen.flowsimulation.FlowSimulation");

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
    //addPort(measuredDataPort_); // Currently ignored.
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

    // Note: measuredDataPort is optional!
    if(measuredDataPort_.hasData() && !measuredDataPort_.isReady()) {
        setNotReadyErrorMessage("Measured Data Port not ready.");
        return false;
    }

    if(!parameterPort_.isReady()) {
        setNotReadyErrorMessage("Parameter Port not ready.");
        return false;
    }

    // Note: ensemblePort is optional!

    return true;
}

void FlowSimulation::adjustPropertiesToInput() {
    const FlowParameterSetEnsemble* flowParameterSetEnsemble = parameterPort_.getData();
    if(!flowParameterSetEnsemble || flowParameterSetEnsemble->empty()) {
        selectedParametrization_.setMinValue(-1);
        selectedParametrization_.setMaxValue(-1);
    }
    else {
        selectedParametrization_.setMinValue(0);
        selectedParametrization_.setMaxValue(static_cast<int>(flowParameterSetEnsemble->size()) - 1);
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
    auto flowParameterSetEnsemble = parameterPort_.getThreadSafeData();
    if(!flowParameterSetEnsemble || flowParameterSetEnsemble->empty()) {
        throw InvalidInputException("No parameterization", InvalidInputException::S_ERROR);
    }

    if(flowParameterSetEnsemble->getFlowFeatures() == FF_NONE) {
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

    std::string simulationPath = simulationResults_.get() + "/" + flowParameterSetEnsemble->getName() + "/";
    if (!tgt::FileSystem::createDirectoryRecursive(simulationPath)) {
        throw InvalidInputException("Output directory could not be created", InvalidInputException::S_ERROR);
    }

    size_t selectedParametrization = FlowParameterSetEnsemble::ALL_PARAMETER_SETS;
    if(!simulateAllParametrizations_.get()) {
        selectedParametrization = static_cast<size_t>(selectedParametrization_.get());
    }

    return FlowSimulationInput{
            geometryPath,
            measuredData,
            flowParameterSetEnsemble,
            selectedParametrization,
            simulationPath,
            deleteOldSimulations_.get()
    };
}

FlowSimulationOutput FlowSimulation::compute(FlowSimulationInput input, ProgressReporter& progressReporter) const {

    // Needs to be initialized in each new thread to be used.
    olb::olbInit(nullptr, nullptr);

    // Run either all or just a single simulation.
    if(input.selectedParametrization == FlowParameterSetEnsemble::ALL_PARAMETER_SETS) {
        progressReporter.setProgress(0.0f);
        size_t numRuns = input.parameterSetEnsemble->size();
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
    const FlowParameterSetEnsemble& parameterSetEnsemble = *input.parameterSetEnsemble;
    const FlowParameterSet& parameters = parameterSetEnsemble.at(input.selectedParametrization);

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

    singleton::directories().setOutputDir(simulationResultPath);

    const int N = parameters.getSpatialResolution();
    UnitConverterFromResolutionAndRelaxationTime<T, DESCRIPTOR> converter(
            (T) N, // Resolution that charPhysLength is resolved by.
            (T) parameters.getRelaxationTime(), // Relaxation time
            (T) parameters.getCharacteristicLength(),         // charPhysLength: reference length of simulation geometry
            (T) parameters.getCharacteristicVelocity(),       // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
            (T) parameters.getViscosity(),                    // physViscosity: physical kinematic viscosity in __m^2 / s__
            (T) parameters.getDensity()                       // physDensity: physical density in __kg / m^3__
    );

    // Prints the converter log as console output
    converter.print();

    // === 2nd Step: Prepare Geometry ===

    // Instantiation of the STLreader class
    // file name, voxel size in meter, stl unit in meter, outer voxel no., inner voxel no.
    STLreader<T> stlReader(input.geometryPath, converter.getConversionFactorLength(), VOREEN_LENGTH_TO_SI, 1);
    IndicatorLayer3D<T> extendedDomain(stlReader, converter.getConversionFactorLength());

    // Instantiation of a cuboidGeometry with weights
    const int noOfCuboids = 1;
    CuboidGeometry3D<T> cuboidGeometry(extendedDomain, converter.getConversionFactorLength(), noOfCuboids);

    // Instantiation of a loadBalancer
    HeuristicLoadBalancer<T> loadBalancer(cuboidGeometry);

    interruptionPoint();

    // Instantiation of a superGeometry
    LINFO("Preparing Geometry ...");
    SuperGeometry3D<T> superGeometry(cuboidGeometry, loadBalancer, 2);

    prepareGeometry(converter, extendedDomain, stlReader, superGeometry,
                    parameterSetEnsemble, input.selectedParametrization);

    interruptionPoint();

    // === 3rd Step: Prepare Lattice ===
    LINFO("Preparing Lattice ...");
    SuperLattice3D<T, DESCRIPTOR> sLattice(superGeometry);

    SmagorinskyBGKdynamics<T, DESCRIPTOR> bulkDynamics(converter.getLatticeRelaxationFrequency(),
                                                       instances::getBulkMomenta<T, DESCRIPTOR>(),
                                                       parameters.getSmagorinskyConstant());

    prepareLattice(sLattice, converter, bulkDynamics,
                   stlReader, superGeometry,
                   measuredData,
                   parameterSetEnsemble,
                   input.selectedParametrization);

    interruptionPoint();

    // === 4th Step: Main Loop  ===
    const int maxIteration = converter.getLatticeTime(parameterSetEnsemble.getSimulationTime());
    util::ValueTracer<T> converge( converter.getLatticeTime(0.5), 1e-5);
    for (int iteration = 0; iteration <= maxIteration; iteration++) {

        // === 5th Step: Definition of Initial and Boundary Conditions ===
        setBoundaryValues(sLattice, converter, iteration, superGeometry, parameterSetEnsemble, input.selectedParametrization);

        // === 6th Step: Collide and Stream Execution ===
        sLattice.collideAndStream();

        // === 7th Step: Computation and Output of the Results ===
        bool success = getResults(sLattice, converter, iteration, maxIteration, bulkDynamics, superGeometry, stlReader,
                                  parameterSetEnsemble, input.selectedParametrization, simulationResultPath);
        if(!success) {
            break;
        }

        // === 8th Step: Check for convergence.
        converge.takeValue(sLattice.getStatistics().getAverageEnergy(), true);
        if(converge.hasConverged()) {
            LINFO("Simulation converged!");
            break;
        }

        float progress = iteration / (maxIteration + 1.0f);
        progressReporter.setProgress(progress);
    }
    progressReporter.setProgress(1.0f);
    LINFO("Finished simulation run: " << parameters.getName());
}

}   // namespace
