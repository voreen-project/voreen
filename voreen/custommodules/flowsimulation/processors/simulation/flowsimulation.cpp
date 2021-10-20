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

#include "voreen/core/datastructures/geometry/geometrysequence.h"
#include "voreen/core/datastructures/geometry/glmeshgeometry.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/datastructures/volume/volumeminmaxmagnitude.h"
#include "voreen/core/ports/conditions/portconditionvolumelist.h"

#include <olb3D.h>
#ifndef OLB_PRECOMPILED
#include "olb3D.hh"
#endif

#include <thread>

using namespace olb;
using namespace olb::descriptors;
typedef double T;
typedef D3Q19<> DESCRIPTOR;

namespace {

using namespace voreen;

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
    MeasuredDataMapper(UnitConverter<T, DESCRIPTOR> const& converter, const VolumeBase* volume, T multiplier = 1.0f)
        : AnalyticalF3D<T, T>(3)
        , converter_(converter)
        , volume_(volume)
        , multiplier_(multiplier)
    {
        tgtAssert(volume_, "No volume");
        tgtAssert(volume_->getNumChannels() == 3, "Num channels != 3");
        bounds_ = volume_->getBoundingBox().getBoundingBox();
        representation_.reset(new VolumeRAMRepresentationLock(volume_));
        worldToVoxelMatrix_ = volume_->getWorldToVoxelMatrix();
        rwm_ = volume->getRealWorldMapping();
    }
    virtual bool operator() (T output[], const T input[]) {
        // Store simulation positions in world coordinates.
        tgt::vec3 pos = tgt::Vector3<T>::fromPointer(input) / VOREEN_LENGTH_TO_SI;
        if (!bounds_.containsPoint(pos)) {
            return false;
        }

        // Convert to voxel coordinates to perform the lookup.
        pos = worldToVoxelMatrix_ * pos;

        for(size_t i=0; i < (**representation_)->getNumChannels(); i++) {
            output[i] = (**representation_)->getVoxelNormalizedLinear(pos, i);
            output[i] = rwm_.normalizedToRealWorld(output[i]);
            output[i] = converter_.getLatticeVelocity(output[i]);
            output[i] = output[i] * multiplier_;
        }

        return true;
    }

private:
    const UnitConverter<T, DESCRIPTOR>& converter_;
    const VolumeBase* volume_;
    const T multiplier_;
    std::unique_ptr<VolumeRAMRepresentationLock> representation_;
    tgt::Bounds bounds_;
    tgt::mat4 worldToVoxelMatrix_;
    RealWorldMapping rwm_;
};

template<typename T>
struct TypedGeometryPair {

    TypedGeometryPair(const Geometry* lhs, const Geometry* rhs)
        : lhs_(dynamic_cast<const T*>(lhs)), rhs_(dynamic_cast<const T*>(rhs)) {}

    operator bool () {
        return lhs_ != nullptr && rhs_ != nullptr;
    }

    const T* lhs_;
    const T* rhs_;
};

template<typename T>
std::unique_ptr<GlMeshGeometryBase> mergeGeometriesTyped(const T* lhs, const T* rhs) {
    tgtAssert(lhs, "lhs null");
    tgtAssert(rhs, "rhs null");

    auto merged = lhs->clone().release();
    auto mergedTyped = static_cast<T*>(merged);

    // Transform all vertex position to world space.
    auto vertices = mergedTyped->getVertices();
    for(auto& vertex : vertices) {
        vertex.pos_ = mergedTyped->getTransformationMatrix() * vertex.pos_;
    }
    mergedTyped->setTransformationMatrix(tgt::mat4::identity);

    // Copy over and transform vertices.
    for(auto vertex : rhs->getVertices()) {
        vertex.pos_ = rhs->getTransformationMatrix() * vertex.pos_;
        vertices.emplace_back(vertex);
    }
    mergedTyped->setVertices(vertices);

    // Add and adjust indices, if required.
    if (lhs->usesIndexedDrawing()) {
        auto indexOffset = lhs->getNumVertices();
        for (auto index : rhs->getIndices()) {
            mergedTyped->addIndex(index + indexOffset);
        }
    }

    return std::unique_ptr<GlMeshGeometryBase>(mergedTyped);
}


std::unique_ptr<GlMeshGeometryBase> mergeGeometries(const GlMeshGeometryBase* lhs, const GlMeshGeometryBase* rhs) {
    if(lhs->getPrimitiveType() != rhs->getPrimitiveType()) {
        return nullptr;
    }

    if(lhs->getVertexLayout() != rhs->getVertexLayout()) {
        return nullptr;
    }

    if(lhs->getIndexType() != rhs->getIndexType()) {
        return nullptr;
    }

    if(lhs->usesIndexedDrawing() != rhs->usesIndexedDrawing()) {
        return nullptr;
    }

    if(auto typedGeometry = TypedGeometryPair<GlMeshGeometryUInt32Normal>(lhs, rhs)) {
        return mergeGeometriesTyped(typedGeometry.lhs_, typedGeometry.rhs_);
    }

    return nullptr;
}



// Stores data from stl file in geometry in form of material numbers
void prepareGeometry(UnitConverter<T,DESCRIPTOR> const& converter,
                     IndicatorF3D<T>& indicator,
                     STLreader<T>& stlReader,
                     SuperGeometry3D<T>& superGeometry,
                     const FlowParameterSetEnsemble& parametrizationList,
                     const FlowParameterSet& parameters)
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

// Set up the geometry of the simulation.
void prepareLattice( SuperLattice3D<T, DESCRIPTOR>& lattice,
                     UnitConverter<T,DESCRIPTOR> const& converter,
                     Dynamics<T, DESCRIPTOR>& bulkDynamics,
                     STLreader<T>& stlReader,
                     SuperGeometry3D<T>& superGeometry,
                     const VolumeList* measuredData,
                     const FlowParameterSetEnsemble& parameterSetEnsemble,
                     const FlowParameterSet& parameters)
{
    FlowBoundaryCondition wallBoundaryCondition = parameters.getWallBoundaryCondition();
    const T omega = converter.getLatticeRelaxationFrequency();

    // material=0 --> do nothing
    lattice.defineDynamics(superGeometry.getMaterialIndicator(MAT_EMPTY), &instances::getNoDynamics<T, DESCRIPTOR>());

    // material=1 --> bulk dynamics
    lattice.defineDynamics(superGeometry.getMaterialIndicator(MAT_FLUID), &bulkDynamics);

    switch(wallBoundaryCondition) {
    case FBC_BOUZIDI:
        // material=2 --> no dynamics + bouzidi zero velocity
        lattice.defineDynamics(superGeometry.getMaterialIndicator(MAT_WALL), &instances::getNoDynamics<T, DESCRIPTOR>());
        setBouzidiZeroVelocityBoundary<T, DESCRIPTOR>(lattice, superGeometry, MAT_WALL, stlReader);
        break;
    case FBC_BOUNCE_BACK:
        // material=2 --> bounceBack dynamics
        lattice.defineDynamics(superGeometry.getMaterialIndicator(MAT_WALL), &instances::getBounceBack<T, DESCRIPTOR>());
        break;
    case FBC_NONE:
    default:
        lattice.defineDynamics(superGeometry.getMaterialIndicator(MAT_WALL), &instances::getNoDynamics<T, DESCRIPTOR>());
        break;
    }

    for(const auto& indicator : parameterSetEnsemble.getFlowIndicators()) {
        if(indicator.type_ == FIT_VELOCITY) {
            switch(wallBoundaryCondition) {
            case FBC_BOUZIDI:
                // no dynamics + bouzidi velocity (inflow)
                lattice.defineDynamics(superGeometry.getMaterialIndicator(indicator.id_), &instances::getNoDynamics<T, DESCRIPTOR>());
                setBouzidiVelocityBoundary<T,DESCRIPTOR>(lattice, superGeometry, indicator.id_, stlReader);
                break;
            case FBC_BOUNCE_BACK:
                // bulk dynamics + velocity (inflow)
                lattice.defineDynamics(superGeometry.getMaterialIndicator(indicator.id_), &bulkDynamics);
                setInterpolatedVelocityBoundary<T,DESCRIPTOR>(lattice, omega, superGeometry.getMaterialIndicator(indicator.id_));
                break;
            case FBC_NONE:
            default:
                lattice.defineDynamics(superGeometry.getMaterialIndicator(MAT_WALL), &instances::getNoDynamics<T, DESCRIPTOR>());
                break;
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
        std::cout << "[fixupLattice] " << errors << " errors have been fixed" << std::endl;
    }

    // Unsteered simulation.
//    if(!measuredData) {

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
//    }
//    // Steered simulation - currently only initializes the first time step!
//    else {
//        MeasuredDataMapper mapper(measuredData->first());
//        lattice.defineU(superGeometry, MAT_FLUID, mapper);
//    }

    // Lattice initialize
    lattice.initialize();
}

// Generates a slowly increasing sinuidal inflow.
void setBoundaryValues( SuperLattice3D<T, DESCRIPTOR>& sLattice,
                        UnitConverter<T,DESCRIPTOR> const& converter,
                        int iteration,
                        SuperGeometry3D<T>& superGeometry,
                        const VolumeList* measuredData,
                        const FlowParameterSetEnsemble& parameterSetEnsemble,
                        const FlowParameterSet& parameters)
{

    FlowBoundaryCondition wallBoundaryCondition = parameters.getWallBoundaryCondition();
    float inletVelocityMultiplier = parameters.getInletVelocityMultiplier();

    for(const auto& indicator : parameterSetEnsemble.getFlowIndicators()) {
        if (indicator.type_ == FIT_VELOCITY) {

            T targetPhysVelocity = indicator.velocityCurve_(converter.getPhysTime(iteration)) * inletVelocityMultiplier;
            T targetLatticeVelocity = converter.getLatticeVelocity(targetPhysVelocity);

            // This function applies the velocity profile to the boundary condition and the lattice.
            auto applyFlowProfile = [&] (AnalyticalF3D<T,T>& profile) {
                if(wallBoundaryCondition == FBC_BOUZIDI) {
                    defineUBouzidi<T, DESCRIPTOR>(sLattice, superGeometry, indicator.id_, profile);
                } else if(wallBoundaryCondition == FBC_BOUNCE_BACK) {
                    sLattice.defineU(superGeometry.getMaterialIndicator(indicator.id_), profile);
                }
            };

            // Create shortcuts.
            const tgt::vec3& center = indicator.center_;
            const tgt::vec3& normal = indicator.normal_;
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
                // TODO: interpolation over temporal domain.
                float time = converter.getPhysTime(iteration);

                size_t idx = 0;
                while(idx < measuredData->size()-1 && measuredData->at(idx)->getTimestep() < time) idx++;

                // For volume indicators, we normalize the velocity curve to [0, 1]
                // because multiplying the measurement by another velocity is confusing.
                // We keep, however, the velocity multiplier, so that we can basically
                // amplify the measurement.
                T multiplier = tgt::clamp<T>(targetPhysVelocity / indicator.velocityCurve_.getMaxVelocity(), 0, 1);
                MeasuredDataMapper mapper(converter, measuredData->at(idx), multiplier);
                applyFlowProfile(mapper);
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


namespace voreen {

void writeVVDFile(STLreader<T>& stlReader,
                  UnitConverter<T,DESCRIPTOR>& converter,
                  int iteration, int maxIteration,
                  const FlowParameterSetEnsemble& parameterSetEnsemble,
                  const FlowParameterSet& parameters,
                  const std::string& simulationOutputPath,
                  const std::string& name,
                  SuperLatticeF3D<T, DESCRIPTOR>& feature) {

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
            << "<MetaItem name=\"" << "ParameterTurbulenceModel" << "\" type=\"IntMetaData\" value=\"" << parameters.getTurbulenceModel() << "\" />"
            << "<MetaItem name=\"" << "ParameterSmagorinskyConstant" << "\" type=\"FloatMetaData\" value=\"" << parameters.getSmagorinskyConstant() << "\" />"
            << "<MetaItem name=\"" << "ParameterWallBoundaryCondition" << "\" type=\"IntMetaData\" value=\"" << parameters.getWallBoundaryCondition() << "\" />"
            << "<MetaItem name=\"" << "ParameterInletVelocityMultiplier" << "\" type=\"FloatMetaData\" value=\"" << parameters.getInletVelocityMultiplier() << "\" />"
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
                 const FlowParameterSet& parameters,
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
                writeVVDFile(stlReader, converter, iteration, maxIteration, parameterSetEnsemble, parameters,
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
                writeVVDFile(stlReader, converter, iteration, maxIteration, parameterSetEnsemble, parameters,
                             simulationOutputPath, "magnitude", magnitude);
            }
            if(writeVTI) {
                vtmWriter.write(magnitude, iteration);
            }
        }

        if(parameterSetEnsemble.getFlowFeatures() & FF_PRESSURE) {
            SuperLatticePhysPressure3D<T, DESCRIPTOR> pressure(sLattice, converter);
            if(writeVVD) {
                writeVVDFile(stlReader, converter, iteration, maxIteration, parameterSetEnsemble, parameters,
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
                writeVVDFile(stlReader, converter, iteration, maxIteration, parameterSetEnsemble, parameters,
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
    addPort(measuredDataPort_);
    measuredDataPort_.addCondition(new PortConditionVolumeListEnsemble());
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

    PortDataPointer<GlMeshGeometryBase> geometry(nullptr, false);

    if(auto geometrySequence = dynamic_cast<const GeometrySequence*>(geometryDataPort_.getData())) {

        std::vector<const GlMeshGeometryBase*> geometries;

        for(size_t i=0; i<geometrySequence->getNumGeometries(); i++) {
            auto geometry = dynamic_cast<const GlMeshGeometryBase*>(geometrySequence->getGeometry(i));
            if(geometry) {
                geometries.push_back(geometry);
            }
            else {
                throw InvalidInputException("GeometrySequence contains non-GlMeshGeometry", InvalidInputException::S_WARNING);
            }
        }

        if(!geometries.empty()) {
            geometry = PortDataPointer<GlMeshGeometryBase>(geometries.front(), false);
            for (size_t i=1; i<geometries.size(); i++) {
                auto combined = mergeGeometries(geometry, geometries.at(i));

                if(!combined) {
                    throw InvalidInputException("Encountered Mesh of unexpected Type", InvalidInputException::S_ERROR);
                }

                geometry = PortDataPointer<GlMeshGeometryBase>(combined.release(), true);
            }
        }
    }
    else if(auto geometryData = dynamic_cast<const GlMeshGeometryBase*>(geometryDataPort_.getData())) {
        geometry = PortDataPointer<GlMeshGeometryBase>(geometryData, false);
    }

    if (!geometry) {
        throw InvalidInputException("No GlMeshGeometry input", InvalidInputException::S_ERROR);
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
        for(const auto& indicator : flowParameterSetEnsemble->getFlowIndicators()) {
            if(indicator.type_ == FIT_VELOCITY && indicator.flowProfile_ == FP_VOLUME) {
                throw InvalidInputException("Volume input required", InvalidInputException::S_ERROR);
            }
        }

        measuredData = nullptr;
        LINFO("Configuring an unsteered simulation");
    }

    std::string geometryPath = VoreenApplication::app()->getUniqueTmpFilePath(".stl");
    try {
        std::ofstream file(geometryPath);
        geometry->exportAsStl(file);
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
    // Nothing to do (yet).
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
            N, // Resolution that charPhysLength is resolved by.
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
    const int noOfCuboids = std::thread::hardware_concurrency();
    CuboidGeometry3D<T> cuboidGeometry(extendedDomain, converter.getConversionFactorLength(), noOfCuboids);

    // Instantiation of a loadBalancer
    HeuristicLoadBalancer<T> loadBalancer(cuboidGeometry);

    interruptionPoint();

    // Instantiation of a superGeometry
    LINFO("Preparing Geometry ...");
    SuperGeometry3D<T> superGeometry(cuboidGeometry, loadBalancer, 2);

    prepareGeometry(converter, extendedDomain, stlReader, superGeometry,
                    parameterSetEnsemble, parameters);

    interruptionPoint();

    // === 3rd Step: Prepare Lattice ===
    LINFO("Preparing Lattice ...");
    SuperLattice3D<T, DESCRIPTOR> sLattice(superGeometry);

    auto createTurbulenceModel = [&] (FlowTurbulenceModel model) -> Dynamics<T, DESCRIPTOR>* {
        switch(model) {
        case FTM_SMAGORINSKY:
            return new SmagorinskyBGKdynamics<T, DESCRIPTOR>(converter.getLatticeRelaxationFrequency(),
                                                             instances::getBulkMomenta<T, DESCRIPTOR>(),
                                                             parameters.getSmagorinskyConstant());
        //case FTM_SMAGORINSKY_SHEAR_IMPROVED: // Does not compile in OpenLB 1.4.
        //    return new ShearSmagorinskyBGKdynamics<T, DESCRIPTOR>(converter.getLatticeRelaxationFrequency(),
        //                                                          instances::getBulkMomenta<T, DESCRIPTOR>(),
        //                                                          parameters.getSmagorinskyConstant());
        case FTM_SMAGORINSKY_CONSISTENT:
            return new ConSmagorinskyBGKdynamics<T, DESCRIPTOR>(converter.getLatticeRelaxationFrequency(),
                                                                instances::getBulkMomenta<T, DESCRIPTOR>(),
                                                                parameters.getSmagorinskyConstant());
        case FTM_SMAGORINSKY_CONSISTENT_STRAIN:
            return new ConStrainSmagorinskyBGKdynamics<T, DESCRIPTOR>(converter.getLatticeRelaxationFrequency(),
                                                                instances::getBulkMomenta<T, DESCRIPTOR>(),
                                                                parameters.getSmagorinskyConstant());
        case FTM_SMAGORINSKY_DYNAMIC:
            return new ConSmagorinskyBGKdynamics<T, DESCRIPTOR>(converter.getLatticeRelaxationFrequency(),
                                                                instances::getBulkMomenta<T, DESCRIPTOR>(),
                                                                parameters.getSmagorinskyConstant());
        case FTM_BGK:
            return new BGKdynamics<T, DESCRIPTOR>(converter.getLatticeRelaxationFrequency(),
                                                  instances::getBulkMomenta<T, DESCRIPTOR>());
        case FTM_NONE:
        default:
            return nullptr;
        };
    };

    std::unique_ptr<Dynamics<T, DESCRIPTOR>> bulkDynamics(createTurbulenceModel(parameters.getTurbulenceModel()));
    if(!bulkDynamics) {
        LERROR("No bulk dynamics");
        return;
    }

    prepareLattice(sLattice, converter, *bulkDynamics,
                   stlReader, superGeometry,
                   measuredData,
                   parameterSetEnsemble,
                   parameters);

    interruptionPoint();

    // === 4th Step: Main Loop  ===
    const int maxIteration = converter.getLatticeTime(parameterSetEnsemble.getSimulationTime());
    util::ValueTracer<T> converge( converter.getLatticeTime(0.5), 1e-5);
    for (int iteration = 0; iteration <= maxIteration; iteration++) {

        // === 5th Step: Definition of Initial and Boundary Conditions ===
        setBoundaryValues(sLattice, converter, iteration, superGeometry, measuredData, parameterSetEnsemble, parameters);

        // === 6th Step: Collide and Stream Execution ===
        sLattice.collideAndStream();

        // === 7th Step: Computation and Output of the Results ===
        bool success = getResults(sLattice, converter, iteration, maxIteration, *bulkDynamics, superGeometry, stlReader,
                                  parameterSetEnsemble, parameters, simulationResultPath);
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
