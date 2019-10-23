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

#include "wallshearstressextractor.h"

#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/datastructures/geometry/glmeshgeometry.h"
#include "voreen/core/ports/conditions/portconditionvolumetype.h"

#include <olb3D.h>
#ifndef OLB_PRECOMPILED
#include "olb3D.hh"
#endif
#define DESCRIPTOR D3Q19Descriptor

using namespace olb;
using namespace olb::descriptors;
typedef double T;

namespace voreen {

const T VOREEN_LENGTH_TO_SI = 0.001;


class MeasuredDataMapper : public AnalyticalF3D<T, T> {
public:

    MeasuredDataMapper(const VolumeBase* volume)
            : AnalyticalF3D<T, T>(3), volume_(volume) {
        tgtAssert(volume_, "No volume");
        tgtAssert(volume_->getNumChannels() == 3, "Num channels != 3");
        bounds_ = volume_->getBoundingBox(false).getBoundingBox(false);
        physicalToVoxelMatrix_ = volume_->getPhysicalToVoxelMatrix();
        representation_.reset(new VolumeRAMRepresentationLock(volume_));
        typedRepresentation_ = dynamic_cast<const VolumeRAM_3xFloat*>(**representation_);
        tgtAssert(typedRepresentation_, "cast failed");
    }

    bool operator()(T output[], const T input[]) {
        tgt::vec3 rwPos = tgt::Vector3<T>::fromPointer(input);
        if (!bounds_.containsPoint(rwPos)) {
            return false;
        }

        tgt::vec3 voxel = typedRepresentation_->getVoxelLinear(physicalToVoxelMatrix_ * rwPos, 0, false);
        for (size_t i = 0; i < typedRepresentation_->getNumChannels(); i++) {
            output[i] = voxel[i] * VOREEN_LENGTH_TO_SI;
        }

        return true;
    }
private:

    const VolumeBase* volume_;
    tgt::Bounds bounds_;
    tgt::mat4 physicalToVoxelMatrix_;
    std::unique_ptr<VolumeRAMRepresentationLock> representation_;
    const VolumeRAM_3xFloat* typedRepresentation_;
};


const std::string WallShearStressExtractor::loggerCat_("voreen.flowsimulation.wallshearstressextractor");

WallShearStressExtractor::WallShearStressExtractor()
    : AsyncComputeProcessor<ComputeInput, ComputeOutput>()
    , inputVolume_(Port::INPORT, "wallshearstressextractor.inputVolume", "Volume Input")
    , inputGeometry_(Port::INPORT, "wallshearstressextractor.inputGeometry", "Geometry Input")
    , outputVolume_(Port::OUTPORT, "wallshearstressextractor.outputVolume", "Volume Output")
    , viscosity_("viscosity", "Dynamic Viscosity (e-3 kg/(m x s))", 3.5, 3, 4)
    , density_("density", "Density (kg/m^3)", 1000.0f, 1000.0f, 1100.0f)
    , padding_("padding", "Add Padding", 1, 0, 5)
{
    addPort(inputVolume_);
    inputVolume_.addCondition(new PortConditionVolumeChannelCount(3));
    addPort(inputGeometry_);
    addPort(outputVolume_);

    addProperty(viscosity_);
    addProperty(density_);
    addProperty(padding_);
}

Processor* WallShearStressExtractor::create() const {
    return new WallShearStressExtractor();
}

WallShearStressExtractorInput WallShearStressExtractor::prepareComputeInput() {

    const VolumeBase* measuredData = inputVolume_.getThreadSafeData();

    float viscosity = measuredData->getMetaDataValue<FloatMetaData>("ParameterViscosity", viscosity_.get());
    float density = measuredData->getMetaDataValue<FloatMetaData>("ParameterDensity", density_.get());

    const GlMeshGeometryBase* geometryData = dynamic_cast<const GlMeshGeometryBase*>(inputGeometry_.getData());
    if (!geometryData) {
        throw InvalidInputException("Invalid geometry", InvalidInputException::S_WARNING);
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

    // Add a voxel in each dimension
    size_t padding = padding_.get();
    tgt::svec3 dimensions = measuredData->getDimensions();
    dimensions += tgt::svec3(padding*2);
    std::unique_ptr<VolumeRAM_Float> output(new VolumeRAM_Float(dimensions));

    return WallShearStressExtractorInput{
            geometryPath,
            measuredData,
            viscosity,
            density,
            padding,
            std::move(output)
    };
}

WallShearStressExtractorOutput WallShearStressExtractor::compute(WallShearStressExtractorInput input, ProgressReporter& progressReporter) const {

    // Needs to be initialized in each new thread to be used.
    olb::olbInit(nullptr, nullptr);

    VolumeRAMRepresentationLock representation(input.measuredData);
    tgt::svec3 dim = input.measuredData->getDimensions();
    tgt::mat4 voxelToPhysicalMatrix = input.measuredData->getVoxelToPhysicalMatrix();
    float spacing = tgt::min(input.measuredData->getSpacing());
    float length = tgt::max(input.measuredData->getCubeSize());

    std::unique_ptr<VolumeRAM_Float> output = std::move(input.output);
    output->clear();

    const int MAT_EMPTY = 0;
    const int MAT_FLUID = 1;
    const int MAT_WALL  = 2;
    const int MAT_COUNT = 3;

    UnitConverter<T, DESCRIPTOR> converter(
            spacing,
            1.0,
            length,
            1.0,
            input.viscosity * 0.001 / input.density,
            input.density
    );

    // Setup geometry.
    STLreader<T> stlReader(input.geometryPath, converter.getConversionFactorLength(), 1.0, 1);
    IndicatorLayer3D<T> extendedDomain(stlReader, converter.getConversionFactorLength());
    CuboidGeometry3D<T> cuboidGeometry(extendedDomain, converter.getConversionFactorLength(), 1);
    HeuristicLoadBalancer<T> loadBalancer(cuboidGeometry);
    SuperGeometry3D<T> superGeometry(cuboidGeometry, loadBalancer, 2);
    superGeometry.rename( MAT_EMPTY, MAT_WALL,  extendedDomain );
    superGeometry.rename( MAT_WALL,  MAT_FLUID, stlReader );
    superGeometry.clean();
    superGeometry.innerClean(MAT_COUNT);
    superGeometry.checkForErrors();

    // Setup lattice.
    SuperLattice3D<T, DESCRIPTOR> sLattice(superGeometry);
    BGKdynamics<T, DESCRIPTOR> bulkDynamics(converter.getLatticeRelaxationFrequency(), instances::getBulkMomenta<T, DESCRIPTOR>());
    sOnLatticeBoundaryCondition3D<T, DESCRIPTOR> sBoundaryCondition(sLattice);
    createInterpBoundaryCondition3D<T, DESCRIPTOR>(sBoundaryCondition);
    sOffLatticeBoundaryCondition3D<T, DESCRIPTOR> sOffBoundaryCondition(sLattice);
    createBouzidiBoundaryCondition3D<T, DESCRIPTOR>(sOffBoundaryCondition);
    sLattice.defineDynamics(superGeometry, MAT_EMPTY, &instances::getNoDynamics<T, DESCRIPTOR>());
    sLattice.defineDynamics(superGeometry, MAT_FLUID, &bulkDynamics);
    sLattice.defineDynamics(superGeometry, MAT_WALL, &instances::getNoDynamics<T, DESCRIPTOR>());
    sOffBoundaryCondition.addZeroVelocityBoundary(superGeometry, MAT_WALL, stlReader);
    MeasuredDataMapper mapper(input.measuredData);
    sLattice.defineU(superGeometry, MAT_FLUID, mapper);
    sLattice.initialize();

    // Retrieve result.
    SuperLatticePhysWallShearStress3D<T, DESCRIPTOR> wallShearStress(sLattice, superGeometry, MAT_WALL, converter, stlReader);
    AnalyticalFfromSuperF3D<T> interpolateFeature(wallShearStress, true);

#ifdef VRN_MODULE_OPENMP
#pragma omp parallel for
#endif
    for (size_t z = 0; z < dim.z; z++) {
        for (size_t y = 0; y < dim.y; y++) {
            for (size_t x = 0; x < dim.x; x++) {

                tgt::Vector3<T> pos = voxelToPhysicalMatrix * tgt::vec3(x, y, z);

                T value[1];
                if(interpolateFeature(&value[0], pos.elem)) {
                    output->voxel(x+input.padding, y+input.padding, z+input.padding) = static_cast<float>(value[0]);
                }
            }
        }
#ifndef VRN_MODULE_OPENMP
        progressReporter.setProgress(1.0f * z / dim.z);
#endif
    }

    std::unique_ptr<Volume> volume(new Volume(output.release(), input.measuredData));
    volume->setOffset(input.measuredData->getOffset() - tgt::vec3(input.padding) * spacing);

    progressReporter.setProgress(1.0f);

    return WallShearStressExtractorOutput{
        std::move(volume)
    };
}

void WallShearStressExtractor::processComputeOutput(WallShearStressExtractorOutput output) {
    outputVolume_.setData(output.volume.release());
}

} // namespace voreen
