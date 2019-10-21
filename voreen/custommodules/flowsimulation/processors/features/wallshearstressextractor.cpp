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

#ifdef VRN_MODULE_OPENMP
#include "omp.h"
#endif

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
        : AnalyticalF3D<T, T>(3)
        , volume_(volume)
    {
        tgtAssert(volume_, "No volume");
        tgtAssert(volume_->getNumChannels() == 3, "Num channels != 3");
        representation_.reset(new VolumeRAMRepresentationLock(volume_));
        typedRepresentation_ = dynamic_cast<const VolumeRAM_3xFloat*>(**representation_);
        tgtAssert(typedRepresentation_, "cast failed");
    }
    virtual bool operator() (T output[], const T input[]) {
        tgt::vec3 voxel = typedRepresentation_->voxel(input[0], input[1], input[2]);
        for(size_t i=0; i < typedRepresentation_->getNumChannels(); i++) {
            output[i] = voxel[i] * VOREEN_LENGTH_TO_SI;
        }
        return true;
    }

private:
    const VolumeBase* volume_;
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
{
    addPort(inputVolume_);
    inputVolume_.addCondition(new PortConditionVolumeChannelCount(3));
    addPort(inputGeometry_);
    addPort(outputVolume_);

    addProperty(viscosity_);
    addProperty(density_);
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

    std::unique_ptr<VolumeRAM_Float> output(new VolumeRAM_Float(measuredData->getDimensions()));

    return WallShearStressExtractorInput{
            geometryPath,
            measuredData,
            viscosity,
            density,
            std::move(output)
    };
}

WallShearStressExtractorOutput WallShearStressExtractor::compute(WallShearStressExtractorInput input, ProgressReporter& progressReporter) const {

    // Needs to be initialized in each new thread to be used.
    //olb::olbInit(nullptr, nullptr);

    VolumeRAMRepresentationLock representation(input.measuredData);
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
            spacing / 40,
            (T) length,
            (T) 1.0,
            (T) input.viscosity * 0.001 / input.density,
            (T) input.density
    );

    STLreader<T> stlReader(input.geometryPath, converter.getConversionFactorLength(), 1.0, 1);
    IndicatorLayer3D<T> extendedDomain(stlReader, converter.getConversionFactorLength());
    Cuboid3D<T> cuboid(stlReader, converter.getConversionFactorLength());
    BlockGeometry3D<T> geometry(cuboid);
    geometry.rename(MAT_EMPTY, MAT_WALL, extendedDomain);
    geometry.rename(MAT_WALL, MAT_FLUID, stlReader);
    geometry.clean();
    geometry.innerClean(MAT_COUNT);
    geometry.checkForErrors();
    BGKdynamics<T, DESCRIPTOR> bulkDynamics(converter.getLatticeRelaxationFrequency(), instances::getBulkMomenta<T, DESCRIPTOR>());
    BlockLattice3D<T, DESCRIPTOR> lattice(geometry.getNx(), geometry.getNy(), geometry.getNz(), geometry);
    lattice.defineDynamics(geometry, MAT_EMPTY, &instances::getNoDynamics<T, DESCRIPTOR>());
    lattice.defineDynamics(geometry, MAT_FLUID, &bulkDynamics);
    lattice.defineDynamics(geometry, MAT_WALL, &instances::getNoDynamics<T, DESCRIPTOR>());
    //lattice.defineDynamics(geometry, MAT_WALL, &instances::getBounceBack<T, DESCRIPTOR>());

    MeasuredDataMapper mapper(input.measuredData);
    //lattice.defineU(geometry, MAT_WALL, mapper);
    lattice.defineU(geometry, MAT_FLUID, mapper);

    lattice.initialize();
    lattice.collideAndStream();

    BlockLatticePhysWallShearStress3D<T, DESCRIPTOR> wallShearStress(lattice, geometry, MAT_WALL, converter, stlReader);

    for(int z = 1; z < lattice.getNz()-1; z++) {
        for(int y = 1; y < lattice.getNy()-1; y++) {
            for(int x = 1; x < lattice.getNx()-1; x++) {
                //if(geometry.get(x, y, z) == MAT_WALL) {
                    // Calculate wall shear stress.
                    T value[1];
                    int position[] = {x, y, z};
                    if(wallShearStress(value, position)) {
                        output->voxel(x, y, z) = value[0];
                    }
                //}
            }
        }
        progressReporter.setProgress(1.0f * z / lattice.getNz());
    }

    std::unique_ptr<VolumeBase> volume(new Volume(output.release(), input.measuredData));
    return WallShearStressExtractorOutput{
        std::move(volume)
    };
}

void WallShearStressExtractor::processComputeOutput(WallShearStressExtractorOutput output) {
    outputVolume_.setData(output.volume.release());
}

} // namespace voreen
