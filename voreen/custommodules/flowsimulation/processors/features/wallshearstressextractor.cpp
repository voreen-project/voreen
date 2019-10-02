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

#ifdef VRN_MODULE_OPENMP
#include "omp.h"
#endif

#include <olb3D.h>
#define DESCRIPTOR D3Q19Descriptor

using namespace olb;
using namespace olb::descriptors;
typedef double T;

namespace voreen {

const T VOREEN_LENGTH_TO_SI = 0.001;

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
    olb::olbInit(nullptr, nullptr);

    VolumeRAMRepresentationLock representation(input.measuredData);
    RealWorldMapping rwm = input.measuredData->getRealWorldMapping();
    tgt::svec3 dimensions = representation->getDimensions();
    std::unique_ptr<VolumeRAM_Float> output = std::move(input.output);
    output->clear();
/*
    // TODO: linker error
    UnitConverterFromResolutionAndLatticeVelocity<T, DESCRIPTOR> converter(
            tgt::min(dimensions),
            1.0,
            (T) 1.0,
            (T) 1.0,
            (T) input.viscosity * 0.001 / input.density,
            (T) input.density
    );

    STLreader<T> stlReader(input.geometryPath, converter.getConversionFactorLength(), VOREEN_LENGTH_TO_SI, 1);
    Cuboid3D<T> cuboid(stlReader, converter.getConversionFactorLength());
    BlockGeometry3D<T> geometry(cuboid);
    BlockLattice3D<T, DESCRIPTOR> lattice(dimensions.x, dimensions.y, dimensions.z, geometry);
    BlockLatticePhysWallShearStress3D<T, DESCRIPTOR> wallShearStress(lattice, geometry, 2,
                                                                     converter, stlReader);

//#pragma omp parallel for
    for(int z = 0; z < lattice.getNz(); z++) {
        for(int y = 0; y < lattice.getNy(); y++) {
            for(int x = 0; x < lattice.getNx(); x++) {

                // Retrieve the actual velocity first.
                T velocity[3];
                for(size_t i=0; i<3; i++) {
                    velocity[i] = representation->getVoxelNormalized(x, y, z, i);
                    velocity[i] = rwm.normalizedToRealWorld(velocity[i]) * VOREEN_LENGTH_TO_SI;
                }

                // Set it to the lattice.
                lattice.get(x, y, z).defineU(velocity);

                // Calculate wall shear stress.
                T value[1];
                int position [] = {x, y, z};
                if(wallShearStress(value, position)) {
                    output->voxel(x, y, z) = value[0];
                }
            }
        }
        progressReporter.setProgress(1.0f * z / lattice.getNz());
    }
*/
    std::unique_ptr<VolumeBase> volume(new Volume(output.release(), input.measuredData));
    return WallShearStressExtractorOutput{
        std::move(volume)
    };
}

void WallShearStressExtractor::processComputeOutput(WallShearStressExtractorOutput output) {
    //Volume* volume = new Volume(output.volume.release(), inputVolume_.getData());
    outputVolume_.setData(output.volume.release());
}

} // namespace voreen
