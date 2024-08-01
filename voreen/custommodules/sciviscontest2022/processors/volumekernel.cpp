/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2024 University of Muenster, Germany,                        *
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

#include "volumekernel.h"

#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/datastructures/volume/volumebase.h"

namespace voreen {

const std::string VolumeKernel::loggerCat_ = "VolumeKernel";

VolumeKernel::VolumeKernel() 
    : Processor()
    , inport_( Port::INPORT, "jacobianInport", "Jacobian of the volume of interest" )
    , outputVolume_( Port::OUTPORT, "outputVolume", "Scalar representation of the pyrogenic vorticity field." )
    , kernelSizeProperty_("kernelSizeProperty", "Kenerl Size: ", 3, 2, 1000, Processor::InvalidationLevel::INVALID_PARAMETERS)
    , sigmaProperty_("sigmaProperty", "Sigma: ", 0.6f, 0.1f, 5.f)
{
    addPort(inport_);
    addPort(outputVolume_);

    addProperty(kernelSizeProperty_);
    addProperty(sigmaProperty_);
}

void VolumeKernel::initializeNormalKernel(float* values, int kernelSize, float sigma) const {
    auto f0 = 1.f / (sigma * std::sqrt(3.1415926f));
    auto f1 = -1.f / (2.f * sigma * sigma);
    auto hSize = kernelSize / 2;
    for (auto x = 0; x < kernelSize; ++x) {
        auto x_ = x - hSize;
        float value = f0 * std::exp(x_ * x_ * f1);
        values[x] = value;
    }
}

void VolumeKernel::process() {
    auto kernelSize = kernelSizeProperty_.get();

    auto kernel = std::unique_ptr<float[]>(new float[kernelSize]);

    initializeNormalKernel(kernel.get(), kernelSize, sigmaProperty_.get());

    auto volume = inport_.getData();
    VolumeRAMRepresentationLock volumeLock(volume);
    auto dimensions = volume->getDimensions();

    const auto* volumeData = dynamic_cast<const VolumeRAM_Float*>(*volumeLock);
    if (volumeData == nullptr)
        throw std::runtime_error("Expected scalar volume as inport!");
    
    auto convolutedVolume = std::unique_ptr<VolumeRAM_Float>(new VolumeRAM_Float(volume->getDimensions()));
    auto xVolume = std::unique_ptr<VolumeRAM_Float>(new VolumeRAM_Float(volume->getDimensions()));
    auto yVolume = std::unique_ptr<VolumeRAM_Float>(new VolumeRAM_Float(volume->getDimensions()));

    auto hSize = kernelSize / 2;

#ifdef VRN_MODULE_OPENMP
#pragma omp parallel for
#endif
    for (auto z = 0; z < dimensions.z; ++z) {
        for (auto y = 0; y < dimensions.y; ++y) {
            
            // apply kernel to x-axis
            for (auto x = 0; x < dimensions.x; ++x) {
                float convoluted = 0.f;
                tgt::svec3 pos(x, y, z);
                for (auto i = 0; i < kernelSize; ++i) {
                    auto x_ = x + (i - hSize);
                    if (x_ < 0 || x_ >= dimensions.x)
                        continue;
                    
                    convoluted += volumeData->voxel(x_, y, z) * kernel[i];
                }
                convoluted /= kernelSize;

                xVolume->voxel(pos) = convoluted;
            }
        }
    }

#ifdef VRN_MODULE_OPENMP
#pragma omp parallel for
#endif
    for (auto z = 0; z < dimensions.z; ++z) {
        for (auto x = 0; x < dimensions.x; ++x) {
            
            // apply kernel to x-axis
            for (auto y = 0; y < dimensions.y; ++y) {
                float convoluted = 0.f;
                tgt::svec3 pos(x, y, z);
                for (auto i = 0; i < kernelSize; ++i) {
                    auto y_ = y + (i - hSize);
                    if (y_ < 0 || y_ >= dimensions.x)
                        continue;
                    
                    convoluted += xVolume->voxel(x, y_, z) * kernel[i];
                }
                convoluted /= kernelSize;

                yVolume->voxel(pos) = convoluted;
            }
        }
    }


#ifdef VRN_MODULE_OPENMP
#pragma omp parallel for
#endif
    for (auto x = 0; x < dimensions.x; ++x) {
        for (auto y = 0; y < dimensions.y; ++y) {
            
            // apply kernel to x-axis
            for (auto z = 0; z < dimensions.z; ++z) {
                float convoluted = 0.f;
                tgt::svec3 pos(x, y, z);
                for (auto i = 0; i < kernelSize; ++i) {
                    auto z_ = z + (i - hSize);
                    if (z_ < 0 || z_ >= dimensions.x)
                        continue;
                    
                    convoluted += yVolume->voxel(x, y, z_) * kernel[i];
                }
                convoluted /= kernelSize;

                convolutedVolume->voxel(pos) = convoluted;
            }
        }
    }

    auto* output = new Volume(convolutedVolume.release(), volume);
    output->setRealWorldMapping(RealWorldMapping()); // Override to default rwm.
    output->setModality(Modality("convoluted_" + volume->getModality().getName()));
    outputVolume_.setData(output);
}

}
