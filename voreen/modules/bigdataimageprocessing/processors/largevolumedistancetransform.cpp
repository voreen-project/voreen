/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2021 University of Muenster, Germany,                        *
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

#include "largevolumedistancetransform.h"
#include "voreen/core/utils/stringutils.h"
#include "voreen/core/datastructures/volume/volumeminmax.h"

#include "tgt/vector.h"

namespace voreen {

const std::string LargeVolumeDistanceTransform::loggerCat_("voreen.bigdataimageprocessing.LargeVolumeDistanceTransform");

LargeVolumeDistanceTransform::LargeVolumeDistanceTransform()
    : AsyncComputeProcessor<ComputeInput, ComputeOutput>()
    , inport_(Port::INPORT, "volumehandle.input", "Volume Input")
    , outport_(Port::OUTPORT, "volumehandle.output", "Volume Output",false)
    , outputVolumeFilePath_("outputVolumeFilePath", "Output volume file path", "Output volume file path", "", LZ4SliceVolumeBase::FILE_EXTENSION)
{
    addPort(inport_);
    addPort(outport_);

    addProperty(outputVolumeFilePath_);
}

LargeVolumeDistanceTransform::~LargeVolumeDistanceTransform() {}

Processor* LargeVolumeDistanceTransform::create() const {
    return new LargeVolumeDistanceTransform();
}

void LargeVolumeDistanceTransform::adjustPropertiesToInput() {
    //
}

template<typename OutputFormat>
void processDispatch(const VolumeBase& input, std::unique_ptr<Volume>& output, const std::string& outputPath, ProgressReporter& progressReporter) {
}

LargeVolumeDistanceTransform::ComputeInput LargeVolumeDistanceTransform::prepareComputeInput() {
    if(outputVolumeFilePath_.get().empty()) {
        throw InvalidInputException("No volume path specified", InvalidInputException::S_ERROR);
    }

    if(!inport_.hasData()) {
        throw InvalidInputException("No input", InvalidInputException::S_WARNING);
    }

    return LargeVolumeDistanceTransform::ComputeInput {
        outputVolumeFilePath_.get(),
        inport_.getThreadSafeData()
    };
}


template<typename T>
static inline T square(T i) {
    return i*i;
}

LargeVolumeDistanceTransform::ComputeOutput LargeVolumeDistanceTransform::compute(LargeVolumeDistanceTransform::ComputeInput input, ProgressReporter& progressReporter) const {
    const VolumeBase& vol = *input.inputVolume_;
    const tgt::svec3 dim = vol.getDimensions();
    const tgt::svec3 sliceDim(dim.x, dim.y, 1);

    LZ4SliceVolumeBuilder<float> gBuilder(input.outputPath_,
            LZ4SliceVolumeMetadata(dim)
            .withOffset(vol.getOffset())
            .withSpacing(vol.getSpacing())
            .withPhysicalToWorldTransformation(vol.getPhysicalToWorldMatrix()));

    const float binarizationThreshold = 0.5f;

    VolumeAtomic<float> gSlice(sliceDim);
    // z-scan 1
    {
        const size_t z = 0;
        std::unique_ptr<VolumeRAM> inputSlice(vol.getSlice(z));
        for(size_t y = 0; y < dim.y; ++y) {
            for(size_t x = 0; x < dim.x; ++x) {
                tgt::svec3 slicePos(x,y,0);
                float val = inputSlice->getVoxelNormalized(slicePos);
                float& g = gSlice.voxel(slicePos);
                if(val < binarizationThreshold) {
                    // Background
                    g = 0;
                } else {
                    // Foreground
                    g = std::numeric_limits<float>::infinity();
                }
            }
        }
        gBuilder.pushSlice(gSlice);
    }
    for(size_t z = 1; z < dim.z; ++z) {
        progressReporter.setProgress(static_cast<float>(z)/dim.z);
        std::unique_ptr<VolumeRAM> inputSlice(vol.getSlice(z));


        for(size_t y = 0; y < dim.y; ++y) {
            for(size_t x = 0; x < dim.x; ++x) {
                tgt::svec3 slicePos(x,y,0);
                float val = inputSlice->getVoxelNormalized(slicePos);
                float& g = gSlice.voxel(slicePos);
                if(val < binarizationThreshold) {
                    // Background
                    g = 0;
                } else {
                    // Foreground
                    g += 1;
                }
            }
        }
        gBuilder.pushSlice(gSlice);
    }
    auto gvol = std::move(gBuilder).finalize();

    // z-scan 2
    {
        auto prevSlice = gvol.getWriteableSlice(dim.z-1);
        for(int z = dim.z-2; z >= 0; --z) {
            progressReporter.setProgress(static_cast<float>(z)/dim.z);

            auto gSlice = gvol.getWriteableSlice(z);
            for(size_t y = 0; y < dim.y; ++y) {
                for(size_t x = 0; x < dim.x; ++x) {
                    tgt::svec3 slicePos(x,y,0);
                    float& gPrev = prevSlice->voxel(slicePos);
                    float& g = gSlice->voxel(slicePos);

                    if(gPrev < g) {
                        g = gPrev + 1;
                    }
                }
            }
            prevSlice = std::move(gSlice);
        }
    }

    LZ4SliceVolumeBuilder<float> builder(input.outputPath_+ "foobar", //TODO
            LZ4SliceVolumeMetadata(dim)
            .withOffset(vol.getOffset())
            .withSpacing(vol.getSpacing())
            .withPhysicalToWorldTransformation(vol.getPhysicalToWorldMatrix()));

    // y and x scans
    auto tmpSlice = VolumeAtomic<float>(sliceDim);
    for(size_t z=0; z<dim.z; ++z) {
        progressReporter.setProgress(static_cast<float>(z)/dim.z);

        auto gSlice = gvol.getWriteableSlice(z);
        auto fSlice = builder.getNextWriteableSlice();

        {
            const int n = dim.y;
            std::vector<int> v(n+1, 0); // Locations of parabolas in lower envelope
            std::vector<float> z(n+1, 0); // Locations of boundaries between parabolas

            for(int x = 0; x < dim.x; ++x) {
                auto f = [&] (int i) {
                    tgt::svec3 slicePos(x,i,0);
                    float g = gSlice->voxel(slicePos);
                    // We have to square here since we have not done it in the z-passes
                    // TODO: maybe we can do that
                    return square(g);
                };

                v[0] = 0;
                z[0] = -std::numeric_limits<float>::infinity();
                z[1] = std::numeric_limits<float>::infinity();
                int k = 0;

                for(int q = 1; q < n; ++q) {
                    float fq = f(q);
                    if(std::isinf(fq)) {
                        continue;
                    }
jmp:
                    int vk = v[k];
                    float s = ((fq - f(vk)) + (square(q) - square(vk))) / (2*(q - vk)); //note: q > vk
                    tgtAssert(!std::isnan(s), "bleh"); // this will happen if fv and fq are inf

                    if(s <= z[k]) {
                        if(k > 0) {
                            k -= 1;
                            goto jmp;
                        } else {
                            v[k] = q;
                            z[k] = s;
                            z[k+1] = std::numeric_limits<float>::infinity();
                        }
                    } else {
                        k += 1;
                        v[k] = q;
                        z[k] = s;
                        z[k+1] = std::numeric_limits<float>::infinity();
                    }
                }

                int maxk = k;
                k = 0;
                for(int q = 0; q < n; ++q) {
                    while(z[k+1] < q) {
                        tgtAssert(k < maxk, "bleh");
                        k += 1;
                    }
                    tgt::svec3 slicePos(x,q,0);
                    int vk = v[k];
                    float val = f(vk) + square(q-vk);
                    tgtAssert(val >= 0, "temporary check");
                    tmpSlice.voxel(slicePos) = val;
                }
            }
        }

        {
            const int n = dim.x;
            std::vector<int> v(n+1, 0); // Locations of parabolas in lower envelope
            std::vector<float> z(n+1, 0); // Locations of boundaries between parabolas

            for(int y = 0; y < dim.y; ++y) {
                auto f = [&] (int i) {
                    tgt::svec3 slicePos(i,y,0);
                    float g = tmpSlice.voxel(slicePos);
                    return g;
                };

                v[0] = 0;
                z[0] = -std::numeric_limits<float>::infinity();
                z[1] = std::numeric_limits<float>::infinity();
                int k = 0;

                for(int q = 1; q < n; ++q) {
                    float fq = f(q);
                    if(std::isinf(fq)) {
                        continue;
                    }
jmp2:
                    int vk = v[k];
                    float s = ((fq - f(vk)) + (square(q) - square(vk))) / (2*(q - vk)); //note: q > vk
                    tgtAssert(!std::isnan(s), "bleh"); // this will happen if fv and fq are inf

                    if(s <= z[k]) {
                        if(k > 0) {
                            k -= 1;
                            goto jmp2;
                        } else {
                            v[k] = q;
                            z[k] = s;
                            z[k+1] = std::numeric_limits<float>::infinity();
                        }
                    } else {
                        k += 1;
                        v[k] = q;
                        z[k] = s;
                        z[k+1] = std::numeric_limits<float>::infinity();
                    }
                }

                int maxk = k;
                k = 0;
                for(int q = 0; q < n; ++q) {
                    while(z[k+1] < q) {
                        tgtAssert(k < maxk, "bleh");
                        k += 1;
                    }
                    tgt::svec3 slicePos(q,y,0);
                    int vk = v[k];
                    float val = f(vk) + square(q-vk);
                    tgtAssert(val >= 0, "temporary check");
                    fSlice->voxel(slicePos) = std::sqrt(val);
                }
            }
        }
    }

    progressReporter.setProgress(1.f);
    return {
        std::move(builder).finalize().toVolume()
        //std::move(gvol).toVolume()
    };
}
void LargeVolumeDistanceTransform::processComputeOutput(LargeVolumeDistanceTransform::ComputeOutput output) {
    outport_.setData(output.outputVolume_.release());
}

}   // namespace
