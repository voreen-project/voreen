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


static const uint32_t INF = std::numeric_limits<uint32_t>::max();
static uint32_t inf_add(uint32_t x, uint32_t y) {
    uint32_t res = x + y;
    res |= -(res < x);
    return res;
}

static uint32_t square(uint32_t i) {
    tgtAssert(i != INF, "i is inf");
    return i*i;
}

LargeVolumeDistanceTransform::ComputeOutput LargeVolumeDistanceTransform::compute(LargeVolumeDistanceTransform::ComputeInput input, ProgressReporter& progressReporter) const {
    const VolumeBase& vol = *input.inputVolume_;
    const tgt::svec3 dim = vol.getDimensions();
    const tgt::svec3 sliceDim(dim.x, dim.y, 1);

    LZ4SliceVolumeBuilder<uint32_t> gBuilder(input.outputPath_,
            LZ4SliceVolumeMetadata(dim)
            .withOffset(vol.getOffset())
            .withSpacing(vol.getSpacing())
            .withPhysicalToWorldTransformation(vol.getPhysicalToWorldMatrix())
            .withRealWorldMapping(RealWorldMapping::createDenormalizingMapping<uint32_t>()));

    //VolumeAtomic<uint32_t> prev;
    const float binarizationThreshold = 0.5f;

    VolumeAtomic<uint32_t> gSlice(sliceDim);
    // z-scan 1
    {
        const size_t z = 0;
        std::unique_ptr<VolumeRAM> inputSlice(vol.getSlice(z));
        for(size_t y = 0; y < dim.y; ++y) {
            for(size_t x = 0; x < dim.x; ++x) {
                tgt::svec3 slicePos(x,y,0);
                float val = inputSlice->getVoxelNormalized(slicePos);
                uint32_t& g = gSlice.voxel(slicePos);
                if(val < binarizationThreshold) {
                    // Background
                    g = 0;
                } else {
                    // Foreground
                    g = INF;
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
                uint32_t& g = gSlice.voxel(slicePos);
                if(val < binarizationThreshold) {
                    // Background
                    g = 0;
                } else {
                    // Foreground
                    g = inf_add(g, 1);
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
                    uint32_t& gPrev = prevSlice->voxel(slicePos);
                    uint32_t& g = gSlice->voxel(slicePos);

                    if(gPrev < g) {
                        g = inf_add(gPrev, 1);
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
    auto tmpSlice = VolumeAtomic<uint32_t>(sliceDim);
    for(size_t z=0; z<dim.z; ++z) {
        progressReporter.setProgress(static_cast<float>(z)/dim.z);

        auto gSlice = gvol.getWriteableSlice(z);
        auto fSlice = builder.getNextWriteableSlice();

        {
            const int m = dim.y;
            std::vector<uint32_t> s(m, 0);
            std::vector<int64_t> t(m, 0);

            for(size_t x = 0; x < dim.x; ++x) {
                auto f = [&] (int64_t z, int64_t i) {
                    tgtAssert(z != INF, "z is inf");
                    tgtAssert(i != INF, "i is inf");
                    tgt::svec3 slicePos(x,i,0);
                    auto diff = z-i;
                    int64_t g = gSlice->voxel(slicePos);
                    if(g == INF) {
                        return INF;
                    }
                    return (uint32_t)(diff*diff + g*g);
                };
                auto sep = [&] (uint32_t i, uint32_t u) {
                    tgt::svec3 uPos(x,u,0);
                    tgt::svec3 iPos(x,i,0);

                    //TODO: see how to resolve INF - INF. Maybe reference
                    // Distance Transforms of Sampled Functions (2012)
                    int64_t gu = gSlice->voxel(uPos);
                    int64_t gi = gSlice->voxel(iPos);
                    tgtAssert(gi != INF, "inf gi");
                    if(gu == INF) {
                        return (int64_t)INF;
                    }

                    return (gu*gu - gi*gi + square(u) - square(i))/(2*(u-i));
                };
                s[0] = 0;
                t[0] = 0;
                int q = 0;

                for(size_t u = 1; u < m; ++u) {
                    while(q >= 0 && f(t[q], s[q]) > f(t[q], u)) {
                        q -= 1;
                    }
                    if(q < 0) {
                        q = 0;
                        s[0] = u;
                    } else {
                        int64_t w = sep(s[q], u) + 1;
                        if(w < (int64_t)m) {
                            q += 1;
                            s[q] = u;
                            t[q] = w;
                        }
                    }
                }
                for(int u = m-1; u >= 0; --u) {
                    tgt::svec3 slicePos(x,u,0);
                    tmpSlice.voxel(slicePos) = f(u, s[q]);
                    if(u == t[q]) {
                        q -= 1;
                    }
                }
            }
        }

        {
            const int m = dim.x;
            std::vector<uint32_t> s(m, 0);
            std::vector<int64_t> t(m, 0);

            for(size_t y = 0; y < dim.y; ++y) {
                auto f = [&] (int64_t z, int64_t i) {
                    tgtAssert(z != INF, "z is inf");
                    tgtAssert(i != INF, "i is inf");
                    tgt::svec3 slicePos(i,y,0);
                    auto diff = z-i;
                    int64_t t = tmpSlice.voxel(slicePos);
                    if(t == INF) {
                        return INF;
                    }
                    return (uint32_t)(diff*diff + t);
                };
                auto sep = [&] (uint32_t i, uint32_t u) {
                    tgt::svec3 uPos(u,y,0);
                    tgt::svec3 iPos(i,y,0);

                    int64_t tu = tmpSlice.voxel(uPos);
                    int64_t ti = tmpSlice.voxel(iPos);
                    tgtAssert(ti != INF, "inf ti");
                    if(tu == INF) {
                        return (int64_t)INF;
                    }

                    return (tu - ti + square(u) - square(i))/(2*(u-i));
                };
                s[0] = 0;
                t[0] = 0;
                int q = 0;

                for(size_t u = 1; u < m; ++u) {
                    while(q >= 0 && f(t[q], s[q]) > f(t[q], u)) {
                        q -= 1;
                    }
                    if(q < 0) {
                        q = 0;
                        s[0] = u;
                    } else {
                        int64_t w = sep(s[q], u) + 1;
                        if(w < (int64_t)m) {
                            q += 1;
                            s[q] = u;
                            t[q] = w;
                        }
                    }
                }
                for(int u = m-1; u >= 0; --u) {
                    tgt::svec3 slicePos(u,y,0);
                    uint32_t dist_squared = f(u, s[q]);
                    float dist = std::sqrt((float)dist_squared);
                    fSlice->voxel(slicePos) = dist;
                    if(u == t[q]) {
                        q -= 1;
                    }
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
