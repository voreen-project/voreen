/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2019 University of Muenster, Germany,                        *
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

#include "flow3d.h"
#include "modules/flowreen/utils/flowmath.h"
#include "streamlinetexture.h"
#include <memory.h>

#include <typeinfo>

namespace voreen {



template<typename T>
StreamlineTexture<T>::StreamlineTexture(const tgt::ivec2& dimensions, const bool autoFree)
    : SimpleTexture<T>(dimensions, autoFree),
    counter_(new size_t[numElements_])
{
    memset(counter_, 0, sizeof(size_t) * numElements_);
}

template<typename T>
StreamlineTexture<T>::StreamlineTexture(const tgt::ivec3& dimensions, const bool autoFree)
    : SimpleTexture<T>(dimensions, autoFree),
    counter_(new size_t[numElements_])
{
    memset(counter_, 0, sizeof(size_t) * numElements_);
}

template<typename T>
StreamlineTexture<T>::~StreamlineTexture() {
    if (autoFree_ == true)
        free();
}

template<typename T>
void StreamlineTexture<T>::free() {
    delete [] counter_;
}

template<typename T>
const typename StreamlineTexture<T>::StreamlineTextureElement StreamlineTexture<T>::operator[](size_t index) const {
    if (index >= numElements_)
        index = numElements_ - 1;
    return StreamlineTextureElement(data_[index], counter_[index]);
}

template<typename T>
const typename StreamlineTexture<T>::StreamlineTextureElement StreamlineTexture<T>::operator[](const tgt::ivec2& position) const {
    size_t index = SimpleTexture<T>::positionToIndex(position);
    return operator[](index);
}

template<typename T>
const typename StreamlineTexture<T>::StreamlineTextureElement StreamlineTexture<T>::operator[](const tgt::ivec3& position) const {
    size_t index = SimpleTexture<T>::positionToIndex(position);
    return operator[](index);
}

template<typename T>
typename StreamlineTexture<T>::StreamlineTextureElement StreamlineTexture<T>::operator[](size_t index) {
    if (index >= numElements_)
        index = numElements_ - 1;
    return StreamlineTextureElement(data_[index], counter_[index]);
}

template<typename T>
typename StreamlineTexture<T>::StreamlineTextureElement StreamlineTexture<T>::operator[](const tgt::ivec2& position) {
    size_t index = SimpleTexture<T>::positionToIndex(position);
    return operator[](index);
}

template<typename T>
typename StreamlineTexture<T>::StreamlineTextureElement StreamlineTexture<T>::operator[](const tgt::ivec3& position) {
    size_t index = SimpleTexture<T>::positionToIndex(position);
    return operator[](index);
}

// ----------------------------------------------------------------------------

template<typename T>
T* StreamlineTexture<T>::integrateDraw(const Flow2D& flow, const size_t textureScaling,
                               const size_t sampling, const tgt::vec2& thresholds)
{
    StreamlineTexture<T> outputTexture(flow.dimensions_ * static_cast<int>(textureScaling));
    const tgt::ivec2& outputTexSize = outputTexture.getDimensions().xy();
    const int delta = static_cast<int>(sampling);
    const float stepSize = (0.5f / textureScaling);
    const float length = static_cast<float>(tgt::max(flow.dimensions_)) * 1.10f;

    size_t numStreamlines = 0;

    for (int y = 0; y < outputTexSize.y; y += delta) {
        for (int x = 0; x < outputTexSize.x; x += delta) {
            tgt::ivec2 ir0(x, y);
            if (outputTexture[ir0].counter_ > 0)
                continue;

            tgt::ivec2 error(0, 0);
            tgt::vec2 r0 = flow.slicePosToFlowPos(ir0, outputTexSize, &error);
            if (flow.lookupFlow(r0) == tgt::vec2::zero)
                continue;

            std::deque<tgt::vec2> streamline =
                FlowMath::computeStreamlineRungeKutta(flow, r0, length, stepSize, 0, thresholds);
            if (streamline.size() <= 1)
                continue;

            ++numStreamlines;
            float uniRand = FlowMath::uniformRandom();
            T gray = T(uniRand);
            if (typeid(T) == typeid(unsigned char))
                gray = T(uniRand * 255.0f);

            std::deque<tgt::ivec2> streamlineVP =
                flow.flowPosToSlicePos(streamline, outputTexSize, error);

            for (std::deque<tgt::ivec2>::const_iterator it = streamlineVP.begin();
                it != streamlineVP.end(); ++it)
            {
                const tgt::ivec2& p = *it;
                outputTexture[p].elem_ += static_cast<unsigned char>((++(outputTexture[p].counter_) * gray));
            }   // for (it
        }   // for (x
    }   // for (y

    size_t unhitPixels = 0;
    const size_t numPixels = outputTexture.getNumElements();
    for (size_t i = 0; i < numPixels; ++i) {
        if (outputTexture[i].counter_ > 1)
            outputTexture[i].elem_ /= T(outputTexture[i].counter_ * (outputTexture[i].counter_ + 1) / 2.0f);
        else if (outputTexture[i].counter_ <= 0)
            ++unhitPixels;
    }
    std::cout << "#streamlines = " << numStreamlines << ", #unhit pixels = " << unhitPixels
        << " (" << static_cast<float>(100 * unhitPixels) / static_cast<float>(numPixels) << " %)\n";

    T* result = new T[numPixels];
    memcpy(result, outputTexture.getData(), numPixels * sizeof(T));
    return result;
}

template<typename T>
T* StreamlineTexture<T>::integrateDraw(const Flow3D& flow, const size_t textureScaling,
                                       const size_t sampling, const tgt::vec2& thresholds)
{
    StreamlineTexture<T> outputTexture(flow.dimensions_ * static_cast<int>(textureScaling));
    const tgt::ivec3& outputTexSize = outputTexture.getDimensions();
    const int delta = static_cast<int>(sampling);
    const float stepSize = (0.5f / textureScaling);
    const float length = static_cast<float>(tgt::max(flow.dimensions_)) * 1.10f;

    size_t numStreamlines = 0;

    for (int z = 0; z < outputTexSize.z; z += delta) {
        for (int y = 0; y < outputTexSize.y; y += delta) {
            for (int x = 0; x < outputTexSize.x; x += delta) {
                tgt::ivec3 ir0(x, y, z);
                if (outputTexture[ir0].counter_ > 0)
                    continue;

                tgt::ivec3 error(0, 0, 0);
                tgt::vec3 r0 = flow.toFlowPosition(ir0, outputTexSize, &error);
                if (flow.lookupFlow(r0) == tgt::vec3::zero)
                    continue;

                std::deque<tgt::vec3> streamline =
                    FlowMath::computeStreamlineRungeKutta(flow, r0, length, stepSize, 0, thresholds);
                if (streamline.size() <= 1)
                    continue;

                ++numStreamlines;
                float uniRand = FlowMath::uniformRandom();
                T gray = T(uniRand);
                if (typeid(T) == typeid(unsigned char))
                    gray = T(uniRand * 255.0f);

                std::deque<tgt::ivec3> streamlineTexture =
                    flow.toTexturePosition(streamline, outputTexSize, error);

                for (std::deque<tgt::ivec3>::const_iterator it = streamlineTexture.begin();
                    it != streamlineTexture.end(); ++it)
                {
                    const tgt::ivec3& p = *it;
                    outputTexture[p].elem_ += static_cast<unsigned char>(++(outputTexture[p].counter_) * gray);
                }   // for (it
            }   // for (x
        }   // for (y
    } // for (z

    size_t unhitPixels = 0;
    const size_t numVoxels = outputTexture.getNumElements();
    for (size_t i = 0; i < numVoxels; ++i) {
        if (outputTexture[i].counter_ > 1)
            outputTexture[i].elem_ /= T(outputTexture[i].counter_ * (outputTexture[i].counter_ + 1) / 2);
        else if (outputTexture[i].counter_ <= 0)
            ++unhitPixels;
    }
    std::cout << "# streamlines = " << numStreamlines << ", # unhit voxels = " << unhitPixels
        << " (" << static_cast<float>(100 * unhitPixels) / static_cast<float>(numVoxels) << " %)\n";

    // as the texture will be destroyed when this methods returns, the
    // data need to be copied and returned.
    //
    T* result = new T[numVoxels];
    memcpy(result, outputTexture.getData(), numVoxels * sizeof(T));
    return result;
}

// ----------------------------------------------------------------------------

template<typename T>
T* StreamlineTexture<T>::fastLIC(const Flow2D& flow, const SimpleTexture<float>& inputTexture,
                                 const size_t textureScaling, const size_t sampling,
                                 const int maxKernelSize, const tgt::vec2& thresholds,
                                 const bool useAdaptiveKernelSize)
{
    const tgt::ivec2 inputTexSize(inputTexture.getDimensions().xy());
    if (flow.dimensions_ != inputTexSize)
        return 0;

    StreamlineTexture<T> outputTexture(inputTexSize * static_cast<int>(textureScaling));
    const tgt::ivec2& outputTexSize = outputTexture.getDimensions().xy();
    const int delta = static_cast<int>(sampling);
    const float stepSize = (0.5f / textureScaling);
    const float length = static_cast<float>(tgt::max(flow.dimensions_)) * 1.10f;

    size_t numStreamlines = 0;
    for (int y = 0; y < outputTexSize.y; y += delta) {
        for (int x = 0; x < outputTexSize.x; x += delta) {
            tgt::ivec2 ir0(x, y);
            if (outputTexture[ir0].counter_ > 0)
                continue;

            // get the coordinates of the pixel in the input texture which corresponds
            // to this position in the output texture and calculate its position within
            // the flow.
            //
            tgt::ivec2 r0Input(ir0 / static_cast<int>(textureScaling));
            tgt::ivec2 errorInput(0, 0);
            tgt::vec2 r0 = flow.slicePosToFlowPos(r0Input, inputTexSize, &errorInput);
            const tgt::vec2 v = flow.lookupFlow(r0);

            if (v == tgt::vec2::zero)
                continue;

            // start streamline computation
            //
            int indexR0 = 0;
            std::deque<tgt::vec2> streamlineD =
                FlowMath::computeStreamlineRungeKutta(flow, r0, length, stepSize, &indexR0, thresholds);
            if (streamlineD.size() <= 1)
                continue;

            // also determine the round-off error which occurs if the flow positions was
            // converted back directly to the coordinates of the output textures.
            //
            tgt::ivec2 errorOutput(0, 0);
            flow.slicePosToFlowPos(ir0, outputTexSize, &errorOutput);

            ++numStreamlines;
            std::vector<tgt::vec2> streamline = FlowMath::dequeToVector(streamlineD);

            // copy the streamline for second coordinate conversion
            //
            std::vector<tgt::vec2> streamlineCopy(streamline);

            // convert the streamline into dimensions of the input texture
            //
            std::vector<tgt::ivec2> streamlineInput =
                flow.flowPosToSlicePos(streamline, inputTexSize, errorInput);

            // also convert the streamline into dimensions of the output texture
            //
            std::vector<tgt::ivec2> streamlineOutput =
                flow.flowPosToSlicePos(streamlineCopy, outputTexSize, errorOutput);

            // calculate initial intensity for the starting pixel
            //
            int L = maxKernelSize;
            if (useAdaptiveKernelSize == true)
                L = static_cast<int>(tgt::round(maxKernelSize * (tgt::length(v) / thresholds.y)));
            const float k = 1.0f / (((2 * L) + 1));
            float intensity0 = k * fastLICIntensity(inputTexture, indexR0, L, streamlineInput);

            // determine the affected pixel in the output texture and add the
            // initial intensity
            //
            tgt::ivec2& outputTexCoord = streamlineOutput[indexR0];
            outputTexture[outputTexCoord].elem_ += T(intensity0);
            ++(outputTexture[outputTexCoord].counter_);

            // trace streamline in forward direction and update intensity
            //
            float intensity = intensity0;
            int left = indexR0 + L + 1;
            int right = indexR0 - L;
            const int numPoints = static_cast<int>(streamlineInput.size());

            for (int i = (indexR0 + 1); i < numPoints; ++i, ++left, ++right) {
                int l = (left >= numPoints) ? (numPoints - 1) : left;
                int r = (right <= 0) ? 0 : right;

                intensity += (inputTexture[streamlineInput[l]] - inputTexture[streamlineInput[r]]) * k;

                outputTexCoord = streamlineOutput[i];
                outputTexture[outputTexCoord].elem_ += T(intensity);
                ++(outputTexture[outputTexCoord].counter_);
            }   // for (i

            // trace streamline in backward direction and update intensity
            //
            intensity = intensity0;
            left = indexR0 - L - 1;
            right = indexR0 + L;
            for (int i = (indexR0 - 1); i >= 0; --i, --left, --right) {
                int l = (left <= 0) ? 0 : left;
                int r = (right >= numPoints) ? (numPoints - 1) : right;

                intensity += (inputTexture[streamlineInput[l]] - inputTexture[streamlineInput[r]]) * k;

                outputTexCoord = streamlineOutput[i];
                outputTexture[outputTexCoord].elem_ += T(intensity);
                ++(outputTexture[outputTexCoord].counter_);
            }   // for (i
        }   // for (x
    }   // for (y

    size_t unhitPixels = 0;
    const size_t numPixels = outputTexture.getNumElements();
    for (size_t i = 0; i < numPixels; ++i) {
        if (outputTexture[i].counter_ > 1)
            outputTexture[i].elem_ /= static_cast<unsigned char>(outputTexture[i].counter_);
        else
            if (outputTexture[i].counter_ <= 0)
                ++unhitPixels;
    }
    std::cout << "# streamlines = " << numStreamlines << ", # unhit pixels = " << unhitPixels
        << " (" << static_cast<float>(100 * unhitPixels) / static_cast<float>(numPixels) << " %)\n";

    T* result = new T[numPixels];
    memcpy(result, outputTexture.getData(), sizeof(T) * numPixels);
    return result;
}

// private methods
//
template<typename T>
float StreamlineTexture<T>::fastLICIntensity(const SimpleTexture<float>& inputTexture,
                                         const int indexR0, const int kernelSize,
                                         const std::vector<tgt::ivec2>& streamline)
{
    float intensity = inputTexture[streamline[indexR0]];
    for (int n = 1; n <= kernelSize; ++n) {
        if (static_cast<size_t>(indexR0 + n) < streamline.size())
            intensity += inputTexture[streamline[indexR0 + n]];

        if ((indexR0 - n) >= 0)
            intensity += inputTexture[streamline[indexR0 - n]];
    }
    return intensity;
}

// ----------------------------------------------------------------------------

template class StreamlineTexture<float>;
template class StreamlineTexture<unsigned char>;

}   // namespace
