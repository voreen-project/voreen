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

#ifndef VRN_STREAMLINETEXTURE_H
#define VRN_STREAMLINETEXTURE_H

#include "modules/flowreen/datastructures/deprecated/simpletexture.h"

namespace voreen {

class Flow2D;
class Flow3D;

template<typename T>
class StreamlineTexture : public SimpleTexture<T> {

using SimpleTexture<T>::data_;
using SimpleTexture<T>::autoFree_;
using SimpleTexture<T>::numElements_;
using SimpleTexture<T>::dimensions_;

public:
    struct StreamlineTextureElement {
        StreamlineTextureElement(T& elem, size_t& counter)
            : elem_(elem), counter_(counter)
        {}

        operator T() { return elem_; }

        T& elem_;
        size_t& counter_;
    };

public:
    StreamlineTexture(const tgt::ivec2& dimensions, const bool autoFree = true);
    StreamlineTexture(const tgt::ivec3& dimensions, const bool autoFree = true);
    virtual ~StreamlineTexture();

    virtual void free();

    const typename StreamlineTexture<T>::StreamlineTextureElement operator[](size_t index) const;
    const typename StreamlineTexture<T>::StreamlineTextureElement operator[](const tgt::ivec2& position) const;
    const typename StreamlineTexture<T>::StreamlineTextureElement operator[](const tgt::ivec3& position) const;

    typename StreamlineTexture<T>::StreamlineTextureElement operator[](size_t index);
    typename StreamlineTexture<T>::StreamlineTextureElement operator[](const tgt::ivec2& position);
    typename StreamlineTexture<T>::StreamlineTextureElement operator[](const tgt::ivec3& position);

    static T* integrateDraw(const Flow2D& flow, const size_t textureScaling,
        const size_t sampling, const tgt::vec2& thresholds = tgt::vec2(0.0f));

    static T* integrateDraw(const Flow3D& flow, const size_t textureScaling,
        const size_t sampling, const tgt::vec2& thresholds = tgt::vec2(0.0f));

    static T* fastLIC(const Flow2D& flow, const SimpleTexture<float>& inputTexture,
        const size_t textureScaling, const size_t sampling, const int maxKernelSize,
        const tgt::vec2& threshold = tgt::vec2(0.0f), const bool useAdaptiveKernelSize = false);

private:
    // prevent objects of this class from being copied
    //
    StreamlineTexture(const StreamlineTexture&);
    StreamlineTexture& operator=(const StreamlineTexture&);

    static float fastLICIntensity(const SimpleTexture<float>& inputTexture,
        const int indexR0, const int kernelSize, const std::vector<tgt::ivec2>& streamline);

private:
    size_t* const counter_;
};

}   // namespace

#endif
