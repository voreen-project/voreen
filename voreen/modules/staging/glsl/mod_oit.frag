/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany.                        *
 * Visualization and Computer Graphics Group <http://viscg.uni-muenster.de>        *
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

#extension GL_ARB_shading_language_420pack : require
#extension GL_ARB_shader_atomic_counters : require
#extension GL_ARB_shader_image_load_store : require
#extension GL_ARB_shader_storage_buffer_object : require
#extension GL_ARB_gpu_shader5 : require

struct LinkedListStruct {
    uint color_; // we pack the color into 4x8 bit / stores normal for proxy geometries
    float depth_;
    uint next_;
    uint proxyGeometryId_;
};

vec4 getColor(in LinkedListStruct lls) {
    return unpackUnorm4x8(lls.color_);
}

void setColor(inout LinkedListStruct lls, in vec4 c) {
    lls.color_ = packUnorm4x8(c);
}

vec3 getNormal(in LinkedListStruct lls) {
    vec3 temp;
    temp = unpackUnorm4x8(lls.color_).xyz;
    temp -= 0.5;
    temp *= 2.0;
    return temp;
}

void setNormal(inout LinkedListStruct lls, in vec3 n) {
    vec3 temp = n;
    temp *= 0.5;
    temp += 0.5;
    lls.color_ = packUnorm4x8(vec4(clamp(temp, vec3(0.0), vec3(1.0)), 1.0));
}

// TODO define bindings
coherent uniform layout (binding=3, r32ui) uimage2D headOffset;  // image containing head pointer
//coherent restrict readonly uniform layout (binding = 3, r32ui) uimage2D headOffset;  // image containing head pointer

coherent layout(std140, binding=4) buffer LinkedList {
    LinkedListStruct linkedList_[];
};

struct ProxyGeometryStruct {
    mat4 textureToWorldMatrix_;
    mat4 worldToTextureMatrix_;
    vec3 cameraPositionPhysical_;
    int volumeId_;
    vec3 lightPositionPhysical_;
};

coherent layout(std140, binding=5) buffer ProxyGeometries {
    ProxyGeometryStruct proxyGeometries_[];
};

#ifdef NUM_PAGES
layout (binding = 0, offset = 0) uniform atomic_uint counter[NUM_PAGES]; // atomic counter
#else
layout (binding = 0, offset = 0) uniform atomic_uint counter; // atomic counter
#endif

#ifdef NUM_PAGES
uint getNewOffset(uint page) {
    uint newOffset = atomicCounterIncrement(counter[page]);
#else
uint getNewOffset() {
    uint newOffset = atomicCounterIncrement(counter);
#endif
    if(newOffset > PAGE_SIZE)
        discard;

#ifdef NUM_PAGES
    newOffset += (page * PAGE_SIZE);
    //newOffset = newOffset + (page * uint(PAGE_SIZE));
    //newOffset = (newOffset * NUM_PAGES) + page;
#endif
    return newOffset;
}

void pushFront(LinkedListStruct lls) {
    ivec2 coord = ivec2(gl_FragCoord.xy - vec2(0.5));

#ifdef NUM_PAGES
    uint page = uint(mod(coord.x, float(NUM_PAGES)));
    //int page = (int(coord.y) % NUM_PAGES);

    uint newOffset = getNewOffset(page);
#else
    uint newOffset = getNewOffset();
#endif

    uint prevOffset = imageAtomicExchange(headOffset, coord, newOffset);

    lls.next_ = prevOffset;
    linkedList_[newOffset] = lls;
}

