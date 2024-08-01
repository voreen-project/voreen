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

#version 330

// At shader compile time we decide whether the target of outputFragment will be the
// framebuffer or OIT buffers:
#ifdef USE_TRANSPARENCY

#extension GL_ARB_shading_language_420pack : enable
#extension GL_ARB_shader_atomic_counters : enable
#extension GL_ARB_shader_image_load_store : enable
#extension GL_ARB_gpu_shader5 : enable

// Append Fragments from an image to linked list
// Taken from OpenGL programming guide p. 609ff

// This is the atomic counter used to allocate items in the
// linked list
layout (binding = 0, offset = 0) uniform atomic_uint indexCounter_;
// Head pointer 2D buffer
layout (binding = 0, r32ui) uniform uimage2D headPointerImage_;
// Linked list 1D buffer
layout (binding = 1, rgba32ui) writeonly uniform uimageBuffer listBuffer_;

// Turn on early fragment testing
layout (early_fragment_tests) in;

// Output a fragment with explicit depth information
void outputFragment(vec4 color, float depth) {

    // Allocate a fragment in headPointerImage_ by reserving an index
    uint index = atomicCounterIncrement(indexCounter_);

    // Insert our new fragment into the list of the output fragment:
    // Replace the current head with the index previously reserved.
    uint old_head = imageAtomicExchange(headPointerImage_, ivec2(gl_FragCoord.xy), index);

    // Now build the fragment list item
    uvec4 item;

    // item.x = index of the previous head of the fragment list
    item.x = old_head;

    // item.yz = the color (packed into two uints => 16 bit per channel)
    item.y = packUnorm2x16(color.rg);
    item.z = packUnorm2x16(color.ba);

    // item.w = depth
    item.w = floatBitsToUint(depth);

    // item.w = unused
    //item.w = 0u;

    // Finally write the item to the previously reserved location in the fragment list
    imageStore(listBuffer_, int(index), item);
}

// Output a fragment using the depth from gl_FragCoord
void outputFragment(vec4 color) {
    outputFragment(color, gl_FragCoord.z);
}
#else

out vec4 output_color;

// Output a fragment with explicit depth information
void outputFragment(vec4 color, float depth) {
    output_color = color;
    gl_FragDepth = depth;
}

// Output a fragment using the depth from gl_FragCoord
void outputFragment(vec4 color) {
    output_color = color;
}
#endif
