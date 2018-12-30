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

#version 330

#extension GL_ARB_shading_language_420pack : require
#extension GL_ARB_shader_atomic_counters : require
#extension GL_ARB_shader_image_load_store : require
#extension GL_ARB_gpu_shader5 : require

// Head pointer 2D buffer
// uniform usampler2D headPointerImage_;
layout (binding = 0, r32ui) readonly uniform uimage2D headPointerImage_;

// Linked list 1D buffer
//uniform usamplerBuffer listBuffer_;
layout (binding = 1, rgba32ui) readonly uniform uimageBuffer listBuffer_;

// Small buffer to hold all of the fragments corresponding to this pixel
uvec4 fragments[MAX_FRAGMENTS];

layout (location = 0) out vec4 output_color;

// Traverse the linked list, place all of the fragments into the
// fragments[] array and return the number of fragments retrieved
// from the list.
int build_local_fragment_list(void)
{
    uint current;
    int frag_count = 0;
    // Get the initial head pointer from the header-pointer image
    current = imageLoad(headPointerImage_, ivec2(gl_FragCoord.xy)).x;
    // While we haven’t reached the end of the list or exhausted
    // the storage available in fragments[]...
    while (current != 0xFFFFFFFFu && frag_count < MAX_FRAGMENTS)
    {
        // Read an item from the linked list
        uvec4 item = imageLoad(listBuffer_, int(current));
        // item.x contains the "next" pointer - update current
        current = item.x;
        // Store the fragment in the array
        fragments[frag_count] = item;
        // Update the fragment count
        frag_count++;
    }
    // Done - return the fragment count
    return frag_count;
}

// Simple bubble-sort for sorting the fragments[] array
// After this fragments[] will be sorted nearest to farthest from front to back
void sort_fragment_list(int frag_count)
{
    for (int i = 0; i < frag_count; i++)
    {
        for (int j = i + 1; j < frag_count; j++)
        {
            // The depth of each fragment is bit-encoded into the
            // .w channel of the fragment array. Unpack it here.
            float depth_i = uintBitsToFloat(fragments[i].w);
            float depth_j = uintBitsToFloat(fragments[j].w);
            // Compare depth and if the comparison fails...
            if (depth_i > depth_j)
            {
                // Swap the fragments in the array
                uvec4 temp = fragments[i];
                fragments[i] = fragments[j];
                fragments[j] = temp;
            }
        }
    }
}

// Alpha blending with correct final alpha
vec4 blend(vec4 current_color, vec4 new_color)
{
    //return mix(current_color, new_color, new_color.a);
    return vec4(mix(current_color.rgb, new_color.rgb, new_color.a),
        1.0f - (1.0f - new_color.a)*(1.0f - current_color.a));
}

// Function for calculating the final output color. Walks the
// fragments[] array and blends each pixel on top of each other
//
// Note: frag_count is expected to be > 0!
vec4 calculate_final_color(int frag_count)
{
    // Initialize the final color output with the farthest fragment
    vec4 final_color = vec4(0,0,0,0);

    // fragments[] are sorted nearest to farthest, so we have to blend from back to front.
    for (int i = frag_count-1; i >= 0; --i)
    {
        // The color is stored packed into the .y channel of the
        // fragment vector. Unpack it here.
        vec4 frag_color;
        frag_color.rg = unpackUnorm2x16(fragments[i].y);
        frag_color.ba = unpackUnorm2x16(fragments[i].z);
        // Now call the blending function.
        final_color = blend(final_color, frag_color);
    }
    // Undo implicit alpha premultiplication
    final_color.rgb /= final_color.a;

    return final_color;
}

void main(void)
{
    int frag_count;
    // Traverse the list and build an array of fragments
    frag_count = build_local_fragment_list();
#if 0
    if(frag_count > 0) {
        output_color = vec4(0,float(frag_count)/MAX_FRAGMENTS,0,1);
        //output_color = vec4(unpackUnorm4x8(fragments[0].y).xyz, fragments[0].z);
    } else {
        output_color = vec4(1,0,0,1);
    }
#else
    if(frag_count > 0) {
        // If we have a fragment to display:

        // Sort the array in depth order
        sort_fragment_list(frag_count);
        // Blend the sorted fragments together to compute the final output color
        output_color = calculate_final_color(frag_count);
        // The depth is the depth of the nearest fragment
        gl_FragDepth = uintBitsToFloat(fragments[0].w);
    } else {
        discard;
    }
#endif
}
