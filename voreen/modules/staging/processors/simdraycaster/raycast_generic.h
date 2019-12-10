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

//#ifndef VRN_RAYCASTER_GENERIC
//#define VRN_RAYCASTER_GENERIC

/* The includeguards are commented out, because this file is included multiple
   times with different defines to generate code for the different modes
   of the raycaster. This happens from raycaster_sse3.cpp und raycaster_sse41.cpp

   The macros that must be defined are:
   - FN_NAME
   - BITS
   - BRICKED_VOLUME
   - MIP
   - PREINTEGRATION
   
   In the scope there must alao be a function _my_floor_ps, which should
   have the same semantics as _mm_floor_ps in newer sse version.
*/
using namespace tgt;

// convenience macros
#define SSE_GET_COMP(vec, comp) ((float*)&(vec))[comp]
#define swizzle4(ssea, sseb, x,y,z,w) _mm_shuffle_ps(ssea, sseb, _MM_SHUFFLE(w,z,y,x) )

#ifndef RAYCASTERPARAMTERS_DEFINED
#define RAYCASTERPARAMTERS_DEFINED
// Configuration of the raycaster
struct RayCasterParameters{
    int offset;
    int linelength;
    int count;
    int stride;
    vec4 *entryBuffer;
    vec4 *exitBuffer;
    int bytesPerVoxel;
    vec4* tf;
    int tfsize;
    float tfscale;
    float tfoffset;
    uint32_t* img;
    ivec3 volDim;
    const void* vol;
    float stepsize;
    const size_t* offsetx;
    const size_t* offsety;
    const size_t* offsetz;
};
#endif

/**
 * Definition of the main raycasting function.
 * This is mostly carefully optimized code
 */
void FN_NAME(RayCasterParameters p){
    // Get settings in local variables, for performance reasons on some compilers
    size_t      offset    = p.offset;
    vec4*       entrybuf  = p.entryBuffer;
    vec4*       exitbuf   = p.exitBuffer;
    vec3        voldim1   = vec3(p.volDim)-vec3(1);
    size_t      linecount = p.linelength;
    const TYPE* vol       = (TYPE*)p.vol;
    ivec3       voldim    = p.volDim;
    uint32_t *  pixbuf    = p.img;
    vec4*       tf        = p.tf;
    int         tfsize    = p.tfsize;
    float       stepsize  = p.stepsize;

    // load volume configuration into sse registers
    const __m128 voldim_f4_p1 = _mm_set_ps(0.0f, voldim1.z+1.0f, voldim1.y+1.0f, voldim1.x+1.0f);
    const __m128 voldim_f4_m1 = _mm_set_ps(0.0f, voldim1.z-1.0f, voldim1.y-1.0f, voldim1.x-1.0f);
    const __m128 voloffset    = _mm_set_ps(0.0f, -0.5f, -0.5f, -0.5f);

    // Initialize constants
    const __m128 _1 = _mm_set_ps1( 1.f );
    const __m128 _0 = _mm_setzero_ps();

    // transfer function configuration
    const __m128 tfsizereg   = _mm_set_ps(0, 0, 0, 1.0f*tfsize-1);
    const __m128 tfoffsetreg = _mm_set_ps(0, 0, 0, p.tfoffset);
#if BITS==8
    __m128 tfnormalizationfactor = _mm_set1_ps(1.0f*(tfsize-1)*p.tfscale/0xff);
#elif BITS==16
    __m128 tfnormalizationfactor = _mm_set1_ps(1.0f*(tfsize-1)*p.tfscale/0xffff);
#else
#error
#endif

    // offset for the current voxel in all buffers
    size_t bufferindex = offset;

    // configuration for loading data from the volume.
#ifdef BRICKED_VOLUME
    const size_t* offsetx = p.offsetx;
    const size_t* offsety = p.offsety;
    const size_t* offsetz = p.offsetz;
#else
    size_t f1 = (size_t)voldim.x*(size_t)voldim.y;
    size_t f2 = voldim.x;
#endif

    // Loop line by line of block
    for(int line = 0; line != p.count; line++){
    for(int row = 0; row != linecount; row++)
    {
        // Load entry and exit points
        __m128 entry = _mm_add_ps(_mm_mul_ps(_mm_load_ps(&entrybuf[bufferindex].x), voldim_f4_p1), voloffset);
        __m128 exit = _mm_add_ps(_mm_mul_ps(_mm_load_ps(&exitbuf[bufferindex].x), voldim_f4_p1), voloffset);

        // test if entry == exit
        __m128 isequal = _mm_cmpeq_ps(exit, entry);
        if ((_mm_movemask_ps(isequal)&7)==7){
            // nothing to to. Write black transparent to the buffer
            pixbuf[bufferindex] = 0x00000000;
        }else{
            // Here begins the initialization of the ray
            __m128 raycastingPosition = entry;
            
            __m128 raycastingRange = _mm_sub_ps(exit, entry);

            // distance between entry and exit point
            __m128 len = _mm_mul_ps(raycastingRange, raycastingRange);
            len = _mm_hadd_ps(len, len);
            len = _mm_hadd_ps(len, len);
            len = _mm_sqrt_ps(len);

            float lens = SSE_GET_COMP(len, 0);

            // calculate number of setps
            int steps = (int)(lens/stepsize);

            // steps should never be 0!
            steps += (steps==0);

            float isteps = 1.0f/steps;
            __m128 stepAdvanceVector = _mm_div_ps(_mm_mul_ps(raycastingRange, _mm_set1_ps(stepsize)), len);

            // initialize results
            #ifdef MIP
            __m128 maxintensity = _mm_set1_ps(0.0f);
            #else
            __m128 result = _0;
            int prev = 0;
            #endif

            float * tffloat = &tf[0].x;

            // Raycasting loop
            #ifndef MIP
            for(int j = 0; j != steps && SSE_GET_COMP(result, 3) < 0.95; j++) {
            #else
            // Don't use early ray termination for MIP-Rendering
            for(int j = 0; j != steps; j++) {
            #endif
                // calculate coordinate for sampling and 
                // make sure it is in range for sampling
                // on advanced in every direction
                __m128 interpolationPos = _my_floor_ps(raycastingPosition);
                interpolationPos = _mm_min_ps(interpolationPos, voldim_f4_m1);
                interpolationPos = _mm_max_ps(interpolationPos, _mm_setzero_ps());

                // calculate linear interpolation weights
                __m128 weights = _mm_sub_ps(raycastingPosition, interpolationPos);
                __m128 _1mweights =_mm_sub_ps(_1, weights);

                __m128i interpolationPosInt = _mm_cvtps_epi32(interpolationPos);
                
                // Get components of the position
                size_t x = _mm_extract_epi32(interpolationPosInt, 0);
                size_t y = _mm_extract_epi32(interpolationPosInt, 1);
                size_t z = _mm_extract_epi32(interpolationPosInt, 2);

                // Generate registers with weights in efficient order
                __m128 z0011  = swizzle4(weights, _1mweights, 2, 2, 2, 2); // z z 1-z 1-z
                __m128 z0101 = swizzle4(z0011, z0011, 0, 2, 0, 2); // z 1-z z 1-z
                __m128 y0011 = swizzle4(weights, _1mweights, 1, 1, 1, 1); // y y 1-y 1-y
                __m128 x1111 = swizzle4(_1mweights, _1mweights, 0, 0, 0, 0); // 1-x 1-x 1-x 1-x 
                __m128 x0000 = swizzle4(weights, weights, 0, 0, 0, 0); // x x x x 

                // load data from volume
#ifndef BRICKED_VOLUME

                // linear volume
                size_t addr = x+y*f2+z*f1;
                //assert(addr < hmul(voldim));
                //assert(addr+f1+f2 < hmul(voldim));
                //assert(addr >= 0);
#if BITS == 8 || BITS == 16
                __m128 v1 = _mm_cvtepi32_ps(_mm_set_epi32(vol[addr], vol[addr+f1],
                    vol[addr+f2], vol[addr+f2+f1]));

                __m128 v2 = _mm_cvtepi32_ps(_mm_set_epi32(vol[addr+1], vol[addr+f1+1],
                    vol[addr+f2+1], vol[addr+f1+f2+1]));
#else
#error
#endif

#else
                // bricked volume
                size_t x0 = offsetx[x];
                size_t x1 = offsetx[x+1];

                size_t y0 = offsety[y];
                size_t y1 = offsety[y+1];

                size_t z0 = offsetz[z];
                size_t z1 = offsetz[z+1];
                __m128 v1 = _mm_cvtepi32_ps(_mm_set_epi32(vol[x0+y0+z0], vol[x0+y0+z1],
                    vol[x0+y1+z0], vol[x0+y1+z1]));

                __m128 v2 = _mm_cvtepi32_ps(_mm_set_epi32(vol[x1+y0+z0], vol[x1+y0+z1],
                    vol[x1+y1+z0], vol[x1+y1+z1]));
#endif
                // Do interpolation
                __m128 yz = _mm_mul_ps(y0011, z0101);
                yz = _mm_mul_ps(yz, tfnormalizationfactor);
                __m128 resultForSample = _mm_mul_ps(yz, _mm_add_ps(_mm_mul_ps(v1, x1111),
                                                     _mm_mul_ps(v2, x0000)));
                resultForSample = _mm_hadd_ps(resultForSample, resultForSample);
                resultForSample = _mm_hadd_ps(resultForSample, resultForSample);

                // apply transfer function
#ifndef MIP
                // apply real world mapping and domain adjustment 
                // real world mapping and domain
                resultForSample = _mm_add_ss(resultForSample, tfoffsetreg);
                resultForSample = _mm_max_ss(resultForSample, _mm_setzero_ps());
                resultForSample = _mm_min_ss(resultForSample, tfsizereg);

                // load color from tf
                float tfcoord = SSE_GET_COMP(resultForSample, 0);
#ifndef PREINTEGRATION
                // without preintegration
                __m128 tfcolor = _mm_load_ps(tffloat+4*(int)tfcoord);
#else
                // with preintegration
                int ptxcoord = (prev*tfsize+(int)tfcoord);
                //assert(ptxcoord < tfsize*tfsize);
                __m128 tfcolor = _mm_load_ps(tffloat+4*ptxcoord);
                prev = (int)tfcoord;
#endif
                // do blending for dvr
                __m128 _1_minus_result_a = _mm_sub_ps(_1, swizzle4(result, result, 3, 3, 3, 3));
                resultForSample = _mm_add_ss(resultForSample, tfoffsetreg);

                result = _mm_add_ps(result, _mm_mul_ps(_1_minus_result_a, tfcolor));

#else
                // update intensity for mip
                maxintensity = _mm_max_ps(maxintensity, resultForSample);
                
#endif
                // advance to the next sampling position
                raycastingPosition = _mm_add_ps(raycastingPosition, stepAdvanceVector);
            }
#ifdef MIP
            // tf lookup for mip
            // real world mapping and domain
            maxintensity = _mm_add_ss(maxintensity, tfoffsetreg);
            maxintensity = _mm_max_ss(maxintensity, _mm_setzero_ps());
            maxintensity = _mm_min_ss(maxintensity, tfsizereg);
            
            // load color from tf
            float tfcoord = SSE_GET_COMP(maxintensity, 0);
            //assert(tfcoord >= 0 && tfcoord < tfsize);
            __m128 result = _mm_load_ps(tffloat+4*(int)tfcoord);
#endif
            // convert color to RGBA uint8
            __m128i color_as_int = _mm_cvtps_epi32(_mm_mul_ps(result, _mm_set1_ps(255.0f)));
            color_as_int = _mm_packs_epi32(color_as_int, color_as_int);
            color_as_int = _mm_packus_epi16(color_as_int, color_as_int);

            // write to buffer
            _mm_store_ss(reinterpret_cast<float*>(&pixbuf[bufferindex]), _mm_castsi128_ps(color_as_int));
        }
        // advance in line of pixelblock
        bufferindex+=1;
    }
    // advance in row of pixelblock
    bufferindex+=p.stride;
    }
}
//#endif

