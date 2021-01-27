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

/*
    This files defines the raycaster functions for sse3.
    It uses a relativly slow implementation of floor.
*/
#include "tgt/vector.h"

#include <stdint.h>
#include <assert.h>

// SSE 3 Header
#include <pmmintrin.h>
#ifndef WIN32
#include <immintrin.h>
#endif

namespace voreen{
namespace simdraycaster{


// http://dss.stephanierct.com/DevBlog/?p=8
inline __m128 _mm_floor_ps2(const __m128& x){
    __m128i v0 = _mm_setzero_si128();
    __m128i v1 = _mm_cmpeq_epi32(v0,v0);
    __m128i ji = _mm_srli_epi32( v1, 25);
    __m128 j = _mm_castsi128_ps(_mm_slli_epi32( ji, 23)); //create vector 1.0f
    __m128i i = _mm_cvttps_epi32(x);
    __m128 fi = _mm_cvtepi32_ps(i);
    __m128 igx = _mm_cmpgt_ps(fi, x);
    j = _mm_and_ps(igx, j);
    return _mm_sub_ps(fi, j);
}


#define _my_floor_ps _mm_floor_ps2

/************************************************************
 *                          DVR                             *
 ************************************************************/
#define BITS 8
#define FN_NAME raycaster_uint8_t_sse3
#define TYPE uint8_t
#include "raycast_generic.h"
#undef FN_NAME
#undef BITS
#undef TYPE

#define BITS 16
#define FN_NAME raycaster_uint16_t_sse3
#define TYPE uint16_t
#include "raycast_generic.h"
#undef FN_NAME
#undef BITS
#undef TYPE

#define BRICKED_VOLUME
#define BITS 8
#define FN_NAME raycaster_uint8_t_bricked_sse3
#define TYPE uint8_t
#include "raycast_generic.h"
#undef FN_NAME
#undef BITS
#undef TYPE
#undef BRICKED_VOLUME

#define BRICKED_VOLUME
#define BITS 16
#define FN_NAME raycaster_uint16_t_bricked_sse3
#define TYPE uint16_t
#include "raycast_generic.h"
#undef FN_NAME
#undef BITS
#undef TYPE
#undef BRICKED_VOLUME


/************************************************************
 *                  DVR Preintegration                      *
 ************************************************************/
#define PREINTEGRATION
#define BITS 8
#define FN_NAME raycaster_uint8_t_preintegration_sse3
#define TYPE uint8_t
#include "raycast_generic.h"
#undef FN_NAME
#undef BITS
#undef TYPE

#define BITS 16
#define FN_NAME raycaster_uint16_t_preintegration_sse3
#define TYPE uint16_t
#include "raycast_generic.h"
#undef FN_NAME
#undef BITS
#undef TYPE

#define BRICKED_VOLUME
#define BITS 8
#define FN_NAME raycaster_uint8_t_bricked_preintegration_sse3
#define TYPE uint8_t
#include "raycast_generic.h"
#undef FN_NAME
#undef BITS
#undef TYPE
#undef BRICKED_VOLUME

#define BRICKED_VOLUME
#define BITS 16
#define FN_NAME raycaster_uint16_t_bricked_preintegration_sse3
#define TYPE uint16_t
#include "raycast_generic.h"
#undef FN_NAME
#undef BITS
#undef TYPE
#undef BRICKED_VOLUME

/************************************************************
 *                          MIP                             *
 ************************************************************/
#define MIP
#define BITS 8
#define FN_NAME raycaster_uint8_t_mip_sse3
#define TYPE uint8_t
#include "raycast_generic.h"
#undef FN_NAME
#undef BITS
#undef TYPE

#define BITS 16
#define FN_NAME raycaster_uint16_t_mip_sse3
#define TYPE uint16_t
#include "raycast_generic.h"
#undef FN_NAME
#undef BITS
#undef TYPE

#define BRICKED_VOLUME
#define BITS 8
#define FN_NAME raycaster_uint8_t_bricked_mip_sse3
#define TYPE uint8_t
#include "raycast_generic.h"
#undef FN_NAME
#undef BITS
#undef TYPE
#undef BRICKED_VOLUME

#define BRICKED_VOLUME
#define BITS 16
#define FN_NAME raycaster_uint16_t_bricked_mip_sse3
#define TYPE uint16_t
#include "raycast_generic.h"
#undef FN_NAME
#undef BITS
#undef TYPE
#undef BRICKED_VOLUME

#undef MIP
}
}
