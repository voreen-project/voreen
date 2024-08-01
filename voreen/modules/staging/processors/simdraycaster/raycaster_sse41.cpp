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

/*
    This files defines the raycaster functions for sse4.1.
    It uses the build-in and hopefully fast implementation of floor.
*/
#include "tgt/vector.h"

#include <stdint.h>
#include <assert.h>

// SSE41 HEADER
#include <smmintrin.h>

namespace voreen{
namespace simdraycaster{

#define _my_floor_ps _mm_floor_ps

/************************************************************
 *                          DVR                             *
 ************************************************************/
#define BITS 8
#define FN_NAME raycaster_uint8_t_sse41
#define TYPE uint8_t
#include "raycast_generic.h"
#undef FN_NAME
#undef BITS
#undef TYPE

#define BITS 16
#define FN_NAME raycaster_uint16_t_sse41
#define TYPE uint16_t
#include "raycast_generic.h"
#undef FN_NAME
#undef BITS
#undef TYPE

#define BRICKED_VOLUME
#define BITS 8
#define FN_NAME raycaster_uint8_t_bricked_sse41
#define TYPE uint8_t
#include "raycast_generic.h"
#undef FN_NAME
#undef BITS
#undef TYPE
#undef BRICKED_VOLUME

#define BRICKED_VOLUME
#define BITS 16
#define FN_NAME raycaster_uint16_t_bricked_sse41
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
#define FN_NAME raycaster_uint8_t_preintegration_sse41
#define TYPE uint8_t
#include "raycast_generic.h"
#undef FN_NAME
#undef BITS
#undef TYPE

#define BITS 16
#define FN_NAME raycaster_uint16_t_preintegration_sse41
#define TYPE uint16_t
#include "raycast_generic.h"
#undef FN_NAME
#undef BITS
#undef TYPE

#define BRICKED_VOLUME
#define BITS 8
#define FN_NAME raycaster_uint8_t_bricked_preintegration_sse41
#define TYPE uint8_t
#include "raycast_generic.h"
#undef FN_NAME
#undef BITS
#undef TYPE
#undef BRICKED_VOLUME

#define BRICKED_VOLUME
#define BITS 16
#define FN_NAME raycaster_uint16_t_bricked_preintegration_sse41
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
#define FN_NAME raycaster_uint8_t_mip_sse41
#define TYPE uint8_t
#include "raycast_generic.h"
#undef FN_NAME
#undef BITS
#undef TYPE

#define BITS 16
#define FN_NAME raycaster_uint16_t_mip_sse41
#define TYPE uint16_t
#include "raycast_generic.h"
#undef FN_NAME
#undef BITS
#undef TYPE

#define BRICKED_VOLUME
#define BITS 8
#define FN_NAME raycaster_uint8_t_bricked_mip_sse41
#define TYPE uint8_t
#include "raycast_generic.h"
#undef FN_NAME
#undef BITS
#undef TYPE
#undef BRICKED_VOLUME

#define BRICKED_VOLUME
#define BITS 16
#define FN_NAME raycaster_uint16_t_bricked_mip_sse41
#define TYPE uint16_t
#include "raycast_generic.h"
#undef FN_NAME
#undef BITS
#undef TYPE
#undef BRICKED_VOLUME

#undef MIP
}
}