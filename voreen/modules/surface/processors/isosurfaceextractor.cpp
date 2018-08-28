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

#include "isosurfaceextractor.h"

#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/datastructures/volume/volumeram.h"
#include "voreen/core/datastructures/volume/volumeminmax.h"
#include "voreen/core/datastructures/volume/volumefactory.h"
#include "voreen/core/datastructures/geometry/glmeshgeometry.h"
#include "voreen/core/datastructures/callback/lambdacallback.h"

#include "tgt/vector.h"
#include "tgt/init.h"
#include "tgt/tgt_gl.h"
#include "tgt/tgt_math.h"

#ifdef VRN_MODULE_OPENMP
    #include "omp.h"
#endif

namespace {
// Lookup tables taken from http://paulbourke.net/geometry/polygonise/marchingsource.cpp (public domain)
const uint32_t EDGE_INDICES[256]={
    0x0  , 0x109, 0x203, 0x30a, 0x406, 0x50f, 0x605, 0x70c,
    0x80c, 0x905, 0xa0f, 0xb06, 0xc0a, 0xd03, 0xe09, 0xf00,
    0x190, 0x99 , 0x393, 0x29a, 0x596, 0x49f, 0x795, 0x69c,
    0x99c, 0x895, 0xb9f, 0xa96, 0xd9a, 0xc93, 0xf99, 0xe90,
    0x230, 0x339, 0x33 , 0x13a, 0x636, 0x73f, 0x435, 0x53c,
    0xa3c, 0xb35, 0x83f, 0x936, 0xe3a, 0xf33, 0xc39, 0xd30,
    0x3a0, 0x2a9, 0x1a3, 0xaa , 0x7a6, 0x6af, 0x5a5, 0x4ac,
    0xbac, 0xaa5, 0x9af, 0x8a6, 0xfaa, 0xea3, 0xda9, 0xca0,
    0x460, 0x569, 0x663, 0x76a, 0x66 , 0x16f, 0x265, 0x36c,
    0xc6c, 0xd65, 0xe6f, 0xf66, 0x86a, 0x963, 0xa69, 0xb60,
    0x5f0, 0x4f9, 0x7f3, 0x6fa, 0x1f6, 0xff , 0x3f5, 0x2fc,
    0xdfc, 0xcf5, 0xfff, 0xef6, 0x9fa, 0x8f3, 0xbf9, 0xaf0,
    0x650, 0x759, 0x453, 0x55a, 0x256, 0x35f, 0x55 , 0x15c,
    0xe5c, 0xf55, 0xc5f, 0xd56, 0xa5a, 0xb53, 0x859, 0x950,
    0x7c0, 0x6c9, 0x5c3, 0x4ca, 0x3c6, 0x2cf, 0x1c5, 0xcc ,
    0xfcc, 0xec5, 0xdcf, 0xcc6, 0xbca, 0xac3, 0x9c9, 0x8c0,
    0x8c0, 0x9c9, 0xac3, 0xbca, 0xcc6, 0xdcf, 0xec5, 0xfcc,
    0xcc , 0x1c5, 0x2cf, 0x3c6, 0x4ca, 0x5c3, 0x6c9, 0x7c0,
    0x950, 0x859, 0xb53, 0xa5a, 0xd56, 0xc5f, 0xf55, 0xe5c,
    0x15c, 0x55 , 0x35f, 0x256, 0x55a, 0x453, 0x759, 0x650,
    0xaf0, 0xbf9, 0x8f3, 0x9fa, 0xef6, 0xfff, 0xcf5, 0xdfc,
    0x2fc, 0x3f5, 0xff , 0x1f6, 0x6fa, 0x7f3, 0x4f9, 0x5f0,
    0xb60, 0xa69, 0x963, 0x86a, 0xf66, 0xe6f, 0xd65, 0xc6c,
    0x36c, 0x265, 0x16f, 0x66 , 0x76a, 0x663, 0x569, 0x460,
    0xca0, 0xda9, 0xea3, 0xfaa, 0x8a6, 0x9af, 0xaa5, 0xbac,
    0x4ac, 0x5a5, 0x6af, 0x7a6, 0xaa , 0x1a3, 0x2a9, 0x3a0,
    0xd30, 0xc39, 0xf33, 0xe3a, 0x936, 0x83f, 0xb35, 0xa3c,
    0x53c, 0x435, 0x73f, 0x636, 0x13a, 0x33 , 0x339, 0x230,
    0xe90, 0xf99, 0xc93, 0xd9a, 0xa96, 0xb9f, 0x895, 0x99c,
    0x69c, 0x795, 0x49f, 0x596, 0x29a, 0x393, 0x99 , 0x190,
    0xf00, 0xe09, 0xd03, 0xc0a, 0xb06, 0xa0f, 0x905, 0x80c,
    0x70c, 0x605, 0x50f, 0x406, 0x30a, 0x203, 0x109, 0x0
};
const int32_t TRIANGLE_VERTEX_INDICES[256][16] =
{{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 1, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 8, 3, 9, 8, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 8, 3, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {9, 2, 10, 0, 2, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {2, 8, 3, 2, 10, 8, 10, 9, 8, -1, -1, -1, -1, -1, -1, -1},
    {3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 11, 2, 8, 11, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 9, 0, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 11, 2, 1, 9, 11, 9, 8, 11, -1, -1, -1, -1, -1, -1, -1},
    {3, 10, 1, 11, 10, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 10, 1, 0, 8, 10, 8, 11, 10, -1, -1, -1, -1, -1, -1, -1},
    {3, 9, 0, 3, 11, 9, 11, 10, 9, -1, -1, -1, -1, -1, -1, -1},
    {9, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {4, 3, 0, 7, 3, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 1, 9, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {4, 1, 9, 4, 7, 1, 7, 3, 1, -1, -1, -1, -1, -1, -1, -1},
    {1, 2, 10, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {3, 4, 7, 3, 0, 4, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1},
    {9, 2, 10, 9, 0, 2, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1},
    {2, 10, 9, 2, 9, 7, 2, 7, 3, 7, 9, 4, -1, -1, -1, -1},
    {8, 4, 7, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {11, 4, 7, 11, 2, 4, 2, 0, 4, -1, -1, -1, -1, -1, -1, -1},
    {9, 0, 1, 8, 4, 7, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1},
    {4, 7, 11, 9, 4, 11, 9, 11, 2, 9, 2, 1, -1, -1, -1, -1},
    {3, 10, 1, 3, 11, 10, 7, 8, 4, -1, -1, -1, -1, -1, -1, -1},
    {1, 11, 10, 1, 4, 11, 1, 0, 4, 7, 11, 4, -1, -1, -1, -1},
    {4, 7, 8, 9, 0, 11, 9, 11, 10, 11, 0, 3, -1, -1, -1, -1},
    {4, 7, 11, 4, 11, 9, 9, 11, 10, -1, -1, -1, -1, -1, -1, -1},
    {9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {9, 5, 4, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 5, 4, 1, 5, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {8, 5, 4, 8, 3, 5, 3, 1, 5, -1, -1, -1, -1, -1, -1, -1},
    {1, 2, 10, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {3, 0, 8, 1, 2, 10, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1},
    {5, 2, 10, 5, 4, 2, 4, 0, 2, -1, -1, -1, -1, -1, -1, -1},
    {2, 10, 5, 3, 2, 5, 3, 5, 4, 3, 4, 8, -1, -1, -1, -1},
    {9, 5, 4, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 11, 2, 0, 8, 11, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1},
    {0, 5, 4, 0, 1, 5, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1},
    {2, 1, 5, 2, 5, 8, 2, 8, 11, 4, 8, 5, -1, -1, -1, -1},
    {10, 3, 11, 10, 1, 3, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1},
    {4, 9, 5, 0, 8, 1, 8, 10, 1, 8, 11, 10, -1, -1, -1, -1},
    {5, 4, 0, 5, 0, 11, 5, 11, 10, 11, 0, 3, -1, -1, -1, -1},
    {5, 4, 8, 5, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1},
    {9, 7, 8, 5, 7, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {9, 3, 0, 9, 5, 3, 5, 7, 3, -1, -1, -1, -1, -1, -1, -1},
    {0, 7, 8, 0, 1, 7, 1, 5, 7, -1, -1, -1, -1, -1, -1, -1},
    {1, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {9, 7, 8, 9, 5, 7, 10, 1, 2, -1, -1, -1, -1, -1, -1, -1},
    {10, 1, 2, 9, 5, 0, 5, 3, 0, 5, 7, 3, -1, -1, -1, -1},
    {8, 0, 2, 8, 2, 5, 8, 5, 7, 10, 5, 2, -1, -1, -1, -1},
    {2, 10, 5, 2, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1},
    {7, 9, 5, 7, 8, 9, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1},
    {9, 5, 7, 9, 7, 2, 9, 2, 0, 2, 7, 11, -1, -1, -1, -1},
    {2, 3, 11, 0, 1, 8, 1, 7, 8, 1, 5, 7, -1, -1, -1, -1},
    {11, 2, 1, 11, 1, 7, 7, 1, 5, -1, -1, -1, -1, -1, -1, -1},
    {9, 5, 8, 8, 5, 7, 10, 1, 3, 10, 3, 11, -1, -1, -1, -1},
    {5, 7, 0, 5, 0, 9, 7, 11, 0, 1, 0, 10, 11, 10, 0, -1},
    {11, 10, 0, 11, 0, 3, 10, 5, 0, 8, 0, 7, 5, 7, 0, -1},
    {11, 10, 5, 7, 11, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 8, 3, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {9, 0, 1, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 8, 3, 1, 9, 8, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1},
    {1, 6, 5, 2, 6, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 6, 5, 1, 2, 6, 3, 0, 8, -1, -1, -1, -1, -1, -1, -1},
    {9, 6, 5, 9, 0, 6, 0, 2, 6, -1, -1, -1, -1, -1, -1, -1},
    {5, 9, 8, 5, 8, 2, 5, 2, 6, 3, 2, 8, -1, -1, -1, -1},
    {2, 3, 11, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {11, 0, 8, 11, 2, 0, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1},
    {0, 1, 9, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1},
    {5, 10, 6, 1, 9, 2, 9, 11, 2, 9, 8, 11, -1, -1, -1, -1},
    {6, 3, 11, 6, 5, 3, 5, 1, 3, -1, -1, -1, -1, -1, -1, -1},
    {0, 8, 11, 0, 11, 5, 0, 5, 1, 5, 11, 6, -1, -1, -1, -1},
    {3, 11, 6, 0, 3, 6, 0, 6, 5, 0, 5, 9, -1, -1, -1, -1},
    {6, 5, 9, 6, 9, 11, 11, 9, 8, -1, -1, -1, -1, -1, -1, -1},
    {5, 10, 6, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {4, 3, 0, 4, 7, 3, 6, 5, 10, -1, -1, -1, -1, -1, -1, -1},
    {1, 9, 0, 5, 10, 6, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1},
    {10, 6, 5, 1, 9, 7, 1, 7, 3, 7, 9, 4, -1, -1, -1, -1},
    {6, 1, 2, 6, 5, 1, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1},
    {1, 2, 5, 5, 2, 6, 3, 0, 4, 3, 4, 7, -1, -1, -1, -1},
    {8, 4, 7, 9, 0, 5, 0, 6, 5, 0, 2, 6, -1, -1, -1, -1},
    {7, 3, 9, 7, 9, 4, 3, 2, 9, 5, 9, 6, 2, 6, 9, -1},
    {3, 11, 2, 7, 8, 4, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1},
    {5, 10, 6, 4, 7, 2, 4, 2, 0, 2, 7, 11, -1, -1, -1, -1},
    {0, 1, 9, 4, 7, 8, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1},
    {9, 2, 1, 9, 11, 2, 9, 4, 11, 7, 11, 4, 5, 10, 6, -1},
    {8, 4, 7, 3, 11, 5, 3, 5, 1, 5, 11, 6, -1, -1, -1, -1},
    {5, 1, 11, 5, 11, 6, 1, 0, 11, 7, 11, 4, 0, 4, 11, -1},
    {0, 5, 9, 0, 6, 5, 0, 3, 6, 11, 6, 3, 8, 4, 7, -1},
    {6, 5, 9, 6, 9, 11, 4, 7, 9, 7, 11, 9, -1, -1, -1, -1},
    {10, 4, 9, 6, 4, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {4, 10, 6, 4, 9, 10, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1},
    {10, 0, 1, 10, 6, 0, 6, 4, 0, -1, -1, -1, -1, -1, -1, -1},
    {8, 3, 1, 8, 1, 6, 8, 6, 4, 6, 1, 10, -1, -1, -1, -1},
    {1, 4, 9, 1, 2, 4, 2, 6, 4, -1, -1, -1, -1, -1, -1, -1},
    {3, 0, 8, 1, 2, 9, 2, 4, 9, 2, 6, 4, -1, -1, -1, -1},
    {0, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {8, 3, 2, 8, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1},
    {10, 4, 9, 10, 6, 4, 11, 2, 3, -1, -1, -1, -1, -1, -1, -1},
    {0, 8, 2, 2, 8, 11, 4, 9, 10, 4, 10, 6, -1, -1, -1, -1},
    {3, 11, 2, 0, 1, 6, 0, 6, 4, 6, 1, 10, -1, -1, -1, -1},
    {6, 4, 1, 6, 1, 10, 4, 8, 1, 2, 1, 11, 8, 11, 1, -1},
    {9, 6, 4, 9, 3, 6, 9, 1, 3, 11, 6, 3, -1, -1, -1, -1},
    {8, 11, 1, 8, 1, 0, 11, 6, 1, 9, 1, 4, 6, 4, 1, -1},
    {3, 11, 6, 3, 6, 0, 0, 6, 4, -1, -1, -1, -1, -1, -1, -1},
    {6, 4, 8, 11, 6, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {7, 10, 6, 7, 8, 10, 8, 9, 10, -1, -1, -1, -1, -1, -1, -1},
    {0, 7, 3, 0, 10, 7, 0, 9, 10, 6, 7, 10, -1, -1, -1, -1},
    {10, 6, 7, 1, 10, 7, 1, 7, 8, 1, 8, 0, -1, -1, -1, -1},
    {10, 6, 7, 10, 7, 1, 1, 7, 3, -1, -1, -1, -1, -1, -1, -1},
    {1, 2, 6, 1, 6, 8, 1, 8, 9, 8, 6, 7, -1, -1, -1, -1},
    {2, 6, 9, 2, 9, 1, 6, 7, 9, 0, 9, 3, 7, 3, 9, -1},
    {7, 8, 0, 7, 0, 6, 6, 0, 2, -1, -1, -1, -1, -1, -1, -1},
    {7, 3, 2, 6, 7, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {2, 3, 11, 10, 6, 8, 10, 8, 9, 8, 6, 7, -1, -1, -1, -1},
    {2, 0, 7, 2, 7, 11, 0, 9, 7, 6, 7, 10, 9, 10, 7, -1},
    {1, 8, 0, 1, 7, 8, 1, 10, 7, 6, 7, 10, 2, 3, 11, -1},
    {11, 2, 1, 11, 1, 7, 10, 6, 1, 6, 7, 1, -1, -1, -1, -1},
    {8, 9, 6, 8, 6, 7, 9, 1, 6, 11, 6, 3, 1, 3, 6, -1},
    {0, 9, 1, 11, 6, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {7, 8, 0, 7, 0, 6, 3, 11, 0, 11, 6, 0, -1, -1, -1, -1},
    {7, 11, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {3, 0, 8, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 1, 9, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {8, 1, 9, 8, 3, 1, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1},
    {10, 1, 2, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 2, 10, 3, 0, 8, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1},
    {2, 9, 0, 2, 10, 9, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1},
    {6, 11, 7, 2, 10, 3, 10, 8, 3, 10, 9, 8, -1, -1, -1, -1},
    {7, 2, 3, 6, 2, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {7, 0, 8, 7, 6, 0, 6, 2, 0, -1, -1, -1, -1, -1, -1, -1},
    {2, 7, 6, 2, 3, 7, 0, 1, 9, -1, -1, -1, -1, -1, -1, -1},
    {1, 6, 2, 1, 8, 6, 1, 9, 8, 8, 7, 6, -1, -1, -1, -1},
    {10, 7, 6, 10, 1, 7, 1, 3, 7, -1, -1, -1, -1, -1, -1, -1},
    {10, 7, 6, 1, 7, 10, 1, 8, 7, 1, 0, 8, -1, -1, -1, -1},
    {0, 3, 7, 0, 7, 10, 0, 10, 9, 6, 10, 7, -1, -1, -1, -1},
    {7, 6, 10, 7, 10, 8, 8, 10, 9, -1, -1, -1, -1, -1, -1, -1},
    {6, 8, 4, 11, 8, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {3, 6, 11, 3, 0, 6, 0, 4, 6, -1, -1, -1, -1, -1, -1, -1},
    {8, 6, 11, 8, 4, 6, 9, 0, 1, -1, -1, -1, -1, -1, -1, -1},
    {9, 4, 6, 9, 6, 3, 9, 3, 1, 11, 3, 6, -1, -1, -1, -1},
    {6, 8, 4, 6, 11, 8, 2, 10, 1, -1, -1, -1, -1, -1, -1, -1},
    {1, 2, 10, 3, 0, 11, 0, 6, 11, 0, 4, 6, -1, -1, -1, -1},
    {4, 11, 8, 4, 6, 11, 0, 2, 9, 2, 10, 9, -1, -1, -1, -1},
    {10, 9, 3, 10, 3, 2, 9, 4, 3, 11, 3, 6, 4, 6, 3, -1},
    {8, 2, 3, 8, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1},
    {0, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 9, 0, 2, 3, 4, 2, 4, 6, 4, 3, 8, -1, -1, -1, -1},
    {1, 9, 4, 1, 4, 2, 2, 4, 6, -1, -1, -1, -1, -1, -1, -1},
    {8, 1, 3, 8, 6, 1, 8, 4, 6, 6, 10, 1, -1, -1, -1, -1},
    {10, 1, 0, 10, 0, 6, 6, 0, 4, -1, -1, -1, -1, -1, -1, -1},
    {4, 6, 3, 4, 3, 8, 6, 10, 3, 0, 3, 9, 10, 9, 3, -1},
    {10, 9, 4, 6, 10, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {4, 9, 5, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 8, 3, 4, 9, 5, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1},
    {5, 0, 1, 5, 4, 0, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1},
    {11, 7, 6, 8, 3, 4, 3, 5, 4, 3, 1, 5, -1, -1, -1, -1},
    {9, 5, 4, 10, 1, 2, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1},
    {6, 11, 7, 1, 2, 10, 0, 8, 3, 4, 9, 5, -1, -1, -1, -1},
    {7, 6, 11, 5, 4, 10, 4, 2, 10, 4, 0, 2, -1, -1, -1, -1},
    {3, 4, 8, 3, 5, 4, 3, 2, 5, 10, 5, 2, 11, 7, 6, -1},
    {7, 2, 3, 7, 6, 2, 5, 4, 9, -1, -1, -1, -1, -1, -1, -1},
    {9, 5, 4, 0, 8, 6, 0, 6, 2, 6, 8, 7, -1, -1, -1, -1},
    {3, 6, 2, 3, 7, 6, 1, 5, 0, 5, 4, 0, -1, -1, -1, -1},
    {6, 2, 8, 6, 8, 7, 2, 1, 8, 4, 8, 5, 1, 5, 8, -1},
    {9, 5, 4, 10, 1, 6, 1, 7, 6, 1, 3, 7, -1, -1, -1, -1},
    {1, 6, 10, 1, 7, 6, 1, 0, 7, 8, 7, 0, 9, 5, 4, -1},
    {4, 0, 10, 4, 10, 5, 0, 3, 10, 6, 10, 7, 3, 7, 10, -1},
    {7, 6, 10, 7, 10, 8, 5, 4, 10, 4, 8, 10, -1, -1, -1, -1},
    {6, 9, 5, 6, 11, 9, 11, 8, 9, -1, -1, -1, -1, -1, -1, -1},
    {3, 6, 11, 0, 6, 3, 0, 5, 6, 0, 9, 5, -1, -1, -1, -1},
    {0, 11, 8, 0, 5, 11, 0, 1, 5, 5, 6, 11, -1, -1, -1, -1},
    {6, 11, 3, 6, 3, 5, 5, 3, 1, -1, -1, -1, -1, -1, -1, -1},
    {1, 2, 10, 9, 5, 11, 9, 11, 8, 11, 5, 6, -1, -1, -1, -1},
    {0, 11, 3, 0, 6, 11, 0, 9, 6, 5, 6, 9, 1, 2, 10, -1},
    {11, 8, 5, 11, 5, 6, 8, 0, 5, 10, 5, 2, 0, 2, 5, -1},
    {6, 11, 3, 6, 3, 5, 2, 10, 3, 10, 5, 3, -1, -1, -1, -1},
    {5, 8, 9, 5, 2, 8, 5, 6, 2, 3, 8, 2, -1, -1, -1, -1},
    {9, 5, 6, 9, 6, 0, 0, 6, 2, -1, -1, -1, -1, -1, -1, -1},
    {1, 5, 8, 1, 8, 0, 5, 6, 8, 3, 8, 2, 6, 2, 8, -1},
    {1, 5, 6, 2, 1, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 3, 6, 1, 6, 10, 3, 8, 6, 5, 6, 9, 8, 9, 6, -1},
    {10, 1, 0, 10, 0, 6, 9, 5, 0, 5, 6, 0, -1, -1, -1, -1},
    {0, 3, 8, 5, 6, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {10, 5, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {11, 5, 10, 7, 5, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {11, 5, 10, 11, 7, 5, 8, 3, 0, -1, -1, -1, -1, -1, -1, -1},
    {5, 11, 7, 5, 10, 11, 1, 9, 0, -1, -1, -1, -1, -1, -1, -1},
    {10, 7, 5, 10, 11, 7, 9, 8, 1, 8, 3, 1, -1, -1, -1, -1},
    {11, 1, 2, 11, 7, 1, 7, 5, 1, -1, -1, -1, -1, -1, -1, -1},
    {0, 8, 3, 1, 2, 7, 1, 7, 5, 7, 2, 11, -1, -1, -1, -1},
    {9, 7, 5, 9, 2, 7, 9, 0, 2, 2, 11, 7, -1, -1, -1, -1},
    {7, 5, 2, 7, 2, 11, 5, 9, 2, 3, 2, 8, 9, 8, 2, -1},
    {2, 5, 10, 2, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1},
    {8, 2, 0, 8, 5, 2, 8, 7, 5, 10, 2, 5, -1, -1, -1, -1},
    {9, 0, 1, 5, 10, 3, 5, 3, 7, 3, 10, 2, -1, -1, -1, -1},
    {9, 8, 2, 9, 2, 1, 8, 7, 2, 10, 2, 5, 7, 5, 2, -1},
    {1, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 8, 7, 0, 7, 1, 1, 7, 5, -1, -1, -1, -1, -1, -1, -1},
    {9, 0, 3, 9, 3, 5, 5, 3, 7, -1, -1, -1, -1, -1, -1, -1},
    {9, 8, 7, 5, 9, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {5, 8, 4, 5, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1},
    {5, 0, 4, 5, 11, 0, 5, 10, 11, 11, 3, 0, -1, -1, -1, -1},
    {0, 1, 9, 8, 4, 10, 8, 10, 11, 10, 4, 5, -1, -1, -1, -1},
    {10, 11, 4, 10, 4, 5, 11, 3, 4, 9, 4, 1, 3, 1, 4, -1},
    {2, 5, 1, 2, 8, 5, 2, 11, 8, 4, 5, 8, -1, -1, -1, -1},
    {0, 4, 11, 0, 11, 3, 4, 5, 11, 2, 11, 1, 5, 1, 11, -1},
    {0, 2, 5, 0, 5, 9, 2, 11, 5, 4, 5, 8, 11, 8, 5, -1},
    {9, 4, 5, 2, 11, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {2, 5, 10, 3, 5, 2, 3, 4, 5, 3, 8, 4, -1, -1, -1, -1},
    {5, 10, 2, 5, 2, 4, 4, 2, 0, -1, -1, -1, -1, -1, -1, -1},
    {3, 10, 2, 3, 5, 10, 3, 8, 5, 4, 5, 8, 0, 1, 9, -1},
    {5, 10, 2, 5, 2, 4, 1, 9, 2, 9, 4, 2, -1, -1, -1, -1},
    {8, 4, 5, 8, 5, 3, 3, 5, 1, -1, -1, -1, -1, -1, -1, -1},
    {0, 4, 5, 1, 0, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {8, 4, 5, 8, 5, 3, 9, 0, 5, 0, 3, 5, -1, -1, -1, -1},
    {9, 4, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {4, 11, 7, 4, 9, 11, 9, 10, 11, -1, -1, -1, -1, -1, -1, -1},
    {0, 8, 3, 4, 9, 7, 9, 11, 7, 9, 10, 11, -1, -1, -1, -1},
    {1, 10, 11, 1, 11, 4, 1, 4, 0, 7, 4, 11, -1, -1, -1, -1},
    {3, 1, 4, 3, 4, 8, 1, 10, 4, 7, 4, 11, 10, 11, 4, -1},
    {4, 11, 7, 9, 11, 4, 9, 2, 11, 9, 1, 2, -1, -1, -1, -1},
    {9, 7, 4, 9, 11, 7, 9, 1, 11, 2, 11, 1, 0, 8, 3, -1},
    {11, 7, 4, 11, 4, 2, 2, 4, 0, -1, -1, -1, -1, -1, -1, -1},
    {11, 7, 4, 11, 4, 2, 8, 3, 4, 3, 2, 4, -1, -1, -1, -1},
    {2, 9, 10, 2, 7, 9, 2, 3, 7, 7, 4, 9, -1, -1, -1, -1},
    {9, 10, 7, 9, 7, 4, 10, 2, 7, 8, 7, 0, 2, 0, 7, -1},
    {3, 7, 10, 3, 10, 2, 7, 4, 10, 1, 10, 0, 4, 0, 10, -1},
    {1, 10, 2, 8, 7, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {4, 9, 1, 4, 1, 7, 7, 1, 3, -1, -1, -1, -1, -1, -1, -1},
    {4, 9, 1, 4, 1, 7, 0, 8, 1, 8, 7, 1, -1, -1, -1, -1},
    {4, 0, 3, 7, 4, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {4, 8, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {9, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {3, 0, 9, 3, 9, 11, 11, 9, 10, -1, -1, -1, -1, -1, -1, -1},
    {0, 1, 10, 0, 10, 8, 8, 10, 11, -1, -1, -1, -1, -1, -1, -1},
    {3, 1, 10, 11, 3, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 2, 11, 1, 11, 9, 9, 11, 8, -1, -1, -1, -1, -1, -1, -1},
    {3, 0, 9, 3, 9, 11, 1, 2, 9, 2, 11, 9, -1, -1, -1, -1},
    {0, 2, 11, 8, 0, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {3, 2, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {2, 3, 8, 2, 8, 10, 10, 8, 9, -1, -1, -1, -1, -1, -1, -1},
    {9, 10, 2, 0, 9, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {2, 3, 8, 2, 8, 10, 0, 1, 8, 1, 10, 8, -1, -1, -1, -1},
    {1, 10, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 3, 8, 9, 1, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 9, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 3, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}
};
const tgt::svec3 VERTEX_OFFSETS[] = {
    tgt::svec3(0, 0, 0),
    tgt::svec3(1, 0, 0),
    tgt::svec3(1, 1, 0),
    tgt::svec3(0, 1, 0),
    tgt::svec3(0, 0, 1),
    tgt::svec3(1, 0, 1),
    tgt::svec3(1, 1, 1),
    tgt::svec3(0, 1, 1),
};
tgt::vec3 normalAt(const voreen::VolumeRAM* v, const tgt::vec3& p, long& numNormalErrors) {
    float val_xp = v->getVoxelNormalizedLinear(p + tgt::vec3( 1, 0, 0));
    float val_xm = v->getVoxelNormalizedLinear(p + tgt::vec3(-1, 0, 0));
    float val_yp = v->getVoxelNormalizedLinear(p + tgt::vec3( 0, 1, 0));
    float val_ym = v->getVoxelNormalizedLinear(p + tgt::vec3( 0,-1, 0));
    float val_zp = v->getVoxelNormalizedLinear(p + tgt::vec3( 0, 0, 1));
    float val_zm = v->getVoxelNormalizedLinear(p + tgt::vec3( 0, 0,-1));

    // Central differences seem to provide better results in edge cases
    float dv_dx = val_xm - val_xp;
    float dv_dy = val_ym - val_yp;
    float dv_dz = val_zm - val_zp;
    tgt::vec3 n = tgt::vec3(dv_dx, dv_dy, dv_dz);
    if(n == tgt::vec3::zero) {
        // Central differences failed. Try forward differences
        float central_val = v->getVoxelNormalizedLinear(p); //Should be about equal to isovalue, but we better get the exact value
        n.x = central_val - val_xp;
        n.y = central_val - val_yp;
        n.z = central_val - val_zp;
        if(n == tgt::vec3::zero) {
            // Forward differences failed, too. Try backward differences
            n.x = val_xm - central_val;
            n.y = val_ym - central_val;
            n.z = val_zm - central_val;
            if(n == tgt::vec3::zero) {
                // Well, none of the above worked. Just set the normal to (1,0,0) and report an error.
                numNormalErrors++;
                n = tgt::vec3(1,0,0);
            }
        }
    }
    return n;
}
void prepareVertex(const voreen::VolumeRAM* v, float isoValue, const tgt::svec3& p1, const tgt::svec3& p2, tgt::vec3& pos /*out*/, tgt::vec3& normal /*out*/, long& numNormalErrors /*out*/) {
    float a = (isoValue - v->getVoxelNormalized(p2)) / (v->getVoxelNormalized(p1) - v->getVoxelNormalized(p2));
    //float a = 0.5f;
    pos = a*tgt::vec3(p1) + (1-a)*tgt::vec3(p2);
    normal = normalAt(v, pos, numNormalErrors);
    normal = tgt::normalize(normal); //normal cannot be zero
}
void processCell(const voreen::VolumeRAM* v, const tgt::svec3& p, float isoValue, voreen::GlMeshGeometryUInt32Normal* mesh, long& numNormalErrors /*out*/) {
    uint8_t index = 0;
    if(v->getVoxelNormalized(p + VERTEX_OFFSETS[0]) < isoValue) index |= (1 << 0);
    if(v->getVoxelNormalized(p + VERTEX_OFFSETS[1]) < isoValue) index |= (1 << 1);
    if(v->getVoxelNormalized(p + VERTEX_OFFSETS[2]) < isoValue) index |= (1 << 2);
    if(v->getVoxelNormalized(p + VERTEX_OFFSETS[3]) < isoValue) index |= (1 << 3);
    if(v->getVoxelNormalized(p + VERTEX_OFFSETS[4]) < isoValue) index |= (1 << 4);
    if(v->getVoxelNormalized(p + VERTEX_OFFSETS[5]) < isoValue) index |= (1 << 5);
    if(v->getVoxelNormalized(p + VERTEX_OFFSETS[6]) < isoValue) index |= (1 << 6);
    if(v->getVoxelNormalized(p + VERTEX_OFFSETS[7]) < isoValue) index |= (1 << 7);
    if(index == 0 || index == 0xff) {
        return;
    }
    tgt::vec3 positions[12];
    tgt::vec3 normals[12];
    const uint16_t edgeIndex = EDGE_INDICES[index];
    if(edgeIndex & (1 << 0x0)) prepareVertex(v, isoValue, p + VERTEX_OFFSETS[0], p + VERTEX_OFFSETS[1], positions[0x0], normals[0x0], numNormalErrors);
    if(edgeIndex & (1 << 0x1)) prepareVertex(v, isoValue, p + VERTEX_OFFSETS[1], p + VERTEX_OFFSETS[2], positions[0x1], normals[0x1], numNormalErrors);
    if(edgeIndex & (1 << 0x2)) prepareVertex(v, isoValue, p + VERTEX_OFFSETS[2], p + VERTEX_OFFSETS[3], positions[0x2], normals[0x2], numNormalErrors);
    if(edgeIndex & (1 << 0x3)) prepareVertex(v, isoValue, p + VERTEX_OFFSETS[3], p + VERTEX_OFFSETS[0], positions[0x3], normals[0x3], numNormalErrors);
    if(edgeIndex & (1 << 0x4)) prepareVertex(v, isoValue, p + VERTEX_OFFSETS[4], p + VERTEX_OFFSETS[5], positions[0x4], normals[0x4], numNormalErrors);
    if(edgeIndex & (1 << 0x5)) prepareVertex(v, isoValue, p + VERTEX_OFFSETS[5], p + VERTEX_OFFSETS[6], positions[0x5], normals[0x5], numNormalErrors);
    if(edgeIndex & (1 << 0x6)) prepareVertex(v, isoValue, p + VERTEX_OFFSETS[6], p + VERTEX_OFFSETS[7], positions[0x6], normals[0x6], numNormalErrors);
    if(edgeIndex & (1 << 0x7)) prepareVertex(v, isoValue, p + VERTEX_OFFSETS[7], p + VERTEX_OFFSETS[4], positions[0x7], normals[0x7], numNormalErrors);
    if(edgeIndex & (1 << 0x8)) prepareVertex(v, isoValue, p + VERTEX_OFFSETS[0], p + VERTEX_OFFSETS[4], positions[0x8], normals[0x8], numNormalErrors);
    if(edgeIndex & (1 << 0x9)) prepareVertex(v, isoValue, p + VERTEX_OFFSETS[1], p + VERTEX_OFFSETS[5], positions[0x9], normals[0x9], numNormalErrors);
    if(edgeIndex & (1 << 0xA)) prepareVertex(v, isoValue, p + VERTEX_OFFSETS[2], p + VERTEX_OFFSETS[6], positions[0xA], normals[0xA], numNormalErrors);
    if(edgeIndex & (1 << 0xB)) prepareVertex(v, isoValue, p + VERTEX_OFFSETS[3], p + VERTEX_OFFSETS[7], positions[0xB], normals[0xB], numNormalErrors);

    //const int8_t* vertices = &TRIANGLE_VERTEX_INDICES[index][0];
    const int32_t* vertices = &TRIANGLE_VERTEX_INDICES[index][0];
    for(size_t vertexNum = 0; vertexNum < 15 && vertices[vertexNum] != -1; vertexNum += 3) {
        int8_t v0 = vertices[vertexNum    ];
        mesh->addVertex(positions[v0], normals[v0]);
        int8_t v1 = vertices[vertexNum + 1];
        mesh->addVertex(positions[v1], normals[v1]);
        int8_t v2 = vertices[vertexNum + 2];
        mesh->addVertex(positions[v2], normals[v2]);
    }
}

// Smoothing

struct SimpleVol {
    float* data_;
    tgt::ivec3 dim_;
    SimpleVol(const tgt::ivec3& dim)
        : dim_(dim)
    {
        data_ = new float[tgt::hmul(dim)];
    }
    ~SimpleVol() {
        delete[] data_;
    }
    inline float& data(int x, int y, int z) {
        return data_[x+dim_.x*(y+dim_.y*z)];
    }
    inline float sample(int x, int y, int z) {
        return data(std::max(std::min(x, (int)dim_.x-1), 0), std::max(std::min(y, (int)dim_.y-1), 0), std::max(std::min(z, (int)dim_.z-1), 0));
    }
};


// VertexLayout for Mesh generation
struct VertexLayout {
    tgt::vec3 pos, normal;
};

} //anonymous namespace

namespace voreen {

static const GLuint EDGE_INDICES_BINDING_POINT = 0;
static const GLuint TRIANGLE_VERTEX_INDICES_BINDING_POINT = 1;

const std::string IsosurfaceExtractor::loggerCat_("voreen.surface.isosurfaceextractor");

IsosurfaceExtractor::IsosurfaceExtractor()
    : VolumeProcessor()
    , inport_(Port::INPORT, "isosurfaceextractor.inport", "Volume Input")
    , outport_(Port::OUTPORT, "isosurfaceextractor.outport", "Geometry Output")
    , enabledProp_("enabledProp","Enabled",true)
    , isoValueProp_("isoValueProp", "Isovalue", 0.0f, std::numeric_limits<float>::lowest(), std::numeric_limits<float>::max())
    , meshBufferSizeProp_("meshBufferSizeProp", "MeshBuffer Size (Bytes)", 0, 0, std::numeric_limits<int>::max(), Processor::INVALID_RESULT, IntProperty::STATIC, Property::LOD_ADVANCED)
    , surfaceTypeProp_("surfaceTypeProp", "Surface Type")
    , marchingCubesProcessingModeProp_("MarchingCubesProcessingMode", "Processing Mode")
    , marchingCubesShaderProp_("marchingCubesShaderProp", "Marching Cubes Shader", ShaderFileList(tgt::ShaderObject::VERTEX_SHADER, "marchingcubes.vert").add(tgt::ShaderObject::GEOMETRY_SHADER, "marchingcubes.geom"), Processor::INVALID_PROGRAM, Property::LOD_DEBUG)
    , blocksShaderProp_("blocksShaderProp", "Blocks Shader", ShaderFileList(tgt::ShaderObject::VERTEX_SHADER, "blocks.vert").add(tgt::ShaderObject::GEOMETRY_SHADER, "blocks.geom"), Processor::INVALID_PROGRAM, Property::LOD_DEBUG)
    , useSmoothingProp_("useSmoothingProp","Apply Smoothing", false, Processor::INVALID_RESULT, Property::LOD_DEFAULT)
    , maxGradientProp_("maxGradientProp", "Max Gradient Magnitude", -3.0f, -10.0f, 0.0f, Processor::INVALID_RESULT, FloatProperty::STATIC, Property::LOD_DEBUG)
    , maxSmoothingIterationsProp_("maxSmoothingIterationsProp", "Max Iterations", 10, 1, 1000, Processor::INVALID_RESULT, IntProperty::STATIC, Property::LOD_ADVANCED)
    , smoothingStepSizeProp_("smoothingStepSizeProp", "Step Size", 0.0f, std::numeric_limits<float>::lowest(), std::numeric_limits<float>::max(), Processor::INVALID_RESULT, FloatProperty::STATIC, Property::LOD_DEBUG)
    , recalculateSmoothingButtonProp_("recalculateSmoothingButtonProp", "Recalculate", Processor::INVALID_RESULT, Property::LOD_DEFAULT)
    , smoothingProgressProp_("smoothingProgressProp", "Progress")
    , meshBuffer_(0)
    , vao_(0)
    , edgeIndicesBuffer_(0)
    , triangleVertexIndicesBuffer_(0)
    , smoothedVolume_(nullptr)
    , smoothingForced_(false)
{
    addPort(inport_);
    addPort(outport_);

        addProperty(enabledProp_);

        addProperty(surfaceTypeProp_);
            surfaceTypeProp_.setGroupID("extraction");
            surfaceTypeProp_.addOption("marchingcubes", "Marching Cubes", MARCHING_CUBES);
            if(tgt::isInitedGL()) {
                surfaceTypeProp_.addOption("blocks", "Blocks", BLOCKS);
            }
            ON_CHANGE(surfaceTypeProp_, IsosurfaceExtractor, adjustExtractionProperties);

        addProperty(marchingCubesProcessingModeProp_);
            marchingCubesProcessingModeProp_.setGroupID("extraction");
            marchingCubesProcessingModeProp_.addOption("cpu", "CPU", MODE_CPU);
            if(tgt::isInitedGL()) {
                marchingCubesProcessingModeProp_.addOption("gpu", "GPU", MODE_GPU);
                marchingCubesProcessingModeProp_.selectByValue(MODE_GPU); //Default to GPU mode (if available)
            }
            ON_CHANGE(marchingCubesProcessingModeProp_, IsosurfaceExtractor, adjustExtractionProperties);

        addProperty(isoValueProp_);
            isoValueProp_.setGroupID("extraction");
            ON_CHANGE_LAMBDA(isoValueProp_, [this] () {
                    resetSmoothedVolume();
                    useSmoothingProp_.set(false);
                    });


        if(tgt::isInitedGL()) {
            addProperty(meshBufferSizeProp_);
                meshBufferSizeProp_.setGroupID("extraction");
                meshBufferSizeProp_.setStepping(sizeof(VertexLayout) * 3);

            addProperty(marchingCubesShaderProp_);
                marchingCubesShaderProp_.setGroupID("extraction");

            addProperty(blocksShaderProp_);
                blocksShaderProp_.setGroupID("extraction");
        }
    setPropertyGroupGuiName("extraction", "Extraction");
    adjustExtractionProperties();

        addProperty(useSmoothingProp_);
            useSmoothingProp_.setGroupID("smoothing");
            ON_CHANGE(useSmoothingProp_, IsosurfaceExtractor, adjustSmoothingProperties);
        addProperty(maxGradientProp_);
            maxGradientProp_.setGroupID("smoothing");
        addProperty(maxSmoothingIterationsProp_);
            maxSmoothingIterationsProp_.setGroupID("smoothing");
        addProperty(smoothingStepSizeProp_);
            smoothingStepSizeProp_.setGroupID("smoothing");
        addProperty(recalculateSmoothingButtonProp_);
            recalculateSmoothingButtonProp_.setGroupID("smoothing");
            ON_CHANGE_LAMBDA(recalculateSmoothingButtonProp_, [this] () {
                        smoothingForced_ = true;
                    });
        addProperty(smoothingProgressProp_);
            smoothingProgressProp_.setGroupID("smoothing");
    setPropertyGroupGuiName("smoothing", "Binary Volume Smoothing");
    adjustSmoothingProperties();
}

IsosurfaceExtractor::~IsosurfaceExtractor() {
}
void IsosurfaceExtractor::initialize() {
    VolumeProcessor::initialize();

    if(tgt::isInitedGL()) {
        glGenVertexArrays(1, &vao_);

        glGenBuffers(1, &edgeIndicesBuffer_);
        glBindBuffer(GL_ARRAY_BUFFER, edgeIndicesBuffer_);
        glBufferData(GL_ARRAY_BUFFER, sizeof(EDGE_INDICES), EDGE_INDICES, GL_STATIC_DRAW);
        glBindBuffer(GL_ARRAY_BUFFER, 0);

        glGenBuffers(1, &triangleVertexIndicesBuffer_);
        glBindBuffer(GL_ARRAY_BUFFER, triangleVertexIndicesBuffer_);
        glBufferData(GL_ARRAY_BUFFER, sizeof(TRIANGLE_VERTEX_INDICES), TRIANGLE_VERTEX_INDICES, GL_STATIC_DRAW);
        glBindBuffer(GL_ARRAY_BUFFER, 0);

        buildMarchingCubesShader();

        buildBlocksShader();
    }
}

void IsosurfaceExtractor::deinitialize(){
    if(tgt::isInitedGL()) {
        if(meshBuffer_) {
            glDeleteBuffers(1, &meshBuffer_);
        }
        if(triangleVertexIndicesBuffer_) {
            glDeleteBuffers(1, &triangleVertexIndicesBuffer_);
        }
        if(edgeIndicesBuffer_) {
            glDeleteBuffers(1, &edgeIndicesBuffer_);
        }
        if(vao_) {
            glDeleteVertexArrays(1, &vao_);
        }
    }

    resetSmoothedVolume();

    VolumeProcessor::deinitialize();
}

void IsosurfaceExtractor::resetSmoothedVolume() {
    smoothedVolume_ = nullptr;
    smoothingProgressProp_.setProgress(0.0f);
}

void IsosurfaceExtractor::adjustExtractionProperties() {
    SurfaceType t = surfaceTypeProp_.getValue();
    MarchingCubesProcessingMode p = marchingCubesProcessingModeProp_.getValue();

    meshBufferSizeProp_.setVisibleFlag(t == BLOCKS || p == MODE_GPU);

    marchingCubesProcessingModeProp_.setVisibleFlag(t == MARCHING_CUBES);
    marchingCubesShaderProp_.setVisibleFlag(t == MARCHING_CUBES);

    blocksShaderProp_.setVisibleFlag(t == BLOCKS);
}

void IsosurfaceExtractor::adjustSmoothingProperties() {
    bool smoothingEnabled = useSmoothingProp_.get();
    binarizationThresholdProp_.setVisibleFlag(smoothingEnabled);
    maxGradientProp_.setVisibleFlag(smoothingEnabled);
    maxSmoothingIterationsProp_.setVisibleFlag(smoothingEnabled);
    smoothingStepSizeProp_.setVisibleFlag(smoothingEnabled);
    recalculateSmoothingButtonProp_.setVisibleFlag(smoothingEnabled);
    smoothingProgressProp_.setVisibleFlag(smoothingEnabled);
}

void IsosurfaceExtractor::process() {
    if(!enabledProp_.get()) {
        outport_.setData(nullptr);
        return;
    }
    const VolumeBase* invol = inport_.getData();
    if(!invol) {
        return;
    }
    if(inport_.hasChanged()) {
        const VolumeMinMax* mm = invol->getDerivedData<VolumeMinMax>();
        tgtAssert(mm, "No VolumeMinMax");

        bool minAndMaxWereEqual = isoValueProp_.getMinValue() == isoValueProp_.getMaxValue();

        // Find proper values for isoValueProp_
        isoValueProp_.setMinValue(mm->getMin());
        isoValueProp_.setMaxValue(mm->getMax());
        isoValueProp_.adaptDecimalsToRange(2);

        // Restore and interesting surface if we e.g. set a zero surface by accident previously
        if(minAndMaxWereEqual) {
            isoValueProp_.set((isoValueProp_.getMinValue() + isoValueProp_.getMaxValue())*0.5f);
        }


        /*
        float targetStepping = range/100;
        int digits = std::max(0, static_cast<int>(-log10(targetStepping)));
        int shiftValue = std::pow(10, digits);
        isoValueProp_.setMinValue(std::floor(mm->getMin()*shiftValue)/shiftValue);
        isoValueProp_.setMaxValue(std::ceil(mm->getMax()*shiftValue)/shiftValue);
        // Set iso value to value mean of min and max. Theres probably an interesting surface (at least more interesting than min or max).
        // However, this messes with restoring the value when loading a workspace
        // isoValueProp_.set((mm->getMax()+mm->getMin())/2);
        isoValueProp_.setNumDecimals(digits);
        isoValueProp_.setStepping(pow(10.0f, -digits));
        */

        // Find a proper value for meshBufferSizeProp_
        size_t maxBufferSize =
              sizeof(VertexLayout) /* size per Vertex */
            * 3 /* We generate Triangles */
            * 5 /* Up to 5 triangles per voxel */
            * tgt::hmul(invol->getDimensions() - tgt::svec3(1,1,1)) /* Number of cuboid spaces between voxels */;
        // TODO find a better way to handle size limitations of int property?
        maxBufferSize = std::min(maxBufferSize, (size_t)std::numeric_limits<int>::max());
        meshBufferSizeProp_.setMaxValue(maxBufferSize);

        // Do not _reduce_ buffersize on new inport data as this can get very annoying when repeatedly changing
        // a rather noisy input geometry (for which the buffer size is too small).
        size_t minInitialBufferSize = maxBufferSize/10;
        if(meshBufferSizeProp_.get() < minInitialBufferSize) {
            meshBufferSizeProp_.set(minInitialBufferSize);
        }

        resetSmoothedVolume();

        // As an indicator, that the volume still has to be smoothed.
        smoothingProgressProp_.setProgress(0.0f);
    }

    if(smoothingForced_ || useSmoothingProp_.get() && !smoothedVolume_) {
        smoothingForced_ = false;
        // Generate a smoothed volume
        float normalizedBinarizationThreshold = invol->getRealWorldMapping().realWorldToNormalized(isoValueProp_.get());
        VolumeRAM* smoothed = smooth(invol->getRepresentation<voreen::VolumeRAM>(), normalizedBinarizationThreshold);
        if(!smoothed) {
            LERROR("Could not smooth volume.");
            outport_.setData(nullptr);
            return;
        }
        smoothedVolume_.reset(new Volume(smoothed, invol));
    }

    const VolumeBase* currentVolume = useSmoothingProp_.get() && smoothedVolume_ ? smoothedVolume_.get() : invol;
    const float normalizedIsoValue = useSmoothingProp_.get() && smoothedVolume_ ? 0.5f : invol->getRealWorldMapping().realWorldToNormalized(isoValueProp_.get());;

    Geometry* surface = nullptr;
    switch(surfaceTypeProp_.getValue()) {
        case MARCHING_CUBES:
            switch(marchingCubesProcessingModeProp_.getValue()) {
                case MODE_CPU:
                    surface = extractSurfaceMCCPU(currentVolume, normalizedIsoValue);
                    break;
                case MODE_GPU:
                    // We extract surfaces in the spaces between the voxels.
                    // Thus we actually only consider vol->getDimensions() - (1,1,1) cuboid spaces.
                    surface = extractSurfaceGPU(currentVolume, marchingCubesShaderProp_.getShader(), currentVolume->getDimensions() - tgt::svec3(1), normalizedIsoValue);
                    break;
                default:
                    tgtAssert(false, "Invalid extraction mode.");
            }
            break;

        case BLOCKS:
            // We extract surfaces around every voxel
            // Thus we only consider all vol->getDimensions() voxels cuboid spaces.
            surface = extractSurfaceGPU(currentVolume, blocksShaderProp_.getShader(), currentVolume->getDimensions(), normalizedIsoValue);
            break;
        default:
            tgtAssert(false, "Invalid surface type.");
    }
    if(surface) {
        surface->setTransformationMatrix(invol->getVoxelToWorldMatrix());
    } else {
        LWARNING("Could not generate isosurface.");
    }
    outport_.setData(surface);
}
VoreenSerializableObject* IsosurfaceExtractor::create() const {
    return new IsosurfaceExtractor();
}
void IsosurfaceExtractor::buildMarchingCubesShader() {
    marchingCubesShaderProp_.rebuild();
    tgt::Shader* shader = marchingCubesShaderProp_.getShader();
    const GLchar* feedbackVaryings[] = { "pos", "normal" };
    glTransformFeedbackVaryings(shader->getID(), sizeof(feedbackVaryings)/sizeof(*feedbackVaryings), feedbackVaryings, GL_INTERLEAVED_ATTRIBS);

    shader->linkProgram();

    const GLuint edgeIndicesBlockIndex = glGetUniformBlockIndex(shader->getID(), "EDGE_INDICES_BLOCK");
    tgtAssert(edgeIndicesBlockIndex != GL_INVALID_INDEX, "invalid block index.");
    LGL_ERROR;
    glUniformBlockBinding(shader->getID(), edgeIndicesBlockIndex, EDGE_INDICES_BINDING_POINT);
    LGL_ERROR;
    const GLuint triangleVertexIndicesBlockIndex = glGetUniformBlockIndex(shader->getID(), "TRIANGLE_VERTEX_INDICES_BLOCK");
    LGL_ERROR;
    tgtAssert(triangleVertexIndicesBlockIndex != GL_INVALID_INDEX, "invalid block index.");
    glUniformBlockBinding(shader->getID(), triangleVertexIndicesBlockIndex, TRIANGLE_VERTEX_INDICES_BINDING_POINT);
    LGL_ERROR;
}

void IsosurfaceExtractor::buildBlocksShader() {
    blocksShaderProp_.rebuild();
    tgt::Shader* shader = blocksShaderProp_.getShader();
    const GLchar* feedbackVaryings[] = { "pos", "normal" };
    glTransformFeedbackVaryings(shader->getID(), sizeof(feedbackVaryings)/sizeof(*feedbackVaryings), feedbackVaryings, GL_INTERLEAVED_ATTRIBS);

    shader->linkProgram();
}

Geometry* IsosurfaceExtractor::extractSurfaceMCCPU(const voreen::VolumeBase* vol, float normalizedIsoValue) const {
    tgtAssert(vol, "No volume");
    const voreen::VolumeRAM* volram = vol->getRepresentation<voreen::VolumeRAM>();
    voreen::GlMeshGeometryUInt32Normal* mesh = new voreen::GlMeshGeometryUInt32Normal();
    long numNormalErrors = 0;
    for(size_t z = 0; z<vol->getDimensions().z-1; ++z) {
        for(size_t y = 0; y<vol->getDimensions().y-1; ++y) {
            for(size_t x = 0; x<vol->getDimensions().x-1; ++x) {
                processCell(volram, tgt::svec3(x,y,z), normalizedIsoValue, mesh, numNormalErrors);
            }
        }
    }
    if(numNormalErrors) {
        LWARNINGC("voreen.surfaceextraction.surfaceectractor", std::to_string(numNormalErrors) + " normals were set to (1,0,0) because they were (0,0,0)");
    }
    return mesh;
}
Geometry* IsosurfaceExtractor::extractSurfaceGPU(const VolumeBase* vol, tgt::Shader* shader, tgt::svec3 gridSize, float normalizedIsoValue) {
    tgtAssert(tgt::isInitedGL(), "OpenGL not initiated");
    tgtAssert(vol, "No volume");

    if(meshBuffer_) {
        glDeleteBuffers(1, &meshBuffer_);
    }

    if(!shader) {
        LERROR("Could not get surface extraction shader.");
        return nullptr;
    }

    const VolumeGL* volumeGL = vol->getRepresentation<VolumeGL>();
    if(!volumeGL) {
        LERROR("Could not get VolumeGL.");
        return nullptr;
    }
    const tgt::Texture* volumeTexture = volumeGL->getTexture();
    tgtAssert(volumeTexture, "No texture");

    glGenBuffers(1, &meshBuffer_);
    glBindBuffer(GL_ARRAY_BUFFER, meshBuffer_);
    glBufferData(GL_ARRAY_BUFFER, meshBufferSizeProp_.get(), nullptr, GL_STATIC_READ);
    LGL_ERROR;

    GLuint query;
    glGenQueries(1, &query);
    LGL_ERROR;

    shader->activate();
    shader->setUniform("isoValue_", normalizedIsoValue);
    shader->setUniform("volumeDimensions_", tgt::ivec3(vol->getDimensions()));
    LGL_ERROR;

    glBindBufferBase(GL_TRANSFORM_FEEDBACK_BUFFER, 0, meshBuffer_);
    LGL_ERROR;

    glEnable(GL_RASTERIZER_DISCARD);
    glBeginQuery(GL_TRANSFORM_FEEDBACK_PRIMITIVES_WRITTEN, query);
    glBeginTransformFeedback(GL_TRIANGLES);
    LGL_ERROR;

    glActiveTexture(GL_TEXTURE0);
    volumeTexture->bind();
    shader->setUniform("colorTex_", 0);
    volumeTexture->enable();

    LGL_ERROR;

    glBindBufferBase(GL_UNIFORM_BUFFER, EDGE_INDICES_BINDING_POINT, edgeIndicesBuffer_);
    glBindBufferBase(GL_UNIFORM_BUFFER, TRIANGLE_VERTEX_INDICES_BINDING_POINT, triangleVertexIndicesBuffer_);
    LGL_ERROR;

    glBindVertexArray(vao_);
    glDrawArrays(GL_POINTS, 0, tgt::hmul(gridSize));
    glBindVertexArray(0);
    LGL_ERROR;

    glFinish();
    LGL_ERROR;
    glBindBufferBase(GL_UNIFORM_BUFFER, EDGE_INDICES_BINDING_POINT, 0);
    glBindBufferBase(GL_UNIFORM_BUFFER, TRIANGLE_VERTEX_INDICES_BINDING_POINT, 0);
    LGL_ERROR;

    volumeTexture->disable();
    LGL_ERROR;

    glEndTransformFeedback();
    glEndQuery(GL_TRANSFORM_FEEDBACK_PRIMITIVES_WRITTEN);
    glDisable(GL_RASTERIZER_DISCARD);
    LGL_ERROR;

    shader->deactivate();

    GLuint primitivesWritten;
    glGetQueryObjectuiv(query, GL_QUERY_RESULT, &primitivesWritten);
    LGL_ERROR;

    std::vector<VertexLayout> vertices(primitivesWritten * 3 /* triangles! */);
    glGetBufferSubData(GL_TRANSFORM_FEEDBACK_BUFFER, 0, vertices.size()*sizeof(VertexLayout), vertices.data());
    LGL_ERROR;
    glBindBufferBase(GL_TRANSFORM_FEEDBACK_BUFFER, 0, 0);

    voreen::GlMeshGeometryUInt32Normal* mesh = new voreen::GlMeshGeometryUInt32Normal();
    for(const auto& vertex : vertices) {
        mesh->addVertex(vertex.pos, vertex.normal);
    }

    const float maxPrimitives = meshBufferSizeProp_.get()/(
              sizeof(VertexLayout) /* size per Vertex */
            * 3 /* We generate Triangles */
         );
    const float bufferUtilization = primitivesWritten/std::max(1.0f, maxPrimitives);
    if(bufferUtilization >= 0.99f || maxPrimitives == 0) {
        LWARNING("MeshBuffer utilization near or above 100%.");
    }

    return mesh;
}

VolumeRAM* IsosurfaceExtractor::smooth(const VolumeRAM* input, float normalizedBinarizationThreshold) {
    const float maxGradient = std::pow(10, maxGradientProp_.get());
    const int maxIt = maxSmoothingIterationsProp_.get();
    const float startStepsize = smoothingStepSizeProp_.get();
    const tgt::ivec3 dim(input->getDimensions());

    smoothingProgressProp_.reset();
    smoothingProgressProp_.setProgressRange(tgt::vec2(0, 1));

    SimpleVol volume(dim);
    SimpleVol mask(dim);
    SimpleVol gradient(dim);

    // Initialize volume:
    const float initVal = 1.0f/std::sqrt(tgt::hmul(dim));
    for(int z=0; z<dim.z; ++z) {
        for(int y=0; y<dim.y; ++y) {
            for(int x=0; x<dim.x; ++x) {
                float val = input->getVoxelNormalized(x,y,z);
                if(val < normalizedBinarizationThreshold) {
                    volume.data(x,y,z) = -initVal;
                    mask.data(x,y,z) = -initVal;
                } else {
                    volume.data(x,y,z) = initVal;
                    mask.data(x,y,z) = initVal;
                }
            }
        }
    }

    float gradientLength = std::numeric_limits<float>::max();
    float roughness = 1;
    float previousRoughness = 1;
    float stepsize = startStepsize;
    for(int it = 0; gradientLength > maxGradient*maxGradient && it < maxIt; ++it) {
        if(roughness > previousRoughness) {
            stepsize *= 0.5;
        }
        previousRoughness = roughness;
        roughness = 0;
        gradientLength = 0;
#ifdef VRN_MODULE_OPENMP
        #pragma omp parallel for
#endif
        for(int z=0; z<dim.z; ++z) {
            for(int y=0; y<dim.y; ++y) {
                for(int x=0; x<dim.x; ++x) {
                    float gradHere = 0;
                    gradHere += volume.sample(x, y, z) * 6 * 3;

                    gradHere += volume.sample(x-2, y, z);
                    gradHere += volume.sample(x-1, y, z) * -4;
                    gradHere += volume.sample(x+1, y, z) * -4;
                    gradHere += volume.sample(x+2, y, z);

                    gradHere += volume.sample(x, y-2, z);
                    gradHere += volume.sample(x, y-1, z) * -4;
                    gradHere += volume.sample(x, y+1, z) * -4;
                    gradHere += volume.sample(x, y+2, z);

                    gradHere += volume.sample(x, y, z-2);
                    gradHere += volume.sample(x, y, z-1) * -4;
                    gradHere += volume.sample(x, y, z+1) * -4;
                    gradHere += volume.sample(x, y, z+2);


                    gradient.data(x,y,z) = gradHere * stepsize;
                    gradientLength += gradHere * gradHere;

                    float term = 0;
                    term += volume.sample(x, y, z) * -2;
                    term += volume.sample(x-1, y, z);
                    term += volume.sample(x+1, y, z);
                    roughness += term*term;

                    term += volume.sample(x, y, z) * -2;
                    term += volume.sample(x, y-1, z);
                    term += volume.sample(x, y+1, z);
                    roughness += term*term;

                    term += volume.sample(x, y, z) * -2;
                    term += volume.sample(x, y, z-1);
                    term += volume.sample(x, y, z+1);
                    roughness += term*term;
                }
            }
        }

        float volumeLengthSq = 0;
#ifdef VRN_MODULE_OPENMP
        #pragma omp parallel for
#endif
        for(int z=0; z<dim.z; ++z) {
            for(int y=0; y<dim.y; ++y) {
                for(int x=0; x<dim.x; ++x) {
                    float newVal = volume.data(x,y,z) - gradient.data(x,y,z);
                    // Constraints
                    if(mask.data(x,y,z) > 0) {
                        newVal = std::max(newVal, 0.f);
                    } else {
                        newVal = std::min(newVal, 0.f);
                    }
                    volume.data(x,y,z) = newVal;
                    volumeLengthSq += newVal * newVal;
                }
            }
        }
        float volumeLength = std::sqrt(volumeLengthSq);

        if(tgt::isNaN(volumeLength)) {
            return nullptr;
        }

#ifdef VRN_MODULE_OPENMP
        #pragma omp parallel for
#endif
        for(int z=0; z<dim.z; ++z) {
            for(int y=0; y<dim.y; ++y) {
                for(int x=0; x<dim.x; ++x) {
                    volume.data(x,y,z) /= volumeLength;
                }
            }
        }

        LDEBUG("Smooting it: " + std::to_string(it) + " | roughness: " + std::to_string(roughness) + "| gradLengthsq: " + std::to_string(gradientLength) +  " | stepsize: " + std::to_string(stepsize));
        smoothingProgressProp_.setProgress(static_cast<float>(it+1)/maxIt);
    }
    smoothingProgressProp_.setProgress(1.0f);

    float min = std::numeric_limits<float>::max();
    float max = std::numeric_limits<float>::min();
    for(int z=0; z<dim.z; ++z) {
        for(int y=0; y<dim.y; ++y) {
            for(int x=0; x<dim.x; ++x) {
                min = std::min(min, volume.data(x,y,z));
                max = std::max(max, volume.data(x,y,z));
            }
        }
    }

    VolumeRAM* output = VolumeFactory().create(input->getBaseType(), dim);
#ifdef VRN_MODULE_OPENMP
    #pragma omp parallel for
#endif
    for(int z=0; z<dim.z; ++z) {
        for(int y=0; y<dim.y; ++y) {
            for(int x=0; x<dim.x; ++x) {
                float val = volume.data(x,y,z);
                if(mask.data(x,y,z) > 0) {
                    val = 0.5+val/(2*max);
                    tgtAssert(val >= 0.5 && val <= 1, "somethings wrong");
                } else {
                    val = 0.5-val/(2*min);
                    tgtAssert(val >= 0 && val <= 0.5, "somethings wrong");
                }
                output->setVoxelNormalized(val, tgt::svec3(x,y,z));
            }
        }
    }

    return output;
}
} // namespace voreen
