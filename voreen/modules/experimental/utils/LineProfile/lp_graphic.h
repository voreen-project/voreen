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

#ifndef VRN_LP_GRAPHIC_H
#define VRN_LP_GRAPHIC_H

//vectors
#include "tgt/vector.h"

namespace voreen {
/**
 * Utility class beeing used in lineprofile.h
 * This class encapsulate the render functions used in lineprofile::render.
 */
class LP_Graphic {
public:
    /**
     * Function for rendering a sphere.
     *
     * @param position  position of the sphere (worldcoordinates)
     * @param size      size of the sphere (radius)
     * @param color     color of the sphere
     */
    static void drawSphere(const tgt::vec3 position, float size, const tgt::vec4 color);

    /**
     * Function for rendering a line (arrow).
     * This function draws an interval defined by [start,end] of an arrow with length "length".
     *
     * @param start     beginning of the arrow interval
     * @param end       end of the arrow interval
     * @param length    length of the entire arrow (used for calculating tail and arrowhead)
     * @param size      width of the arrow (radius)
     * @param color     color of the arrow
     */
    static void drawLineInterval(float start, float end, float length, float size, const tgt::vec4 color);
};

}   //namespace

#endif // VRN_LP_GRAPHIC_H
