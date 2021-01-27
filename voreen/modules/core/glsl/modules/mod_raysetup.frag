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

/**
 * This module contains all functions which are used for the raysetup
 * as well as the ray traversal within a raycaster.
 */

/**
 * This parameter defines the minimum opacity saturation
 * a ray has to reach before it is terminated.
 * Setting this value to 1.0 disables early ray termination.
 */
const float EARLY_RAY_TERMINATION_OPACITY = 0.95;

/***
 * Calculates the direction of the ray and returns the number
 * of steps and the direction.
 ***/
void raySetup(in vec3 first, in vec3 last, in float samplingStepSize, out vec3 rayDirection, out float tIncr, out float tEnd) {
    // calculate ray parameters
    tIncr = samplingStepSize;

    rayDirection = last - first;
    tEnd = length(rayDirection);
    rayDirection = normalize(rayDirection);
}

/***
 * Applies early ray termination. The current opacity is compared to
 * the maximum opacity. In case it is greater, the opacity is set to
 * 1.0 and true is returned, otherwise false is returned.
 ***/
bool earlyRayTermination(inout float opacity, in float maxOpacity) {
    if (opacity >= maxOpacity) {
        opacity = 1.0;
        return true;
    } else {
        return false;
    }
}

bool earlyRayTermination(inout float opacity1, inout float opacity2, inout float opacity3, in float maxOpacity) {
    if (min(opacity1, min(opacity2, opacity3)) >= maxOpacity) {
        opacity1 = 1.0;
        opacity2 = 1.0;
        opacity3 = 1.0;
        return true;
    } else {
        return false;
    }
}

// Does the hardware support a shader program length that allows us to use a single loop or do
// we need two nested loops?
#if VRN_MAX_PROGRAM_LOOP_COUNT < 256*256
    // Use two nested loops, should be supported everywhere
    #define WHILE(keepGoing) for (int loop0=0; keepGoing && loop0<255; loop0++) { for (int loop1=0; keepGoing && loop1<255; loop1++) {

    #define END_WHILE } }
#else
    #define WHILE(keepGoing) for (int loop=0; keepGoing && loop<VRN_MAX_PROGRAM_LOOP_COUNT; loop++) {

    #define END_WHILE }
#endif

