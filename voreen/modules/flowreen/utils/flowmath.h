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

#ifndef FLOWMATH_H
#define FLOWMATH_H

#include "tgt/vector.h"
#include "tgt/matrix.h"

#include <deque>
#include <vector>
#include <time.h>

namespace voreen {

class Flow2D;
class Flow3D;

class FlowMath {
public:
    template<class Flow, class Vector>
    static std::vector<Vector> computePathline(const std::vector<Flow>& flows, const Vector& r0,
        const float deltaT = 0.5f, float* const lineLength = 0);

    template<class Flow, class Vector>
    static std::deque<Vector> computeStreamlineEuler(const Flow& flow,
        const Vector& r0, const float length = 150.0f, const float stepwidth = 0.5f,
        int* const startIndex = 0, const tgt::vec2& thresholds = tgt::vec2(0.0f));

    template<class Flow, class Vector>
    static std::deque<Vector> computeStreamlineRungeKutta(const Flow& flow,
        const Vector& r0, const float length = 150.0f, const float stepwidth = 0.5f,
        int* const startIndex = 0, const tgt::vec2& thresholds = tgt::vec2(0.0f));

    static tgt::mat4 getTransformationMatrix(const std::vector<tgt::vec3>& streamline,
        const size_t& index, const float scaling = 1.0f);

    template<typename T>
    static std::vector<T> dequeToVector(std::deque<T>& deq);

    static tgt::ivec2 getNearestIntegral(const tgt::vec2& p);
    static tgt::ivec3 getNearestIntegral(const tgt::vec3& p);

    /**
     * Normalizes the given vector if it is not the null-vector. If it
     * is the null-vector, it will not be normalized.
     * tgt::normalize() does not do this.
     */
    template<class Vector>
    static Vector normalize(const Vector& v);

    static float uniformRandom() {
        return (rand() / static_cast<float>(RAND_MAX));
    }

    static tgt::vec2 uniformRandomVec2() {
        return tgt::vec2(uniformRandom(), uniformRandom());
    }

    static tgt::vec3 uniformRandomVec3() {
        return tgt::vec3(uniformRandom(), uniformRandom(), uniformRandom());
    }

    static void resetSeed(unsigned int seed = 0) {
#ifdef VRN_NO_RANDOM
        return;
#else
        if (seed == 0)
            srand(static_cast<unsigned int>(time(0)));
        else
            srand(seed);
#endif
    }

private:
    // prevent instantiation and copying of objects from this class
    //
    FlowMath();
    FlowMath(const FlowMath&);
    FlowMath& operator=(const FlowMath&);

    template<class Flow, class Vector>
    static Vector lintTime(const std::vector<Flow>& flows, const Vector& r, const float time);
};

}   // namespace

#endif
