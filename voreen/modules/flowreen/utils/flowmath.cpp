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

#ifndef FLOWMATH_H
#include "modules/flowreen/datastructures/deprecated/flow3d.h"
#include "modules/flowreen/utils/flowmath.h"
//#include "modules/flowreen/include/streamlinetexture.h"

using tgt::vec3;

namespace voreen {

// public static methods
//

template<class Flow, class Vector>
std::vector<Vector> FlowMath::computePathline(const std::vector<Flow>& flows, const Vector& r0,
    const float deltaT, float* const lineLength)
{
    std::vector<Vector> pathline;

    if (flows.empty() == true)
        return pathline;

    float length = 0.0f;

    Vector r(r0);

    // At first interpolate between the starting position r0 and the position r1 = r0 + v0
    // where the first flow vector v0 is located.
    //
    const size_t aux = static_cast<size_t>(1.0f / deltaT);
    const Vector v1 = flows[0]->lookupFlow(r0);
    for (size_t i = 1; i <= aux; ++i) {
        Vector delta_r = (v1 * (i * deltaT)) * deltaT;
        r += delta_r;
        length += tgt::length(delta_r);
        pathline.push_back(r);
    }

    // Now interpolate between the flow vectors
    //
    const size_t numSteps = static_cast<size_t>((flows.size() - 1) / deltaT);
    for (size_t i = 1; i <= numSteps; ++i) {
        Vector v = lintTime(flows, r, (i * deltaT)) * deltaT;
        r += v;
        length += tgt::length(v);
        pathline.push_back(r);
    }

    if (lineLength != 0)
        *lineLength = length;
    return pathline;
}

template std::vector<tgt::vec3> FlowMath::computePathline(const std::vector<const Flow3D*>& flows,
                                                          const tgt::vec3& r0, const float deltaT,
                                                          float* const lineLength);

// ----------------------------------------------------------------------------

template<class Flow, class Vector>
std::deque<Vector> FlowMath::computeStreamlineEuler(const Flow& flow,
                                           const Vector& r0,
                                           const float length,
                                           const float stepwidth,
                                           int* const startIndex,
                                           const tgt::vec2& thresholds)
{
    // N: number of partitions
    //
    const float h = fabsf(stepwidth);
    float fn = ceilf(fabsf(length) / h);
    const unsigned int N = (length != 0.0f) ? static_cast<unsigned int>(fn) : 0;

    Vector r(r0), r_(r0);

    std::deque<Vector> points;
    points.push_back(r0);
    int indexR0 = 0;

    bool lookupPos = true;  // integrate along the streamline in positive direction?
    bool lookupNeg = true;  // integrate along the streamline in negative direction?
    bool useThresholds = (thresholds != tgt::vec2::zero);

    for (unsigned int i = 0; ((N == 0) || (i < N)) ; ++i) {
        if (lookupPos == true) {
            const Vector& v = flow.lookupFlow(r);
            if (useThresholds == true) {
                float magnitude = tgt::length(v);
                if ((magnitude < thresholds.x) || (magnitude > thresholds.y)) {
                    lookupPos = false;
                    break;
                }
            }

            if (v != Vector::zero) {
                r += normalize(v) * h;

                lookupPos = flow.isInsideBoundings(r);
                if (r == points.back()) // in case of no progress on streamline in this direction...
                    lookupPos = false;
                else if (lookupPos == true)
                    points.push_back(r);
            } else
                lookupPos = false;
        }

        if (lookupNeg == true) {
            const Vector& v = flow.lookupFlow(r_);
            if (useThresholds == true) {
                float magnitude = tgt::length(v);
                if ((magnitude < thresholds.x) || (magnitude > thresholds.y)) {
                    lookupNeg = false;
                    break;
                }
            }

            if (v != Vector::zero) {
                r_ -= normalize(v) * h;

                lookupNeg = flow.isInsideBoundings(r_);
                if (r_ == points.front()) // in case of no progress on streamline in this direction...
                    lookupNeg = false;
                else if (lookupNeg == true) {
                    points.push_front(r_);
                    ++indexR0;
                }
            } else
                lookupNeg = false;
        }

        if ((lookupPos == false) && (lookupNeg == false))
            break;
    }   // for ( ; ; ++i)

    if (startIndex != 0)
        *startIndex = indexR0;

    return points;
}

template std::deque<tgt::vec2> FlowMath::computeStreamlineEuler(
    const Flow2D& flow, const tgt::vec2& r, const float length,
    const float stepwidth, int* const startIndex, const tgt::vec2& thresholds);

template std::deque<tgt::vec3> FlowMath::computeStreamlineEuler(
    const Flow3D& flow, const tgt::vec3& r, const float length,
    const float stepwidth, int* const startIndex, const tgt::vec2& thresholds);

// ----------------------------------------------------------------------------

template<class Flow, class Vector>
std::deque<Vector> FlowMath::computeStreamlineRungeKutta(const Flow& flow,
                                                const Vector& r0,
                                                const float length,
                                                const float stepwidth,
                                                int* const startIndex,
                                                const tgt::vec2& thresholds)
{
    // N: number of partitions
    //
    const float h = fabsf(stepwidth);
    float fn = ceilf(fabsf(length) / h);
    const unsigned int N = (length != 0.0f) ? static_cast<unsigned int>(fn) : 0;

    Vector r(r0), r_(r0);
    Vector k1(0.0f), k2(0.0f), k3(0.0f), k4(0.0f);
    Vector k1_(0.0f), k2_(0.0f), k3_(0.0f), k4_(0.0f);

    std::deque<Vector> points;    // points on streamline in positive direction
    points.push_back(r0);       // avoid that one of the deque is empty
    int indexR0 = 0;            // position of r0 within deque

    bool lookupPos = true;  // integrate along the streamline in positive direction?
    bool lookupNeg = true;  // integrate along the streamline in negative direction?
    bool useThresholds = (thresholds != tgt::vec2::zero);

    for (unsigned int i = 0; ((N == 0) || (i < N)); ++i) {
        if (lookupPos == true) {
            const Vector& v = flow.lookupFlow(r);
            if (useThresholds == true) {
                float magnitude = tgt::length(v);
                if ((magnitude < thresholds.x) || (magnitude > thresholds.y)) {
                    lookupPos = false;
                    break;
                }
            }

            if (v != Vector::zero) {
                k1 = normalize(v) * h;
                k2 = normalize( flow.lookupFlow(r + (k1 / 2.0f)) ) * h;
                k3 = normalize( flow.lookupFlow(r + (k2 / 2.0f)) ) * h;
                k4 = normalize( flow.lookupFlow(r + k3) ) * h;

                r += ((k1 / 6.0f) + (k2 / 3.0f) + (k3 / 3.0f) + (k4 / 6.0f));

                lookupPos = flow.isInsideBoundings(r);
                if (r == points.back()) // in case of no progress on streamline in this direction...
                    lookupPos = false;
                else if (lookupPos == true)
                    points.push_back(r);
            } else
                lookupPos = false;
        }

        if (lookupNeg == true) {
            const Vector& v = flow.lookupFlow(r_);
            if (useThresholds == true) {
                float magnitude = tgt::length(v);
                if ((magnitude < thresholds.x) || (magnitude > thresholds.y)) {
                    lookupNeg = false;
                    break;
                }
            }

            if (v != Vector::zero) {
                k1_ = normalize(v) * h;
                k2_ = normalize( flow.lookupFlow(r_ - (k1_ / 2.0f)) ) * h;
                k3_ = normalize( flow.lookupFlow(r_ - (k2_ / 2.0f)) ) * h;
                k4_ = normalize( flow.lookupFlow(r_ - k3_) ) * h;

                r_ -= ((k1_ / 6.0f) + (k2_ / 3.0f) + (k3_ / 3.0f) + (k4_ / 6.0f));

                lookupNeg = flow.isInsideBoundings(r_);
                if (r_ == points.front()) // in case of no progress on streamline in this direction...
                    lookupNeg = false;
                else if (lookupNeg == true) {
                    points.push_front(r_);
                    ++indexR0;
                }
            } else
                lookupNeg = false;
        }

        if ((lookupPos == false) && (lookupNeg == false))
            break;
    }   // for (; ; ++i)

    if (startIndex != 0)
        *startIndex = indexR0;

    return points;
}

template std::deque<tgt::vec2> FlowMath::computeStreamlineRungeKutta(
    const Flow2D& flow, const tgt::vec2& r, const float length,
    const float stepwidth, int* const startIndex, const tgt::vec2& thresholds);

template std::deque<tgt::vec3> FlowMath::computeStreamlineRungeKutta(
    const Flow3D& flow, const tgt::vec3& r, const float length,
    const float stepwidth, int* const startIndex, const tgt::vec2& thresholds);

// ----------------------------------------------------------------------------

tgt::mat4 FlowMath::getTransformationMatrix(const std::vector<tgt::vec3>& streamline,
                                                        const size_t& index, const float scaling)
{
    if (streamline.empty())
        return tgt::mat4::identity;

    // Create a local tripod for the current element on the streamline:
    // tangent:     points along the streamline's direction
    // normal:      points to origin from the streamlines position
    // binormal:    perpendicular to tangent and normal
    //
    tgt::vec3 tangent(0.0f);

    const tgt::vec3& r = streamline[index];

    // If there is an element left "further on" the streamline take it to
    // to calculate the tangent. Otherwise use the element preceding the
    // current one.
    //
    if (index < (streamline.size() - 1) )
        tangent = normalize(streamline[index + 1] - r);
    else if (index > 0)
        tangent = normalize(r - streamline[index - 1]);
    else
        tangent = normalize(r);

    tgt::vec3 temp(0.0f, 1.0f, 0.0f);
    if (fabsf(tgt::dot(temp, tangent)) >= 0.99f)
        temp = tgt::vec3(0.0f, 0.0f, 1.0f);

    tgt::vec3 binormal(normalize(cross(temp, tangent)));
    tgt::vec3 normal(normalize(cross(tangent, binormal)));

    if (scaling != 1.0f) {
        normal *= scaling;
        binormal *= scaling;
        tangent *= scaling;
    }

    // Create a transformation matrix based on that tripod (Rigid-body or
    // direction-of-flight transformation).
    // The circle given in the plane will be rotated by multiplying
    // that matrix so that is oriented along the stramline afterwards.
    //
    return tgt::mat4(normal.x, binormal.x, tangent.x, r.x,
        normal.y, binormal.y, tangent.y, r.y,
        normal.z, binormal.z, tangent.z, r.z,
        0.0f, 0.0f, 0.0f, 1.0f);
}

// ----------------------------------------------------------------------------

template<typename T>
std::vector<T> FlowMath::dequeToVector(std::deque<T>& deq) {
    // convert the std::deque into a std::vector
    //
    std::vector<T> vec;
    typename std::deque<T>::iterator it = deq.begin();
    while (deq.empty() == false) {
        vec.push_back(*it);
        deq.erase(it);
        it = deq.begin();
    }
    return vec;
}

template std::vector<tgt::vec2> FlowMath::dequeToVector(std::deque<tgt::vec2>& deq);
template std::vector<tgt::vec3> FlowMath::dequeToVector(std::deque<tgt::vec3>& deq);

// ----------------------------------------------------------------------------

tgt::ivec2 FlowMath::getNearestIntegral(const tgt::vec2& p) {
    return tgt::ivec2(static_cast<int>(tgt::round(p.x)),
        static_cast<int>(tgt::round(p.y)));
}

tgt::ivec3 FlowMath::getNearestIntegral(const tgt::vec3& p) {
    return tgt::ivec3(static_cast<int>(tgt::round(p.x)),
        static_cast<int>(tgt::round(p.y)),
        static_cast<int>(tgt::round(p.z)));
}

// ----------------------------------------------------------------------------

template<class Vector>
Vector FlowMath::normalize(const Vector& v) {
    if (v == Vector(0.0))
        return Vector(0.0);
    float len = sqrtf(tgt::dot(v, v));
    return Vector(v / len);
}

template tgt::vec2 FlowMath::normalize(const tgt::vec2& v);
template tgt::vec3 FlowMath::normalize(const tgt::vec3& v);
template tgt::vec4 FlowMath::normalize(const tgt::vec4& v);

// ----------------------------------------------------------------------------

template<class Flow, class Vector>
Vector FlowMath::lintTime(const std::vector<Flow>& flows, const Vector& r, const float time) {
    if (flows.empty() == true)
        return Vector();

    if (time <= 0.0f)
        return flows[0]->lookupFlow(r);

    float fintegral = 0.0f;
    float fract = modff(time, &fintegral);
    size_t t1 = tgt::clamp(static_cast<int>(fintegral), 0, static_cast<int>(flows.size() - 1));
    size_t t2 = tgt::clamp(static_cast<int>(fintegral) + 1, 0, static_cast<int>(flows.size() - 1));

    Vector v1 = (flows[t1]->lookupFlow(r) * (1.0f - fract));
    Vector v2 = (flows[t2]->lookupFlow(r) * fract);
    return (v1 + v2);
}

template tgt::vec3 FlowMath::lintTime(const std::vector<const Flow3D*>& flows,
                                      const tgt::vec3& r, const float time);

// ----------------------------------------------------------------------------

}   // namespace

#endif
