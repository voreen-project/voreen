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

#include "curve.h"

#include "Eigen/Core"
#include "tgt/logmanager.h"

namespace voreen {

Curve::Curve(std::vector<tgt::vec3>&& inputPoints) : 
    inputPoints_(std::move(inputPoints)),
    totalLength_(0.f)
{
}

std::vector<tgt::vec3> Curve::resamplePoints(const std::vector<tgt::vec3>& segment) {
    std::vector<tgt::vec3> resampled;
    if (segment.empty())
        return resampled;

    // Append the first element separately.
    resampled.push_back(segment.front());

    // In case we have more than two samples, we achieve the result by linear interpolation.
    if(segment.size() > 2) {

        std::vector<float> distances(segment.size());

        float totalLength = 0.0f;
        distances[0] = 0.0f;
        for(size_t i = 1; i < segment.size(); i++) {
            totalLength += tgt::distance(segment[i-1], segment[i]);
            distances[i] = totalLength;
        }

        const float segmentLength = totalLength / (segment.size() - 1);
        for (size_t i = 1, pos = 0; i < segment.size() - 1; i++) {

            while (distances[pos+1] < segmentLength * i) pos++;

            const float t = tgt::clamp((segmentLength * i - distances[pos]) / (distances[pos + 1] - distances[pos]), 0.0f, 1.0f);

            tgt::vec3 position = segment[pos] * (1.0f - t) + segment[pos + 1] * t;

            resampled.push_back(position);
        }
    }

    // Append the last element separately.
    resampled.push_back(segment.back());

    return resampled;
}

CubicCurve::CubicCurve(std::vector<tgt::vec3>&& points) :
    Curve(std::move(points))
{
    if (inputPoints_.size() == 0) {
        meta_.emplace_back(Meta{ 0.f, tgt::vec3::zero, tgt::vec3::zero, tgt::vec3::zero, tgt::vec3::zero });
        return;
    }
    if (inputPoints_.size() == 1) {
        meta_.emplace_back(Meta{ 0.f, tgt::vec3::zero, tgt::vec3::zero, tgt::vec3::zero, inputPoints_[0] });
        return;
    }

    tgt::vec3 m = tgt::normalize(inputPoints_[0] - inputPoints_[inputPoints_.size() - 1]);
    for (auto i = 1u; i < inputPoints_.size(); ++i) {
        auto p0 = inputPoints_[i - 1];
        auto p1 = inputPoints_[i];

        // Generate cubic interpolation between p0 and p1
        // p(t)=at³+bt²+ct+d
        // p'(t)=3at²+2bt+c
        // p(0)=p0 => d = p0
        // p'(0)=m => c = m
        // p(1)=p1 => a+b+c+d = p1 => a+b+m+p0 = p1
        // least-squares-method => a = b => a = (p1 - m - p0) / 2

        auto a = (p1 - m - p0) * 0.5f;
        // TODO: better approximation for the length of the segment
        auto l = tgt::length(p1 - p0);
        totalLength_ += l;

        meta_.emplace_back(Meta{ l, a, a, m, p0 });

        // calculate next m = p'(1) [same as p'(0) of next polynomial]
        m = tgt::normalize(a * 5.f + m);
    }
}

tgt::vec3 CubicCurve::evaluate(size_t index, float t) const
{
    const auto& meta = meta_.at(index);
    auto tt = t * t;
    auto ttt = tt * t;
    return meta.a * ttt + meta.b * tt + meta.c * t + meta.d;
}

tgt::vec3 CubicCurve::sample(float t) const
{
    auto requestedLength = t;
    auto cachedLength = 0.f;
    for (auto i = 0u; i < meta_.size(); ++i) {
        if (cachedLength + meta_.at(i).length >= requestedLength)
        {
            // modify, so that the sampled points are approximately equidistant
            auto t_ = (requestedLength - cachedLength) / meta_.at(i).length;
            return evaluate(i, t_);
        }
        cachedLength += meta_.at(i).length;
    }
    return evaluate(meta_.size() - 1, 1.f);
}

}
