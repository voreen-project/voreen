/**********************************************************************
 *                                                                    *
 * tgt - Tiny Graphics Toolbox                                        *
 *                                                                    *
 * Copyright (C) 2005-2018 University of Muenster, Germany,           *
 * Department of Computer Science.                                    *
 *                                                                    *
 * This file is part of the tgt library. This library is free         *
 * software; you can redistribute it and/or modify it under the terms *
 * of the GNU Lesser General Public License version 2.1 as published  *
 * by the Free Software Foundation.                                   *
 *                                                                    *
 * This library is distributed in the hope that it will be useful,    *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of     *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the       *
 * GNU Lesser General Public License for more details.                *
 *                                                                    *
 * You should have received a copy of the GNU Lesser General Public   *
 * License in the file "LICENSE.txt" along with this library.         *
 * If not, see <http://www.gnu.org/licenses/>.                        *
 *                                                                    *
 **********************************************************************/

#include "tgt/naturalcubicspline.h"

namespace tgt {

bool vec2CompareX(const vec2& i, const vec2& j) {
    return (i.x < j.x);
}

bool vec2CompareY(const vec2& i, const vec2& j) {
    return (i.y < j.y);
}

bool vec2Equal(const vec2& i, const vec2& j) {
    return (i.x == j.x);
}

bool vec3Compare(const vec3& i, const vec3& j) {
    return (i.x < j.x);
}

bool vec3Equal(const vec3& i, const vec3& j) {
    return (i.x == j.x);
}

using std::vector;

NaturalCubicSpline::NaturalCubicSpline() {}

NaturalCubicSpline::NaturalCubicSpline(const vector<vec2>& points) {
    setPoints(points);
}

float NaturalCubicSpline::getPoint(float t) const {
    float delta = points_[points_.size()-1].x - points_[0].x;
    float x = t * delta;
    int i = lookup(x);

    float a = coefficients_[i].x;
    float b = coefficients_[i].y;
    float c = coefficients_[i].z;
    float d = coefficients_[i].w;
    float xi = points_[i].x;
    float deltaX = x - xi;
    return a + deltaX*(b + deltaX*(c + deltaX*d));
}

void NaturalCubicSpline::setPoints(const vector<vec2>& points) {
    points_ = points;
    coefficients_.clear();
    coefficients_.resize(points.size() - 1);

    computeCoefficients();
}

void NaturalCubicSpline::computeCoefficients() {
    // implement method presented in "Algorithms for Cubic Spline Interpolation"
    size_t n = points_.size() - 1;
    vector<float> h(n);
    vector<float> b(n);
    for (size_t i = 0; i < n; ++i) {
        h[i] = points_[i+1].x - points_[i].x;
        b[i] = (points_[i+1].y - points_[i].y) / h[i];
    }

    vector<float> u(n);
    vector<float> v(n);
    u[1] = 2 * (h[0] + h[1]);
    v[1] = 6 * (b[1] - b[0]);
    for (size_t i = 2; i < n; ++i) {
        u[i] = 2 * (h[i-1] + h[i]) - (h[i-1] * h[i-1]) / u[i-1];
        v[i] = 6 * (b[i] - b[i-1]) - (h[i-1] * v[i-1] / u[i-1]);
    }

    vector<float> z(n + 1);
    z[n] = 0;
    for (size_t i = n - 1; i > 0; --i) {
        z[i] = (v[i] - h[i]*z[i+1]) / u[i];
    }
    z[0] = 0;

    for (size_t i = 0; i < n; ++i) {
        float aCoeff = points_[i].y;
        float bCoeff = -((h[i]/6.f) * z[i+1]) - ((h[i]/3.f)*z[i]) + ((points_[i+1].y - points_[i].y) / h[i]);
        float cCoeff = z[i]/2.f;
        float dCoeff = ((z[i+1] - z[i]) / (6.f * h[i]));
        coefficients_[i] = vec4(aCoeff, bCoeff, cCoeff, dCoeff);
    }
}

int NaturalCubicSpline::lookup(float x) const {
    const int size = static_cast<int>(points_.size());

    if (x <= points_[0].x)
        return 0;
    else if (x >= points_[size - 2].x)
        return size - 2;
    else {
        int i1 = 0;
        int i2 = size - 2;
        int i3 = 0;

        while (i2 - i1 > 1) {
            i3 = i1 + ((i2 - i1) >> 1);
            if (points_[i3].x > x)
                i2 = i3;
            else
                i1 = i3;
        }
        return i1;
    }
}

vector<vec3> generateKnotPoints(const vector<vec2>& ctrlPoints) {
    vector<vec3> result(ctrlPoints.size());
    //vector<vec2> px(ctrlPoints.size());
    result[0].x = 0.f;
    for (size_t i = 1; i < ctrlPoints.size(); ++i) {
        vec2 delta = ctrlPoints[i] - ctrlPoints[i-1];
        float val = sqrtf(std::abs(delta.x + delta.y));
        result[i].x = val;
        //px[i].x = val + px[i-1].x;
    }

    for (size_t i = 1; i < result.size(); ++i)
        result[i].x = result[i].x + result[i-1].x;

    //vector<vec2> py(px);

    for (size_t i = 0; i < ctrlPoints.size(); ++i) {
        result[i].y = ctrlPoints[i].x;
        result[i].z = ctrlPoints[i].y;
    }

    //vector<vector<vec2> > result(2);
    //result[0] = px;
    //result[1] = py;
    return result;
}

vector<vec4> generateKnotPoints(const vector<vec3>& ctrlPoints) {
    vector<vec4> result(ctrlPoints.size());
    result[0].x = 0.f;
    for (size_t i = 1; i < ctrlPoints.size(); ++i) {
        vec3 delta = ctrlPoints[i] - ctrlPoints[i-1];
        //float val = sqrtf(std::abs(delta.x + delta.y + delta.z));
        float val = sqrtf(std::abs(delta.x) + std::abs(delta.y) + std::abs(delta.z));
        result[i].x = val;
        //px[i].x = val + px[i-1].x;
    }

    for (size_t i = 1; i < result.size(); ++i)
        result[i].x = result[i].x + result[i-1].x;

    //vector<vec2> py(px);
    //vector<vec2> pz(px);

    for (size_t i = 0; i < ctrlPoints.size(); ++i) {
        result[i].y = ctrlPoints[i].x;
        result[i].z = ctrlPoints[i].y;
        result[i].w = ctrlPoints[i].z;
    }

    //vector<vector<vec2> > result(3);
    //result[0] = px;
    //result[1] = py;
    //result[2] = pz;
    return result;
}

} // namespace
