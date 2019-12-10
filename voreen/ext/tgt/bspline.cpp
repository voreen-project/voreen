/**********************************************************************
 *                                                                    *
 * tgt - Tiny Graphics Toolbox                                        *
 *                                                                    *
 * Copyright (C) 2005-2019 University of Muenster, Germany,           *
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

#include "tgt/bspline.h"

namespace tgt {

BSpline::BSpline(const std::vector<vec3>& ctrlPoints, int degree, GLuint stepCount)
    : Curve(stepCount, Bounds(), true, true)
    , degree_(degree)
    , ctrlPoints_(ctrlPoints)
{
    tgtAssert(degree_ > 0, "Degree of a B-Spline must be greater zero.");
    tgtAssert(static_cast<int>(ctrlPoints_.size()) > degree_,
        "The number of control points must exceed the spline's degree.");

    generateKnotVector();

}

BSpline::BSpline(const std::vector<vec3>& ctrlPoints, const std::vector<float>& knotValues,
    int degree, GLuint stepCount)
    : Curve(stepCount, Bounds(), true, true)
    , degree_(degree)
    , ctrlPoints_(ctrlPoints)
{
    tgtAssert(degree_ > 0, "Degree of a B-Spline must be greater zero.");
    tgtAssert(static_cast<int>(ctrlPoints_.size()) > degree_, "The number of control points must exceed the spline's degree.");
    tgtAssert(knotValues.size() == ctrlPoints_.size(), "The number of knot values must match the number of control points.");
    //    tgtAssert( std::adjacent_find(knotValues.begin(), knotValues.end(), std::greater<float>()) == knotValues.end(),
    //        "Knot values must be in non-decreasing order");

    knots_ = knotValues;
    generateKnotVector();
}


vec3 BSpline::getPoint(GLfloat t) {
    // evaluate spline by Cox-deBoor recursion (also known as deCasteljau recursion)
    vec3 point(0.f);
    for (size_t i=0; i<ctrlPoints_.size(); ++i)
        point += evalBasisFunction(static_cast<int>(i), degree_, t)*ctrlPoints_.at(i);

    return point;
}

vec3 BSpline::getDerivative(GLfloat /*t*/) {
    // TODO
    return vec3(0.f);
}

// private
float BSpline::evalBasisFunction(int knotID, int deg, float u) {
    // Cox-deBoor recursion

    if (deg == 0) {
        if ( u >= knots_.at(knotID) && u < knots_.at(knotID+1) )
            return 1.f;
        else
            return 0.f;
    }
    else {
        float factor1, factor2;

        if (knots_.at(knotID+deg) == knots_.at(knotID))
            factor1 = 0.f;
        else
            factor1 = (u - knots_.at(knotID)) / (knots_.at(knotID+deg) - knots_.at(knotID));

        if (knots_.at(knotID+deg+1) == knots_.at(knotID+1))
            factor2 = 0.f;
        else
            factor2 = (knots_.at(knotID+deg+1) - u) / (knots_.at(knotID+deg+1) - knots_.at(knotID+1));

        return factor1*evalBasisFunction(knotID, deg-1, u) + factor2*evalBasisFunction(knotID+1, deg-1, u);
    }

}

// private
void BSpline::generateKnotVector() {
    if (knots_.empty()) {
        // No knot values have been passed => generate uniform distribution.
        // In order to force the spline to interpolate the start and end control points,
        // (degree_+1) knots of value 0.0 or 1.0 are inserted at the beginning and the end
        // of the knot vector, respectively.

        // total number of knot values and number of inner knots
        size_t numKnots = ctrlPoints_.size() + degree_ + 1;
        size_t innerKnots = numKnots - 2*(degree_+1);

        // start knots (for interpolation of the first control point)
        for (int i=0; i<=degree_; ++i) {
            knots_.push_back(0.f);
        }

        // inner knots
        if (innerKnots > 0) {
            float step = 1.f / static_cast<float>(innerKnots+1);
            for (size_t i=1; i<=innerKnots; ++i) {
                knots_.push_back(i*step);
            }
        }

        // end knots (for interpolation of the last control point)
        for (int i=0; i<=degree_; ++i) {
            knots_.push_back(1.f + 1e-6f);
        }

    }
    else {
        // scale and shift passed knot values to interval [0.0, 1.0] and
        // add (degree_+1)/2 knots at the beginning and end of the knot vector,
        // respectively.

        std::vector<float> tempKnots = knots_;
        knots_.clear();

        int knotsToAdd = degree_;

        int firstKnot = tgt::iceil((degree_+1)/2.f) - 1;
        int lastKnot = static_cast<int>(tempKnots.size()) - tgt::ifloor((degree_+1)/2.f);

        float minKnot = tempKnots.at(firstKnot);
        float maxKnot = tempKnots.at(lastKnot);

        float shift = -minKnot;
        float scale = (1.f + 1e-6f) / (maxKnot - minKnot);

        for (int i=0; i<knotsToAdd; ++i) {
            knots_.push_back(0.f);
        }

        for (int i=firstKnot; i<=lastKnot; ++i) {
            knots_.push_back( (tempKnots.at(i)+shift)*scale );
        }

        for (int i=0; i<knotsToAdd; ++i) {
            knots_.push_back(1.f + 1e-6f);
        }
    }
}

} // namespace
