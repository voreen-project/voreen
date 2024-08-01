/**********************************************************************
 *                                                                    *
 * tgt - Tiny Graphics Toolbox                                        *
 *                                                                    *
 * Copyright (C) 2005-2024 University of Muenster, Germany,           *
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

#ifndef TGT_SPLINE_H_
#define TGT_SPLINE_H_

#include <vector>

#include "tgt/curve.h"
#include "tgt/bspline.h"
#include "tgt/types.h"
#include "tgt/vector.h"

namespace tgt {

/**
    A Spline-Class, that implements Curve-Class. It has a user-defined number of controlpoints that build the
    B-Spline curve.  This curve can, for example, be used to run a camera or a light source along itself.
*/
class TGT_API Spline : public Curve {
public:
    //Spline() : Curve (100, Bounds(), true, true) {}
    /// The constructor, by default, sets an empty controlpoint-vector and a tau-value of .5
    Spline(const std::vector<vec3>& points = std::vector<vec3>(), float t = .5f, int stepCount = 100);

    /// This function yields the interpolated position between the current- and the next-position,
    /// depending on the interpolation-variable t
    vec3 interpolateCurrent(float t, const vec3& previous,
        const vec3& current, const vec3& next, const vec3& afternext) const;

    /// This function yields the interpolated first derivation between the derivation at the current-
    /// and that at the next-position, depending on the interpolation-variable t
    vec3 firstDeriveCurrent(float t, const vec3& previous,
        const vec3& current, const vec3& next, const vec3& afternext) const;

    /// This can be used to simply interpolate between all control points that are currently in the
    /// control point-vector
    vec3 interpolate(float t, float total = 1.f) const;

    /// This can be used to simply interpolate between all first derivations at the control points
    vec3 firstDerive(float t, float total = 1.f) const;

    /// setter / getter
    void setTau(float t);
    float getTau() const;

    std::vector<vec3>& getControlPoints();
    const std::vector<vec3>& getControlPoints() const;

    /// Use this to add controlpoints to the spline.
    void addControlPoint(const vec3& p);


    /* These methods implement the virtuals of the Curve-Class */

    //! Returns the point on the spline at the parameter value t.
    vec3 getPoint(GLfloat t);

    //! Returns the first derivation (tangent) of the spline at the parameter value t.
    vec3 getDerivative(GLfloat t);

private:
    /// This std::vector contains the controlpoints
    std::vector<vec3> controlPoints_;

    /// tau is a measure for how tightly or loosely the curve follows the controlpoints.
    float tau_;
};


} // namespace tgt

#endif //#ifndef _SPLINE_H_

