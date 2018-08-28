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

#ifndef TGT_BSPLINE_H
#define TGT_BSPLINE_H

#include <vector>
#include "tgt/curve.h"
#include "tgt/vector.h"

namespace tgt {

/**
 * Uniform B-Spline of arbitrary degree, which is evaluated by Cox-deBoor recursion.
 *
 * In order to force the spline to interpolate the first and last control point,
 * degree+1 knot values are repeated at the beginning and the end of
 * the knot vector, respectively.
 *
 */
class TGT_API BSpline : public Curve {
public:
    /**
     * Constructor expecting the spline's control points. The inner knot values are
     * distributed uniformly and range from 0.0 till 1.0.
     *
     * @param ctrlPoints The spline's control points. Note: The number of control points must
     *          exceed the spline's degree.
     * @param degree The spline's degree, i.e. the degree of its basis functions. By default,
     *          a cubic B-Spline is generated. Note: The degree must be greater zero.
     * @param stepCount Only relevant, if the spline is to be rendered: It specifies the number
     *          of line segments the spline is divided into.
     */
    BSpline(const std::vector<vec3>& ctrlPoints, int degree = 3, GLuint stepCount = 50);

    /**
     * Constructor expecting the spline's control points and the corresponding knot values.
     * The passed knot-vector is rescaled uniformly and shifted to match the interval [0.0, 1.0].
     *
     * @param ctrlPoints The spline's control points. Note: The number of control points must
     *          exceed the spline's degree.
     * @param knotValues The spline's knot values. Note: The number of knot values must match
     *          the number of control points and must be passed in non-decreasing order.
     * @param degree The spline's degree, i.e. the degree of its basis functions. By default,
     *          a cubic B-Spline is generated. Note: The degree must be greater zero.
     * @param stepCount Only relevant, if the spline is to be rendered: It specifies the number
     *          of line segments the spline is divided into.
     */
    BSpline(const std::vector<vec3>& ctrlPoints, const std::vector<float>& knotValues, int degree = 3, GLuint stepCount = 50);

    /**
     * Evaluates the B-Spline at an arbitrary position t within the interval [0.0,1.0].
     */
    vec3 getPoint(GLfloat t);

    /**
     * Calculates the B-Spline's first derivative (tangent) at an arbitrary position t
     * within the interval [0.0,1.0]. Not implemented, yet.
     */
    vec3 getDerivative(GLfloat t);

private:
    // Cox-deBoor recursion
    float evalBasisFunction(int knotID, int deg, float u);

    // generates the spline's knot values
    void generateKnotVector();

    // degree and control points
    int degree_;
    std::vector<vec3> ctrlPoints_;

    // knot vector (internally generated)
    std::vector<float> knots_;
};

} // namespace

#endif // TGT_BSPLINE_H
