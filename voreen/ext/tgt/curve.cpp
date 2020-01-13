/**********************************************************************
 *                                                                    *
 * tgt - Tiny Graphics Toolbox                                        *
 *                                                                    *
 * Copyright (C) 2005-2020 University of Muenster, Germany,           *
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

#include "tgt/curve.h"
#include "tgt/types.h"

namespace tgt {


Curve::Curve(GLuint stepCount, const Bounds &bounds, bool _static, bool visible)
    : Renderable (bounds, _static, visible)
{
    setStepCount(stepCount);
    setDrawStyle(LINE);
}

void Curve::render() {
    if (visible_)
        render(0.0, 1.0);
}

void Curve::setStepCount(GLuint stepCount) {
    tgtAssert(stepCount > 0, "Step count is expected to be greater zero.");
    stepCount_ = stepCount;
}

GLuint Curve::getStepCount() {
    return stepCount_;
}

void Curve::setDrawStyle(DrawStyle drawStyle) {
    drawStyle_ = drawStyle;
}

tgt::Curve::DrawStyle Curve::getDrawStyle() {
    return drawStyle_;
}

void Curve::render(GLfloat startParam, GLfloat endParam) {
    if (!visible_)
        return ;

    tgtAssert(startParam >= 0 && endParam <= 1.0 && startParam <= endParam,
              "Expected: 0 <= startParam <= endParam <= 1.0 ");
    GLfloat step = 1.f / stepCount_;
    tgtAssert(step > 0, "Step must be positive.");
    if (step > 0) {
        if (drawStyle_ == LINE) {
            glBegin(GL_LINE_STRIP);
        } else if (drawStyle_ == POINT) {
            glBegin(GL_POINTS);
        }
        for (GLfloat t = startParam; t < endParam; t += step) {
            glVertex3fv(getPoint(t).elem);
        }
        glVertex3fv(getPoint(endParam).elem);
        glEnd();
    }
}


Circle::Circle()
    : Curve(100, Bounds(), false, true)
{
    setRadius(1.0);
    setCenter(vec3(0.0, 0.0, 0.0));
    setNormal(vec3(0.0, 0.0, 1.0));
    setRotationMatrix();
}

Circle::Circle(GLfloat radius, const vec3& center, const vec3& normal, GLuint stepCount /*= 100*/, bool _static /*= false*/, bool visible /*= true*/ )
    : Curve(stepCount, Bounds(), _static, visible)
{
    setRadius(radius);
    setCenter(center);
    setNormal(normal);
    setRotationMatrix();
}

// The circle point is first calculated in the x-y-plane, then translated by its center
// and then rotated by multiplying it with the prior calculated rotation matrix.
vec3 Circle::getPoint(GLfloat t) {
    tgtAssert(t >= 0.0f && t <= 1.0f, "Parameter t must be inside interval [0, 1].");
    t = t * 2.0f * PIf;
    vec3 position = radius_ * vec3(cosf(t), sinf(t), 0.0f);
    position = rotationMatrix_ * position;
    position += center_;
    return position;
}

// see getPoint
vec3 Circle::getDerivative(GLfloat t) {
    tgtAssert(t >= 0.0f && t <= 1.0f, "Parameter t must be inside interval [0, 1.0].");
    t = t * 2.0f * PIf;
    vec3 derive = radius_ * vec3( -sinf(t), cosf(t), 0.0f);
    derive = rotationMatrix_ * derive;
    return derive;
}

void Circle::setNormal(const vec3& normal) {
    tgtAssert((normal.elem[0] != 0) || (normal.elem[1] != 0) || (normal.elem[2] != 0),
              "Normal must not be the nullvector!");
    normal_ = normalize(normal);
    setRotationMatrix();
}

// calculates two perpendicular vectors that are also perpendicular the normal vector
// and creates an rotation matrix out of them.
void Circle::setRotationMatrix() {
    vec3 rx;
    vec3 ry;
    if (dot(normal_, vec3(1.0, 0.0, 0.0)) > 0.99 ) {
        rx = normalize(cross(normal_, vec3(0.0, 1.0, 0.0))) * static_cast<float>(isign(dot(normal_, vec3(0,0,1))));
        ry = normalize(cross(normal_, rx));
    } else {
        ry = normalize(cross(normal_, vec3(1.0, 0.0, 0.0))) * static_cast<float>(isign(dot(normal_, vec3(0,0,1))));
        rx = normalize(cross(ry, normal_));
    }
// warning: statement has no effect
//     sign(1);
    rotationMatrix_ = transpose(mat3(rx, ry, normal_));
}

void Circle::setRadius(GLfloat radius) {
    tgtAssert(radius > 0.f, "Radius must be positive." );
    radius_ = radius;
}

GLfloat Circle::getRadius() {
    return radius_;
}

void Circle::setCenter(const vec3& center) {
    center_ = center;
}

vec3 Circle::getCenter() {
    return center_;
}

vec3 Circle::getNormal() {
    return normal_;
}


Ellipse::Ellipse()
    : Curve(100, Bounds(), false, true)
{
    setCenter(vec3(0.0, 0.0, 0.0));
    setSemiMajorAxisLength(1.0);
    setSemiMinorAxisLength(1.0);
    setNormal(vec3(0.0, 0.0, 1.0));
    setSemiMajorAxisDirection(vec3(1.0, 0.0, 0.0));
}

Ellipse::Ellipse(const vec3& center, GLfloat semiMajorAxisLength, GLfloat semiMinorAxisLength, vec3 normal, vec3 semiMajorAxisDirection, GLuint stepCount /*= 100*/, bool _static /*= false*/, bool visible /*= true*/)
    : Curve(stepCount, Bounds(), _static, visible)
{
    setCenter(center);
    setSemiMajorAxisLength(semiMajorAxisLength);
    setSemiMinorAxisLength(semiMinorAxisLength);
    setNormal(normal);
    setSemiMajorAxisDirection(semiMajorAxisDirection);
}

// The ellipse point is first calculated in the x-y-plane, then translated by its center
// and then rotated by multiplying it with the prior calculated rotation matrix.
vec3 Ellipse::getPoint(GLfloat t) {
    tgtAssert(t >= 0.0f && t <= 1.0f, "Parameter t must be inside interval [0, 1.0].");
    t = t * 2.0f * PIf;
    vec3 position = vec3(a_ * cosf(t), b_ * sinf(t), 0.0f);
    position = rotationMatrix_ * position;
    position += center_;
    return position;
}

// see getPoint
vec3 Ellipse::getDerivative(GLfloat t) {
    tgtAssert(t >= 0.0f && t <= 1.0f, "Parameter t must be inside interval [0, 1.0].");
    t = t * 2 * PIf;
    vec3 derive = vec3( -a_ * sinf(t), b_ * cosf(t), 0.f);
    derive = rotationMatrix_ * derive;
    return derive;
}


// Creates an rotation matrix from the normal vector, the semi major and semi minor
// axis direction
void Ellipse::setRotationMatrix() {
    rotationMatrix_ = transpose(mat3(semiMajor_, semiMinor_, normal_));
}

// sets normal an recalculates semi major and semi minor axis direction
void Ellipse::setNormal(const vec3& normal) {
    tgtAssert((normal.elem[0] != 0) || (normal.elem[1] != 0) || (normal.elem[2] != 0),
              "Normal must not be the nullvector!");
    normal_ = normalize(normal);
    semiMinor_ = cross(normal_, semiMajor_);
    semiMajor_ = cross(semiMinor_, normal_);
    setRotationMatrix();
}

// semi major axis direction is projected into the plane, which is perpendicular to the normal
// semi minor axis direction is calculated from normal and semi major axis direction
void Ellipse::setSemiMajorAxisDirection(const vec3& semiMajor) {
    tgtAssert((semiMajor.elem[0] != 0) || (semiMajor.elem[1] != 0) || (semiMajor.elem[2] != 0),
              "Direction of the semi major axis must not be the nullvector!");
    semiMinor_ = normalize(cross(normal_, semiMajor));
    semiMajor_ = cross(semiMinor_, normal_);
    setRotationMatrix();
}

void Ellipse::setCenter(const vec3& center) {
    center_ = center;
}

vec3 Ellipse::getCenter() {
    return center_;
}

void Ellipse::setSemiMajorAxisLength(GLfloat a) {
    tgtAssert(a > 0, "Semi major axis length must be positive.");
    a_ = a;
}

GLfloat Ellipse::getSemiMajorAxisLength() {
    return a_;
}

void Ellipse::setSemiMinorAxisLength(GLfloat b) {
    tgtAssert(b > 0, "Semi minor axis length must be positive.");
    b_ = b;
}

GLfloat Ellipse::getSemiMinorAxisLength() {
    return b_;
}

vec3 Ellipse::getNormal() {
    return normal_;
}

vec3 Ellipse::getSemiMajorAxisDirection() {
    return semiMajor_;
}

vec3 Ellipse::getSemiMinorAxisDirection() {
    return semiMinor_;
}

} // namespace tgt
