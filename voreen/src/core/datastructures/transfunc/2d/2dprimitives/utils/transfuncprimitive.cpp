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

#include "voreen/core/datastructures/transfunc/2d/2dprimitives/utils/transfuncprimitive.h"

#include "tgt/glmath.h"
#include "tgt/matrixstack.h"

#include "tgt/immediatemode/immediatemode.h"

namespace voreen {

TransFuncPrimitive::TransFuncPrimitive()
    : fuzziness_(1.f)
{}

size_t TransFuncPrimitive::getNumControlPoints() const {
    return controlPoints_.size();
}

void TransFuncPrimitive::setColor(const tgt::col4& color) {
    for (std::vector<TransFuncPrimitiveControlPoint>::iterator i = controlPoints_.begin(); i != controlPoints_.end(); ++i)
        i->color_ = color;
}

const TransFuncPrimitiveControlPoint& TransFuncPrimitive::getControlPoint(size_t index) const {
    return controlPoints_.at(index);
}
TransFuncPrimitiveControlPoint& TransFuncPrimitive::getControlPoint(size_t index) {
    return controlPoints_.at(index);
}

void TransFuncPrimitive::setFuzziness(float f) {
    fuzziness_ = f;
}

float TransFuncPrimitive::getFuzziness() const {
    return fuzziness_;
}

void TransFuncPrimitive::move(const tgt::vec2& offset) {
    for (std::vector<TransFuncPrimitiveControlPoint>::iterator i = controlPoints_.begin(); i != controlPoints_.end(); ++i)
        i->position_ += offset;
}

void TransFuncPrimitive::serialize(Serializer& s) const {
    s.serialize("fuzzy", fuzziness_);
    s.serialize("controlpoints", controlPoints_);
}

void TransFuncPrimitive::deserialize(Deserializer& s) {
    s.deserialize("fuzzy", fuzziness_);
    s.deserialize("controlpoints", controlPoints_);
}

//-----------------------------------------------------------------------------

TransFuncTriangle::TransFuncTriangle()
    : TransFuncPrimitive()
{
    // set some default values...
    controlPoints_.push_back(TransFuncPrimitiveControlPoint(tgt::vec2(0.3f, 0.3f), tgt::col4(255)));
    controlPoints_.push_back(TransFuncPrimitiveControlPoint(tgt::vec2(0.6f, 0.3f), tgt::col4(255)));
    controlPoints_.push_back(TransFuncPrimitiveControlPoint(tgt::vec2(0.45f, 0.6f), tgt::col4(255)));
}

TransFuncTriangle::TransFuncTriangle(const tgt::vec2& a, const tgt::vec2& b, const tgt::vec2& c, const tgt::col4& col)
    : TransFuncPrimitive()
{
    // add the control points
    controlPoints_.push_back(TransFuncPrimitiveControlPoint(a, col));
    controlPoints_.push_back(TransFuncPrimitiveControlPoint(b, col));
    controlPoints_.push_back(TransFuncPrimitiveControlPoint(c, col));
}

void TransFuncTriangle::paint() {

    // compute the center of the triangle and its color as the average of the control points
    tgt::vec2 center(0.f);
    tgt::vec4 centerColor(0.f);
    for (size_t i = 0; i < 3; ++i) {
        center += controlPoints_.at(i).position_;
        centerColor += controlPoints_.at(i).getColorNormalized();
    }
    center /= 3.f;
    centerColor /= 3.f;

    // compute middle points for each edge between adjacent faded out control points at the edge of the triangle
    tgt::vec2 middlePositions[3];
    tgt::vec4 middleColors[3];
    for (size_t i = 0; i < 3; ++i) {
        middlePositions[i] = 0.5f * controlPoints_.at(i).position_ + 0.5f * controlPoints_.at((i+1) % 3).position_;
        middleColors[i] = tgt::vec4(0.5f * controlPoints_.at(i).getColorNormalized().xyz() + 0.5f * controlPoints_.at((i+1) % 3).getColorNormalized().xyz(), 0.f);
    }

    // compute inner triangle coordinates
    std::vector<TransFuncPrimitiveControlPoint> innerPoints(3);
    for (size_t i = 0; i < 3; ++i) {
        innerPoints.at(i).position_ = fuzziness_ * controlPoints_.at(i).position_ + (1.f - fuzziness_) * center;
        innerPoints.at(i).color_ = controlPoints_.at(i).color_;
    }

    // set transformation
    MatStack.translate(0.f, 0.f, -0.5f);

    // paint outer triangles using triangle strips for creating the fuzziness
    // the middle points are needed to create a reasonable interpolation for different vertex colors
    IMode.begin(tgt::ImmediateMode::TRIANGLE_STRIP);
        IMode.color(tgt::vec4(controlPoints_.at(0).getColorNormalized().xyz(), 0.f));
        IMode.vertex(controlPoints_.at(0).position_);

        IMode.color(innerPoints.at(0).getColorNormalized());
        IMode.vertex(innerPoints.at(0).position_);

        IMode.color(middleColors[2]);
        IMode.vertex(middlePositions[2]);

        IMode.color(innerPoints.at(2).getColorNormalized());
        IMode.vertex(innerPoints.at(2).position_);

        IMode.color(tgt::vec4(controlPoints_.at(2).getColorNormalized().xyz(), 0.f));
        IMode.vertex(controlPoints_.at(2).position_);

        IMode.color(middleColors[1]);
        IMode.vertex(middlePositions[1]);
    IMode.end();
    IMode.begin(tgt::ImmediateMode::TRIANGLE_STRIP);
        IMode.color(innerPoints.at(2).getColorNormalized());
        IMode.vertex(innerPoints.at(2).position_);

        IMode.color(innerPoints.at(1).getColorNormalized());
        IMode.vertex(innerPoints.at(1).position_);

        IMode.color(middleColors[1]);
        IMode.vertex(middlePositions[1]);

        IMode.color(tgt::vec4(controlPoints_.at(1).getColorNormalized().xyz(), 0.f));
        IMode.vertex(controlPoints_.at(1).position_);
    IMode.end();
    IMode.begin(tgt::ImmediateMode::TRIANGLE_STRIP);
        IMode.color(tgt::vec4(controlPoints_.at(1).getColorNormalized().xyz(), 0.f));
        IMode.vertex(controlPoints_.at(1).position_);

        IMode.color(innerPoints.at(1).getColorNormalized());
        IMode.vertex(innerPoints.at(1).position_);

        IMode.color(middleColors[0]);
        IMode.vertex(middlePositions[0]);

        IMode.color(innerPoints.at(0).getColorNormalized());
        IMode.vertex(innerPoints.at(0).position_);

        IMode.color(tgt::vec4(controlPoints_.at(0).getColorNormalized().xyz(), 0.f));
        IMode.vertex(controlPoints_.at(0).position_);
    IMode.end();

    // paint the inner triangle which contains the full unfaded colors
    IMode.begin(tgt::ImmediateMode::TRIANGLES);
        for (size_t i = 0; i < 3; ++i) {
            IMode.color(innerPoints.at(i).getColorNormalized());
            IMode.vertex(innerPoints.at(i).position_);
        }
    IMode.end();

    //clean up
    MatStack.translate(0.f, 0.f, 0.5f);
    IMode.color(tgt::vec4::one);
}

float TransFuncTriangle::getClosestControlPointDist(const tgt::vec2& pos) {
    float min = distance(pos, controlPoints_.at(0).position_);
    float d;
    for (size_t i = 1; i < 3; ++i) {
        d = distance(pos, controlPoints_.at(i).position_);
        if (d < min)
            min = d;
    }
    return min;
}

TransFuncPrimitive* TransFuncTriangle::clone() const {
    TransFuncTriangle* prim = new TransFuncTriangle();

    prim->fuzziness_ = fuzziness_;
    prim->controlPoints_ = controlPoints_;

    return prim;
}

//-----------------------------------------------------------------------------

TransFuncQuad::TransFuncQuad()
    : TransFuncPrimitive()
{
    // set some default values...
    for (size_t i = 0; i < 4; ++i) {
        TransFuncPrimitiveControlPoint a(tgt::vec2(static_cast<float>(((i+1) / 2) % 2) / 2.f, static_cast<float>(i > 1 ? 0 : 1) / 2.f), tgt::col4(255));
        controlPoints_.push_back(a);
    }
}

TransFuncQuad::TransFuncQuad(const tgt::vec2& center, float size, const tgt::col4& col)
    : TransFuncPrimitive()
{
    //for (std::vector<TransFuncPrimitiveControlPoint>::iterator i = controlPoints_.begin(); i != controlPoints_.end(); ++i) {
    for (size_t i = 0; i < 4; i++) {
        TransFuncPrimitiveControlPoint a;
        a.color_ = col;
        controlPoints_.push_back(a);
    }
    size *= 0.5f;
    controlPoints_.at(0).position_ = center + tgt::vec2(-size, -size);
    controlPoints_.at(1).position_ = center + tgt::vec2( size, -size);
    controlPoints_.at(2).position_ = center + tgt::vec2( size,  size);
    controlPoints_.at(3).position_ = center + tgt::vec2(-size,  size);
}

void TransFuncQuad::paint() {

    // compute the center of the quad and its color as the average of the control points
    tgt::vec2 center(0.f);
    tgt::vec4 centerColor(0.f);
    for (size_t i = 0; i < 4; ++i) {
        center += controlPoints_.at(i).position_;
        centerColor += controlPoints_.at(i).getColorNormalized();
    }
    center /= 4.f;
    centerColor /= 4.f;

    // compute middle points for each edge between adjacent faded out control points at the edge of the quad
    tgt::vec2 middlePositions[4];
    tgt::vec4 middleColors[4];
    for (size_t i = 0; i < 4; ++i) {
        middlePositions[i] = 0.5f * controlPoints_.at(i).position_ + 0.5f * controlPoints_.at((i+1) % 4).position_;
        middleColors[i] = tgt::vec4(0.5f * controlPoints_.at(i).getColorNormalized().xyz() + 0.5f * controlPoints_.at((i+1) % 4).getColorNormalized().xyz(), 0.f);
    }

    // compute inner quad coordinates
    std::vector<TransFuncPrimitiveControlPoint> innerPoints(4);
    for (size_t i = 0; i < 4; ++i) {
        innerPoints.at(i).position_ = fuzziness_ * controlPoints_.at(i).position_ + (1.f - fuzziness_) * center;
        innerPoints.at(i).color_ = controlPoints_.at(i).color_;
    }

    // set transformation
    MatStack.translate(0.f, 0.f, -0.5f);

    // paint outer triangles using triangle strips for creating the fuzziness
    // the middle points are needed to create a reasonable interpolation for different vertex colors
    IMode.begin(tgt::ImmediateMode::TRIANGLE_STRIP);
        IMode.color(tgt::vec4(controlPoints_.at(0).getColorNormalized().xyz(), 0.f));
        IMode.vertex(controlPoints_.at(0).position_);

        IMode.color(innerPoints.at(0).getColorNormalized());
        IMode.vertex(innerPoints.at(0).position_);

        IMode.color(middleColors[3]);
        IMode.vertex(middlePositions[3]);

        IMode.color(innerPoints.at(3).getColorNormalized());
        IMode.vertex(innerPoints.at(3).position_);

        IMode.color(tgt::vec4(controlPoints_.at(3).getColorNormalized().xyz(), 0.f));
        IMode.vertex(controlPoints_.at(3).position_);

        IMode.color(middleColors[2]);
        IMode.vertex(middlePositions[2]);
    IMode.end();
    IMode.begin(tgt::ImmediateMode::TRIANGLE_STRIP);
        IMode.color(innerPoints.at(3).getColorNormalized());
        IMode.vertex(innerPoints.at(3).position_);

        IMode.color(innerPoints.at(2).getColorNormalized());
        IMode.vertex(innerPoints.at(2).position_);

        IMode.color(middleColors[2]);
        IMode.vertex(middlePositions[2]);

        IMode.color(tgt::vec4(controlPoints_.at(2).getColorNormalized().xyz(), 0.f));
        IMode.vertex(controlPoints_.at(2).position_);
    IMode.end();
    IMode.begin(tgt::ImmediateMode::TRIANGLE_STRIP);
        IMode.color(tgt::vec4(controlPoints_.at(2).getColorNormalized().xyz(), 0.f));
        IMode.vertex(controlPoints_.at(2).position_);

        IMode.color(innerPoints.at(2).getColorNormalized());
        IMode.vertex(innerPoints.at(2).position_);

        IMode.color(middleColors[1]);
        IMode.vertex(middlePositions[1]);

        IMode.color(innerPoints.at(1).getColorNormalized());
        IMode.vertex(innerPoints.at(1).position_);

        IMode.color(tgt::vec4(controlPoints_.at(1).getColorNormalized().xyz(), 0.f));
        IMode.vertex(controlPoints_.at(1).position_);

        IMode.color(middleColors[0]);
        IMode.vertex(middlePositions[0]);
    IMode.end();
    IMode.begin(tgt::ImmediateMode::TRIANGLE_STRIP);
        IMode.color(innerPoints.at(1).getColorNormalized());
        IMode.vertex(innerPoints.at(1).position_);

        IMode.color(innerPoints.at(0).getColorNormalized());
        IMode.vertex(innerPoints.at(0).position_);

        IMode.color(middleColors[0]);
        IMode.vertex(middlePositions[0]);

        IMode.color(tgt::vec4(controlPoints_.at(0).getColorNormalized().xyz(), 0.f));
        IMode.vertex(controlPoints_.at(0).position_);
    IMode.end();

    // paint the inner quad which contains the full unfaded colors
    IMode.begin(tgt::ImmediateMode::TRIANGLE_FAN);
        IMode.color(centerColor);
        IMode.vertex(center);
        for (size_t i = 0; i < 4; ++i) {
            IMode.color(innerPoints.at(i).getColorNormalized());
            IMode.vertex(innerPoints.at(i).position_);
        }
        IMode.color(innerPoints.at(0).getColorNormalized());
        IMode.vertex(innerPoints.at(0).position_);
    IMode.end();

    //clean up
    MatStack.translate(0.f, 0.f, 0.5f);
    IMode.color(tgt::vec4::one);
}

float TransFuncQuad::getClosestControlPointDist(const tgt::vec2& pos) {
    float min = distance(pos, controlPoints_.at(0).position_);
    float d;
    for (size_t i = 1; i < 4; ++i) {
        d = distance(pos, controlPoints_.at(i).position_);
        if (d < min)
            min = d;
    }
    return min;
}

TransFuncPrimitive* TransFuncQuad::clone() const {
    TransFuncQuad* prim = new TransFuncQuad();

    prim->fuzziness_ = fuzziness_;
    prim->controlPoints_ = controlPoints_;

    return prim;
}

//-----------------------------------------------------------------------------

TransFuncBanana::TransFuncBanana()
    : TransFuncPrimitive()
    , steps_(20)
{
    // set some default values
    controlPoints_.push_back(TransFuncPrimitiveControlPoint(tgt::vec2(0.25f), tgt::col4(255)));
    controlPoints_.push_back(TransFuncPrimitiveControlPoint(tgt::vec2(0.5f, 0.5f), tgt::col4(255)));
    controlPoints_.push_back(TransFuncPrimitiveControlPoint(tgt::vec2(0.5f, 0.33f), tgt::col4(255)));
    controlPoints_.push_back(TransFuncPrimitiveControlPoint(tgt::vec2(0.74f, 0.25f), tgt::col4(255)));
}

TransFuncBanana::TransFuncBanana(const tgt::vec2& a, const tgt::vec2& b1, const tgt::vec2& b2, const tgt::vec2& c, const tgt::col4& col)
    : TransFuncPrimitive()
    , steps_(20)
{
    // add the control points
    controlPoints_.push_back(TransFuncPrimitiveControlPoint(a, col));
    controlPoints_.push_back(TransFuncPrimitiveControlPoint(b1, col));
    controlPoints_.push_back(TransFuncPrimitiveControlPoint(b2, col));
    controlPoints_.push_back(TransFuncPrimitiveControlPoint(c, col));
}

void TransFuncBanana::paint() {
    MatStack.translate(0.f, 0.f, -0.5f);
    //glColor4ubv(color_.elem);
    paintInner();
    MatStack.translate(0.f, 0.f, 0.5f);
}

void TransFuncBanana::paintInner() {
    float t;
    tgt::vec2 v1, v2, t1, t2, t3, t4, tc;
    tgt::vec4 v1Color, v2Color, t1Color, t2Color, t3Color, t4Color, tcColor;

    // compute a middle point for the upper bezier curve (t1) and the lower bezier curve (t2)
    t1 = (2.f * controlPoints_.at(1).position_) - (0.5f * controlPoints_.at(0).position_) - (0.5f * controlPoints_.at(3).position_);
    t2 = (2.f * controlPoints_.at(2).position_) - (0.5f * controlPoints_.at(0).position_) - (0.5f * controlPoints_.at(3).position_);
    t1Color = (2.f * controlPoints_.at(1).getColorNormalized()) - (0.5f * controlPoints_.at(0).getColorNormalized()) - (0.5f * controlPoints_.at(3).getColorNormalized());
    t2Color = (2.f * controlPoints_.at(2).getColorNormalized()) - (0.5f * controlPoints_.at(0).getColorNormalized()) - (0.5f * controlPoints_.at(3).getColorNormalized());

    tc = (t1 + t2) / 2.f;   // tc is the center between the middle points of the two curves
    tcColor = (t1Color + t2Color) / 2.f;

    // compute inner middle points which are the middle points of inner curves that have full color (while the outer curves are faded out using the fuzziness parameter)
    t3 = fuzziness_ * t1 + (1.f - fuzziness_) * tc;
    t4 = fuzziness_ * t2 + (1.f - fuzziness_) * tc;
    t3Color = fuzziness_ * t1Color + (1.f - fuzziness_) * tcColor;
    t4Color = fuzziness_ * t2Color + (1.f - fuzziness_) * tcColor;

    // fill the space between the two bezier curves:
    IMode.begin(tgt::ImmediateMode::TRIANGLE_STRIP);
        // start with the space that fades from the inner fully colored space to the upper curve (therefore using t1 for the upper part and t3 for the lower part)
        IMode.color(controlPoints_.at(0).getColorNormalized());
        IMode.vertex(controlPoints_.at(0).position_);
        for (size_t i = 0; i < steps_; ++i) {
            t = static_cast<float>(i) / static_cast<float>(steps_ - 1);
            // find vertices (and colors) using bezier interpolation
            v1 = (((1 - t) * (1 - t)) * controlPoints_.at(0).position_) + ((2 * (1 - t) * t) * t1) + ((t * t) * controlPoints_.at(3).position_);
            v2 = (((1 - t) * (1 - t)) * controlPoints_.at(0).position_) + ((2 * (1 - t) * t) * t3) + ((t * t) * controlPoints_.at(3).position_);

            v1Color = (((1 - t) * (1 - t)) * controlPoints_.at(0).getColorNormalized()) + ((2 * (1 - t) * t) * t1Color) + ((t * t) * controlPoints_.at(3).getColorNormalized());
            v2Color = (((1 - t) * (1 - t)) * controlPoints_.at(0).getColorNormalized()) + ((2 * (1 - t) * t) * t3Color) + ((t * t) * controlPoints_.at(3).getColorNormalized());

            // render the vertices with the corresponding colors
            IMode.color(tgt::vec4(v1Color.xyz(), 0.f));
            IMode.vertex(v1);
            IMode.color(v2Color);
            IMode.vertex(v2);
        }
        IMode.color(controlPoints_.at(3).getColorNormalized());
        IMode.vertex(controlPoints_.at(3).position_);

        // now render the inner part with full color which uses t3 and t4 as the inner middle control points
        IMode.color(controlPoints_.at(0).getColorNormalized());
        IMode.vertex(controlPoints_.at(0).position_);
        for (size_t i = 0; i < steps_; ++i) {
            t = static_cast<float>(i) / static_cast<float>(steps_ - 1);
            v1 = (((1 - t) * (1 - t)) * controlPoints_.at(0).position_) + ((2 * (1 - t) * t) * t3) + ((t * t) * controlPoints_.at(3).position_);
            v2 = (((1 - t) * (1 - t)) * controlPoints_.at(0).position_) + ((2 * (1 - t) * t) * t4) + ((t * t) * controlPoints_.at(3).position_);

            v1Color = (((1 - t) * (1 - t)) * controlPoints_.at(0).getColorNormalized()) + ((2 * (1 - t) * t) * t3Color) + ((t * t) * controlPoints_.at(3).getColorNormalized());
            v2Color = (((1 - t) * (1 - t)) * controlPoints_.at(0).getColorNormalized()) + ((2 * (1 - t) * t) * t4Color) + ((t * t) * controlPoints_.at(3).getColorNormalized());

            IMode.color(v1Color);
            IMode.vertex(v1);
            IMode.color(v2Color);
            IMode.vertex(v2);
        }
        IMode.color(controlPoints_.at(3).getColorNormalized());
        IMode.vertex(controlPoints_.at(3).position_);

        // now render the lower part which fades out from the curve using t4 to the curve using t2
        IMode.color(controlPoints_.at(0).getColorNormalized());
        IMode.vertex(controlPoints_.at(0).position_);
        for (size_t i = 0; i < steps_; ++i) {
            t = static_cast<float>(i) / static_cast<float>(steps_ - 1);
            v1 = (((1 - t) * (1 - t)) * controlPoints_.at(0).position_) + ((2 * (1 - t) * t) * t4) + ((t * t) * controlPoints_.at(3).position_);
            v2 = (((1 - t) * (1 - t)) * controlPoints_.at(0).position_) + ((2 * (1 - t) * t) * t2) + ((t * t) * controlPoints_.at(3).position_);

            v1Color = (((1 - t) * (1 - t)) * controlPoints_.at(0).getColorNormalized()) + ((2 * (1 - t) * t) * t4Color) + ((t * t) * controlPoints_.at(3).getColorNormalized());
            v2Color = (((1 - t) * (1 - t)) * controlPoints_.at(0).getColorNormalized()) + ((2 * (1 - t) * t) * t2Color) + ((t * t) * controlPoints_.at(3).getColorNormalized());

            IMode.color(v1Color);
            IMode.vertex(v1);
            IMode.color(tgt::vec4(v2Color.xyz(), 0.f));
            IMode.vertex(v2);
        }
        IMode.color(controlPoints_.at(3).getColorNormalized());
        IMode.vertex(controlPoints_.at(3).position_);

    IMode.end();
    //clean up
     IMode.color(tgt::vec4::one);

}

float TransFuncBanana::getClosestControlPointDist(const tgt::vec2& pos) {
    float min = distance(pos, controlPoints_.at(0).position_);
    float d;
    for (size_t i = 1; i < 4; ++i) {
        d = distance(pos, controlPoints_.at(i).position_);
        if (d < min)
            min = d;
    }
    return min;
}

TransFuncPrimitive* TransFuncBanana::clone() const {
    TransFuncBanana* prim = new TransFuncBanana();

    prim->fuzziness_ = fuzziness_;
    prim->controlPoints_ = controlPoints_;
    prim->steps_ = steps_;

    return prim;
}

void TransFuncBanana::serialize(Serializer& s) const {
    TransFuncPrimitive::serialize(s);
    s.serialize("steps", steps_);
}

void TransFuncBanana::deserialize(Deserializer& s) {
    TransFuncPrimitive::deserialize(s);
    s.deserialize("steps", steps_);
}

} // namespace voreen
