/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2019 University of Muenster, Germany,                        *
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

#include "flow2d.h"
#include "modules/flowreen/utils/flowmath.h"

namespace voreen {

Flow2D::Flow2D(const tgt::vec2* const flow2D, const tgt::ivec2& dimensions, const DataOrientation& orientation,
               const float minValue, const float maxValue, const float maxMagnitude, bool freeOnDestruction)
    : flow2D_(flow2D),
    dimensions_(dimensions),
    orientation_(orientation),
    axisPermutation_(getAxisPermutation(orientation)),
    minValue_(minValue),
    maxValue_(maxValue),
    maxMagnitude_(maxMagnitude),
    freeOnDestruction_(freeOnDestruction)
{
}

Flow2D::~Flow2D() {
    if (freeOnDestruction_ == true)
        delete [] flow2D_;
}

tgt::ivec2 Flow2D::cellNumberToPos(const size_t n) const {
    tgt::ivec2 pos(0, 0);
    const int& i = axisPermutation_[0];
    const int& j = axisPermutation_[1];

    pos[i] = (n % dimensions_[i]);
    pos[j] = (static_cast<int>(n) / dimensions_[j]);
    return pos;
}

tgt::ivec2 Flow2D::flowPosToSlicePos(const tgt::vec2& fp, const tgt::ivec2& sliceSize,
                                        const tgt::ivec2& offset) const
{
    tgt::vec2 aux(fp / static_cast<tgt::vec2>(dimensions_));
    return tgt::ivec2(
        (static_cast<int>(tgt::round(aux.x * sliceSize.x)) + offset.x) % sliceSize.x,
        (static_cast<int>(tgt::round(aux.y * sliceSize.y)) + offset.y) % sliceSize.y);
}

std::vector<tgt::ivec2> Flow2D::flowPosToSlicePos(std::vector<tgt::vec2>& fps,
        const tgt::ivec2& sliceSize, const tgt::ivec2& offset) const
{
    std::vector<tgt::ivec2> sps;    // slice positions
    while (fps.empty() == false) {
        std::vector<tgt::vec2>::iterator it = fps.begin();

        // no call to overloaded flowPosToViewportPos() for performance
        //
        tgt::vec2 aux(*it / static_cast<tgt::vec2>(dimensions_));
        tgt::ivec2 sp(
            (static_cast<int>(tgt::round(aux.x * sliceSize.x)) + offset.x) % sliceSize.x,
            (static_cast<int>(tgt::round(aux.y * sliceSize.y)) + offset.y) % sliceSize.y);

        sps.push_back(sp);
        fps.erase(it);
    }
    return sps;
}

std::deque<tgt::ivec2> Flow2D::flowPosToSlicePos(std::deque<tgt::vec2>& fps,
        const tgt::ivec2& sliceSize, const tgt::ivec2& offset) const
{
    std::deque<tgt::ivec2> sps;    // slice positions
    while (fps.empty() == false) {
        std::deque<tgt::vec2>::iterator it = fps.begin();

        // no call to overloaded flowPosToViewportPos() for performance
        //
        tgt::vec2 aux(*it / static_cast<tgt::vec2>(dimensions_));
        tgt::ivec2 sp(
            (static_cast<int>(tgt::round(aux.x * sliceSize.x)) + offset.x) % sliceSize.x,
            (static_cast<int>(tgt::round(aux.y * sliceSize.y)) + offset.y) % sliceSize.y);

        sps.push_back(sp);
        fps.erase(it);
    }
    return sps;
}

tgt::ivec2 Flow2D::getNearestCellPos(const tgt::vec2& r) const {
    return tgt::ivec2(tgt::clamp(static_cast<int>(r.x + 0.5f), 0, dimensions_.x),
        tgt::clamp(static_cast<int>(r.y + 0.5f), 0, dimensions_.y));
    /*return tgt::ivec2(tgt::clamp(static_cast<int>(r.x), 0, dimensions_.x),
        tgt::clamp(static_cast<int>(r.y), 0, dimensions_.y));*/
}

tgt::vec3* Flow2D::getNormalizedFlow() const {
    const size_t textureSize = dimensions_.x * dimensions_.y;
    tgt::vec3* texture = new tgt::vec3[textureSize];
    memset(texture, 0, textureSize * sizeof(tgt::vec3));

    const float range = maxValue_ - minValue_;
    for (size_t i = 0; i < textureSize; ++i) {
        if (flow2D_[i] != tgt::vec2(0.0f))
            texture[i] = tgt::vec3((flow2D_[i].x - minValue_) / range, (flow2D_[i].y - minValue_) / range, 0.0);
    }
    return texture;
}

tgt::mat3 Flow2D::getTransformationMatrix(const tgt::vec2& r, const float scaling) const {
    tgt::vec2 tangent = FlowMath::normalize(lookupFlow(r));
    if (tangent == tgt::vec2::zero)
        return tgt::mat3(scaling, 0.0f, 0.0f,
                         0.0f, scaling, 0.0f,
                         0.0f, 0.0f, 1.0f);

    tgt::vec2 normal(tangent.y, -tangent.x);

    if (scaling != 1.0f) {
        normal *= scaling;
        tangent *= scaling;
    }

    return tgt::mat3(normal.x, tangent.x, r.x,
        normal.y, tangent.y, r.y,
        0.0f, 0.0f, 1.0f);
}

bool Flow2D::isInsideBoundings(const tgt::vec2& r) const {
    for (size_t i = 0; i < 2; ++i) {
        if (r[i] < 0.0f)
            return false;
        else if (static_cast<int>(ceilf(r[i])) >= dimensions_[i])
            return false;
    }
    return true;
}

size_t Flow2D::posToCellNumber(const tgt::ivec2& pos) const {
    const int& i = axisPermutation_[0];
    const int& j = axisPermutation_[1];
    return (pos[i] + pos[j] * dimensions_[i]);
}

tgt::ivec2 Flow2D::slicePosToFlowPos(const tgt::ivec2& vp, const tgt::ivec2& sliceSize,
                                        tgt::ivec2* const error) const
{
    tgt::vec2 p(static_cast<float>(vp.x % sliceSize.x) / sliceSize.x,
        static_cast<float>(vp.y % sliceSize.y) / sliceSize.y);

    tgt::ivec2 r = getNearestCellPos(p * static_cast<tgt::vec2>(dimensions_));
    if (error != 0)
        *error = (vp - flowPosToSlicePos(r, sliceSize));
    return r;
}

// static public methods
//

tgt::ivec2 Flow2D::getAxisPermutation(const DataOrientation& orientation) {
    if (orientation == YX)
        return tgt::ivec2(1, 0);
    return tgt::ivec2(0, 1);
}

// private methods
//

tgt::vec2 Flow2D::lookupFlowBilinear(const tgt::vec2& r) const {
    tgt::ivec2 irll = tgt::clamp(static_cast<tgt::ivec2>(r), tgt::ivec2::zero, dimensions_);
    tgt::ivec2 irlr = tgt::clamp(tgt::ivec2(irll.x + 1, irll.y), tgt::ivec2::zero, dimensions_);
    tgt::ivec2 irul = tgt::clamp(tgt::ivec2(irll.x, irll.y + 1), tgt::ivec2::zero, dimensions_);
    tgt::ivec2 irur = tgt::clamp(tgt::ivec2(irlr.x, irul.y), tgt::ivec2::zero, dimensions_);

    const tgt::vec2& vll = flow2D_[posToCellNumber(irll)];
    const tgt::vec2& vlr = flow2D_[posToCellNumber(irlr)];
    const tgt::vec2& vul = flow2D_[posToCellNumber(irul)];
    const tgt::vec2& vur = flow2D_[posToCellNumber(irur)];

    // interpolation in x-direction
    //
    float fintegral = 0.0f;
    float b = modff(r.x, &fintegral);
    float a = 1.0f - b;
    tgt::vec2 vl = vll * a + vlr * b;
    tgt::vec2 vu = vul * a + vur * b;

    // interpolation in y-direction
    //
    b = modff(r.y, &fintegral);
    a = 1.0f - b;
    return (vl * a + vu * b);
}

}   // namespace
