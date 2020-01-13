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

#include "catmullromspline.h"

//using tgt::vec2;
using std::vector;

namespace tgt {

CatmullRomSpline::CatmullRomSpline() {}

CatmullRomSpline::CatmullRomSpline(const vector<vec2>& points) {
    setPoints(points);
}

void CatmullRomSpline::setPoints(const vector<vec2>& points) {
    points_ = points;

    if (points_.size() > 1) {
        const vec2& first = points[0];
        const vec2& second = points[1];

        float deltaX = second.x - first.x;
        float deltaY = second.y - first.y;

        if (deltaX == 0.f)
            deltaX += 0.001f;

        if (deltaY == 0.f)
            deltaY += 0.001f;

        vec2 virtualPoint(first.x - deltaX, first.y - deltaY);
        points_.insert(points_.begin(), virtualPoint);
    }
}

float CatmullRomSpline::getPoint(float t) const {
    if (points_.size() == 1) {
        return points_[0].y;
    }
    else {
        float delta = points_[points_.size()-1].x - points_[1].x;
        float x = t * delta;

        size_t i = 0;
        while (points_[i].x < x)
            ++i;

        float x2 = points_[i].x;
        float p2 = points_[i].y;
        float x1 = points_[i-1].x;
        float p1 = points_[i-1].y;
        float m1;
        if (i != 1) {
            float p0 = points_[i-2].y;
            m1 = 0.5f * (p2 - p0);
        }
        else
            m1 = p2 - p1;

        i++;
        float m2;
        if (i < points_.size()) {
            float p3 = points_[i].y;
            m2 = 0.5f * (p3 - p1);
        }
        else
            m2 = p2 - p1;

        const float t2 = (x - x1)/(x2 - x1);
        const float h00 = (1 + 2*t2)*(1 - t2)*(1 - t2);
        const float h10 = t2*(1 - t2)*(1 - t2);
        const float h01 = t2*t2*(3-2*t2);
        const float h11 = t2*t2*(t2 - 1);
        return h00*p1 + h10*m1 + h01*p2 + h11*m2;
    }
}

} // namespace
