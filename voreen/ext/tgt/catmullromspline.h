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

#ifndef VRN_CATMULLROMSPLINE_H
#define VRN_CATMULLROMSPLINE_H

#include <vector>
#include "tgt/types.h"
#include "tgt/vector.h"

namespace tgt {

class TGT_API SplineTmp {
public:
    virtual float getPoint(float t) const = 0;
}; // temporary replacement for "real" spline class

class TGT_API CatmullRomSpline : public SplineTmp {
public:
    CatmullRomSpline();
    CatmullRomSpline(const std::vector<tgt::vec2>& points);

    void setPoints(const std::vector<tgt::vec2>& points);
    float getPoint(float t) const;

private:
    std::vector<tgt::vec2> points_;
};

} // namespace

#endif // VRN_CATMULLROMSPLINE_H
