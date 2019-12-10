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

#include "voreen/core/animation/interpolation/boundingboxinterpolationfunctions.h"
#include "voreen/core/animation/interpolation/basicintinterpolation.h"
#include "voreen/core/animation/interpolation/basicfloatinterpolation.h"
#include <cmath>
namespace voreen {


// MSVC11 does not support round!?!
static float roundToNearest(float num) {
    return (num > 0.0) ? std::floor(num + 0.5f) : std::ceil(num - 0.5f);
}

static int interpolate(int s, int e, float t){
    t = tgt::clamp(t, 0.0f, 1.0f);
    float i = static_cast<float>(s)*(1.0f-t)+static_cast<float>(e)*t;
    return static_cast<int>(roundToNearest(i));

}

tgt::IntBounds LinearInterpolate(tgt::IntBounds startvalue,
                                 tgt::IntBounds endvalue,
                                 float time){
    tgt::ivec3 sllf = startvalue.getLLF();
    tgt::ivec3 surb = startvalue.getURB();

    tgt::ivec3 ellf = endvalue.getLLF();
    tgt::ivec3 eurb = endvalue.getURB();

    tgt::ivec3 llf = tgt::ivec3(interpolate(sllf.x, ellf.x, time),
                                interpolate(sllf.y, ellf.y, time),
                                interpolate(sllf.z, ellf.z, time));

    tgt::ivec3 urb = tgt::ivec3(interpolate(surb.x, eurb.x, time),
                                interpolate(surb.y, eurb.y, time),
                                interpolate(surb.z, eurb.z, time));
    return tgt::IntBounds(llf, urb);
}


tgt::IntBounds IntBoundingboxLinearInterpolationFunction::interpolate(tgt::IntBounds startvalue,
                                                                   tgt::IntBounds endvalue,
                                                                   float time) const{
    return LinearInterpolate(startvalue, endvalue, time);
}

tgt::IntBounds IntBoundingboxStartInterpolationFunction::interpolate(tgt::IntBounds startvalue,
                                                                  tgt::IntBounds endvalue,
                                                                  float time) const{
    return LinearInterpolate(startvalue, endvalue, 0);
}

tgt::IntBounds IntBoundingboxEndInterpolationFunction::interpolate(tgt::IntBounds startvalue,
                                                                tgt::IntBounds endvalue,
                                                                float time) const{
    return LinearInterpolate(startvalue, endvalue, 1);
}

tgt::IntBounds IntBoundingboxStartEndInterpolationFunction::interpolate(tgt::IntBounds startvalue,
                                                                     tgt::IntBounds endvalue,
                                                                     float time) const{
    if (time < 0.5f)
        return LinearInterpolate(startvalue, endvalue, 0);
    else
        return LinearInterpolate(startvalue, endvalue, 1);
}

tgt::Bounds LinearInterpolate(tgt::Bounds startvalue,
                                 tgt::Bounds endvalue,
                                 float time){
    tgt::vec3 sllf = startvalue.getLLF();
    tgt::vec3 surb = startvalue.getURB();

    tgt::vec3 ellf = endvalue.getLLF();
    tgt::vec3 eurb = endvalue.getURB();

    tgt::vec3 llf = tgt::vec3(BasicFloatInterpolation::linearInterpolation(sllf.x, ellf.x, time),
                                BasicFloatInterpolation::linearInterpolation(sllf.y, ellf.y, time),
                                BasicFloatInterpolation::linearInterpolation(sllf.z, ellf.z, time));

    tgt::vec3 urb = tgt::vec3(BasicFloatInterpolation::linearInterpolation(surb.x, eurb.x, time),
                                BasicFloatInterpolation::linearInterpolation(surb.y, eurb.y, time),
                                BasicFloatInterpolation::linearInterpolation(surb.z, eurb.z, time));
    return tgt::Bounds(llf, urb);
}


tgt::Bounds FloatBoundingboxLinearInterpolationFunction::interpolate(tgt::Bounds startvalue,
                                                                   tgt::Bounds endvalue,
                                                                   float time) const{
    return LinearInterpolate(startvalue, endvalue, time);
}

tgt::Bounds FloatBoundingboxStartInterpolationFunction::interpolate(tgt::Bounds startvalue,
                                                                  tgt::Bounds endvalue,
                                                                  float time) const{
    return LinearInterpolate(startvalue, endvalue, 0);
}

tgt::Bounds FloatBoundingboxEndInterpolationFunction::interpolate(tgt::Bounds startvalue,
                                                                tgt::Bounds endvalue,
                                                                float time) const{
    return LinearInterpolate(startvalue, endvalue, 1);
}

tgt::Bounds FloatBoundingboxStartEndInterpolationFunction::interpolate(tgt::Bounds startvalue,
                                                                     tgt::Bounds endvalue,
                                                                     float time) const{
    if (time < 0.5f)
        return LinearInterpolate(startvalue, endvalue, 0);
    else
        return LinearInterpolate(startvalue, endvalue, 1);
}

} // namespace voreen
