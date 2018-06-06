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

#include "arbitraryvolumeclipping.h"

#include "tgt/glmath.h"

using tgt::ivec2;
using tgt::vec3;
using tgt::mat4;
using tgt::Texture;

namespace voreen {

const std::string ArbitraryVolumeClipping::loggerCat_("voreen.ArbitraryVolumeClipping");

ArbitraryVolumeClipping::ArbitraryVolumeClipping()
    : Processor()
    , depth_("depth", "Depth", 0.0f, 0.0f, 100.0f)
    , azimuth_("azimuth", "Azimuth", 0.0f, 0.0f, 360.0f)
    , elevation_("elevation", "Elevation", 0.0f, -90.0f, 90.0f)
    , plane_("planeNormal", "Plane Normal", vec3(1.0f, 0.0f, 0.0f), vec3(-5.0f), vec3(5.0f))
    , planeDist_("planePosition", "Plane Distance", 0.0f, -1000.0f, 1000.0f)
    , inport_(Port::INPORT, "volume", "Volume Input")
{
    addProperty(depth_);
    addProperty(azimuth_);
    addProperty(elevation_);

    addProperty(plane_);
    plane_.setInvalidationLevel(VALID); //output only properties
    addProperty(planeDist_);
    planeDist_.setInvalidationLevel(VALID);

    addPort(inport_);
}

ArbitraryVolumeClipping::~ArbitraryVolumeClipping() {
}

void ArbitraryVolumeClipping::process() {
    const VolumeBase* vol = inport_.getData();
    if(!vol)
        return;

    float az = (azimuth_.get() / 360.0f) * tgt::PIf * 2.0f;
    float el = ((elevation_.get() + 90.0f) / 180.0f) * tgt::PIf;
    vec3 n;
    n.x = sinf(el) * cosf(az);
    n.y = sinf(el) * sinf(az);
    n.z = cosf(el);

    vec3 center = (vol->getLLF() + vol->getURB()) * 0.5f;
    float lengthDiagonal = length(vol->getCubeSize());
    float distCenter = ((depth_.get() / 100.0f) - 0.5f) * lengthDiagonal;

    tgt::plane pl(n, -distCenter);
    pl = pl.transform(tgt::mat4::createTranslation(center));
    pl = pl.transform(vol->getPhysicalToWorldMatrix());

    plane_.set(pl.n);
    planeDist_.set(pl.d);
}


} // namespace voreen
