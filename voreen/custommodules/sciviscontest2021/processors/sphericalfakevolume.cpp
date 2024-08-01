/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2024 University of Muenster, Germany,                        *
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

#include "sphericalfakevolume.h"

#include "voreen/core/datastructures/volume/volumedecorator.h"
#include "voreen/core/datastructures/meta/templatemetadata.h"

namespace voreen {

const std::string SphericalFakeVolume::loggerCat_("voreen.SphericalFakeVolume");

SphericalFakeVolume::SphericalFakeVolume()
        : VolumeCreateBase()
        , outport_(Port::OUTPORT, "volumehandle.output", "Volume Output", false)
        , radiusMax_("radiusMax", "Set Radius Max", 1000.0f, 0.0f, 10000.0f)
        , shiftFactor_("shiftFactor", "Set shift factor", 100.0f, 0.0f, 1000.0f)
        , regenerate_("regenerate", "Regenerate Volume", Processor::INVALID_RESULT)
{
    addPort(outport_);

    //Radius
    addProperty(radiusMax_);
    addProperty(shiftFactor_);
    radiusMax_.setGroupID("radius");
    shiftFactor_.setGroupID("radius");
    setPropertyGroupGuiName("radius","Radius");
}

Processor* SphericalFakeVolume::create() const {
    return new SphericalFakeVolume();
}

void SphericalFakeVolume::process() {
    //set offset to half radius divided shiftfactor
    float offset = -radiusMax_.get()/shiftFactor_.get()/2;

    //set spacing to radius divided by shiftfactor
    float spacing = radiusMax_.get()/shiftFactor_.get();

    auto* target = new VolumeRAM_UInt8(tgt::svec3(2, 2, 2));
    target->fill(std::numeric_limits<uint8_t>::max());

    auto volume = new Volume(target, tgt::vec3(spacing), tgt::vec3(offset));
    outport_.setData(volume);
}

}   // namespace
