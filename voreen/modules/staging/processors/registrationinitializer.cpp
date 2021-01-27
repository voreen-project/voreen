/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2021 University of Muenster, Germany,                        *
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

#include "registrationinitializer.h"
#include "voreen/core/datastructures/volume/volume.h"

namespace voreen {

using tgt::ivec2;
using tgt::vec3;
using tgt::vec4;
using tgt::mat4;

const std::string RegistrationInitializer::loggerCat_("voreen.RegistrationInitializer");

RegistrationInitializer::RegistrationInitializer()
    : Processor()
    , staticPort_(Port::INPORT, "static", "static Volume", Processor::VALID)
    , movingPort_(Port::INPORT, "moving", "Moving Volume", Processor::VALID)
    , transformMatrix_("transformMatrix", "Transformation Matrix", tgt::mat4::identity, tgt::mat4(-2000.0), tgt::mat4(2000.0), Processor::VALID)
    , autoCompute_("autoCompute", "Auto compute", false)
    , initialize_("initialize", "Initialize Transformation", Processor::VALID)
{
    addPort(staticPort_);
    addPort(movingPort_);

    addProperty(transformMatrix_);
    addProperty(autoCompute_);
    addProperty(initialize_);
    initialize_.onChange(MemberFunctionCallback<RegistrationInitializer>(this, &RegistrationInitializer::update));
}

RegistrationInitializer::~RegistrationInitializer() {}

Processor* RegistrationInitializer::create() const {
    return new RegistrationInitializer();
}

void RegistrationInitializer::process() {
    if(autoCompute_.get()) {
        update();
    }
}

// generates a transformation from a coordinate system centered in the volume to world space
mat4 getCenteredTransform(const VolumeBase* vol) {
    tgt::ivec3 dims = vol->getDimensions();
    vec3 c = vec3(dims) * vec3(0.5f);
    vec3 e1 = c; e1.x = static_cast<float>(dims.x);
    vec3 e2 = c; e2.y = static_cast<float>(dims.y);
    vec3 e3 = c; e3.z = static_cast<float>(dims.z);

    c = vol->getVoxelToWorldMatrix() * c;
    e1 = vol->getVoxelToWorldMatrix() * e1; e1 = e1 - c; e1 = normalize(e1);
    e2 = vol->getVoxelToWorldMatrix() * e2; e2 = e2 - c; e2 = normalize(e2);
    e3 = vol->getVoxelToWorldMatrix() * e3; e3 = e3 - c; e3 = normalize(e3);

    return mat4(e1.x, e2.x, e3.x, -c.x,
                e1.y, e2.y, e3.y, -c.y,
                e1.z, e2.z, e3.z, -c.z,
                0.0f, 0.0f, 0.0f, 1.0f);
}

void RegistrationInitializer::update() {
    const VolumeBase* stat = staticPort_.getData();
    const VolumeBase* move = movingPort_.getData();

    if(!stat || !move)
        return;

    mat4 volToWorldStatic = getCenteredTransform(stat);
    mat4 worldToVolStatic = mat4::identity;
    volToWorldStatic.invert(worldToVolStatic);

    mat4 volToWorldMoving = getCenteredTransform(move);

    transformMatrix_.set(worldToVolStatic * volToWorldMoving);
}

} // namespace
