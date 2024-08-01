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

#include "idraycaster.h"
#include "voreen/core/utils/voreenpainter.h"
#include "voreen/core/utils/voreenqualitymode.h"

#include "voreen/core/ports/conditions/portconditionvolumetype.h"
#include "voreen/core/properties/cameraproperty.h"

#include "tgt/textureunit.h"

using tgt::vec3;
using tgt::TextureUnit;

namespace voreen {

IDRaycaster::IDRaycaster()
    : VolumeRaycaster()
    , firstHitpointsPort_(Port::INPORT, "image.firsthitpoints", "First-Hit-Points Input")
    , raycastPrg_(0)
    , camera_("camera", "Camera", tgt::Camera(vec3(0.f, 0.f, 3.5f), vec3(0.f, 0.f, 0.f), vec3(0.f, 1.f, 0.f)))
    , penetrationDepth_("penetrationDepth", "Penetration Depth", 0.05f, 0.f, 0.5f)
{
    addPort(firstHitpointsPort_);

    addProperty(penetrationDepth_ );
    addProperty(camera_);
}

Processor* IDRaycaster::create() const {
    return new IDRaycaster();
}

void IDRaycaster::initialize() {
    VolumeRaycaster::initialize();
    raycastPrg_ = ShdrMgr.loadSeparate("passthrough.vert", "rc_id.frag", generateHeader(), false);
}

void IDRaycaster::deinitialize() {
    ShdrMgr.dispose(raycastPrg_);
    raycastPrg_ = 0;

    VolumeRaycaster::deinitialize();
}

/**
 * Compile and link the shader program
 */
void IDRaycaster::compile() {
    // evaluate state and configure shader appropriately
    raycastPrg_->setHeaders(generateHeader());
    raycastPrg_->rebuild();
}

/**
 * Performs the raycasting.
 *
 * Initialize two texture units with the entry and exit params and renders
 * a screen aligned quad.
 */
void IDRaycaster::process() {

    if (!volumeInport_.isReady())
        return;

    if(getInvalidationLevel() >= Processor::INVALID_PROGRAM)
        compile();
    LGL_ERROR;

    outport1_.activateTarget();
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // don't render when in interaction mode
    if (QualityMode.isInteractionMode()) {
        outport1_.deactivateTarget();
        TextureUnit::setZeroUnit();
        return;
    }

    TextureUnit entryUnit, entryDepthUnit, exitUnit, exitDepthUnit, firstHitPointUnit;
    // bind entry params
    entryPort_.bindTextures(entryUnit.getEnum(), entryDepthUnit.getEnum());

    // bind first hit points
    firstHitpointsPort_.bindColorTexture(firstHitPointUnit.getEnum());

    // bind exit params
    exitPort_.bindTextures(exitUnit.getEnum(), exitDepthUnit.getEnum());

    // vector containing the volumes to bind
    std::vector<VolumeStruct> volumes;

    TextureUnit volUnit;
    volumes.push_back(VolumeStruct(
        volumeInport_.getData(),
        &volUnit,
        "segmentation_","segmentationParameters_",
        GL_CLAMP_TO_EDGE,
        tgt::vec4(0.f),
        GL_NEAREST)
    );

    // initialize shader
    raycastPrg_->activate();
    tgt::Camera cam = camera_.get();
    setGlobalShaderParameters(raycastPrg_, &cam);
    bindVolumes(raycastPrg_, volumes, &cam, lightPosition_.get());
    raycastPrg_->setUniform("entryPoints_", entryUnit.getUnitNumber());
    raycastPrg_->setUniform("entryPointsDepth_", entryDepthUnit.getUnitNumber());
    entryPort_.setTextureParameters(raycastPrg_, "entryParameters_");
    raycastPrg_->setUniform("firstHitPoints_", firstHitPointUnit.getUnitNumber());
    firstHitpointsPort_.setTextureParameters(raycastPrg_, "firstHitParameters_");
    raycastPrg_->setUniform("exitPoints_", exitUnit.getUnitNumber());
    raycastPrg_->setUniform("exitPointsDepth_", exitDepthUnit.getUnitNumber());
    entryPort_.setTextureParameters(raycastPrg_, "exitParameters_");

    raycastPrg_->setUniform("penetrationDepth_", penetrationDepth_.get());
    renderQuad();

    raycastPrg_->deactivate();
    outport1_.deactivateTarget();
    TextureUnit::setZeroUnit();
    LGL_ERROR;
}


} // namespace voreen
