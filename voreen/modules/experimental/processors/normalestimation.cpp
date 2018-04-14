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

#include "normalestimation.h"

#include "voreen/core/properties/cameraproperty.h"
#include "voreen/core/utils/voreenpainter.h"

#include "tgt/assert.h"
#include "tgt/glmath.h"
#include "tgt/vector.h"
#include "tgt/textureunit.h"

using tgt::vec3;
using tgt::TextureUnit;

namespace voreen {

NormalEstimation::NormalEstimation()
    : ImageProcessor("pp_normalestimation")
    , camera_("camera", "Camera", tgt::Camera(vec3(0.f, 0.f, 3.5f), vec3(0.f, 0.f, 0.f), vec3(0.f, 1.f, 0.f)))
    , inport_(Port::INPORT, "image.inport")
    , outport_(Port::OUTPORT, "image.outport", "image.outport", true)
{
    addProperty(camera_);

    addPort(inport_);
    addPort(outport_);
}

Processor* NormalEstimation::create() const {
    return new NormalEstimation();
}

void NormalEstimation::process() {
    outport_.activateTarget();

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    TextureUnit colorUnit;
    inport_.bindColorTexture(colorUnit.getEnum());

    // initialize shader
    program_->activate();
    tgt::Camera cam = camera_.get();
    setGlobalShaderParameters(program_, &cam);
    program_->setUniform("firstHitPointsNormalEstimation_", colorUnit.getUnitNumber());

    renderQuad();

    program_->deactivate();
    LGL_ERROR;
}

} // namespace
