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

#include "voreen/core/processors/volumerenderer.h"
#include "voreen/core/utils/glsl.h"

#include "tgt/gpucapabilities.h"
#include "tgt/textureunit.h"

using tgt::vec3;
using tgt::vec4;
using tgt::Color;
using tgt::TextureUnit;

namespace voreen {

const std::string VolumeRenderer::loggerCat_("voreen.VolumeRenderer");

VolumeRenderer::VolumeRenderer()
    : RenderProcessor()
    , lightPosition_("lightPosition", "Light Source Position", tgt::vec4(2.3f, 1.5f, 1.5f, 1.f),
                     tgt::vec4(-10000.f), tgt::vec4(10000.f))
    , lightAmbient_("lightAmbient", "Ambient Light", tgt::vec4(0.4f, 0.4f, 0.4f, 1.f))
    , lightDiffuse_("lightDiffuse", "Diffuse Light", tgt::vec4(0.8f, 0.8f, 0.8f, 1.f))
    , lightSpecular_("lightSpecular", "Specular Light", tgt::vec4(0.6f, 0.6f, 0.6f, 1.f))
    , lightAttenuation_("lightAttenuation", "Attenuation", tgt::vec3(1.f, 0.f, 0.f))
    , applyLightAttenuation_("applyLightAttenuation", "Apply Attenuation", false, Processor::INVALID_PROGRAM)
    , materialAmbient_("materialAmbient", "Ambient material color", tgt::vec4(1.f, 1.f, 1.f, 1.f))
    , materialDiffuse_("materialDiffuse", "Diffuse material color", tgt::vec4(1.f, 1.f, 1.f, 1.f))
    , materialSpecular_("materialSpecular", "Specular material color", tgt::vec4(1.f, 1.f, 1.f, 1.f))
    , materialShininess_("materialShininess", "Shininess", 60.f, 0.1f, 128.f)
{}

void VolumeRenderer::setGlobalShaderParameters(tgt::Shader* shader, const tgt::Camera* camera, tgt::ivec2 screenDim) {
    RenderProcessor::setGlobalShaderParameters(shader, camera, screenDim);

    shader->setIgnoreUniformLocationError(true);

    // light source position in world coordinates
    shader->setUniform("lightSource_.position_", lightPosition_.get().xyz());
    shader->setUniform("lightSource_.attenuation_", lightAttenuation_.get());
    shader->setUniform("lightSource_.ambientColor_", lightAmbient_.get().xyz());
    shader->setUniform("lightSource_.diffuseColor_", lightDiffuse_.get().xyz());
    shader->setUniform("lightSource_.specularColor_", lightSpecular_.get().xyz());

    shader->setUniform("shininess_", materialShininess_.get());

    shader->setIgnoreUniformLocationError(false);
}

bool VolumeRenderer::bindVolumes(tgt::Shader* shader, const std::vector<VolumeStruct>& volumes,
                                 const tgt::Camera* camera, const tgt::vec4& lightPosition) {
    shader->setIgnoreUniformLocationError(true);

    bool success = true;
    for (size_t i=0; i < volumes.size(); ++i) {
        const VolumeStruct& volumeStruct = volumes[i];
        const VolumeGL* volumeGL = volumeStruct.volume_->getRepresentation<VolumeGL>();
        if (!volumeGL || !volumeGL->getTexture()) {
            LERROR("OpenGL volume not available");
            success = false;
            continue;
        }

        const TextureUnit* texUnit = volumeStruct.texUnit_;
        if (!texUnit) {
            LERROR("No texture unit while binding volumes");
            success = false;
            continue;
        }

        success &= GLSL::bindVolumeTexture(volumeStruct.volume_, texUnit, volumeStruct.filterMode_, volumeStruct.wrapMode_, volumeStruct.borderColor_);

        // set volume meta-data
        GLSL::setUniform(shader, volumeStruct.volumeIdentifier_, volumeStruct.volumeStructIdentifier_, volumeStruct.volume_, texUnit, camera, lightPosition);

        LGL_ERROR;
    }

    shader->setIgnoreUniformLocationError(false);
    LGL_ERROR;

    return success;
}

VolumeRenderer::VolumeStruct::VolumeStruct()
    : volume_(0)
{}

VolumeRenderer::VolumeStruct::VolumeStruct(const VolumeBase* volume, const tgt::TextureUnit* texUnit,
                                           const std::string& volumeIdentifier,
                                           const std::string& volumeStructIdentifier,
                                           GLint wrapMode,
                                           tgt::vec4 borderColor, GLint filterMode)
                                           : volume_(volume),
                                           texUnit_(texUnit),
                                           wrapMode_(wrapMode),
                                           borderColor_(borderColor),
                                           filterMode_(filterMode),
                                           volumeIdentifier_(volumeIdentifier),
                                           volumeStructIdentifier_(volumeStructIdentifier)
{}

} // namespace voreen
