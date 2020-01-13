/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2020 University of Muenster, Germany,                        *
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

#ifndef VRN_GLSL_H
#define VRN_GLSL_H

#include "voreen/core/voreencoreapi.h"
#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/datastructures/volume/slice/slicehelper.h"
#include "voreen/core/properties/optionproperty.h"
#include "tgt/shadermanager.h"
#include "tgt/camera.h"
#include "tgt/textureunit.h"

namespace voreen {

    /**
     * Utility class for glsl utility functions.
     */
class VRN_CORE_API GLSL {
public:

static void setUniform(tgt::Shader* shader, const std::string& volumeUniform, const std::string& structUniform, const VolumeBase* vh, const tgt::TextureUnit* texUnit, const tgt::Camera* camera = 0, const tgt::vec4& lightPosition = tgt::vec4(0.f));

static bool bindVolumeTexture(const VolumeBase* vh, const tgt::TextureUnit* texUnit, GLint filterMode = GL_LINEAR, GLint wrapMode = GL_CLAMP_TO_EDGE, tgt::vec4 borderColor = tgt::vec4(0.f));

static void setUniform(tgt::Shader* shader, const std::string& textureUniform, const std::string& structUniform, const SliceTexture* sl, const tgt::TextureUnit* texUnit);

static bool bindSliceTexture(const SliceTexture* slice, const tgt::TextureUnit* texUnit, GLint filterMode = GL_LINEAR, GLint wrapMode = GL_CLAMP_TO_EDGE, tgt::vec4 borderColor = tgt::vec4(0.f));

static tgt::Texture* createArrayTexture(const std::vector<std::string>& texNames, tgt::Texture::Filter filterMode = tgt::Texture::LINEAR,
                                              tgt::Texture::Wrapping wrapMode = tgt::Texture::CLAMP_TO_EDGE, tgt::vec4 borderColor = tgt::vec4(0.f));

static std::string generateStandardShaderHeader(const tgt::GpuCapabilities::GlVersion* version = 0);

// Helper functions for shading mode:

/// Fill a given property with available shading modes
static void fillShadingModesProperty(StringOptionProperty& shadeMode);

/// Generate a shader define based on the chosen shading mode
static std::string getShaderDefine(std::string shadeMode, std::string functionName, std::string n = "n", std::string pos = "pos", std::string lPos = "lPos", std::string cPos = "cPos", std::string ka = "ka", std::string kd = "kd", std::string ks = "ks");

}; // end class

}  // namespace voreen

#endif
