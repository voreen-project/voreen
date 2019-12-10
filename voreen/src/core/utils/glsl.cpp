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

#include "voreen/core/utils/glsl.h"

#include "tgt/framebufferobject.h"
#include "tgt/immediatemode/immediatemode.h"
#include "voreen/core/datastructures/volume/volumegl.h"
#include "tgt/texture2darray.h"

namespace voreen {

using tgt::vec3;
using std::string;

void GLSL::setUniform(tgt::Shader* shader, const std::string& volumeUniform, const std::string& structUniform, const VolumeBase* vh, const tgt::TextureUnit* texUnit, const tgt::Camera* camera, const tgt::vec4& lightPosition) {
    if(texUnit)
        shader->setUniform(volumeUniform, texUnit->getUnitNumber());

    // volume size, i.e. dimensions of the proxy geometry in world coordinates
    shader->setUniform(structUniform + ".datasetDimensions_", tgt::vec3(vh->getDimensions()));
    shader->setUniform(structUniform + ".datasetDimensionsRCP_", vec3(1.f) / tgt::vec3(vh->getDimensions()));

    // volume spacing, i.e. voxel size
    shader->setUniform(structUniform + ".datasetSpacing_", vh->getSpacing());
    shader->setUniform(structUniform + ".datasetSpacingRCP_", vec3(1.f) / vh->getSpacing());

    // volume's size in its physical coordinates
    shader->setUniform(structUniform + ".volumeCubeSize_", vh->getCubeSize());
    shader->setUniform(structUniform + ".volumeCubeSizeRCP_", vec3(1.f) / vh->getCubeSize());

    shader->setUniform(structUniform + ".volumeOffset_", vh->getOffset());

    shader->setUniform(structUniform + ".numChannels_", static_cast<GLint>(vh->getNumChannels()));

    // volume's transformation matrix
    shader->setUniform(structUniform + ".physicalToWorldMatrix_", vh->getPhysicalToWorldMatrix());

    tgt::mat4 invTm = vh->getWorldToPhysicalMatrix();
    shader->setUniform(structUniform + ".worldToPhysicalMatrix_", invTm);

    shader->setUniform(structUniform + ".worldToTextureMatrix_", vh->getWorldToTextureMatrix());
    shader->setUniform(structUniform + ".textureToWorldMatrix_", vh->getTextureToWorldMatrix());

    // camera position in volume object coords
    if (camera)
        shader->setUniform(structUniform + ".cameraPositionPhysical_", invTm*camera->getPositionWithOffsets());

    // light position in volume object coords
    shader->setUniform(structUniform + ".lightPositionPhysical_", (invTm*lightPosition).xyz());

    LGL_ERROR;

    // bit depth of the volume, TODO: check correctness of use for multi-channel data sets
    shader->setUniform(structUniform + ".bitDepth_", (GLint)(vh->getBytesPerVoxel() * 8));

    // construct shader real-world mapping by combining volume rwm and pixel transfer mapping
    RealWorldMapping rwm = vh->getRealWorldMapping();
    shader->setUniform(structUniform + ".rwmScale_", rwm.getScale());
    shader->setUniform(structUniform + ".rwmOffset_", rwm.getOffset());
}

bool GLSL::bindVolumeTexture(const VolumeBase* vh, const tgt::TextureUnit* texUnit, GLint filterMode, GLint wrapMode, tgt::vec4 borderColor) {
    const VolumeGL* volumeGL = vh->getRepresentation<VolumeGL>();
    if (!volumeGL || !volumeGL->getTexture()) {
        LWARNINGC("voreen.glsl", "No volume texture while binding volumes");
        return false;
    }
    const VolumeTexture* volumeTex = volumeGL->getTexture();

    texUnit->activate();

    volumeTex->bind();

    // texture filtering
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, filterMode);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, filterMode);
    LGL_ERROR;

    // texture wrapping
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, wrapMode);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, wrapMode);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, wrapMode);
    glTexParameterfv(GL_TEXTURE_3D, GL_TEXTURE_BORDER_COLOR, static_cast<tgt::Vector4<GLfloat> >(borderColor).elem);
    LGL_ERROR;

    texUnit->setZeroUnit();

    return true;
}

void GLSL::setUniform(tgt::Shader* shader, const std::string& textureUniform, const std::string& structUniform, const SliceTexture* slice, const tgt::TextureUnit* texUnit) {
    tgtAssert(shader, "null pointer passed");
    tgtAssert(slice, "null pointer passed");
    tgtAssert(texUnit, "null pointer passed");

    // check state of slice shader
    if (!shader->isActivated()) {
        LERRORC("voreen.glsl", "bindSliceTexture() slice shader not activated");
        return;
    }

    // pass tex unit
    if (texUnit)
        shader->setUniform(textureUniform, texUnit->getUnitNumber());
    LGL_ERROR;

    // determine real-world-mapping
    RealWorldMapping rwm = slice->getRealWorldMapping();

    // pass TextureParameters struct values to shader
    //shader->setIgnoreUniformLocationError(true);
    shader->setUniform(structUniform + ".dimensions_", tgt::vec2(slice->getSliceDimensions()));
    shader->setUniform(structUniform + ".dimensionsRCP_", tgt::vec2(1.0f) / tgt::vec2(slice->getSliceDimensions()));
    shader->setUniform(structUniform + ".matrix_", tgt::mat4::identity);
    shader->setUniform(structUniform + ".numChannels_", (int)slice->getNumChannels());
    shader->setUniform(structUniform + ".rwmScale_", rwm.getScale());
    shader->setUniform(structUniform + ".rwmOffset_", rwm.getOffset());
    //shader->setIgnoreUniformLocationError(false);
    LGL_ERROR;
}

bool GLSL::bindSliceTexture(const SliceTexture* slice, const tgt::TextureUnit* texUnit,
    GLint filterMode /*= GL_LINEAR*/, GLint wrapMode /*= GL_CLAMP_TO_EDGE*/, tgt::vec4 borderColor /*= tgt::vec4(0.f)*/)
{
    tgtAssert(slice, "null pointer passed");
    tgtAssert(texUnit, "null pointer passed");

    /// bind slicetexture to unit
    texUnit->activate();
    slice->bind();
    LGL_ERROR;

    // set texture filtering
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, filterMode);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, filterMode);
    LGL_ERROR;

    // set texture wrapping
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, wrapMode);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, wrapMode);
    glTexParameterfv(GL_TEXTURE_2D, GL_TEXTURE_BORDER_COLOR, borderColor.elem);
    LGL_ERROR;

    texUnit->setZeroUnit();

    return true;
}

tgt::Texture* GLSL::createArrayTexture(const std::vector<std::string>& texNames, tgt::Texture::Filter filterMode, tgt::Texture::Wrapping wrapMode, tgt::vec4 borderColor) {
    LINFOC("glsl", "Loading textures...");
    std::vector<tgt::Texture*> singleTexes;
    for(size_t i = 0; i < texNames.size(); i++) {
        if(tgt::FileSystem::fileExists(texNames.at(i)))
            singleTexes.push_back(TexMgr.load(texNames.at(i), filterMode, false, false, true));
        else {
            LERRORC("glsl", "Could not find texture " << texNames.at(i));
            for(size_t i = 0; i < singleTexes.size(); i++)
                TexMgr.dispose(singleTexes.at(i));
            return 0;
        }
    }

    std::string header = generateStandardShaderHeader();
    header += "#define NO_DEPTH_TEX\n";
    tgt::Shader* textureProgram = ShdrMgr.loadSeparate("passthrough.vert", "fullscreenquad.geom", "copyimage.frag", header, false);
    if(!textureProgram) {
        LERRORC("glsl", "Failed to load texture to texture array shader");
        return 0;
    }

    GLuint vao;
    glGenVertexArrays(1, &vao);
    glBindVertexArray(vao);
    GLuint vbo;
    glGenBuffers(1, &vbo);
    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    tgt::vec3 tmp(0.f);
    glBufferData(GL_ARRAY_BUFFER, sizeof(tgt::vec3), tmp.elem, GL_STATIC_DRAW);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(tgt::vec3), 0);

    tgt::FramebufferObject fbo;
    tgt::svec2 maxSize = tgt::svec2(0, 0);

    for(size_t i = 0; i < singleTexes.size(); i++) {
        tgt::svec2 curSize = singleTexes.at(i)->getDimensions().xy();
        if(curSize.x > maxSize.x)
            maxSize.x = curSize.x;
        if(curSize.y > maxSize.y)
            maxSize.y = curSize.y;
    }

    tgt::Texture* firstTex = singleTexes.at(0);
    tgt::Texture2DArray* arrayTex = new tgt::Texture2DArray(tgt::svec3(maxSize, singleTexes.size()),firstTex->getGLFormat(),firstTex->getGLInternalFormat(),
                                              firstTex->getGLDataType(),filterMode,wrapMode);
    arrayTex->uploadTexture();
    LGL_ERROR;

    MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
    MatStack.pushMatrix();
    MatStack.loadIdentity();

    MatStack.matrixMode(tgt::MatrixStack::PROJECTION);
    MatStack.pushMatrix();
    MatStack.loadIdentity();

    GLint viewport[4];
    glGetIntegerv(GL_VIEWPORT, viewport);
    glViewport(0, 0, maxSize.x, maxSize.y);

    textureProgram->activate();
    fbo.activate();
    LGL_ERROR;

    glDisable(GL_DEPTH_TEST);
    glDisable(GL_BLEND);

    GLenum colAtt = GL_COLOR_ATTACHMENT0;
    glDrawBuffers(1, &colAtt);

    // pass texture parameters to the shader
    textureProgram->setUniform("screenDimRCP_", 1.f / tgt::vec2(maxSize));
    textureProgram->setUniform("texParams_.matrix_", tgt::mat4::identity);

    tgt::TextureUnit texUnit;
    texUnit.activate();
    textureProgram->setUniform("colorTex_", texUnit.getUnitNumber());

    for(size_t curLayer = 0; curLayer < singleTexes.size(); curLayer++) {
        singleTexes.at(curLayer)->bind();
        fbo.attachTexture(arrayTex, colAtt, 0, (int) curLayer);
        LGL_ERROR;
        glDrawArrays(GL_POINTS, 0, 1);
        LGL_ERROR;
    }

    glEnable(GL_DEPTH_TEST);
    fbo.deactivate();

    glViewport(viewport[0], viewport[1], viewport[2], viewport[3]);
    MatStack.matrixMode(tgt::MatrixStack::PROJECTION);
    MatStack.popMatrix();
    MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
    MatStack.popMatrix();

    textureProgram->deactivate();
    tgt::TextureUnit::setZeroUnit();
    ShdrMgr.dispose(textureProgram);
    for(size_t i = 0; i < singleTexes.size(); i++)
        TexMgr.dispose(singleTexes.at(i));

    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindVertexArray(0);
    glDeleteBuffers(1,&vbo);
    glDeleteVertexArrays(1,&vao);

    LGL_ERROR;
    return arrayTex;
}

std::string GLSL::generateStandardShaderHeader(const tgt::GpuCapabilities::GlVersion* version) {
    if (!tgt::GpuCapabilities::isInited()) {
        return "";
    }

    using tgt::GpuCapabilities;

    tgt::GpuCapabilities::GlVersion useVersion;

    //use supplied version if available, else use highest available.
    //if no version is supplied, use up tp 1.30 as default.
    if (version && GpuCaps.getShaderVersion() >= *version)
        useVersion = *version;
    else if(GpuCaps.getShaderVersion() > GpuCapabilities::GlVersion::SHADER_VERSION_410)
        useVersion = GpuCapabilities::GlVersion::SHADER_VERSION_410;
    else
        useVersion = GpuCaps.getShaderVersion();

    std::stringstream versionHeader;
    versionHeader << (useVersion.major)() << (useVersion.minor)();

    std::string header = "#version " + versionHeader.str();

    if(header.length() < 12)
        header += "0";

    //Run in compability mode to use deprecated functionality
#ifdef VRN_OPENGL_COMPATIBILITY_PROFILE
    if (useVersion >= GpuCapabilities::GlVersion::SHADER_VERSION_150)
        header += " compatibility";
    else if(useVersion >= GpuCapabilities::GlVersion::SHADER_VERSION_140)
        header += "\n#extension GL_ARB_compatibility : enable";
#else
    header += " core";
#endif

    header += "\n";

    if (useVersion >= GpuCapabilities::GlVersion::SHADER_VERSION_450)
        header += "#define GLSL_VERSION_450\n";
    if (useVersion >= GpuCapabilities::GlVersion::SHADER_VERSION_440)
        header += "#define GLSL_VERSION_440\n";
    if (useVersion >= GpuCapabilities::GlVersion::SHADER_VERSION_430)
        header += "#define GLSL_VERSION_430\n";
    if (useVersion >= GpuCapabilities::GlVersion::SHADER_VERSION_420)
        header += "#define GLSL_VERSION_420\n";
    if (useVersion >= GpuCapabilities::GlVersion::SHADER_VERSION_410)
        header += "#define GLSL_VERSION_410\n";
    if (useVersion >= GpuCapabilities::GlVersion::SHADER_VERSION_400)
        header += "#define GLSL_VERSION_400\n";
    if (useVersion >= GpuCapabilities::GlVersion::SHADER_VERSION_330)
        header += "#define GLSL_VERSION_330\n";
    if (useVersion >= GpuCapabilities::GlVersion::SHADER_VERSION_150) {
        header += "#define GLSL_VERSION_150\n";
        header += "#extension GL_ARB_explicit_attrib_location : enable\n";
    }
    if (useVersion >= GpuCapabilities::GlVersion::SHADER_VERSION_140)
        header += "#define GLSL_VERSION_140\n";
    if (useVersion >= GpuCapabilities::GlVersion::SHADER_VERSION_130) {
        header += "#define GLSL_VERSION_130\n";
        header += "precision highp float;\n";
    }
    if (useVersion >= GpuCapabilities::GlVersion::SHADER_VERSION_120)
        header += "#define GLSL_VERSION_120\n";
    if (useVersion >= GpuCapabilities::GlVersion::SHADER_VERSION_110)
        header += "#define GLSL_VERSION_110\n";

    bool couldReadMaxLoopCount = false;
#ifdef VRN_OPENGL_COMPATIBILITY_PROFILE
    if (GLEW_NV_fragment_program2) {
        GLint i = -1;
        glGetProgramivARB(GL_FRAGMENT_PROGRAM_ARB, GL_MAX_PROGRAM_LOOP_COUNT_NV, &i);
        if (i > 0) {
            std::ostringstream o;
            o << i;
            header += "#define VRN_MAX_PROGRAM_LOOP_COUNT " + o.str() + "\n";
            couldReadMaxLoopCount = true;
        }
    }
#endif
    if(!couldReadMaxLoopCount) {
        header += "#define VRN_MAX_PROGRAM_LOOP_COUNT 256*256\n";
    }

    //
    // add some defines needed for workarounds in the shader code
    //
// COREPORT_TODO
#ifdef VRN_OPENGL_COMPATIBILITY_PROFILE
    if (GLEW_ARB_draw_buffers)
#endif
        header += "#define VRN_GLEW_ARB_draw_buffers\n";

    #ifdef __APPLE__
        header += "#define VRN_OS_APPLE\n";
        if (GpuCaps.getVendor() == GpuCaps.GPU_VENDOR_ATI)
            header += "#define VRN_VENDOR_ATI\n";
    #endif

    if (GpuCaps.getShaderVersion() >= GpuCapabilities::GlVersion::SHADER_VERSION_130) {
        // define output for single render target
        header += "//$ @name = \"outport0\"\n";
        header += "out vec4 FragData0;\n";
    }
    else {
        header += "#define FragData0 gl_FragData[0]\n";
    }

    header += IMode.generateMatstackUniformHeader();

    return header;
}

void GLSL::fillShadingModesProperty(StringOptionProperty& shadeMode) {
    shadeMode.addOption("none",                   "none"                   );
    shadeMode.addOption("phong-diffuse",          "Phong (Diffuse)"        );
    shadeMode.addOption("phong-specular",         "Phong (Specular)"       );
    shadeMode.addOption("phong-diffuse-ambient",  "Phong (Diffuse+Amb.)"   );
    shadeMode.addOption("phong-diffuse-specular", "Phong (Diffuse+Spec.)"  );
    shadeMode.addOption("phong",                  "Phong (Full)"           );
    shadeMode.addOption("toon",                   "Toon"                   );
    shadeMode.addOption("cook-torrance",          "Cook-Torrance"          );
    shadeMode.addOption("oren-nayar",             "Oren-Nayar"             );
    shadeMode.addOption("ward",                   "Ward (Isotropic)"       );
    shadeMode.select("phong");
}

string GLSL::getShaderDefine(string shadeMode, string functionName, string n, string pos, string lPos, string cPos, string ka, string kd, string ks) {
    string headerSource = "#define " + functionName + "(" + n + ", " + pos + ", " + lPos + ", " + cPos + ", " + ka + ", " + kd + ", " + ks + ") ";
    if (shadeMode == "none")
        headerSource += "" + ka + ";\n";
    else if (shadeMode == "phong-diffuse")
        headerSource += "phongShadingD(" + n + ", " + pos + ", " + lPos + ", " + cPos + ", " + kd + ");\n";
    else if (shadeMode == "phong-specular")
        headerSource += "phongShadingS(" + n + ", " + pos + ", " + lPos + ", " + cPos + ", " + ks + ");\n";
    else if (shadeMode == "phong-diffuse-ambient")
        headerSource += "phongShadingDA(" + n + ", " + pos + ", " + lPos + ", " + cPos + ", " + kd + ", " + ka + ");\n";
    else if (shadeMode == "phong-diffuse-specular")
        headerSource += "phongShadingDS(" + n + ", " + pos + ", " + lPos + ", " + cPos + ", " + kd + ", " + ks + ");\n";
    else if (shadeMode == "phong")
        headerSource += "phongShading(" + n + ", " + pos + ", " + lPos + ", " + cPos + ", " + ka + ", " + kd + ", " + ks + ");\n";
    else if (shadeMode == "toon")
        headerSource += "toonShading(" + n + ", " + pos + ", " + lPos + ", " + cPos + ", " + kd + ", 3);\n";
    else if (shadeMode == "cook-torrance")
        headerSource += "cookTorranceShading(" + n + ", " + pos + ", " + lPos + ", " + cPos + ", " + ka + ", " + kd + ", " + ks + ");\n";
    else if (shadeMode == "oren-nayar")
        headerSource += "orenNayarShading(" + n + ", " + pos + ", " + lPos + ", " + cPos + ", " + ka + ", " + kd + ");\n";
    else if (shadeMode == "lafortune")
        headerSource += "lafortuneShading(" + n + ", " + pos + ", " + lPos + ", " + cPos + ", " + ka + ", " + kd + ", " + ks + ");\n";
    else if (shadeMode == "ward")
        headerSource += "wardShading(" + n + ", " + pos + ", " + lPos + ", " + cPos + ", " + ka + ", " + kd + ", " + ks + ");\n";

    return headerSource;
}

} // namespace
