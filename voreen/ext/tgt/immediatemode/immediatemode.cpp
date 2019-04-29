/**********************************************************************
 *                                                                    *
 * tgt - Tiny Graphics Toolbox                                        *
 *                                                                    *
 * Copyright (C) 2005-2018 University of Muenster, Germany,           *
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

#include "immediatemode.h"

#include "tgt/tgt_gl.h"
#include "tgt/exception.h"
#include "tgt/matrixstack.h"
#include "tgt/shadermanager.h"

namespace tgt {

    const std::string ImmediateMode::loggerCat_("tgt.Intermediate.Mode");

    ImmediateMode::ImmediateMode() :
    started_(false),
    lightingEnabled_(false),
    texMode_(TEXNONE),
    materialColor_(MATCOLNONE),
    currentColor_(1,1,1,1),
    currentTexcoord_(0,0,0,0),
    currentNormal_(0,0,0),
    shaders_(),
    clipPlanes_()
    {
        // Initialize Shaders as not loaded
        std::fill(shaders_.begin(), shaders_.end(), nullptr);
        // Default Initialization for ClipPlanes
        std::fill(clipPlanes_.begin(), clipPlanes_.end(), ClipPlane());

        glGenBuffers(1, &vertexBuffer_);

        samplerNames_.push_back("none");
        samplerNames_.push_back("texture1D_");
        samplerNames_.push_back("texture2D_");
        samplerNames_.push_back("texture3D_");
    }

    ImmediateMode::~ImmediateMode() {
        glDeleteBuffers(1, &vertexBuffer_);
    }

    void ImmediateMode::setTextureMode(ImmediateMode::TexturMode texmode) {
        texMode_ = texmode;
    }

    ImmediateMode::TexturMode ImmediateMode::textureMode() const {
        return texMode_;
    }

    void ImmediateMode::enableLighting() {
        lightingEnabled_ = true;
    }

    void ImmediateMode::disableLighting() {
        lightingEnabled_ = false;
    }

    void ImmediateMode::setLightSource(const tgt::ImmediateMode::LightSource &lightSource) {
        lightSource_ = lightSource;
    }

    ImmediateMode::LightSource ImmediateMode::getLightSource() const {
        return lightSource_;
    }

    void ImmediateMode::setMaterial(const tgt::ImmediateMode::Material &material) {
        material_ = material;
    }

    ImmediateMode::Material ImmediateMode::getMaterial() const {
        return material_;
    }

    void ImmediateMode::setMaterialColor(MaterialColor matColor) {
        materialColor_ = matColor;
    }

    void ImmediateMode::begin(PolygonMode polygonMode) {
        if (started_) {
            throw Exception("Called ImmediateMode::begin multiple times.");
        }
        started_ = true;
        polygonMode_ = polygonMode;
    }

    void ImmediateMode::end() {
        if (!started_) {
            throw Exception("No corresponding begin call for ImmediateMode::end.");
        }
        started_ = false;

        if (vertices_.size() == 0) {
            return;
        }

        //TODO,FIXME: AMD HACK: REMOVE IT
        if(GpuCaps.getVendor() == GpuCapabilities::GPU_VENDOR_ATI) {
            //loading the standard workspace twice causes a freeze in glBufferData.
            //most likely caused by context problems. This behavior should be tested with new driver versions.
            glDeleteBuffers(1, &vertexBuffer_);
            glGenBuffers(1, &vertexBuffer_);
        }

        // a new vao every frame because vaos can't be shared
        // across contexts
        GLuint vao;
        glGenVertexArrays(1, &vao);
        glBindVertexArray(vao);
        glBindBuffer(GL_ARRAY_BUFFER, vertexBuffer_);

        // position
        glEnableVertexAttribArray(0);
        glVertexAttribPointer(0, 4, GL_FLOAT, false, sizeof(Vertex), (GLvoid*) 0);
        // color
        glEnableVertexAttribArray(1);
        glVertexAttribPointer(1, 4, GL_FLOAT, false, sizeof(Vertex), (GLvoid*) offsetof(Vertex, color));
        // normal
        glEnableVertexAttribArray(2);
        glVertexAttribPointer(2, 3, GL_FLOAT, false, sizeof(Vertex), (GLvoid*) offsetof(Vertex, normal));
        // texcoord
        glEnableVertexAttribArray(3);
        glVertexAttribPointer(3, 4, GL_FLOAT, false, sizeof(Vertex), (GLvoid*) offsetof(Vertex, texcoord));


        // this is the location where where OpenGL will read the vertice data from
        GLvoid* finalVertices = vertices_.data();
        GLsizei finalSize = (GLsizei)vertices_.size();

        // we need to convert quad geometry to triangle geometry
        std::vector<Vertex> convertedVertices;
        if (polygonMode_ == QUADS) {
            for (size_t i = 0; i < vertices_.size()-3; i+=4) {
                convertedVertices.push_back(vertices_[i+0]);
                convertedVertices.push_back(vertices_[i+1]);
                convertedVertices.push_back(vertices_[i+2]);

                convertedVertices.push_back(vertices_[i+2]);
                convertedVertices.push_back(vertices_[i+3]);
                convertedVertices.push_back(vertices_[i+0]);
            }

            polygonMode_ = TRIANGLES;

            finalVertices = convertedVertices.data();
            finalSize = (GLsizei)convertedVertices.size();
        } else if (polygonMode_ == FAKE_LINES || polygonMode_ == FAKE_LINE_STRIP) {
            float lineWidth = 1;
            glGetFloatv(GL_LINE_WIDTH, &lineWidth);

            tgt::vec4 viewport = tgt::vec4::zero;
            glGetFloatv(GL_VIEWPORT, viewport.elem);
            tgt::mat4 viewportMat = tgt::mat4::createTranslation(tgt::vec3(viewport.xy(),0)) // [0,w|h] -> [x|y,x+w|y+h]
                                  * tgt::mat4::createScale(tgt::vec3(viewport.zw(), 1))      // [0,1] -> [0,w|h]
                                  * tgt::mat4::createScale(tgt::vec3(0.5,0.5,1))             // [0,2] -> [0,1]
                                  * tgt::mat4::createTranslation(tgt::vec3(1,1,0));          // [-1,1] -> [0,2]

            tgt::mat4 toPixelMat = viewportMat * MatStack.getProjectionMatrix() * MatStack.getModelViewMatrix();
            tgt::mat4 fromPixelMat;
            bool inverted = toPixelMat.invert(fromPixelMat);
            tgtAssert(inverted, "Failed to invert toPixelMat");
            LGL_ERROR;

            auto offsetVertices2 = [&] (const Vertex& v1, const Vertex& v2) {
                    Vertex v11 = v1;
                    Vertex v12 = v1;
                    Vertex v21 = v2;
                    Vertex v22 = v2;

                    tgt::vec4 v1Pix = toPixelMat*v1.position;
                    tgt::vec4 v2Pix = toPixelMat*v2.position;

                    tgt::vec2 pixFromTo = v2Pix.xy() - v1Pix.xy();
                    if(pixFromTo != tgt::vec2::zero) {
                        tgt::vec2 pixTangent = tgt::normalize(pixFromTo);
                        tgt::vec2 pixOrthogonal(pixTangent.y, -pixTangent.x);
                        tgt::vec4 pixOffset(pixOrthogonal*(lineWidth*0.5f), 0, 0);

                        v11.position = fromPixelMat*(v1Pix+pixOffset);
                        v12.position = fromPixelMat*(v1Pix-pixOffset);
                        v21.position = fromPixelMat*(v2Pix+pixOffset);
                        v22.position = fromPixelMat*(v2Pix-pixOffset);
                    }

                    return std::make_tuple(v11, v12, v21, v22);
            };

            if (polygonMode_ == FAKE_LINES) {
                for (size_t i = 0; i < vertices_.size()-1; i+=2) {
                    auto t = offsetVertices2(vertices_[i+0], vertices_[i+1]);

                    convertedVertices.push_back(std::get<0>(t));
                    convertedVertices.push_back(std::get<2>(t));
                    convertedVertices.push_back(std::get<1>(t));

                    convertedVertices.push_back(std::get<1>(t));
                    convertedVertices.push_back(std::get<2>(t));
                    convertedVertices.push_back(std::get<3>(t));
                }
                polygonMode_ = TRIANGLES;
            } else if (polygonMode_ == FAKE_LINE_STRIP) {
                for (size_t i = 0; i < vertices_.size()-1; i+=1) {
                    auto t = offsetVertices2(vertices_[i+0], vertices_[i+1]);

                    convertedVertices.push_back(std::get<1>(t));
                    convertedVertices.push_back(std::get<0>(t));
                    convertedVertices.push_back(std::get<3>(t));
                    convertedVertices.push_back(std::get<2>(t));
                }
                polygonMode_ = TRIANGLE_STRIP;
            }

            finalVertices = convertedVertices.data();
            finalSize = (GLsizei)convertedVertices.size();
        }

        // load standardshader if none is set
        bool noShader = tgt::Shader::getCurrentProgram() == 0;

        // Only used if no shader is bound and clipping is enabled,
        // but needs to be declared here to have access before and
        // after rendering, to save restore GL_CLIP_DISTANCEX
        std::array<bool, MAX_CLIP_PLANES> savedGlClipDistanceStatus;

        tgt::Shader* shader = NULL;
        if (noShader) {
            ShaderConfiguration shaderConfig = 0;
            // determine which shader to use based on the current textur target
            switch (texMode_) {
                case TEXNONE:
                    shaderConfig = NO_TEX;
                    break;
                case TEX1D:
                    shaderConfig = TEX_1D;
                    break;
                case TEX2D:
                    shaderConfig = TEX_2D;
                    break;
                case TEX3D:
                    shaderConfig = TEX_3D;
                    break;
                default:
                    shaderConfig = NO_TEX;
            }
            bool clipping = clippingEnabled();
            if(clipping) {
                shaderConfig |= CLIPPING;
            }
            shader = getShader(shaderConfig);

            shader->activate();

            if(clipping) {
                for(unsigned int i = 0; i < MAX_CLIP_PLANES; ++i) {
                    GLuint clipDistanceID = GL_CLIP_DISTANCE0 + i;
                    savedGlClipDistanceStatus[i] = glIsEnabled(clipDistanceID);
                    if(getClipPlane(i).enabled) {
                        glEnable(clipDistanceID);
                    } else {
                        glDisable(clipDistanceID);
                    }
                }
                setClippingUniforms(*shader);
                LGL_ERROR;
            }

            shader->setUniform("lightingEnabled_", lightingEnabled_);
            setLightingUniforms(*shader);
            LGL_ERROR;
            if (texMode_ != TEXNONE) {
                shader->setUniform("sampler_", 0);
            }
        }

        setMatstackUniforms();

        //LGL_ERROR;
        glBufferData(GL_ARRAY_BUFFER, finalSize*sizeof(Vertex), finalVertices, GL_DYNAMIC_DRAW);
        glDrawArrays(static_cast<GLenum>(polygonMode_), 0, finalSize);
        glBufferData(GL_ARRAY_BUFFER, 0, 0, GL_DYNAMIC_DRAW);
        vertices_.clear();
        glBindVertexArray(0);
        glDeleteVertexArrays(1, &vao);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
        if (noShader) {
            shader->deactivate();

            if(clippingEnabled()) {
                for(unsigned int i = 0; i < MAX_CLIP_PLANES; ++i) {
                    GLuint clipDistanceID = GL_CLIP_DISTANCE0 + i;
                    if(savedGlClipDistanceStatus[i]) {
                        glEnable(clipDistanceID);
                    } else {
                        glDisable(clipDistanceID);
                    }
                }
            }
        }
    }

    void ImmediateMode::vertex(const vec2& vector) {
        uploadVertex(vec4(vector, 0, 1));
    }

    void ImmediateMode::vertex(const vec3& vector) {
        uploadVertex(vec4(vector, 1));
    }

    void ImmediateMode::vertex(const vec4& vector) {
        uploadVertex(vector);
    }

    void ImmediateMode::texcoord(float vector) {
        currentTexcoord_ = vec4(vector, 0, 0, 1);
    }

    void ImmediateMode::texcoord(const vec2& vector) {
        currentTexcoord_ = vec4(vector, 0, 1);
    }

    void ImmediateMode::texcoord(const vec3& vector) {
        currentTexcoord_ = vec4(vector, 1);
    }

    void ImmediateMode::texcoord(const vec4& vector) {
        currentTexcoord_ = vector;
    }

    void ImmediateMode::color(const vec3 &color) {
        currentColor_ = vec4(color, 1);
    }

    void ImmediateMode::color(const vec4 &color) {
        currentColor_ = color;
    }

    tgt::vec4 ImmediateMode::getCurrentColor() const {
        return currentColor_;
    }

    void ImmediateMode::normal(const tgt::vec3 &normal) {
        currentNormal_ = normal;
    }

    void ImmediateMode::enableClipPlane(unsigned int numClipPlane) {
        getClipPlane(numClipPlane).enabled = true;

    }

    void ImmediateMode::disableClipPlane(unsigned int numClipPlane) {
        getClipPlane(numClipPlane).enabled = false;
    }

    void ImmediateMode::setClipPlaneEquation(unsigned int numClipPlane, const tgt::plane& plane) {
        getClipPlane(numClipPlane).plane = plane;
    }

    void ImmediateMode::setModelSpaceClipPlaneEquation(unsigned int numClipPlane, const tgt::plane& plane) {
        tgt::mat4 inverseMV;
        tgt::vec4 modelSpaceEquation = plane.toVec4();
        bool invertible = MatStack.getModelViewMatrix().invert(inverseMV);
        tgtAssert(invertible, "Singular MV matrix.");
        tgt::vec4 eyeSpaceEquation = tgt::transpose(inverseMV)*modelSpaceEquation;
        setClipPlaneEquation(numClipPlane, tgt::plane(eyeSpaceEquation.xyz(), eyeSpaceEquation.w));
    }

    tgt::plane ImmediateMode::getClipPlaneEquation(unsigned int numClipPlane) const {
        return getClipPlane(numClipPlane).plane;
    }

    tgt::ImmediateMode::ClipPlane& ImmediateMode::getClipPlane(unsigned int numClipPlane) {
        tgtAssert(numClipPlane < MAX_CLIP_PLANES, "Invalid clip plane number.");

        return clipPlanes_[numClipPlane];
    }

    const tgt::ImmediateMode::ClipPlane& ImmediateMode::getClipPlane(unsigned int numClipPlane) const {
        tgtAssert(numClipPlane < MAX_CLIP_PLANES, "Invalid clip plane number.");

        return clipPlanes_[numClipPlane];
    }

    void ImmediateMode::uploadVertex(const vec4 &position) {
        if (!started_) {
            throw Exception("vertex can only be specified in a begin/end block");
        }
        Vertex vertex;
        vertex.position = position;
        vertex.texcoord = currentTexcoord_;
        vertex.color = currentColor_;
        vertex.normal = currentNormal_;
        vertices_.push_back(vertex);
    }

    void ImmediateMode::setMatstackUniforms() const {
        Shader* shader = ShdrMgr.getActiveShader();
        if(shader)
            setMatstackUniforms(shader);
    }

    void ImmediateMode::setMatstackUniforms(Shader *shader) const {
        shader->setIgnoreUniformLocationError(true);
        shader->setUniform("modelViewMatrixStack_", MatStack.getModelViewMatrix());
        shader->setUniform("projectionMatrixStack_", MatStack.getProjectionMatrix());
        shader->setUniform("modelViewProjectionMatrixStack_", MatStack.getProjectionMatrix() * MatStack.getModelViewMatrix());

        mat4 viewInverse;
        MatStack.getModelViewMatrix().invert(viewInverse);
        mat4 projInverse;
        MatStack.getProjectionMatrix().invert(projInverse);
        shader->setUniform("modelViewMatrixInverseStack_", viewInverse);
        shader->setUniform("projectionMatrixInverseStack_", projInverse);
        shader->setUniform("modelViewProjectionMatrixInverseStack_", viewInverse * projInverse);
        shader->setIgnoreUniformLocationError(false);
    }

    void ImmediateMode::setLightingUniforms(tgt::Shader& shader) const {
        shader.setIgnoreUniformLocationError(true);

        // set light source
        std::string prefix = "lightSource_.";
        shader.setUniform(prefix+"position_", lightSource_.position);
        shader.setUniform(prefix+"ambientColor_", lightSource_.ambientColor);
        shader.setUniform(prefix+"diffuseColor_", lightSource_.diffuseColor);
        shader.setUniform(prefix+"specularColor_", lightSource_.specularColor);
        shader.setUniform(prefix+"attenuation_", lightSource_.attenuation);

        // material
        prefix = "material_.";
        shader.setUniform(prefix+"ambientColor_", material_.ambientColor);
        shader.setUniform(prefix+"diffuseColor_", material_.diffuseColor);
        shader.setUniform(prefix+"specularColor_", material_.specularColor);
        shader.setUniform(prefix+"shininess_", material_.shininess);

        // materialColor
        shader.setUniform("materialColorAmbient_", (bool)(materialColor_ & MATCOLAMBIENT));
        shader.setUniform("materialColorDiffuse_", (bool)(materialColor_ & MATCOLDIFFUSE));
        shader.setUniform("materialColorSpecular_", (bool)(materialColor_ & MATCOLSPECULAR));

        shader.setIgnoreUniformLocationError(false);
    }

    void ImmediateMode::setClippingUniforms(tgt::Shader& shader) const {
        std::vector<tgt::vec4> clipPlaneEquations;
        for(const ClipPlane& clipPlane : clipPlanes_) {
            if(clipPlane.enabled) {
                clipPlaneEquations.push_back(clipPlane.plane.toVec4());
            } else {
                clipPlaneEquations.push_back(tgt::vec4::zero);
            }
        }
        tgtAssert(clipPlaneEquations.size() == clipPlanes_.size(), "Clip planes size mismatch");
        shader.setUniform("clip_planes", clipPlaneEquations.data(), clipPlaneEquations.size());
    }

    std::string ImmediateMode::generateMatstackUniformHeader() const {
        std::string header;

        header += "uniform mat4 projectionMatrixStack_;\n";
        header += "uniform mat4 modelViewMatrixStack_;\n";
        header += "uniform mat4 modelViewProjectionMatrixStack_;\n";
        header += "uniform mat4 projectionMatrixInverseStack_;\n";
        header += "uniform mat4 modelViewMatrixInverseStack_;\n";
        header += "uniform mat4 modelViewProjectionMatrixInverseStack_;\n";

        return header;
    }

    std::string ImmediateMode::generateHeader(ShaderConfiguration config) const {
        unsigned int texdim = config & 3;
        unsigned int clippingMask = static_cast<unsigned int>(CLIPPING);
        bool clippingEnabled = (config & clippingMask) == clippingMask;

        std::stringstream header;
        header << "#version 150 core\n";
        header << "#extension GL_ARB_explicit_attrib_location : enable\n";
        header << "#define TEXDIM " << texdim << "\n";
        if(clippingEnabled) {
            header << "#define CLIPPING_ENABLED\n";
            header << "#define NUM_CLIP_PLANES " << MAX_CLIP_PLANES << "\n";
        }
        return header.str() + generateMatstackUniformHeader();
    }

    tgt::Shader* ImmediateMode::getShader(ShaderConfiguration config) {
        tgtAssert(config < NUM_SHADER_CONFIGURATIONS, "Invalid ShaderConfiguration");
        tgt::Shader*& shader = shaders_[config];
        if(!shader) {
            shader = ShdrMgr.load("immediatemode/immediatemode", generateHeader(config), false);
            ShdrMgr.registerStaticShader(shader);
        }
        return shader;
    }

    bool ImmediateMode::clippingEnabled() const {
        for(ClipPlane plane : clipPlanes_) {
            if(plane.enabled) {
                return true;
            }
        }
        return false;
    }

} // namespace tgt
