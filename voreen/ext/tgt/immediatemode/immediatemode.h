/**********************************************************************
 *                                                                    *
 * tgt - Tiny Graphics Toolbox                                        *
 *                                                                    *
 * Copyright (C) 2005-2024 University of Muenster, Germany,           *
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

#ifndef TGT_IMMEDIATEMODE_H
#define TGT_IMMEDIATEMODE_H

#include "tgt/vector.h"
#include "tgt/plane.h"
#include "tgt/tgt_gl.h"

#include <vector>
#include <array>
#include <limits>
#include <boost/static_assert.hpp>

namespace tgt {

    class Shader;
    class ImmediateMode;
    #ifdef DLL_TEMPLATE_INST
        template class TGT_API Singleton<ImmediateMode>;
    #endif

    /**
     *  Singleton class helping you transition to OpenGL Core profile by providing
     *  immediatemode style rendering functionality.
     */
    class TGT_API ImmediateMode : public Singleton<ImmediateMode> {
    public:
        /**
         *  Defines the different primitives which can be used with begin.
         */
        enum PolygonMode {
            POINTS = GL_POINTS,
            LINES = GL_LINES,
            LINE_STRIP = GL_LINE_STRIP,
            LINE_LOOP = GL_LINE_LOOP,
            TRIANGLES = GL_TRIANGLES,
            TRIANGLE_STRIP = GL_TRIANGLE_STRIP,
            TRIANGLE_FAN = GL_TRIANGLE_FAN,
            QUADS = GL_QUADS,
            QUAD_STRIP = GL_TRIANGLE_STRIP, // it is basically the same
            POLYGON = GL_TRIANGLE_FAN,      // it is the same for convex polygons
            FAKE_LINES = 0xfeedbeee,        // Emulate lines with triangles
            FAKE_LINE_STRIP = 0xfeedbeef,   // Emulate line strip with triangle strip
        };

        /**
         *  Defines the texture targets which can be used when not using
         *  a custom shader.
         */
        enum TexturMode {
            TEXNONE,
            TEX1D = GL_TEXTURE_1D,
            TEX2D = GL_TEXTURE_2D,
            TEX3D = GL_TEXTURE_3D
        };

        /**
         *  Defines all properties of the light source when lighting is enabled.
         *  Default values are chosen in a way that no light is applied and attenuation is disabled.
         */
        struct LightSource {
            LightSource():
            position(0,0,1,0),
            ambientColor(0,0,0),
            diffuseColor(0,0,0),
            specularColor(0,0,0),
            attenuation(1,0,0) {}

            /**
             *  The position of the light source.
             *  Default value: (0,0,1,0)
             */
            vec4 position;

            /**
             *  Ambient color of the light source.
             *  Default value: (0,0,0)
             */
            vec3 ambientColor;

            /**
             *  Diffuse color of the light source.
             *  Default value: (0,0,0)
             */
            vec3 diffuseColor;

            /**
             *  Specular Color of the light source.
             *  Default value: (0,0,0)
             */
            vec3 specularColor;

            /**
             *  The attenuation.
             *  The components stand for constant, linear and quadric attenuation.
             *  (1,0,0) disables the attenuation.
             *  Default value: (0,0,0)
             */
            vec3 attenuation;
        };

        /**
         *  Defines the properties of the currently active material.
         *  Either these or the attribute/texture color values are used for lighting
         *  depending on the value set for materialColor.
         */
        struct Material {
            Material():
            ambientColor(1,1,1),
            diffuseColor(1,1,1),
            specularColor(1,1,1),
            shininess(5) {}

            /**
             *  Ambient color of the material.
             *  Default value: (1,1,1)
             */
            vec3 ambientColor;
            /**
             *  Diffuse color of the material.
             *  Default value: (1,1,1)
             */
            vec3 diffuseColor;
            /**
             *  Specular Color of the material.
             *  Default value: (1,1,1)
             */
            vec3 specularColor;

            /**
             *  Shininess of the material used for the specular term.
             *  Default value: 5
             */
            float shininess;
        };

        /**
         *  Bitfield to choose which material colors should been
         *  taken from the attribut color or texture value.
         *  Components which are not enabled are taken from the material.
         */
        enum MaterialColor {
            MATCOLNONE = 0,
            MATCOLAMBIENT = 2,
            MATCOLDIFFUSE = 4,
            MATCOLSPECULAR = 8
        };

        /**
         *  A Clip plane that defines the plane and whether it is enabled
         *  or not. Note that the plane is stored as in eye coordinates!
         */
        struct ClipPlane {
            ClipPlane():
            plane(0,0,0,0),
            enabled(false) {}
            /**
             * Plane (in eye coordinates)!
             */
            tgt::plane plane;
            /**
             * Wheter or not this plane is currently enabled
             */
            bool enabled;
        };

        /**
         * Maximum number of clip planes,
         * i.e. one can define and enable/disable
         * clip planes 0 to MAX_CLIP_PLANES-1.
         *
         * More or less equivalent to GL_MAX_CLIP_PLANES
         */
        static const unsigned int MAX_CLIP_PLANES = 6;

        /**
         *  Constructor.
         */
        ImmediateMode();
        /**
         *  Destructor.
         */
        ~ImmediateMode();

        /**
         *  Normalizes (unsigned) integer inputs.
         *
         *  @param iVal An integer value which will be normalized to floating point based on
         *              it's value range.
         *
         *  @return [0,1] for unsigned inputs and [-1,1] for signed.
         */
        template<typename T>
        float normalize(T integerValue) {
            BOOST_STATIC_ASSERT(std::numeric_limits<T>::is_integer);

            size_t maxVal = powf(2, sizeof(T) * 8);
            if (std::numeric_limits<T>::is_signed) {
                return (2.0f * integerValue + 1.0f) / (maxVal - 1.0f);
            } else {
                return integerValue / (maxVal - 1.0f);
            }
        }

        /**
         *  Normalizes (unsigned) integer vector inputs.
         *
         *  @param iVal An integer value which will be normalized to floating point based on
         *              it's value range.
         *
         *  @return [0,1] for unsigned inputs and [-1,1] for signed.
         */
        template<typename T>
        vec2 normalize(Vector2<T> integerVec) {
            return tgt::vec2(normalize(integerVec.x), normalize(integerVec.y));
        }

        /**
         *  Normalizes (unsigned) integer vector inputs.
         *
         *  @param iVal An integer value which will be normalized to floating point based on
         *              it's value range.
         *
         *  @return [0,1] for unsigned inputs and [-1,1] for signed.
         */
        template<typename T>
        vec4 normalize(Vector4<T> integerVec) {
            return tgt::vec4(normalize(integerVec.x), normalize(integerVec.y), normalize(integerVec.z), normalize(integerVec.w));
        }

        /**
         *  Normalizes (unsigned) integer vector inputs.
         *
         *  @param iVal An integer value which will be normalized to floating point based on
         *              it's value range.
         *
         *  @return [0,1] for unsigned inputs and [-1,1] for signed.
         */
        template<typename T>
        vec3 normalize(Vector3<T> integerVec) {
            return tgt::vec3(normalize(integerVec.x), normalize(integerVec.y), normalize(integerVec.z));
        }

        /**
         *  Specifiy to which texture target you bound the texture you want to use.
         *  This implementation does not support multi texturing. Therefore your texture
         *  must be bound to unit 0. To use multi texturing you have to bind a custom
         *  shader (then this setting does not have any effect).
         *  Is called automatically by tgt::Texture::enable(). Therefore you don't have to
         *  call this method in most cases.
         *
         *  @param texmode The texture target of texture unit 0 to which you have bound your texture.
         */
        void setTextureMode(TexturMode texmode);

        /**
         *  Get the current texture mode.
         *
         *  @return The current texture mode.
         */
        TexturMode textureMode() const;

        /**
         *  Enables basic phong lighting with a single light.
         *  Don't forget to set the light colors and material.
         *  Disabled by default.
         */
        void enableLighting();

        /**
         *  Disables lighting (use solid coloring).
         */
        void disableLighting();

        /**
         *  Set the position, color and attenuation of the light.
         *  The position must be supplied in eye space.
         *  Has no effect when lighting is disabled.
         *
         *  @param lightSource  The light source.
         *                      Pass a standard constructed LightSource to set it to default values.
         */
        void setLightSource(const LightSource& lightSource);

        /**
         *  Get the current position, color and attenuation of the light.
         */
        LightSource getLightSource() const;

        /**
         *  Sets the material which is used when lighting is activated.
         *  Is only used for light components which are not set by setColorMaterial.
         *  Can not be used between glBegin() and glEnd().
         *
         *  @param material The material.
         *                  Pass a standard constructed LightSource to set it to default values.
         */
        void setMaterial(const Material& material);

        /**
         *  Gets the current material which is used when lighting is activated.
         */
        Material getMaterial() const;

        /**
         *  Sets for which light component whe attribute/texture color should be used.
         *  When set for a light component it replaces the value set in the material for this component.
         *
         *  @param matColor This is a bitfield. You can pass multiple values by using binary OR.
         *                  Initial value is MATCOLNONE.
         */
        void setMaterialColor(MaterialColor matColor);

        /**
         *  Generates matrix uniform declarations for the matrices definded in tgt::MatrixStack.
         *  This is called automatically by GLSL::generateStandardShaderHeader. Therefore you
         *  only need to call this method if you don't use a generated header.
         *
         *  @return The textual representation of the uniform declaration.
         */
        std::string generateMatstackUniformHeader() const;

        /**
         *  Sets the current values of the tgt::MatrixStack to the currently bound shader.
         *  This is called automatically by the ImmediateMode. You don't have to call this when
         *  using ImmediateMode.
         */
        void setMatstackUniforms() const;

        /**
         *  Sets the current values of the tgt::MatrixStack to the specified shader.
         *  This is called automatically by the ImmediateMode. You don't have to call this when
         *  using ImmediateMode.
         *
         *  @param shader The shader program which should be modified (get it's uniforms set)
         */
        void setMatstackUniforms(tgt::Shader *shader) const;

        /**
         *  Start rendering primitives. Equivalent to glBegin(GLenum).
         *
         *  @param polygonMode Sets which primitive the specified vertices form.
         */
        void begin(PolygonMode polygonMode);

        /**
         *  Depicts that we are done with specifying vertices.
         *  Equivalent to glEnd().
         */
        void end();

        /**
         *  Specifiy a new vertex.
         *  Equivalent to glVertex2f(GLfloat, GLfloat)
         *
         *  @param vector The vertex of question.
         */
        void vertex(const tgt::vec2 &vector);

        inline void vertex(float x, float y) {vertex(tgt::vec2(x, y));}

        /**
         *  Specifiy a new vertex.
         *  Equivalent to glVertex3f(GLfloat, GLfloat, GLfloat)
         *
         *  @param vector The vertex of question.
         */
        void vertex(const tgt::vec3 &vector);

        inline void vertex(float x, float y, float z) {vertex(tgt::vec3(x, y, z));}

        /**
         *  Specifiy a new vertex.
         *  Equivalent to glVertex4f(GLfloat, GLfloat, GLfloat)
         *
         *  @param vector The vertex of question.
         */
        void vertex(const tgt::vec4 &vector);

        inline void vertex(float x, float y, float z, float w) {vertex(tgt::vec4(x, y, z, w));}

        /**
         *  Specify the current texcoord.
         *  Equivalent to glTexCoord1f(GLfloat)
         *
         *  @param vector The texcoord of question.
         */
        void texcoord(float vector);

        /**
         *  Specify the current texcoord.
         *  Equivalent to glTexCoord2f(GLfloat, GLfloat)
         *
         *  @param vector The texcoord vertex of question.
         */
        void texcoord(const tgt::vec2 &vector);

        inline void texcoord(float x, float y) {texcoord(tgt::vec2(x, y));}

        /**
         *  Specify the current texcoord.
         *  Equivalent to glTexCoord2f(GLfloat, GLfloat, GLfloat)
         *
         *  @param vector The texcoord vertex of question.
         */
        void texcoord(const tgt::vec3 &vector);

        inline void texcoord(float x, float y, float z) {texcoord(tgt::vec3(x, y, z));}

        /**
         *  Specify the current texcoord.
         *  Equivalent to glTexCoord4f(GLfloat, GLfloat, GLfloat, GLfloat)
         *
         *  @param vector The texcoord vertex of question.
         */
        void texcoord(const tgt::vec4 &vector);

        inline void texcoord(float x, float y, float z, float w) {texcoord(tgt::vec4(x, y, z, w));}

        /**
         *  Specify the current color.
         *  Equivalent to glColor3f(GLfloat, GLfloat, GLfloat)
         *
         *  @param vector The color of question.
         */
        void color(const tgt::vec3 &color);

        inline void color(float r, float g, float b) {color(tgt::vec3(r, g, b));}

        /**
         *  Specify the current color.
         *  Equivalent to glColor4f(GLfloat, GLfloat, GLfloat, GLfloat)
         *
         *  @param vector The color of question.
         */
        void color(const tgt::vec4 &color);

        inline void color(float r, float g, float b, float a) {color(tgt::vec4(r, g, b, a));}

        /**
         * Returns the current color.
         */
        tgt::vec4 getCurrentColor() const;

        /**
         *  Specify the current normal. Equivalent to glNormal3f(GLfloat, GLfloat, GLfloat)
         *
         *  @param vector The normal of question.
         */
        void normal(const tgt::vec3 &normal);

        void normal(float x, float y, float z) {normal(tgt::vec3(x, y, z));}

        /**
         *  Enable the clip plane specified by numClipPlane. Equivalent to glEnable(GL_CLIP_PLANE?).
         *
         *  @param numClipPlane The number of the clip plane to be enabled. Has to be < MAX_CLIP_PLANES!
         */
        void enableClipPlane(unsigned int numClipPlane);

        /**
         *  Disable the clip plane specified by numClipPlane. Equivalent to glDisable(GL_CLIP_PLANE?).
         *
         *  @param numClipPlane The number of the clip plane to be disabled. Has to be < MAX_CLIP_PLANES!
         */
        void disableClipPlane(unsigned int numClipPlane);

        /**
         *  Set the parameters of the clip plane in _eye coordinate space_.
         *  \note Therefore this is NOT equivalent to glClipPlane! See setModelSpaceClipPlane!
         *
         *  @param numClipPlane The number of the clip plane to be changed. Has to be < MAX_CLIP_PLANES!
         *  @param plane The new clip plane in eye coordinate space
         */
        void setClipPlaneEquation(unsigned int numClipPlane, const tgt::plane& plane);

        /**
         *  Set the parameters of the clip plane in model space, according to the current
         *  modelview transformation (defined by MatStack).
         *  This is equivalent to glClipPlane.
         *
         *  @param numClipPlane The number of the clip plane to be changed. Has to be < MAX_CLIP_PLANES!
         *  @param plane The new clip plane in model space
         */
        void setModelSpaceClipPlaneEquation(unsigned int numClipPlane, const tgt::plane& plane);

        /**
         *  Get the parameters of the clip plane in _eye coordinate space_.
         *
         *  @param numClipPlane The number of the clip plane to be gotten. Has to be < MAX_CLIP_PLANES!
         *
         *  @return The clip plane in eye coordinate space
         */
        tgt::plane getClipPlaneEquation(unsigned int numClipPlane) const;

        /**
         *  Get direct access to the specified clip plane.
         *  \note the clip plane is specified in eye coordinate space!
         *
         *  @param numClipPlane The number of the clip plane to be returned. Has to be < MAX_CLIP_PLANES!
         *
         *  @return Reference to the specified clip plane.
         */
        ClipPlane& getClipPlane(unsigned int numClipPlane);
        const ClipPlane& getClipPlane(unsigned int numClipPlane) const;

    private:

        /**
         *  Configuration of a default ImmediateMode shader.
         *  Bits 0 and 1 are used to specify the texture mode.
         *  Bit 2 is used to specify whether clipping is enabled or not.
         */
        typedef unsigned int ShaderConfiguration;
        enum ShaderConfigurationProperty {
            NO_TEX = 0,
            TEX_1D = 1,
            TEX_2D = 2,
            TEX_3D = 3,
            CLIPPING = 4,
        };
        const static unsigned int NUM_SHADER_CONFIGURATIONS = 8;

        /**
         *  The vertex structure which is fed to OpenGL.
         */
        struct Vertex {
            tgt::vec4 position;
            tgt::vec4 color;
            tgt::vec3 normal;
            tgt::vec4 texcoord;
        };
        BOOST_STATIC_ASSERT(sizeof(Vertex) ==
                            3 * 4 * sizeof(GLfloat) + 3 * sizeof(GLfloat));

        /**
         *  Sets the lighting parameters of the given shader according to
         *  the current state.
         *
         *  @param shader The shader the uniforms will be set for.
         */
        void setLightingUniforms(tgt::Shader& shader) const;

        /**
         *  Helper function which adds a vertex including
         *  all it's attributes to the verte buffer.
         *
         *  @param position The vertex which is added to the vertex buffer.
         */
        void uploadVertex(const tgt::vec4& position);

        /**
         *  Generates the header for the default shader.
         *  Call generateMatstackUniformHeader().
         *  We can not use GLSL::generateStandardShaderHeader because tgt must not depend on voreen.
         *
         *  @return The header which can be prepended to the shader code.
         */
        std::string generateHeader(ShaderConfiguration config) const;


        /**
         *  Get a shader according to the given ShaderConfiguration
         *  to be used when no other shader is bound.
         *  The shader may be compiled first.
         *
         *  @return A fitting shader.
         */
        tgt::Shader* getShader(ShaderConfiguration config);

        /**
         *  Determines wheter clipping is enabled, i.e. whether any
         *  of the MAX_CLIP_PLANES clip planes are enabled.
         *
         *  @return True if any clip plane is enabled, false otherwise.
         */
        bool clippingEnabled() const;

        /**
         *  Sets the clipping parameters of the given shader according to
         *  the current state.
         *
         *  @param shader The shader the uniforms will be set for.
         */
        void setClippingUniforms(tgt::Shader& shader) const;

        // ----------------------------------------------------------------------------


        /**
         *  Obligatory variable for the logging facility.
         */
        static const std::string loggerCat_;

        /**
         *  Holds the names of the sampler uniform variables.
         */
        std::vector<std::string> samplerNames_;

        /**
         *  The fallback shaders to simulate fixed function behaviour.
         */
        std::array<tgt::Shader*, NUM_SHADER_CONFIGURATIONS> shaders_;


        /*
         *  --------------------------------------------------------
         *  Various variables contributing to the internal
         *  state of the immediate mode simulator
         *  --------------------------------------------------------
         */

        bool started_;
        bool lightingEnabled_;
        PolygonMode polygonMode_;
        TexturMode texMode_;
        LightSource lightSource_;
        Material material_;
        MaterialColor materialColor_;
        std::array<ClipPlane, MAX_CLIP_PLANES> clipPlanes_;

        tgt::vec4 currentTexcoord_;
        tgt::vec4 currentColor_;
        tgt::vec3 currentNormal_;

        std::vector<Vertex> vertices_;
        GLuint vertexBuffer_;
    };

} // namespace tgt

#define IMode tgt::ImmediateMode::getRef()

#endif // TGT_IMMEDIATEMODE_H
