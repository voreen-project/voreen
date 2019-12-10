/**********************************************************************
 *                                                                    *
 * tgt - Tiny Graphics Toolbox                                        *
 *                                                                    *
 * Copyright (C) 2005-2019 University of Muenster, Germany,           *
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

#ifndef TGT_MATERIAL_H
#define TGT_MATERIAL_H

#include "tgt/vector.h"
#include "tgt/tgt_gl.h"
#include "tgt/texturemanager.h"

namespace tgt {

class Texture;

/**
    This class is used to specify material parameters for
    the Phong lighting model and read the currently used
    material parameters from the OpenGL status machine.
*/
class Material {

public:

    /// This struct is used to set which face the material parameters are specified for.
    enum Face {
        FRONT   = GL_FRONT,                    /*!< Front face (direction of normals).                    */
        BACK    = GL_BACK,                    /*!< Back face (perpendicular to direction of normals).    */
        FRONT_AND_BACK  = GL_FRONT_AND_BACK    /*!< Material parameters are set for both face            */
    };

    // Lighting types.
    enum LightingType {
        CONSTANT,
        LAMBERT,
        PHONG,
        BLINN,
        UNKNOWN
    };

    /// Transparency modes supported by the COMMON profile of COLLADA.
    enum TransparencyMode {
        /** Takes the transparency information from the color's alpha channel, where the
            value 1.0 is opaque. */
        A_ONE,

        /** Takes the transparency information from the color's red, green, and blue channels,
            where the value 0.0 is opaque, with each channel modulated independently. */
        RGB_ZERO
    };


    /* Constructors */

    /// Standard constructor. Reflectance parameters and emission are set to OpenGL defaults.
    Material(Face face = FRONT) {
        setValues(vec4(0.2f,0.2f,0.2f,1.0f), vec4(0.8f,0.8f,0.8f,1.0f), vec4(0.0f,0.0f,0.0f,1.0f),
                  0.0, vec4(0.0f,0.0f,0.0f,1.0f), face);
    };

    //! Constructor. Ambient and diffuse reflectance are set to ambient_and_diffuse.
    /// No specular reflection and emission.
    Material(const vec4& ambient_and_diffuse, Face face = FRONT){
        setValues(ambient_and_diffuse, ambient_and_diffuse, vec4(0.f,0.f,0.f,1.f),
                  0.0, vec4(0.f,0.f,0.f,1.f), face);
    };

    //! Constructor. Ambient and diffuse reflectance are specified.
    /// No specular reflection and emission.
    Material(const vec4& ambient, const vec4& diffuse, Face face = FRONT){
        setValues(ambient, diffuse, vec4(0.f,0.f,0.f,1.f), 0.0, vec4(0.f,0.f,0.f,1.f), face);
    };

    //! Constructor. Ambient, diffuse, specular reflectance are specified. No emission.
    /// shininess must be inside the range [0,128].
    Material(const vec4& ambient, const vec4& diffuse, const vec4& specular,
             GLfloat shininess, Face face = FRONT){
        setValues(ambient, diffuse, specular, shininess, vec4(0.f,0.f,0.f,1.f), face);
    };

    //! Constructor. All material parameters are specified.
    /// shininess must be inside the range [0,128].
    Material(const vec4& ambient, const vec4& diffuse, const vec4& specular,
             GLfloat shininess, const vec4& emission, Face face = FRONT){
        setValues(ambient, diffuse, specular, shininess, emission, face);
    };

    //! Destructor.
    ~Material() {
        if ( textureDiffuse_ ) TexMgr.dispose( textureDiffuse_ );
        if ( textureReflective_ ) TexMgr.dispose( textureReflective_ );
        if ( textureTransparent_ ) TexMgr.dispose( textureTransparent_ );
    }

    //! Applies the material parameters for the prior specified face.
    /// \note{This is the only method of this class that changes the OpenGL state.}
    void activate(){
#ifdef VRN_OPENGL_COMPATIBILITY_PROFILE
        glMaterialfv(face_, GL_AMBIENT, ambient_.elem);
        glMaterialfv(face_, GL_DIFFUSE, diffuse_.elem);
        glMaterialfv(face_, GL_SPECULAR, specular_.elem);
        glMaterialf(face_, GL_SHININESS, shininess_);
        glMaterialfv(face_, GL_EMISSION, emission_.elem);
#else
        LWARNINGC("tgt.Material", "apply: not supported in Core-Profile");
#endif
    };

    /// Reads the material parameters from OpenGL and saves them to the object.
    void read(){
#ifdef VRN_OPENGL_COMPATIBILITY_PROFILE
        if (face_ == FRONT_AND_BACK) face_ = FRONT;
        glGetMaterialfv(face_, GL_AMBIENT, ambient_.elem);
        glGetMaterialfv(face_, GL_DIFFUSE, diffuse_.elem);
        glGetMaterialfv(face_, GL_SPECULAR, specular_.elem);
        glGetMaterialfv(face_, GL_SHININESS, &shininess_);
        glGetMaterialfv(face_, GL_EMISSION, emission_.elem);
#else
        LWARNINGC("tgt.Material", "read: not supported in Core-Profile");
#endif
    };


    /* Getters and Setters */

    //! Set all material parameters at once.
    /// shininess must be inside the range [0,128].
    void setValues(const vec4& ambient, const vec4& diffuse, const vec4& specular, GLfloat shininess, const vec4& emission, Face face = FRONT){
        tgtAssert(shininess >= 0 && shininess <= 128,
            "Invalid Shininess value. Must be inside the range [0,128].");
        ambient_ = ambient;
        diffuse_ = diffuse;
        specular_ = specular;
        shininess_ = shininess;
        emission_ = emission;
        face_ = face;

        textureDiffuse_ = NULL;
        textureReflective_ = NULL;
        textureTransparent_ = NULL;

        isReflective_ = false;
        reflectivity_ = vec4 ( 0.f, 0.f, 0.f, 0.f );

        transparencyMode_ = RGB_ZERO;
        isRefractive_ = false;
        indexOfRefraction_ = 1.0f;

        lightingType_ = UNKNOWN;
    };

    /// Set ambient reflectance.
    void setAmbient(const vec4& ambient) { ambient_ = ambient; };

    /// Get ambient reflectance.
    vec4 getAmbient(){ return ambient_;    };

    /// Set diffuse reflectance.
    void setDiffuse(const vec4& diffuse) { diffuse_ = diffuse; };

    /// Get diffuse reflectance.
    vec4 getDiffuse(){ return diffuse_; };

    /// Set specular reflectance.
    void setSpecular(const vec4& specular) { specular_ = specular; };

    /// Get specular reflectance.
    vec4 getSpecular() { return specular_; };

    //! Set specular exponent.
    /// Must be inside the range [0,128]
    void setShininess(GLfloat shininess){
        tgtAssert(shininess >= 0 && shininess <= 128,
            "Invalid Shininess value. Must be inside the range [0,128].");
        shininess_ = shininess;
    };

    /// Get specular exponent.
    GLfloat getShininess() { return shininess_; };

    /// Set emission.
    void setEmission(const vec4& emission) { emission_ = emission; };

    /// Get emission.
    vec4 getEmission() { return emission_; };

    /// Set face which the material properties are specified for:
    void setFace(Face face) { face_ = face; };

    /// Get face.
    Face getFace() { return face_; };



    /// Set material translucency.
    void setTranslucency( vec4 t ) { translucency_ = t;  };
    /// Get material translucency.
    vec4 getTranslucency() { return translucency_; };

    /// Set transparency mode.
    void setTransparencyMode( TransparencyMode tm ) { transparencyMode_ = tm;  };
    /// Get transparency mode.
    TransparencyMode getTransparencyMode() { return transparencyMode_; };

    /// Set diffuse Texture.
    void setTextureDiffuse( Texture* t ) { textureDiffuse_ = t;  };
    /// Get diffuse Texture.
    Texture* getTextureDiffuse() { return textureDiffuse_; };

    /// Set reflective Texture.
    void setTextureReflective( Texture* t ) { textureReflective_ = t;  };
    /// Get reflective Texture.
    Texture* getTextureReflective() { return textureReflective_; };

    /// Set transparent Texture.
    void setTextureTransparent( Texture* t ) { textureTransparent_ = t;  };
    /// Get transparent Texture.
    Texture* getTextureTransparent() { return textureTransparent_; };

    /// Set refractive state.
    void setRefractive( bool r ) { isRefractive_ = r; };
    /// Get refractive state.
    bool getRefractive() { return isRefractive_; };

    /// Set Index of refraction.
    void setIndexOfRefraction( GLfloat ior ) { indexOfRefraction_ = ior; };
    /// Get Index of refraction.
    GLfloat getIndexOfRefraction() { return indexOfRefraction_; };

    /// Set reflectivity color.
    void setReflectivity( vec4 r ) { reflectivity_ = r; };
    /// Get reflectitiy color.
    vec4 getReflectivity() { return reflectivity_; };

    /// Set reflective state.
    void setReflective( bool r ) { isReflective_ = r; };
    /// Get reflective state.
    bool getReflective() { return isReflective_; };

    /// Set Lighting Type.
    void setLightingType( LightingType l ) { lightingType_ = l; };
    /// Get Lighting Type.
    LightingType getLightingType() { return lightingType_; };

private:

    // ambient, diffuse, specular reflectance and light emission
    vec4 ambient_;
    vec4 diffuse_;
    vec4 specular_;
    GLfloat shininess_;
    vec4 emission_;

    // face which the material properties are specified for.
    Face face_;

    // FIXME: Are these attributes needed?
    // The below attributes were added to have full collada support.
    // Check if really needed and change functions accordingly, especially the state changing ones.
    // Getter & setter already implemented.
    bool isReflective_;
    vec4 reflectivity_;
    vec4 translucency_;

    TransparencyMode transparencyMode_;

    bool isRefractive_;
    GLfloat indexOfRefraction_;

    LightingType lightingType_;

    Texture* textureDiffuse_;
    Texture* textureReflective_;
    Texture* textureTransparent_;
};

} // namespace tgt

#endif //TGT_MATERIAL_H
