/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany,                        *
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

#ifndef VRN_GEOMETRYRENDERER_H
#define VRN_GEOMETRYRENDERER_H

#include "voreen/core/processors/geometryrendererbase.h"

#include "voreen/core/properties/optionproperty.h"
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/floatproperty.h"
#include "voreen/core/properties/shaderproperty.h"
#include "voreen/core/properties/lightsourceproperty.h"
#include "voreen/core/properties/colorproperty.h"

namespace voreen {

/**
 * Basic processor for rendering arbitrary geometry, simply taking
 * a Geometry object through its inport und calling render() on it.
 */
class VRN_CORE_API GeometryRenderer : public GeometryRendererBase {
public:
    GeometryRenderer();
    ~GeometryRenderer();
    virtual Processor* create() const;

    virtual std::string getCategory() const   { return "Geometry";         }
    virtual std::string getClassName() const  { return "GeometryRenderer"; }
    virtual CodeState getCodeState() const    { return CODE_STATE_STABLE;  }

    virtual bool isReady() const;

    virtual bool usesTransparency() const;

    virtual void render();
    virtual void renderTransparent();

    virtual tgt::Bounds getBoundingBox() const;

    /**
     * Calls render() on the Geometry object.
     */
    virtual void render(ShaderProperty& shaderProp);

    virtual void process();

    void render(tgt::Shader* shader);
protected:
    virtual void setDescriptions() {
        setDescription("Renders the passed geometry object with adjustable polygon mode and lighting settings. Additionally, a texture passed trough the render inport can be mapped onto the geometry (texture coordinates have to be set by the geometry object).");
    }

    virtual void adjustPropertiesToInput();
    virtual void updatePropertyVisibilities();

    bool isUsingSolidColor() const;

    void invalidateShader();
    void setupShaders();

    std::string generateShaderDefines();

    GeometryPort inport_;
    RenderPort texPort_;
    GLEnumOptionProperty polygonMode_;

    FloatProperty pointSize_;
    FloatProperty lineWidth_;

    ColorProperty solidColor_;          ///< color for rendering meshes without color information in the vertices

    BoolProperty mapTexture_;
    StringOptionProperty textureMode_;

    BoolProperty enableLighting_;
    LightSourceProperty lightPosition_; ///< The position of the light source in world coordinates
    ColorProperty lightAmbient_;        ///< The light source's ambient color
    ColorProperty lightDiffuse_;        ///< The light source's diffuse color
    ColorProperty lightSpecular_;       ///< The light source's specular color
    FloatProperty materialShininess_;   ///< The material's specular exponent

    ShaderProperty shaderOpaqueProp_;
    ShaderProperty shaderTransparentProp_;
    BoolProperty enableClipping_;
    FloatVec3Property planeNormal_;
    FloatProperty planeDistance_;

    int lastShaderConfiguration_; ///< -1, if no shader has been used
    bool shaderNeedsRebuild_;

    tgt::FramebufferObject* fbo_;
    tgt::Texture* arrayTex_;
};

}

#endif // VRN_GEOMETRYRENDERER_H

