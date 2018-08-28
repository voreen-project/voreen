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

#ifndef VRN_ORIENTATIONOVERLAY_H
#define VRN_ORIENTATIONOVERLAY_H

#include "voreen/core/processors/imageprocessor.h"

#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/filedialogproperty.h"
#include "voreen/core/properties/cameraproperty.h"
#include "voreen/core/properties/floatproperty.h"
#include "voreen/core/properties/optionproperty.h"
#include "voreen/core/properties/shaderproperty.h"

#include "voreen/core/datastructures/volume/slice/slicehelper.h"

#include "tgt/texture.h"

namespace voreen {

    class TriangleMeshGeometryBase;
    class TriangleMeshGeometryUInt16IndexedColorNormal;
/**
 * Renders an orientation cube over the input rendering and
 * optionally draws a tripod.
 */
class VRN_CORE_API OrientationOverlay : public ImageProcessor {
public:
    /** Enums representing all supported overlays */
    enum OrientationType{
        OT_AXES_2D,             ///< renders 2 axes
        OT_AXES_3D,             ///< renders 3 axes
        OT_COLOR_CUBE,          ///< renders a colored cube
        OT_TEXTURE_CUBE,        ///< renders a textured cube
        OT_COLOR_TEXTURE_CUBE   ///< renders a colored and textured cube
    };

    OrientationOverlay();
    ~OrientationOverlay();
    virtual Processor* create() const;

    virtual std::string getClassName() const    { return "OrientationOverlay"; }
    virtual std::string getCategory() const     { return "Image Processing"; }
    virtual CodeState getCodeState() const      { return CODE_STATE_STABLE; }

    virtual bool isReady() const;

protected:
    virtual void setDescriptions() {
        setDescription("Adds an orientation overlay, e.g., a cube or a tripod, to the input rendering indicating the current orientation of the camera.");
    }

    virtual void beforeProcess();
    virtual void process();
    virtual void initialize();
    virtual void deinitialize();

private:
    /** Renders the axis overlay */
    void renderGeometry();
    void createAxesGeometry(bool createXAxis = true, bool createYAxis = true, bool createZAxis = true);
    void createCubeGeometry(bool colored, bool textured);
    void createOverlayGeometry();

    TriangleMeshGeometryUInt16IndexedColorNormal* createArrowGeometry(tgt::vec4 col);

    template<class T>
    void createCubeGeometry();
    template<class T>
    T* createQuadGeometry(tgt::vec4 color);

    /** Adjusts the visibility of all properties according to the current orientation type */
    void adjustPropertyVisibility();
    /** Forces the axis geometry to be recreated. */
    void invalidateGeometry();

    //----------------------------
    //  Texture handling
    //----------------------------
protected:
    /** Loads (and create) needed textures. */
    void loadTextures();
    /** Sets reloadTextures_ flag to true. */
    void reloadTextures();

    // texture handles
    tgt::Texture* frontTex_;
    tgt::Texture* backTex_;
    tgt::Texture* topTex_;
    tgt::Texture* leftTex_;
    tgt::Texture* bottomTex_;
    tgt::Texture* rightTex_;
    std::string textureNames_[6];

    bool reloadTextures_;
    bool loadingTextures_;          ///< Set to true during texture loading for preventing
                                    ///  multiple/cyclic execution of loadTextures().

private:
    //ports
    RenderPort inport_;             ///< Input rendering the orientation overlay is drawn onto.
    RenderPort outport_;            ///< Output: input + orientation overlay
    RenderPort privatePort_;        ///< stored the rendered orientation overlay.
    //basic
    BoolProperty enableProp_;                               ///< used to disable overlay
    OptionProperty<OrientationType> orientationTypeProp_;   ///< determining the current used overlay
    CameraProperty cameraProp_;
    //position
    FloatProperty shiftXProp_;              ///< Distance to shift cube and axis horizontally.
    FloatProperty shiftYProp_;              ///< Distance to shift cube and axis vertically.
    //axis
    FloatProperty axisSizeProp_;            ///< Length of axes indicating orientation.
    OptionProperty<SliceAlignment> axesAlignmentProp_; ///< alignment in 2D mode
    BoolProperty renderAxesLabelsProp_;    ///< if true, the axes are labeled
    //cube
    FloatProperty cubeSizeProp_;            ///< Size of cube indicating orientation.
    //textures
    FileDialogProperty filenameFrontProp_;  ///< Filename of front texture.
    FileDialogProperty filenameBackProp_;   ///< Filename of back texture.
    FileDialogProperty filenameTopProp_;    ///< Filename of top texture.
    FileDialogProperty filenameBottomProp_; ///< Filename of bottom texture.
    FileDialogProperty filenameLeftProp_;   ///< Filename of left texture.
    FileDialogProperty filenameRightProp_;  ///< Filename of Right texture.
    //shader
    ShaderProperty shaderProp_;             ///< shader property used to render geometry

    TriangleMeshGeometryBase* currentGeometry_; ///< stores the current geometry
    TriangleMeshGeometryBase* letterXGeometry_;  ///< geometry of the axes labels
    TriangleMeshGeometryBase* letterYGeometry_;  ///< geometry of the axes labels
    TriangleMeshGeometryBase* letterZGeometry_;  ///< geometry of the axes labels
    bool geometryMustBeRecreated_;          ///< true at the beginning and if the axis size has been changed
    const float overlayBaseLength_;         ///< the base length of the overlay

    static const std::string loggerCat_; ///< category used in logging
};

template<class T>
T* OrientationOverlay::createQuadGeometry(tgt::vec4 color) {
    float size = overlayBaseLength_*cubeSizeProp_.get();
    T* geom = new T();
    geom->addQuadGeometry(size,size,tgt::vec3::zero,color);
    return geom;
}

template<class T>
void OrientationOverlay::createCubeGeometry() {
    float shift = overlayBaseLength_*cubeSizeProp_.get()/2.f;

    T* xp_Quad = createQuadGeometry<T>(tgt::vec4(1.f,0.f,0.f,1.f));
    xp_Quad->setTransformationMatrix(tgt::Matrix::createTranslation(tgt::vec3(shift,0.f,0.f)) * tgt::Matrix::createRotationY(1.570796327f) * tgt::Matrix::createRotationZ(1.570796327f));
    T* xn_Quad = createQuadGeometry<T>(tgt::vec4(1.f,0.f,0.f,1.f));
    xn_Quad->setTransformationMatrix(tgt::Matrix::createTranslation(tgt::vec3(-shift,0.f,0.f)) * tgt::Matrix::createRotationY(-1.570796327f) * tgt::Matrix::createRotationZ(-1.570796327f));
    T* yp_Quad = createQuadGeometry<T>(tgt::vec4(0.f,1.f,0.f,1.f));
    yp_Quad->setTransformationMatrix(tgt::Matrix::createTranslation(tgt::vec3(0.f,shift,0.f)) * tgt::Matrix::createRotationX(-1.570796327f) * tgt::Matrix::createRotationZ(3.141692654f));
    T* yn_Quad = createQuadGeometry<T>(tgt::vec4(0.f,1.f,0.f,1.f));
    yn_Quad->setTransformationMatrix(tgt::Matrix::createTranslation(tgt::vec3(0.f,-shift,0.f)) * tgt::Matrix::createRotationX(1.570796327f));
    T* zp_Quad = createQuadGeometry<T>(tgt::vec4(0.f,0.f,1.f,1.f));
    zp_Quad->setTransformationMatrix(tgt::Matrix::createTranslation(tgt::vec3(0.f,0.f,shift)) );
    T* zn_Quad = createQuadGeometry<T>(tgt::vec4(0.f,0.f,1.f,1.f));
    zn_Quad->setTransformationMatrix(tgt::Matrix::createTranslation(tgt::vec3(0.f,0.f,-shift)) * tgt::Matrix::createRotationX(3.141692654f)); //pi

    if(xp_Quad->supportsTextureData()) {
        xp_Quad->loadTextureData(filenameLeftProp_.get());
        xn_Quad->loadTextureData(filenameRightProp_.get());
        yp_Quad->loadTextureData(filenameBackProp_.get());
        yn_Quad->loadTextureData(filenameFrontProp_.get());
        zp_Quad->loadTextureData(filenameTopProp_.get());
        zn_Quad->loadTextureData(filenameBottomProp_.get());
    }

    xp_Quad->addMesh(xn_Quad);
    xp_Quad->addMesh(yp_Quad);
    xp_Quad->addMesh(yn_Quad);
    xp_Quad->addMesh(zp_Quad);
    xp_Quad->addMesh(zn_Quad);

    delete xn_Quad;
    delete yp_Quad;
    delete yn_Quad;
    delete zp_Quad;
    delete zn_Quad;

    currentGeometry_ =  xp_Quad;
}


} // namespace

#endif
