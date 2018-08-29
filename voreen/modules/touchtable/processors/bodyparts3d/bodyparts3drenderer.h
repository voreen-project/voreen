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

#ifndef VRN_BODYPARTS3DRENDERER_H
#define VRN_BODYPARTS3DRENDERER_H

#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/eventproperty.h"
#include "voreen/core/properties/floatproperty.h"
#include "voreen/core/properties/shaderproperty.h"
#include "voreen/core/properties/optionproperty.h"
#include "voreen/core/properties/cameraproperty.h"
#include "voreen/core/properties/colorproperty.h"
#include "voreen/core/properties/stringproperty.h"
#include "voreen/core/properties/buttonproperty.h"
#include "voreen/core/properties/lightsourceproperty.h"

#include "voreen/core/ports/textport.h"

#include "voreen/core/interaction/idmanager.h"
#include "voreen/core/processors/renderprocessor.h"
#include "voreen/core/datastructures/geometry/trianglemeshgeometry.h"

#include "modules/touchtable/ports/bodyparts3d/bodypartsport.h"
#include "modules/touchtable/ports/bodyparts3d/bodypartstextureport.h"

#include "tgt/shadermanager.h"

namespace voreen{

class CameraInteractionHandler;

class VRN_CORE_API BodyParts3DRenderer : public RenderProcessor {
public:
    BodyParts3DRenderer();
    virtual ~BodyParts3DRenderer();

    virtual Processor* create() const {return new BodyParts3DRenderer();}
    virtual std::string getClassName() const { return "BodyParts3DRenderer";}
    virtual std::string getCategory() const  { return "Renderer";}
    virtual CodeState getCodeState() const   { return CODE_STATE_EXPERIMENTAL;}

    virtual bool isReady() const;

    virtual void adjustRenderOutportSizes();

    std::string generateHeader(const tgt::GpuCapabilities::GlVersion* version = 0);

    void compile();

protected:

    virtual void setDescriptions(){
        setDescription("This is the Processor responsible for rendering BodyParts3D");
    };

    virtual void initialize();
    virtual void beforeProcess();
    virtual void process();
    virtual void deinitialize();

    /**
     * Renders the Geometries to the outport
     */
    void renderGeometry();

    /**
     * Renders the Geometries to the picking Port with the color converted from the index they have in portInput_
     */
    void renderPickingPort();

    /**
     * Renders a Texture of the frontal view of the geometry
     *
     * @param geometry Geometry for which a texture is generated
     * @return pointer to an uploaded texture
     */
    tgt::Texture* renderGeometryTexture(TriangleMeshGeometrySingleColor* geometry);

    virtual void updatePropertyVisibilities();

    /**
     * Manages Events for the processor or passes them to the inport
     *
     * @param e Event to be handled by the processor
     */
    void onEvent(tgt::Event* e);

    /**
     * Checks wether or not a geometry was hit by the mouseevent and highlights it
     *
     * In case onEvent caught a mouse event and processor is in highlight-mode (i.e. action_ is false),
     * this function will highlight the geometry at mouseposition (if there is one)
     *
     * @param e Mouseevent which is checked for a hit with a geometry
     */
    void highlightGeometry(tgt::MouseEvent* e);

    /**
     * Checks wether or not a geometry was hit by the mouseevent and removes it
     *
     * In case onEvent caught a mouse event and processor is in remove-mode (i.e. action_ is true),
     * this function will remove the geometry at mouseposition (if there is one) and put in at the back of
     * removedParts_.
     *
     * @param e Mouseevent which is checked for a hit with a geometry
     */
    void removeGeometry(tgt::MouseEvent* e);

    /**
     * Reverts the removal of the last geometry
     */
    void addGeometry();

    /**
     * restores the geometry with the id restoreSignal_ from removedParts_ to portInput_
     *
     * Called when restoreSignal is changed. Intended to be invoked when the BodyPartsWidget sends a Signal
     */
    void restore();

    /**
     * Resets the current highlight
     */
    void resetHighlight();
private:
    /**
     * Calculates the Boundingbox of all input Geometries
     */
    virtual void calculateBoundingBox();

    std::vector<std::pair<TriangleMeshGeometrySingleColor*,std::string> > portInput_;
    std::vector<std::pair<TriangleMeshGeometrySingleColor*,std::string> > removedParts_;
    std::vector<std::pair<tgt::Texture*,std::string> > textures_;

    //ports
    BodyPartsPort inPort_;
    GeometryPort boundingBoxPort_;
    BodyPartsTexturePort texturePort_;
    RenderPort renderOutPort_;
    TextPort textOutPort_;

    RenderPort pickingPort_; ///< Port for picking
    RenderPort privateRenderPort_; ///< port for texture generation

    tgt::Shader* pickingShaderProgram_; ///< picking shader
    tgt::Shader* textureShaderProgram_;
    ShaderProperty shaderProgram_; ///< render
    StringOptionProperty shadeMode_; ///< dropdown menu for shader selection

    CameraProperty camera_;

    LightSourceProperty lightPosition_; ///< Relative Light Position in space
    ColorProperty lightAmbient_; ///< ambient light color
    ColorProperty lightDiffuse_; ///< diffuse light color
    ColorProperty lightSpecular_; ///< specular light color
    FloatProperty shininess_; ///< material shininess

    BoolProperty enableClipping_; ///< enables clipping
    FloatProperty planeManipulatorScale_; ///< scales the size of the planemanipulator down to fit non cubic bounding boxes
    FloatProperty planeDistance_; ///< clipping plane distance
    FloatVec3Property planeNormal_; ///< clipping plane normal

    BoolProperty enablePicking_; ///< disables the picking rendering to increase performance if wished
    BoolProperty action_; ///< switches between highlighting and removing on mouse/touchevent
    BoolProperty invertPlane_; ///< while true the uniform plane given to the shaders will be inverted

    CameraInteractionHandler* cameraHandler_;

    IDManager idManager_; ///< idManager mainly being used for color generation and rendertarget management

    tgt::mat4 transformationMatrix_; ///< describes a scale and translation so that the geometry is alway in [-1,1]x[-1,1]x[-1,1]
    tgt::Bounds boundingBox_; ///< BoundingBox containing all geometrys displayed

    IntProperty highlight_; ///< vector index of the highlighted geometry
    ButtonProperty revert_; ///< resotres the last removed geometry
    IntProperty restoreSignal_; ///< onChange: restores the geometry with index restoreSignal_.get() in removedParts_

    static const std::string loggerCat_;
};

}

#endif //VRN_BODYPARTS3DRENDERER_H
