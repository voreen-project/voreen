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

#ifndef VRN_CAMERAPOSITIONRENDERER_H
#define VRN_CAMERAPOSITIONRENDERER_H

#include "voreen/core/processors/geometryrendererbase.h"

#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/cameraproperty.h"
#include "voreen/core/properties/shaderproperty.h"
#include "voreen/core/properties/colorproperty.h"
#include "voreen/core/properties/floatproperty.h"
#include "voreen/core/properties/optionproperty.h"

namespace voreen {

/**
 * Depicts the position of a camera in a 3D scene.
 * The camera position to display is retrieved from the processor's camera property. To use this processor,
 * the automatic camera linking should be disabled.
 */
class VRN_CORE_API CameraPositionRenderer : public GeometryRendererBase {
public:
    CameraPositionRenderer();
    ~CameraPositionRenderer();
    virtual Processor* create() const;

    virtual std::string getClassName() const { return "CameraPositionRenderer"; }
    virtual std::string getCategory() const  { return "Geometry";               }
    virtual CodeState getCodeState() const   { return CODE_STATE_BROKEN;        }

    virtual std::string generateHeader(const tgt::GpuCapabilities::GlVersion* version = 0);
    virtual void initialize();
    virtual void deinitialize();
    virtual void render();

private:
    virtual void setDescriptions() {
        setDescription("Depicts the position of a camera in a 3D scene. The camera position to display is retrieved from the processor's camera property. "\
                       "To use this processor, the automatic camera linking should be disabled.");
        sizeReferencePort_.setDescription("The port is used to get the current screen size. Connect the outport of the processor, "\
                                          "which camera property has been linked.");
        displayCameraProp_.setDescription("Camera, which will be visualized and which has to be linked. Also the outport of the linked camera property "\
                                          "has to be connected to the render inport of this processor.");
        myShaderProp_.setDescription("Handles interaction with the custom 'camposition' shaders.");
        renderFrustumProp_.setDescription("Enables the rendering of the frustum borders.");
        frustLineColorProp_.setDescription("Sets the color of frustum outline (and near and far plane, if active).");
        renderPlanesProp_.setDescription("Enables the rendering of near plane and far plane, on the condition that the frustum is rendered.");
        opacPlanesProp_.setDescription("Sets the opacity of near and far plane.");
        optionModelProp_.setDescription("Enables the rendering of the camera model.");
        opacModelProp_.setDescription("Sets the opacity of the camera model.");
        colorModelProp_.setDescription("Sets the color of the camera model.");
        scalingModelProp_.setDescription("Sets the scale of the camera model.");
        optionArrowsProp_.setDescription("Handles if and how the directions the linked camera's relative directions are rendered.");
        xLengthProp_.setDescription("Sets the length of the rendered coordinate system's x-axis.");
        yLengthProp_.setDescription("Sets the length of the rendered coordinate system's y-axis.");
        zLengthProp_.setDescription("Sets the length of the rendered coordinate system's z-axis.");
        scalingArrowsProp_.setDescription("Sets the scale of the direction arrows.");
    }

private:
    //------------------
    //  Helper Functions
    //------------------

    /** Creates the frustum VBO. It is stored as vbo_[0] and uses vao_[0] for the draw call. */
    void createFrustumBuffer();
    /** Draws Frustum*/
    void renderFrustum();
    /** Manages visibility of frustum rendering options*/
    void onRenderFrustumChanged();
    /** Manages visibility of planes rendering options*/
    void onRenderPlanesChanged();
    /** Creates the camera model VBO. It is stored as vbo_[1] and uses vao_[1] for the draw call. */
    void createCameraModelBuffer();
    /** Draws Camera Model*/
    void renderCameraModel();
    /** Manages visibility of camera rendering options*/
    void onRenderModelChanged();
    /** Creates the Arrows VBO. It is stored as vbo_[2] and uses vao_[2] for the draw call. */
    void createArrowsBuffer();
    /** Draws Arrows Model*/
    void renderArrows();
    /** Manages visibility of arrows rendering options*/
    void onRenderArrowsChanged();

    //------------------
    //  Members
    //------------------
    //Ports
    RenderPort sizeReferencePort_;              ///< Used to get the current screen size
    //Properties
    CameraProperty displayCameraProp_;          ///< camera, which has to be linked 
    ShaderProperty myShaderProp_;               ///< shader used to display frustum and camera model
    BoolProperty renderFrustumProp_;            ///< enables the rendering of the frustum
    ColorProperty frustLineColorProp_;          ///< color of the lines describing the frustum
    FloatProperty frustLineWidthProp_;          ///< width of the lines describing the frustum
    BoolProperty renderPlanesProp_;             ///< enables rendering the near and far plane
    FloatProperty opacPlanesProp_;              ///< opacity of near and far plane
    StringOptionProperty optionModelProp_;      ///< enables rendering of camera model
    FloatProperty opacModelProp_;               ///< opacity of camera model
    ColorProperty colorModelProp_;              ///< color of camera model
    FloatProperty lineModelProp_;               ///< line width of camera model
    FloatProperty scalingModelProp_;            ///< scale of camera model
    StringOptionProperty optionArrowsProp_;     ///< enables the rendering of the direction arrows
    FloatProperty scalingArrowsProp_;           ///< scale of direction arrows
    FloatProperty xLengthProp_;                 ///< length value of x axis
    FloatProperty yLengthProp_;                 ///< length value of y axis
    FloatProperty zLengthProp_;                 ///< length value of z axis

    //Variables
    GLuint vbo_[3];                             ///< OpenGL Vertex Buffer Object
    GLuint vao_[3];                             ///< OpenGL Vertex Array Object
    tgt::mat4 vInv_;                            ///< inverted View Matrix of secondary camera, passed to shader
    tgt::mat4 pvInv_;                           ///< inverted (Projection Matrix * View Matrix) of secondary camera, passed to shader
};

}

#endif
