/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2024 University of Muenster, Germany,                        *
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

#ifndef VRN_PLANEMANIPULATION_H
#define VRN_PLANEMANIPULATION_H

#include "voreen/core/processors/geometryrendererbase.h"
#include "voreen/core/properties/eventproperty.h"
#include "voreen/core/properties/vectorproperty.h"
#include "voreen/core/properties/floatproperty.h"
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/cameraproperty.h"
#include "voreen/core/properties/buttonproperty.h"
#include "voreen/core/properties/colorproperty.h"
#include "voreen/core/properties/shaderproperty.h"

#include "voreen/core/ports/volumeport.h"
#include "voreen/core/ports/geometryport.h"

#include "voreen/core/datastructures/geometry/glmeshgeometry.h"

#include <vector>

namespace voreen {

/**
 * Renders 3D-control elements to manipulate an arbitrary oriented clipping plane
 * using linking to propagate property changes to @c GeometryClipping processors.
 */
class PlaneManipulation : public GeometryRendererBase {
private:
    /** Definition of @c Anchor data structure. */
    typedef tgt::vec3 Anchor;
    /** Definition of @c Line data structure. */
    typedef std::pair<Anchor*, Anchor*> Line;

    //----------------------//
    // Override Functions   //
    //----------------------//
public:
    /** Constructor. */
    PlaneManipulation();
    /** Destructor */
    ~PlaneManipulation();
    virtual Processor* create() const;
    virtual std::string getClassName() const override { return "PlaneManipulation"; }
    virtual std::string getCategory() const override  { return "Geometry";          }
    virtual CodeState getCodeState() const override   { return CODE_STATE_TESTING;  }

protected:
    virtual void setDescriptions() override {
        setDescription("3D control elements to manipulate an arbitrarily oriented geometry clipping plane: \
                <ul> \
                <li>Use linking to connect the properties to a GeometryClipping processor.</li>\
                <li>Use the <b>left mouse button</b> to change the plane orientation by rotating the manipulator handle through \
                grabbing one of its tips or shift the plane along its normal by grabbing the cylinder.</li>\
                <li>Use the <b>right mouse button</b> to shift the manipulator on the plane without changing the plane itself \
                (only if the plane, ie. its intersection with the input bounding box, is visible).</li>\
                </ul>");
        //linking
        enable_.setDescription("Disables the rendering and can be linked to GeometryClipping.");
        invert_.setDescription("Link to GeometryClipping processor to invert the plane's orientation.");
        planeNormal_.setDescription("Normal of the plane, link with GeometryClipping processor.");
        planePosition_.setDescription("Distance of the plane from the origin (in world coordinates). Link with GeometryClipping processor.");
        //general
        generalColor_.setDescription("Color of the plane and manipulator handle rendering.");
        resetPlane_.setDescription("This resets the plane so that the plane normal is (0,0,-1) and the plane goes through the center of the input bounding box. \
                                    The manipulator handle is set to its default position on the plane.");
        //plane
        renderPlane_.setDescription("Selects if the geometry of the plane (ie. its intersection with the input bounding box) is rendered.");
        lineWidth_.setDescription("Width of the lines that are rendered for the plane visualization.");
        grabbedPlaneOpacity_.setDescription("Opacity of the plane, which is rendered during interaction with the manipulator handle.");
        alwaysUseOpacity_.setDescription("Always visualize the plane with the defined opacity.");
        renderAnchors_.setDescription("Render anchors at bounding box intersections with the clipping plane. This is only an optical feature.");
        //manipulator
        renderManipulator_.setDescription("Selects if the manipulator handle is rendered. If disabled this also prevents any interaction with the handle.");
        manipulatorScale_.setDescription("Scaling factor to adjust the size of the manipulator handle, which is computed from the bounding box of the input volume or geometry.");
        grabbedElementColor_.setDescription("Color of the manipulator element that is currently grabbed.");
        blockedElementColor_.setDescription("Color of an element that is grabbed, but may currently not be moved (e.g. moving the manipulator outside the scene");
        resetManipulatorPos_.setDescription("This resets the manipulator handle to the default position (without changing the plane), \
                                            which is the point on the plane nearest to the center of the input bounding box.");
    }

    /** Used to initialize the shader and geometries. */
    virtual void initialize() override;
    /** Used to deinitialize the shader and geometries. */
    virtual void deinitialize() override;
    /** Invalidates (rebuilds) the shader. */
    virtual void invalidate(int inv = INVALID_RESULT) override;
    /** Porcessor should be ready, if at least one port is connected. */
    virtual bool isReady() const override;
    /** Initialize ID manager and also register objects for picking. */
    virtual void setIDManager(IDManager* idm) override;
    /** Renders the manipulator.
     * @see GeometryRendererBase::render
     */
    virtual void render() override;
    /** Renders the manipulator using transparency.
     * @see GeometryRendererBase::render
     */
    virtual void renderTransparent() override;
    /** Renders the manipulator using transparency depending on parameter.
     * @see GeometryRendererBase::renderTransparent
     */
    void render(bool useTransparency);
    /** Determines, whether he processesor uses transparency in the next rendering run
     * @see GeometryRendererBase::usesTransparency
     */
    virtual bool usesTransparency() const override;
    /** Renders the manipulator for picking.
     * @see GeometryRendererBase::renderPicking
     */
    virtual void renderPicking() override;
    /** Returns the Enclosing Bounding Box in world coordinates.
     * @see GeometryRendererBase::getBoundingBox
     */
    virtual tgt::Bounds getBoundingBox() const override;

private:
    //----------------------//
    // Render Functions     //
    //----------------------//
    /** Renders the plane and anchors if enabled. */
    void renderPlane(bool useTransparency);
    /**
     * Paints the manipulator, ie. cylinder (shaft) and tips.
     * Takes into account which parts of the handle are on the other side of the plane and renders those transparent (depending on the angle to the camera).
     */
    void paintManipulator(bool useTransparency);

    //----------------------//
    // Anchors and Linesns  //
    //----------------------//
    std::vector<Anchor> anchors_;   ///< List of all anchors (clipping plane intersection points with the volume's / geometry's bounding box)
    std::vector<Line> lines_;       ///< List of all lines (connection of two anchors to a line, which is in one of the bounding box surface planes)
    GlMeshGeometryUInt16Normal* anchorGeometry_;    ///< geometry used for anchor rendering
    tgt::vec4 actualPlane_;                         ///< Actually visualized clipping plane (normal.xyz,distPosition). (Used to trigger updateAnchorsAndLines())

    /** Updates anchors and lines. */
    void updateAnchorsAndLines();
     /**
     * Add the intersection point of the given edge and plane equation,
     * if it is within the edge (ie. between start and end vertex).
     *
     * @param startVertex start of the edge
     * @param endVertex end of the edge
     * @param plane the plane equation
     */
    void addAnchor(const tgt::vec3& startVertex, const tgt::vec3& endVertex, const tgt::vec4& plane);
    /**
     * Add line connecting two anchors lying on the given plane.
     *
     * @note: All anchors must exist before the @c addLine function is called,
     *        since the anchors are checked against the given plane.
     *
     * @param plane plane equation
     */
    void addLine(const tgt::vec4& plane);

    //----------------------//
    // Manipulator          //
    //----------------------//
    GlMeshGeometryUInt16Normal* tipGeometry_;             ///< Geometry used to render the manipulator tips
    GlMeshGeometryUInt16Normal* topCylinderGeometry_;     ///< Geometry used to render the top manipulator cylinders
    GlMeshGeometryUInt16Normal* bottomCylinderGeometry_;  ///< Geometry used to render the bottom manipulator cylinders
    float manipulatorLength_;       ///< length of the manipulator (computed from scene bounding box)
    float manipulatorDiameter_;     ///< diameter of the manipulator cylinder (computed from scene bounding box)
    float manipulatorTipRadius_;    ///< radius of the manipulator tip (computed from scene bounding box)
    FloatVec3Property manipulatorCenter_;   ///< center position of the manipulator (for linking)
    tgt::vec3 topManipulatorPos_;           ///< position of the top manipulator tip
    tgt::vec3 bottomManipulatorPos_;        ///< position of the bottom manipulator tip
    tgt::vec3 manipulatorOrientation_;      ///< orientation of the manipulator (ie. normalized plane normal)
    /** Renders a manipulator tip, ie. a sphere. The shader must be activated and general uniforms must be set first. */
    void paintManipulatorTip(const tgt::vec3& center, tgt::Shader* shader);
    /** Renders the cylinder of the manipulator, ie. the shaft of the handle. The shader must be activated and general uniforms must be set first. */
    void paintManipulatorCylinder(const tgt::vec3& center, const tgt::vec3& direction, tgt::Shader* shader);
    /** Renders the top half of the manipulator cylinder. The shader must be activated and general uniforms must be set first. */
    void paintManipulatorCylinderTopHalf(const tgt::vec3& center, const tgt::vec3& direction, tgt::Shader* shader);
    /** Renders the bottom half of the manipulator cylinder. The shader must be activated and general uniforms must be set first. */
    void paintManipulatorCylinderBottomHalf(const tgt::vec3& center, const tgt::vec3& direction, tgt::Shader* shader);

    //----------------------//
    // Events and Callbacks //
    //----------------------//
    bool manipulatorTopTipIsGrabbed_;       ///< Flag used in planeManipulation() for event handling
    bool manipulatorBottomTipIsGrabbed_;    ///< Flag used in planeManipulation() for event handling
    bool manipulatorCylinderIsGrabbed_;     ///< Flag used in planeManipulation() for event handling
    bool wholeManipulatorIsGrabbed_;        ///< Flag used in planeManipulation() for event handling
    bool invalidMovement_; ///< Flag that determines if current movement is invalid (e.g. moving the handle outside of the scene)
    /**
     * Callback function for arbitrary clipping plane manipulation.
     * @param e mouse event which is translated into arbitrary clipping plane manipulation(s)
     */
    void planeManipulation(tgt::MouseEvent* e);
    /**
     * Checks if the scene bounding box changed.
     * Computes the new plane position range from the scene bounds and sets the position if necessary.
     */
    void checkInportsForBoundingBox();
    /** Updates the manipulator length, diameter, and tip radius. */
    void updateManipulatorScale();
    /** Computes the default position of the manipulator on the plane if no manipulator element is grabbed. */
    void setManipulatorToStdPos();
    /** Sets the normal to (0,0,-1) and the position to the center of the bounding box. */
    void resetPlane();
    /** Updates the tip positions using center position, and orientation of the plane manipulator.
     * @param updatePosition if true, the planePosition_ is updated. Only used if manipulaterCenter_ is linked.
     */
    void updateManipulatorElements(bool updatePosition = false);
    /** Update propertie visibility. */
    void updatePropertyVisibility();
    /** Determines, whether the plane should be rendered for the current input and property states **/
    bool shouldRenderPlane() const;

    //----------------------//
    // Helper               //
    //----------------------//
    float shiftOffset_;                     ///< offset between the mouse position and the manipulator center when shifting the handle
    float scale_;                           ///< scale factor for shifting the manipulator to avoid numerical problems with very small / large scenes

    tgt::Bounds physicalSceneBounds_;           ///< bounds of input geometry in physical coordinates
    tgt::Bounds worldSceneBounds_;              ///< bounds of input geometry in world coordinates
    tgt::mat4 physicalToWorld_;                 ///< physical-to-world matrix for the input geometry

    tgt::Bounds manipulatorBounds_;             ///< bounds for the manipulator handle center position

    //----------------------//
    // Ports and Properties //
    //----------------------//
        //Ports
    VolumePort volumeInport_;               ///< Volume inport to determine AABB dimension.
    GeometryPort geometryInport_;           ///< alternative for using volumeInport_
        //Properties
            //linking
    BoolProperty enable_;                   ///< Enable flag to be linked with the clipping processor
    BoolProperty invert_;                   ///< for linking with GeometryClipping processor, has no functionality on its own
    FloatVec3Property planeNormal_;         ///< Clipping plane normal
    FloatProperty planePosition_;           ///< Clipping plane position (distance to the world origin)
            //general
    ColorProperty generalColor_;            ///< Color of lines, anchors and the actual plane (if it is rendered during manipulation)
    ButtonProperty resetPlane_;             ///< reset the plane to normal = (0,0,-1) and the position to the center of the bounding box
            //plane
    BoolProperty renderPlane_;              ///< Determines whether the plane is shown
    FloatProperty lineWidth_;               ///< Clipping plane boundary lines width.
    FloatProperty grabbedPlaneOpacity_;     ///< opacity of the plane that is rendered while the manipulator is grabbed
    BoolProperty alwaysUseOpacity_;         ///< Always render the plane
    BoolProperty renderAnchors_;            ///< Render Anchors at clipping plane
            //manipulator
    BoolProperty renderManipulator_;        ///< determines wether the manipulator is shown
    FloatProperty manipulatorScale_;        ///< scale to resize the manipulator
    ColorProperty grabbedElementColor_;     ///< Grabbed color for manipulator parts
    ColorProperty blockedElementColor_;     ///< color of a grabbed element that may currently not be moved
    ButtonProperty resetManipulatorPos_;    ///< reset position of the manipulator to normal * position
    BoolProperty autoReset_;                ///< always reset manipulator after interaction (only if plane visible)?
            //debug
    FloatVec2Property positionRange_;                 ///< determines the range for the plane position, depending on the scene bounds
    ShaderProperty glMeshShaderOpaqueProp_;           ///< ShaderProperty for rendering the GlMeshGeometries
    ShaderProperty glMeshShaderTransparentProp_;      ///< ShaderProperty for oit-rendering the GlMeshGeometries
    ShaderProperty iModeShaderTransparentProp_;       ///< ShaderProperty for oit-rendering things via immediate mode
    ShaderProperty iModeShaderOpaqueProp_;            ///< ShaderProperty for rendering things via immediate mode
    EventProperty<PlaneManipulation>* moveEventProp_; ///< Mouse event property.

    static const std::string loggerCat_;        ///< MEOW
};

} // namespace voreen

#endif // VRN_PLANEMANIPULATION_H
