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

#ifndef VRN_GEOMETRYCLIPPINGWIDGET_H
#define VRN_GEOMETRYCLIPPINGWIDGET_H

#include "voreen/core/processors/geometryrendererbase.h"
#include "voreen/core/properties/eventproperty.h"
#include "voreen/core/properties/vectorproperty.h"
#include "voreen/core/properties/floatproperty.h"
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/colorproperty.h"

#include "voreen/core/ports/geometryport.h"

#include <vector>

namespace voreen {

/**
 * Renders 3D-control elements to manipulate an arbitrary oriented clipping plane
 * using linking to propagate property changes to @c GeometryClipping processors.
 */
class VRN_CORE_API GeometryClippingWidget : public GeometryRendererBase {
public:
    GeometryClippingWidget();
    ~GeometryClippingWidget();
    virtual Processor* create() const;

    virtual std::string getClassName() const  { return "GeometryClippingWidget"; }
    virtual std::string getCategory() const   { return "Geometry";               }
    virtual CodeState getCodeState() const    { return CODE_STATE_OBSOLETE;      }

protected:
    virtual void setDescriptions() {
        setDescription("3D control elements to manipulate an arbitrarily oriented geometry clipping plane. "
                       "<p><strong>DEPRECATED</strong>: Use PlaneManipulation instead!");
    }

    virtual void initialize();

    /// @see GeometryRendererBase::render
    virtual void render();

    /// @see GeometryRendererBase::renderPicking
    virtual void renderPicking();

    /**
     * Initialize ID manager to work with this @c MeshClippingWidget instance.
     *
     * @param idm ID manager
     */
    void setIDManager(IDManager* idm);

private:
    /**
     * Accuracy for floating point operations within this class.
     */
    static const float EPS;

    /**
     * Callback function for arbitrary clipping plane manipulation.
     *
     * @param e mouse event which is translated into arbitrary clipping plane manipulation(s)
     */
    void planeManipulation(tgt::MouseEvent* e);

    /**
     * Definition of @c Anchor data structure.
     */
    typedef tgt::vec3 Anchor;

    /**
     * Definition of @c Line data structure.
     */
    typedef std::pair<Anchor*, Anchor*> Line;

    /**
     * List of all anchors (clipping plane intersection points with the volume's
     * axis aligned bounding box (AABB)).
     */
    std::vector<Anchor> anchors_;

    /**
     * List of all lines (connection of two anchors to a line, which is in at least one
     * of the volume's axis aligned bounding box (AABB) surface planes)
     */
    std::vector<Line> lines_;

    /**
     * Renders the given anchor.
     *
     * @param anchor Anchor
     */
    void paintAnchor(const Anchor& anchor);

    /**
     * Renders the given line.
     *
     * @param line Line
     */
    void paintLine(const Line& line);

    /**
     * Returns the distance of the given point to the given line in vector equation form.
     *
     * @param point the Point
     * @param linePositionVector position vector of line vector equation
     * @param lineDirectionVector direction vector of line vector equation
     *
     * @return distance between point and line
     */
    float getPointLineDistance(const tgt::vec3& point, const tgt::vec3& linePositionVector, const tgt::vec3& lineDirectionVector) const;

    /**
     * Return the given line vector equation scalar for the projection of the given point.
     *
     * @param point the Point
     * @param linePositionVector position vector of line vector equation
     * @param lineDirectionVector direction vector of line vector equation
     *
     * @return distance between point and line
     */
    float getPointLineProjectionScalar(const tgt::vec2& point, const tgt::vec2& linePositionVector, const tgt::vec2& lineDirectionVector) const;

    /**
     * Determines the intersection point of the given line vector equation and plane equation.
     *
     * If the intersection point is outside of the volume's axis aligned bounding box (AABB),
     * the intersection point is omitted.
     *
     * @param intersectionPoint the variable in which the intersection point should be stored
     * @param linePositionVector position vector of line vector equation
     * @param lineDirectionVector direction vector of line vector equation
     * @param plane the plane equation
     *
     * @return @c true if an intersection point within the AABB is found, otherwise @c false.
     */
    bool getIntersectionPointInAABB(tgt::vec3& intersectionPoint, const tgt::vec3& linePositionVector, const tgt::vec3& lineDirectionVector, const tgt::vec4& plane);

    /**
     * Add the intersection point of the given line vector equation and plane equation,
     * if it is within the axis aligned bounding box (AABB).
     *
     * @param linePositionVector position vector of line vector equation
     * @param lineDirectionVector direction vector of line vector equation
     * @param plane the plane equation
     */
    void addAnchor(const tgt::vec3& linePositionVector, const tgt::vec3& lineDirectionVector, const tgt::vec4& plane);

    /**
     * Add line connecting two anchors lying on the given axis aligned bounding box (AABB) plane.
     *
     * @note: All anchors must exist before the @c addLine function is called,
     *        since the anchors are checked against the given AABB plane.
     *
     * @param planeAABB AABB plane equation
     */
    void addLine(const tgt::vec4& planeAABB);

    /**
     * Updates anchors and lines.
     *
     * Synchronize anchor and line state with ID manager, recreate anchors and lines,
     * and order them for clipping plane visualization on clipping plane manipulation.
     */
    void updateAnchorsAndLines();

    GeometryPort inport_;                 ///< Volume inport to determine AABB dimension.

    BoolProperty renderGeometry_;       ///< Determines whether the widget is shown.

    FloatVec3Property planeNormal_;     ///< Clipping plane normal
    FloatProperty planePosition_;       ///< Clipping plane position (distance to the world origin)

    FloatVec4Property lightPosition_;       ///< Light property.
    ColorProperty planeColor_;          ///< Clipping plane, which is visible on clipping plane manipulation, color.
    ColorProperty anchorColor_;         ///< Anchors color.
    ColorProperty grabbedAnchorColor_;  ///< Grabbed anchor color.
    ColorProperty lineColor_;           ///< Clipping plane boundary lines color.
    FloatProperty width_;               ///< Clipping plane boundary lines width.

    EventProperty<GeometryClippingWidget>* moveEventProp_;  ///< Mouse event property.

    tgt::vec4 actualPlane_;         ///< Actually visualized clipping plane equation.
    bool isAnchorGrabbed_;          ///< Flag to determine whether an anchor is grabbed or not.
    Anchor grabbedAnchor_;          ///< Grabbed anchor.
    Anchor grabbedSecondAnchor_;    ///< The first anchor in the anchor list which is not the grabbed anchor.
    Anchor grabbedThirdAnchor_;     ///< The second anchor in the anchor list which is not the grabbed anchor.

    /**
     * Flag to determine whether the new anchor plane direction vector is the first
     * parameter in cross product calculation call.
     */
    bool crossNewAnchorFirst_;

    bool isLineGrabbed_;                ///< Flag to determine whether a line is grabbed or not.
    Anchor grabbedLineBegin_;           ///< Start anchor/point of the grabbed line.
    Anchor grabbedLineEnd_;             ///< End anchor/point of the grabbed line.
    tgt::vec3 grabbedLinePoint_;        ///< Position vector of the clipping plane movement line vector equation.
    tgt::vec3 grabbedLineDirection_;    ///< Direction vector of the clipping plane movement line vector equation.
    float grabbedPlaneOriginOffset_;    ///< Clipping plane offset at the moment the actual grabbed line was grabbed.

};

} // namespace voreen

#endif // VRN_GEOMETRYCLIPPINGWIDGET_H
