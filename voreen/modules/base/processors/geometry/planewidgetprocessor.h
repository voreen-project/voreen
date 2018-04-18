/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany.                        *
 * Visualization and Computer Graphics Group <http://viscg.uni-muenster.de>        *
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

#ifndef VRN_CLIPPINGPLANEWIDGET_H
#define VRN_CLIPPINGPLANEWIDGET_H

#include "voreen/core/processors/geometryrendererbase.h"
#include "voreen/core/properties/eventproperty.h"
#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/floatproperty.h"
#include "voreen/core/properties/vectorproperty.h"
#include "voreen/core/properties/colorproperty.h"
#include "voreen/core/properties/boundingboxproperty.h"
#include "voreen/core/datastructures/geometry/trianglemeshgeometryindexed.h"
#include "voreen/core/datastructures/geometry/glmeshgeometry.h"
#include "voreen/core/ports/volumeport.h"


#include <vector>

namespace voreen {

/**
 * Renders 3D control elements to manipulate axis-aligned clipping planes.
 *
 * Use linking to connect the properties to a CubeProxyGeometry.
 *
 * @see MeshClipping
 */
class VRN_CORE_API PlaneWidgetProcessor : public GeometryRendererBase {
public:
    PlaneWidgetProcessor();
    ~PlaneWidgetProcessor();
    virtual Processor* create() const;

    virtual std::string getClassName() const    { return "PlaneWidgetProcessor"; }
    virtual std::string getCategory() const     { return "Geometry"; }
    virtual CodeState getCodeState() const      { return CODE_STATE_STABLE; }

    ///Max values of the clipping properties are set here to prevent out of range warnings.
    virtual void invalidate(int inv = INVALID_RESULT);
    virtual void render();
    virtual void renderPicking();
    virtual tgt::Bounds getBoundingBox() const;
    void setIDManager(IDManager* idm);

protected:
    virtual void setDescriptions() {
        setDescription("3D control elements to manipulate axis-aligned clipping planes. Use linking to connect the properties to a CubeProxyGeometry.\
<p>See MeshClipping and MeshClippingWidget for clipping against an arbitrarily oriented plane.</p>");
    }

    virtual std::string generateHeader(const tgt::GpuCapabilities::GlVersion* version = 0);
    virtual void initialize();
    virtual void deinitialize();

    //onChange
    void propertyVisibilityOnChange();
private:
    //onChange
    void inportNewData();

    struct Manipulator {
        Manipulator(int axis, int corner, IntBoundingBoxProperty* prop, bool inverted)
            : axis_(axis)
            , corner_(corner)
            , prop_(prop)
            , inverted_(inverted)
    {
    }

    int axis_;                      //0=x, 1=y, 2=z
    int corner_;                    //0-3
    IntBoundingBoxProperty* prop_;  //The property this arrow is associated with
    bool inverted_;                 //invert the arrows?
    };

    void planeMovement(tgt::MouseEvent* e);
    void planeMovementSync(tgt::MouseEvent* e);

    tgt::vec3 getLlf();
    tgt::vec3 getUrb();

    /**
     * @brief Helper method to get position on the edges of cube.
     *
     * @param axis Along which axis? (x=0, y=1, z=2)
     * @param num Which of the 4 edges (0-3)
     * @param t At which position (0.0 - 1.0)
     *
     * @return The Point on the edge.
     */
    tgt::vec3 getCorner(int axis, int num, float t);

    void paintManipulator(const Manipulator& a, bool useSphereManipulator, tgt::vec4 color, bool lighting);
    void paintManipulator(tgt::vec3 translation, tgt::vec3 rotationAxis, float rotationAngle, float scaleFactor, bool useSphereManipulator, tgt::vec4 color, bool lighting);

    std::vector<Manipulator> manipulators_;
    int grabbed_;                       //-1 if no arrow is grabbed, <indexOfManipulator> otherwise
    bool syncMovement_;                 //true if modifier key was pressed when grabbing an arrow


    tgt::Shader* shader_;
    TriangleMeshGeometryUInt16IndexedNormal arrowGeometry_;
    GlMeshGeometryUInt16Normal sphereGeometry_;

    BoolProperty enable_;
    BoolProperty showHandles_;
    FloatProperty manipulatorScaling_;
    ColorProperty xColor_;
    ColorProperty yColor_;
    ColorProperty zColor_;
    BoolProperty clipEnabledLeftX_;
    BoolProperty clipUseSphereManipulatorLeftX_;
    BoolProperty clipEnabledRightX_;
    BoolProperty clipUseSphereManipulatorRightX_;
    BoolProperty clipEnabledFrontY_;
    BoolProperty clipUseSphereManipulatorFrontY_;
    BoolProperty clipEnabledBackY_;
    BoolProperty clipUseSphereManipulatorBackY_;
    BoolProperty clipEnabledBottomZ_;
    BoolProperty clipUseSphereManipulatorBottomZ_;
    BoolProperty clipEnabledTopZ_;
    BoolProperty clipUseSphereManipulatorTopZ_;

    IntBoundingBoxProperty clipRegion_;

    FloatProperty width_;

    BoolProperty showInnerBB_;
    ColorProperty innerColor_;

    EventProperty<PlaneWidgetProcessor>* moveEventProp_;
    EventProperty<PlaneWidgetProcessor>* syncMoveEventProp_;

    VolumePort inport_;
};

} // namespace voreen

#endif // VRN_CLIPPINGPLANEWIDGET_H
