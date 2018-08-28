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

#ifndef VRN_FLOWDIRECTIONOVERLAY_H
#define VRN_FLOWDIRECTIONOVERLAY_H

#include "voreen/core/processors/imageprocessor.h"

#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/filedialogproperty.h"
#include "voreen/core/properties/cameraproperty.h"
#include "voreen/core/properties/floatproperty.h"
#include "voreen/core/properties/optionproperty.h"
#include "voreen/core/properties/shaderproperty.h"
#include "voreen/core/properties/matrixproperty.h"

#include "voreen/core/datastructures/geometry/glmeshgeometry.h"

namespace voreen {

    class TriangleMeshGeometryBase;
    class TriangleMeshGeometryUInt16IndexedColorNormal;
/**
 * Renders an orientation cube or the axes to visualize streamline orientation color coding.
 */
class VRN_CORE_API FlowDirectionOverlay : public ImageProcessor {
public:
    FlowDirectionOverlay();
    ~FlowDirectionOverlay();
    virtual Processor* create() const           { return new FlowDirectionOverlay(); }
    virtual std::string getClassName() const    { return "FlowDirectionOverlay"; }
    virtual std::string getCategory() const     { return "Overlay"; }
    virtual CodeState getCodeState() const      { return CODE_STATE_STABLE; }
    /** IsReady, if the outport is connected.*/
    virtual bool isReady() const;

protected:
    virtual void setDescriptions() {
        setDescription("Adds a direction overlay for the color coding of the streamlines with respect to the camera. " \
            "Processors supporting a direction encoding like the <b>StreamlineRenderer3D</b> needed to be linked with the rotation matrix property.");
    }

    /** Nothing to do here */
    virtual void initialize();
    /** Delete current Geometry, since it is a gl object. */
    virtual void deinitialize();

    virtual void beforeProcess();
    virtual void process();


    //----------------------
    //  Geometry Handling
    //----------------------
private:

    /**
     * Sets some properties visible / hidden depending on the mode
     */
    void adjustPropertyVisibility();

    /** Renders the overlay */
    void renderGeometry();
    void createDiagonalsGeometry();
    void createDiagonalsAxesGeometry();
    void createAxesGeometry();
    void createCubeGeometry();
    void createSphereGeometry();
    void createOverlayGeometry();
    /** Helper to create an arrow. */
    TriangleMeshGeometryUInt16IndexedColorNormal* createArrowGeometry(tgt::vec4 col);

    //--------------
    //  Callbacks
    //--------------
    /** Forces the geometry to be recreated. */
    void invalidateGeometry();

    //--------------
    //  Members
    //--------------
private:
    /** Enums representing all supported overlays */
    enum FlowDirectionType{
        FD_DIAGONALS_EDGES,
        FD_DIAGONALS,
        FD_DIAGONALS_AXES,
        FD_AXES,                    ///< renders axes
        FD_COLOR_CUBE,              ///< renders a colored cube
        FD_COLOR_SPHERE,             ///< renders a colored sphere
        FD_COLOR_SPHERE_DIAGONALS,   ///< renders a colored sphere and the diagonals
        FD_COLOR_SPHERE_DIAGONALS_EDGES   ///< renders a colored sphere and the diagonals and the edges
    };

    //ports
    RenderPort inport_;             ///< Input rendering the direction overlay is drawn onto.
    RenderPort outport_;            ///< Output: input + direction overlay
    RenderPort privatePort_;        ///< stored the rendered direction overlay.
    //basic
    BoolProperty enableProp_;                               ///< used to disable overlay
    OptionProperty<FlowDirectionType> directionTypeProp_;   ///< determining the current used overlay
    CameraProperty cameraProp_;                             ///< the camera
    ShaderProperty shaderProp_;                             ///< shader property used to render geometry
    //position
    FloatProperty shiftXProp_;              ///< Distance to shift cube and axis horizontally.
    FloatProperty shiftYProp_;              ///< Distance to shift cube and axis vertically.
    FloatProperty overlaySizeProp_;         ///< Size of the overlay.

    FloatProperty sphereScaleFactor_;
    FloatProperty edgeLineWidth_;

    // the rotation matrix for the direction encoding color (not editable, only to be linked, e.g. to the StreamlineRenderer3D)
    FloatMat4Property colorRotationMatrix_;

    //geometry
    TriangleMeshGeometryBase* currentGeometry_; ///< stores the current geometry
    GlMeshGeometryBase* sphereGeometry_;        ///< sphere is handled separately

    bool geometryMustBeRecreated_;              ///< true at the beginning and if the axis size has been changed
    const float overlayBaseSize_;               ///< the base size of the overlay

    static const std::string loggerCat_; ///< category used in logging

};

} // namespace

#endif
