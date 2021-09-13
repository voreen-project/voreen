/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2014 University of Muenster, Germany.                        *
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

#ifndef VRN_COSMOLOGYORIENTATIONOVERLAY_H
#define VRN_COSMOLOGYORIENTATIONOVERLAY_H

#include "voreen/core/processors/imageprocessor.h"

#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/filedialogproperty.h"
#include "voreen/core/properties/cameraproperty.h"
#include "voreen/core/properties/floatproperty.h"
#include "voreen/core/properties/optionproperty.h"
#include "voreen/core/properties/shaderproperty.h"
#include "voreen/core/properties/matrixproperty.h"

#include "tgt/texture.h"

namespace voreen {

    class TriangleMeshGeometryBase;
    class TriangleMeshGeometryUInt16IndexedColorNormal;
/**
 * Renders an orientation cube or the axes to visualize streamline orientation color coding.
 */
class VRN_CORE_API CosmologyOrientationOverlay : public ImageProcessor {
public:
    CosmologyOrientationOverlay();
    ~CosmologyOrientationOverlay();
    virtual Processor* create() const;

    virtual std::string getClassName() const    { return "CosmologyOrientationOverlay"; }
    virtual std::string getCategory() const     { return "Image Processing"; }
    virtual CodeState getCodeState() const      { return CODE_STATE_EXPERIMENTAL; }

    virtual bool isReady() const;

protected:
    virtual void setDescriptions() {
        setDescription("Adds an orientation overlay for the color coding of the streamlines with respect to the camera.");
    }

    virtual void beforeProcess();
    virtual void process();
    virtual void initialize();
    virtual void deinitialize();

private:
    /** Renders the axis overlay */
    void renderGeometry();
    void createDiagonalsGeometry();
    void createDiagonalsAxesGeometry();
    void createAxesGeometry();
    void createCubeGeometry();
    void createOverlayGeometry();

    TriangleMeshGeometryUInt16IndexedColorNormal* createArrowGeometry(tgt::vec4 col);

    /** Adjusts the visibility of all properties according to the current orientation type */
    void adjustPropertyVisibility();
    /** Forces the axis geometry to be recreated. */
    void invalidateGeometry();
    
    /// Invalidates the geometry (and the processor)
    void rotationChanged();
private:
    /** Enums representing all supported overlays */
    enum CosmolpogyOrientationType{
        OT_DIAGONALS_EDGES,
        OT_DIAGONALS,
        OT_DIAGONALS_AXES,
        OT_AXES,                ///< renders axes
        OT_COLOR_CUBE,          ///< renders a colored cube 
    };
    //ports
    RenderPort inport_;             ///< Input rendering the orientation overlay is drawn onto.
    RenderPort outport_;            ///< Output: input + orientation overlay
    RenderPort privatePort_;        ///< stored the rendered orientation overlay.
    //basic
    BoolProperty enableProp_;                               ///< used to disable overlay
    BoolProperty enableProp2_;                               ///< used to disable overlay
    OptionProperty<CosmolpogyOrientationType> orientationTypeProp_;   ///< determining the current used overlay
    CameraProperty cameraProp_;
    //position
    FloatProperty shiftXProp_;              ///< Distance to shift cube and axis horizontally.
    FloatProperty shiftYProp_;              ///< Distance to shift cube and axis vertically.
    //axis
    FloatProperty axisSizeProp_;            ///< Length of axes indicating orientation.
    //cube
    FloatProperty cubeSizeProp_;            ///< Size of cube indicating orientation.
    //shader
    ShaderProperty shaderProp_;             ///< shader property used to render geometry

    // the rotation matrix for the direction encoding color cube (not editable, only to be linked, e.g. from StreamlineRenderer3D)
    FloatMat4Property colorCubeRotation_;

    TriangleMeshGeometryBase* currentGeometry_; ///< stores the current geometry
    bool geometryMustBeRecreated_;          ///< true at the beginning and if the axis size has been changed
    const float axisBaseLength_;        ///< the base length of the axis
    const float cubeBaseLength_;        ///< the base length of the cube

    static const std::string loggerCat_; ///< category used in logging

};

} // namespace

#endif
