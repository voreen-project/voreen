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

#ifndef VRN_BOUNDINGBOXRENDERER_H
#define VRN_BOUNDINGBOXRENDERER_H

#include "voreen/core/processors/geometryrendererbase.h"

#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/floatproperty.h"
#include "voreen/core/properties/color/colorswitchproperty.h"
#include "voreen/core/properties/fontproperty.h"

#include "voreen/core/ports/volumeport.h"
#include "voreen/core/ports/geometryport.h"

namespace voreen {

///Draws bounding box around the data set
class VRN_CORE_API BoundingBoxRenderer : public GeometryRendererBase {
public:
    BoundingBoxRenderer();
    virtual Processor* create() const;

    virtual std::string getClassName() const { return "BoundingBoxRenderer"; }
    virtual std::string getCategory() const  { return "Geometry"; }
    virtual CodeState getCodeState() const   { return CODE_STATE_STABLE; }

    /// Returns true, if one of the inports is connected.
    virtual bool isReady() const;

protected:
    virtual void setDescriptions() {
        setDescription("Draws a bounding box around the input volume or geometry, depending on which inport is connected, and allows to show a grid behind the volume. For volumes, the axis-aligned bounding box is computed in physical coordinates and then transformed into world-coordinates. For geometry, the axis-aligned bounding is computed directly in world coordinates.");
    }

    virtual tgt::Bounds getBoundingBox() const;
    virtual void render();

    void renderLegend();

    /**
     * Initializes all private RenderPorts.
     */
    virtual void initialize();

    /**
     * Deinitializes all private RenderPorts.
     */
    virtual void deinitialize();

private:

    // Layout Options for
    enum AlignmentOptions {
        OPTION_N,
        OPTION_NE,
        OPTION_E,
        OPTION_SE,
        OPTION_S,
        OPTION_SW,
        OPTION_W,
        OPTION_NW,
        OPTION_CENTER,
    };
    friend class OptionProperty<AlignmentOptions>;

    //------------------
    // Callbacks
    //------------------
    void onInputDataChange();       ///< Invoked by volume and geometry inport
    void onTilesSectionChange();    ///< Invoked by tilesSections_

    /// ports
    VolumePort volumeInport_;
    GeometryPort geometryInport_;
    /// enable
    BoolProperty enable_;
    /// line
    FloatProperty width_;
    ColorSwitchProperty bboxColor_;
    /// grid
    BoolProperty showGrid_;
    IntProperty tilesSections_;
    FloatProperty tilesOpacity_;
    /// legend
    BoolProperty showLegend_;
    FontProperty fontProp_;
    FloatProperty legendScaleFactor_;
    OptionProperty<AlignmentOptions> alignment_;
    IntVec2Property legendDisplacement_;

    tgt::mat4 inputToWorldTransformation_;
    GLuint occlusionQueries_[6];    ///< Occlusion query objects to determine visible tiles

    float tileLength_;              ///< Actual tile length calculated in order to achieve round numers
    float maxBoxLength_;            ///< Length of the longest side of the box

    tgt::Bounds boundingBox_;

    static const std::string loggerCat_;
};

}

#endif

