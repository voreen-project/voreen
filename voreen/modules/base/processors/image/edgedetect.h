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

#ifndef VRN_EDGEDETECT_H
#define VRN_EDGEDETECT_H

#include "voreen/core/processors/processor.h"
#include "voreen/core/processors/imageprocessorbypassable.h"

#include "voreen/core/properties/floatproperty.h"
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/optionproperty.h"
#include "voreen/core/properties/colorproperty.h"

namespace voreen {

/**
 * Performs an edge detection on the input image.
 *
 * The detected edge is then colored.
 */
class VRN_CORE_API EdgeDetect : public ImageProcessorBypassable {

public:
    EdgeDetect();
    virtual ~EdgeDetect();
    virtual Processor* create() const;

    virtual std::string getClassName() const { return "EdgeDetect";       }
    virtual std::string getCategory() const  { return "Image Processing"; }
    virtual CodeState getCodeState() const   { return CODE_STATE_STABLE;  }

protected:
    virtual void setDescriptions() {
        setDescription("Performs an edge detection based on the Sobel operator. The detected edge is then drawn in selectable colors, styles, blending modi etc.");
    }

    void process();
    void adjustPropVisibility();

    FloatProperty edgeThreshold_; ///< Threshold value for the gradient length to be considered as an edge
    ColorProperty backgroundColor_; ///< the color of the background
    BoolProperty showImage_; ///< Whether to show the original image or just the edge
    IntOptionProperty blendMode_; ///< Specifies how the edge is colored
    IntOptionProperty edgeStyle_; ///< The edge style can be: contour, silhouette, depth-dependent contour, color-based contour
    IntOptionProperty colorChannel_; ///< In case of color-based contours, which color channel should be used
    ColorProperty edgeColor_; ///< The color used for drawing the edges
    FloatProperty edgeThickness_; ///< the offset used to detect edges, in pixels (affects thickness of contours and silhouettes)

    RenderPort inport_;
    RenderPort outport_;
};

} // namespace voreen

#endif //VRN_EDGEDETECT_H

