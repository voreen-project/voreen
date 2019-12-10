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

#ifndef VRN_CONNECTEDCOMPONENTS2D_H
#define VRN_CONNECTEDCOMPONENTS2D_H

#include "voreen/core/processors/imageprocessor.h"

#include "voreen/core/properties/optionproperty.h"
#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/boolproperty.h"

namespace voreen {

/**
 * Detects connected components in a 2D image and luminance-codes the assigned labels
 * in a 16 bit output image.
 *
 * @see ConnectedComponents3D
 */
class ConnectedComponents2D : public ImageProcessor {
public:
    ConnectedComponents2D();
    virtual Processor* create() const;

    virtual std::string getClassName() const    { return "ConnectedComponents2D"; }
    virtual std::string getCategory() const     { return "Image Processing"; }
    virtual CodeState getCodeState() const      { return CODE_STATE_STABLE; }

protected:
    virtual void setDescriptions() {
        setDescription("Detects connected components in a 2D image and luminance-codes the assigned labels in a 16 bit output image.\
<p><strong>Properties</strong>:\
<ul>\
 <li>Channel: channel(s) to analyze, pixels with a value greater zero in any of the selected channels are regarded as foreground.</li>\
 <li>Connectivity: pixel neighborhood to consider.</li>\
 <li>Min Component Size: components consisting of less pixels are discarded.</li>\
 <li>Max Components: maximal number of connected components that are put out.</li>\
 <li>Component Sorting: determines in which order the component labels are assigned.</li>\
 <li>Binarize Output: all connected components are assigned the same label (max intensity).</li>\
 <li>Stretch Labels: component labels are stretched over the whole intensity range of the output image in order to maximize visual contrast.</li>\
</ul></p>\
<p>See: ConnectedComponents3D</p>");
    }

    virtual void process();

    virtual std::string generateHeader(const tgt::GpuCapabilities::GlVersion* version = 0);
    static const std::string loggerCat_;

private:
    RenderPort inport_;     ///< Image to analyze.
    RenderPort outport_;    ///< Output image storing the components' labels.

    IntOptionProperty channel_;             ///< Channel(s) to analyze: pixels with a value greater zero in any
                                            ///  of the selected channels are regarded as foreground.
    IntOptionProperty connectivity_;        ///< Pixel neighborhood to consider for the analysis.
    IntProperty minComponentSize_;          ///< Components consisting of less pixels are discarded.
    IntProperty maxComponents_;             ///< The maximal number of connected components that are put out.
    StringOptionProperty componentSorting_; ///< Determines in which order the component labels are assigned
    BoolProperty binarizeOutput_;           ///< If enabled, all connected components are assigned the same label (max intensity).
    BoolProperty stretchLabels_;            ///< If enabled, component labels are stretched over the whole intensity range
                                            ///  of the output volume in order to maximize contrast.
};

} // namespace

#endif // VRN_CONNECTEDCOMPONENTS2D_H
