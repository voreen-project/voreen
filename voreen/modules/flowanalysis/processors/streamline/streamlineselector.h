/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2020 University of Muenster, Germany,                        *
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

#ifndef VRN_STREAMLINESELECTOR_H
#define VRN_STREAMLINESELECTOR_H

#include "voreen/core/processors/asynccomputeprocessor.h"

#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/buttonproperty.h"
#include "voreen/core/properties/optionproperty.h"
#include "voreen/core/properties/boundingboxproperty.h"
#include "voreen/core/properties/color/colorswitchproperty.h"

#include "../../ports/streamlinelistport.h"
#include "voreen/core/ports/geometryport.h"

namespace voreen {

struct StreamlineSelectorComputeInput {
    std::unique_ptr<StreamlineListBase> streamlines;
    std::function<bool(const Streamline&)> filter;
};

struct StreamlineSelectorComputeOutput {
    std::unique_ptr<StreamlineListBase> streamlines;
};

/**
 * Used to select(clip) streamlines according to a region of interest.
 * A Box can be defined, in which streamlines must/must not begin/end/intersect.
 *
 * @Note: It uses a background thread to handle changed parameters during calculation.
 */
class StreamlineSelector : public AsyncComputeProcessor<StreamlineSelectorComputeInput, StreamlineSelectorComputeOutput> {
public:
    StreamlineSelector();
    virtual ~StreamlineSelector();

    virtual Processor* create() const { return new StreamlineSelector(); }

    virtual std::string getCategory() const { return "Streamline Processing"; }
    virtual std::string getClassName() const { return "StreamlineSelector"; }
    virtual Processor::CodeState getCodeState() const { return CODE_STATE_TESTING; }

    virtual ComputeInput prepareComputeInput();
    virtual ComputeOutput compute(ComputeInput input, ProgressReporter& progressReporter) const;
    virtual void processComputeOutput(ComputeOutput output);
    virtual bool isReady() const;

protected:
    virtual void setDescriptions() {
        setDescription("Used to select(clip) streamlines according to a region of interest.  Box can be defined, " \
                       "in which streamlines must (not) begin/end/intersect.");
        //ports
        streamlineInport_.setDescription("Streamlines, which should be filtered.");
        streamlineOutport_.setDescription("Filtered streamlines.");
        geometryOutport_.setDescription("LLF and URB of the ROI. Can be used with the BoundingBoxRenderer to visualize the ROI.");
        //properties
        enabled_.setDescription("If disabled, the input is not modified.");
        inside_.setDescription("Defines, if streamlines must (not) be in the region of interest.");
        selectionMode_.setDescription("Defines, if streamlines must begin/end in the region of interest or must intersect it.");
        roi_.setDescription("The region of interest where streamlines must (not) begin/end or intersect.");
        color_.setDescription("The color of the ROI.");
    }

    virtual void adjustPropertiesToInput();
    virtual void afterProcess();
    virtual void dataWillChange(const Port* source);

private:

    /** Enum used to define intersection mode. */
    enum StreamlineSelectionMode {
        STREAM_BEGIN,
        STREAM_INTERSECT,
        STREAM_END
    };

    // ports
    StreamlineListPort streamlineInport_;
    StreamlineListPort streamlineOutport_;
    GeometryPort geometryOutport_;
    // properties
        // enable
    BoolProperty enabled_;                                  ///< toggles the processor on and off
        // config
    OptionProperty<bool> inside_;                           ///< streamlines must be (not) inside the ROI
    OptionProperty<StreamlineSelectionMode> selectionMode_; ///< streamlines must start/end/intersect ROI
    IntBoundingBoxProperty roi_;                            ///< ROI for selecting/clipping
        //roi representation (advanced)
    ColorSwitchProperty color_;                             ///< color for handling de bounding box color

    tgt::IntBounds lastUsedGeometry_;                       ///< last used geometry (in voxel space)

    static const std::string loggerCat_;
};

}   // namespace

#endif  // VRN_STREAMLINESELCETOR_H
