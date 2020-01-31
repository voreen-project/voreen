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

#ifndef VRN_STREAMLINEFILTER_H
#define VRN_STREAMLINEFILTER_H

#include "voreen/core/processors/asynccomputeprocessor.h"

#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/numeric/intervalproperty.h"

#include "../../ports/streamlinelistport.h"
#include "voreen/core/ports/geometryport.h"

namespace voreen {

struct StreamlineFilterComputeInput {
    std::unique_ptr<StreamlineListBase> streamlines;
    tgt::vec2 physicalLengthRange;
    tgt::vec2 curvatureRange;
};

struct StreamlineFilterComputeOutput {
    std::unique_ptr<StreamlineListBase> streamlines;
};

/**
 * Used to filter streamlines according to different conditions.
 *
 * @Note: It uses a background thread to handle changed parameters during calculation.
 */
class StreamlineFilter : public AsyncComputeProcessor<StreamlineFilterComputeInput, StreamlineFilterComputeOutput> {
public:
    StreamlineFilter();
    virtual ~StreamlineFilter();

    virtual Processor* create() const { return new StreamlineFilter(); }

    virtual std::string getCategory() const { return "Streamline Processing"; }
    virtual std::string getClassName() const { return "StreamlineFilter"; }
    virtual Processor::CodeState getCodeState() const { return CODE_STATE_EXPERIMENTAL; }

    virtual ComputeInput prepareComputeInput();
    virtual ComputeOutput compute(ComputeInput input, ProgressReporter& progressReporter) const;
    virtual void processComputeOutput(ComputeOutput output);

protected:
    virtual void setDescriptions() {
        setDescription("Used to select(clip) streamlines according to a region of interest.  Box can be defined, " \
                   "in which streamlines must (not) begin/end/intersect.");
        //ports
        streamlineInport_.setDescription("Streamlines, which should be filtered.");
        streamlineOutport_.setDescription("Filtered streamlines.");
    }

    virtual void adjustPropertiesToInput();

private:

    // ports
    StreamlineListPort streamlineInport_;
    StreamlineListPort streamlineOutport_;

    // properties
    // enable
    BoolProperty enabled_;                       ///< toggles the processor on and off

    FloatIntervalProperty physicalLengthRange_;
    IntIntervalProperty curvatureRange_;

    static const std::string loggerCat_;
};

}   // namespace

#endif
