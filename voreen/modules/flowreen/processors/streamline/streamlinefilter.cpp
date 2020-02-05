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

#include "streamlinefilter.h"

#include "../../datastructures/streamlinelist.h"

namespace voreen {

const std::string StreamlineFilter::loggerCat_("flowreen.StreamlineFilter");

StreamlineFilter::StreamlineFilter()
    : AsyncComputeProcessor()
    //ports
    , streamlineInport_(Port::INPORT, "streamlineInport", "Streamline Input")
    , streamlineOutport_(Port::OUTPORT, "streamlineOutport", "Streamline Output")
    //properties
    , enabled_("enableProp", "Enable", false)
    , physicalLengthRange_("physicalLengthRange", "Physical Length Range (mm)", 0.0f, 0.0f, 1.0f)
    , curvatureRange_("curvatureRange", "Avg. Angle (degrees)", tgt::ivec2(0, 180), 0, 180)
{
    // ports
    addPort(streamlineInport_);
    addPort(streamlineOutport_);

    //properties
    addProperty(enabled_);
    //general
    addProperty(physicalLengthRange_);
    physicalLengthRange_.setGroupID("filter");
    addProperty(curvatureRange_);
    curvatureRange_.setGroupID("filter");
    setPropertyGroupGuiName("filter", "Filter");
}

StreamlineFilter::~StreamlineFilter() {
}

StreamlineFilterComputeInput StreamlineFilter::prepareComputeInput() {
    if(!enabled_.get()) {
        // HACK: passing through the original data causes crashes.
        StreamlineListBase* clone = streamlineInport_.hasData() ? streamlineInport_.getData()->clone() : nullptr;
        streamlineOutport_.setData(clone, true);
        //streamlineOutport_.setData(streamlineInport_.getData(), false);
        throw InvalidInputException("", InvalidInputException::S_IGNORE);
    }

//    // Pipe input through as long as calculation is not finished.
//    if(!streamlineOutport_.hasData()) {
//        streamlineOutport_.setData(streamlineInport_.getData(), false);
//    }

    if(!streamlineInport_.hasData()) {
        throw InvalidInputException("No input", InvalidInputException::S_WARNING);
    }

    std::unique_ptr<StreamlineListBase> streamlines(streamlineInport_.getData()->clone());

    return StreamlineFilterComputeInput {
            std::move(streamlines),
            physicalLengthRange_.get(),
            tgt::vec2(curvatureRange_.get()) * tgt::PIf / 180.0f
    };
}

StreamlineFilterComputeOutput StreamlineFilter::compute(ComputeInput input, ProgressReporter& progressReporter) const {

    std::unique_ptr<StreamlineListBase> streamlines = std::move(input.streamlines);

    size_t numStreamlines = streamlines->getStreamlines().size();
    for (size_t i = 0; i < numStreamlines; i++) {
        size_t index = numStreamlines - i - 1;
        const Streamline& streamline = streamlines->getStreamlines()[index];

        bool remove = false;

        if(streamline.getPhysicalLength() < input.physicalLengthRange.x ||
            streamline.getPhysicalLength() > input.physicalLengthRange.y) {
            remove = true;
        }

        if(!remove && (streamline.getCurvatureStatistics().getMean() < input.curvatureRange.x ||
            streamline.getCurvatureStatistics().getMean() > input.curvatureRange.y)) {
            remove = true;
        }

        if(remove) {
            streamlines->removeStreamline(index);
        }
    }

    return StreamlineFilterComputeOutput {
            std::move(streamlines)
    };
}

void StreamlineFilter::processComputeOutput(ComputeOutput output) {
    streamlineOutport_.setData(output.streamlines.release());
}

void StreamlineFilter::adjustPropertiesToInput() {
    const StreamlineListBase* streamlines = streamlineInport_.getData();
    if(streamlines) {

        float maxPhysicalLength = 0.0f;
        for(const Streamline& streamline : streamlines->getStreamlines()) {
            maxPhysicalLength = std::max(maxPhysicalLength, streamline.getPhysicalLength());
        }

        physicalLengthRange_.setMaxValue(maxPhysicalLength);

        // Do not update currently set values, because we might
        // just want to recalculate pathlines on the same dataset.
        //physicalLengthRange_.set(maxPhysicalLength);
        //curvatureRange_.set(tgt::ivec2(0, 180));
    }
}

}   // namespace
