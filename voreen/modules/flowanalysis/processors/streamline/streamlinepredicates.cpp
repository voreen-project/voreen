/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2021 University of Muenster, Germany,                        *
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

#include "streamlinepredicates.h"

#include "../../datastructures/streamlinelist.h"

namespace voreen {

const std::string StreamlinePredicates::loggerCat_("flowanalysis.StreamlineFilter");

StreamlinePredicates::StreamlinePredicates()
    : AsyncComputeProcessor()
    //ports
    , streamlineInport_(Port::INPORT, "streamlineInport", "Streamline Input")
    , streamlinePredicateVolumeInport_(Port::INPORT, "predicatePort", "Volume Predicate")
    , streamlineOutport_(Port::OUTPORT, "streamlineOutport", "Streamline Output")
    //properties
    , enabled_("enableProp", "Enable", false)
    , predicateVolumeFiltering_("predicateFiltering", "Predicate Filtering")
    , physicalLengthRange_("physicalLengthRange", "Physical Length Range (mm)", 0.0f, 0.0f, 1.0f)
    , curvatureRange_("curvatureRange", "Avg. Angle (degrees)", tgt::ivec2(0, 180), 0, 180)
{
    // ports
    addPort(streamlineInport_);
    addPort(streamlinePredicateVolumeInport_);
    addPort(streamlineOutport_);

    //properties
    addProperty(enabled_);
    //general
    addProperty(predicateVolumeFiltering_);
        predicateVolumeFiltering_.addOption("disable", "Disabled");
        predicateVolumeFiltering_.addOption("intersection", "Any Element");
        predicateVolumeFiltering_.addOption("none", "No Element");
        predicateVolumeFiltering_.addOption("all", "All Elements");
        predicateVolumeFiltering_.setGroupID("filter");
    addProperty(physicalLengthRange_);
        physicalLengthRange_.setTracking(false);
        physicalLengthRange_.setGroupID("filter");
    // Experimental:
    //addProperty(curvatureRange_);
    //    curvatureRange_.setTracking(false);
    //    curvatureRange_.setGroupID("filter");
    setPropertyGroupGuiName("filter", "Filter");
}

StreamlinePredicates::~StreamlinePredicates() {
}

bool StreamlinePredicates::isReady() const {
    if(!streamlineInport_.isReady()) {
        setNotReadyErrorMessage("No streamline input");
        return false;
    }

    // Note: predicate volume is optional.

    if(!streamlineOutport_.isReady()) {
        setNotReadyErrorMessage("No output connected");
        return false;
    }

    return true;
}

StreamlinePredicatesComputeInput StreamlinePredicates::prepareComputeInput() {
    if(!enabled_.get()) {
        streamlineOutport_.setData(streamlineInport_.getData(), false);
        throw InvalidInputException("", InvalidInputException::S_IGNORE);
    }

//    // Pipe input through as long as calculation is not finished.
//    if(!streamlineOutport_.hasData()) {
//        streamlineOutport_.setData(streamlineInport_.getData(), false);
//    }

    if(!streamlineInport_.hasData()) {
        throw InvalidInputException("No input", InvalidInputException::S_WARNING);
    }

    auto streamlinePredicateVolume = streamlinePredicateVolumeInport_.getThreadSafeData();

    using PredicateFunc = std::function<bool(const tgt::vec3&)>;

    std::function<bool(const Streamline&, const PredicateFunc&)> predicateFilter;
    if(!streamlinePredicateVolume || predicateVolumeFiltering_.get() == "disable") {
        predicateFilter = [] (const Streamline&, const PredicateFunc&) -> bool { return false; };
    }
    else if(predicateVolumeFiltering_.get() == "intersection") {
        predicateFilter = [] (const Streamline& streamline, const PredicateFunc& predicate) -> bool {
            for(size_t i=0; i<streamline.getNumElements(); i++) {
                if(predicate(streamline.getElementAt(i).position_)) {
                    return false;
                }
            }
            return true;
        };
    }
    else if(predicateVolumeFiltering_.get() == "none") {
        predicateFilter = [] (const Streamline& streamline, const PredicateFunc& predicate) -> bool {
            for(size_t i=0; i<streamline.getNumElements(); i++) {
                if(predicate(streamline.getElementAt(i).position_)) {
                    return true;
                }
            }
            return false;
        };
    }
    else if(predicateVolumeFiltering_.get() == "all") {
        predicateFilter = [] (const Streamline& streamline, const PredicateFunc& predicate) -> bool {
            for(size_t i=0; i<streamline.getNumElements(); i++) {
                if(!predicate(streamline.getElementAt(i).position_)) {
                    return true;
                }
            }
            return false;
        };
    }
    else {
        tgtAssert(false, "Unimplemented predicate");
    }

    std::unique_ptr<StreamlineListBase> streamlines(streamlineInport_.getData()->clone());

    return StreamlinePredicatesComputeInput {
            std::move(streamlines),
            std::move(streamlinePredicateVolume),
            std::move(predicateFilter),
            physicalLengthRange_.get(),
            tgt::vec2(curvatureRange_.get()) * tgt::PIf / 180.0f
    };
}

StreamlinePredicatesComputeOutput StreamlinePredicates::compute(ComputeInput input, ProgressReporter& progressReporter) const {

    std::unique_ptr<StreamlineListBase> streamlines = std::move(input.streamlines);

    auto streamlinePredicateVolume = std::move(input.streamlinePredicateVolume);

    std::function<bool(const tgt::vec3&)> predicate = [] (const tgt::vec3&) -> bool { return true; };
    if(streamlinePredicateVolume) {
        tgt::mat4 worldToVoxel = streamlinePredicateVolume->getWorldToVoxelMatrix();
        tgt::Bounds bounds = streamlinePredicateVolume->getBoundingBox().getBoundingBox();
        VolumeRAMRepresentationLock lock(streamlinePredicateVolume);
        predicate = [=] (const tgt::vec3& position) -> bool {
            if(!bounds.containsPoint(position)) {
                return false;
            }
            return lock->getVoxelNormalized(worldToVoxel * position) != 0.0f;
        };
    }

    size_t numStreamlines = streamlines->getStreamlines().size();
    for (size_t i = 0; i < numStreamlines; i++) {
        size_t index = numStreamlines - i - 1;
        const Streamline& streamline = streamlines->getStreamlines()[index];

        // Check predicate first.
        bool remove = input.predicateFilter(streamline, predicate);

        if(!remove && (streamline.getPhysicalLength() < input.physicalLengthRange.x ||
            streamline.getPhysicalLength() > input.physicalLengthRange.y)) {
            remove = true;
        }
/*
        if(!remove && (streamline.getCurvatureStatistics().getMean() < input.curvatureRange.x ||
            streamline.getCurvatureStatistics().getMean() > input.curvatureRange.y)) {
            remove = true;
        }
*/
        if(remove) {
            streamlines->removeStreamline(index);
        }
    }

    return StreamlinePredicatesComputeOutput {
            std::move(streamlines)
    };
}

void StreamlinePredicates::processComputeOutput(ComputeOutput output) {
    streamlineOutport_.setData(output.streamlines.release());
}

void StreamlinePredicates::adjustPropertiesToInput() {
    const StreamlineListBase* streamlines = streamlineInport_.getData();
    if(streamlines) {

        predicateVolumeFiltering_.setReadOnlyFlag(!streamlinePredicateVolumeInport_.hasData());

        float maxPhysicalLength = 0.0f;
        for(const Streamline& streamline : streamlines->getStreamlines()) {
            maxPhysicalLength = std::max(maxPhysicalLength, streamline.getPhysicalLength());
        }

        physicalLengthRange_.setMaxValue(maxPhysicalLength);
    }
}

void StreamlinePredicates::dataWillChange(const Port* source) {
    streamlineOutport_.clear();
    AsyncComputeProcessor::dataWillChange(source);
}

}   // namespace
