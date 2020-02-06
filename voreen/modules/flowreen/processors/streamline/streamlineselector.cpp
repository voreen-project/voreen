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

#include "streamlineselector.h"

#include "../../datastructures/streamlinelist.h"

#include "voreen/core/datastructures/geometry/pointlistgeometry.h"

namespace voreen {

const std::string StreamlineSelector::loggerCat_("flowreen.StreamlineSelector");

StreamlineSelector::StreamlineSelector()
    : AsyncComputeProcessor()
    //ports
    , streamlineInport_(Port::INPORT, "streamlineInport", "Streamline Input")
    , streamlineOutport_(Port::OUTPORT, "streamlineOutport", "Streamline Output")
    , geometryOutport_(Port::OUTPORT, "geometryOutport", "ROI Output", Processor::VALID)
    //properties
    , enabled_("enableProp", "Enable", false)
        //generate
    , clearSelection_("clearSelectionProp", "Clear Selection")
        //config
    , inside_("insideProp", "Inside?")
    , selectionMode_("selectionModeProp", "Selection Mode")
    , roi_("roiProp", "Region of Interest")
        //roi settings
    , color_("usedGeometryColorProp", "ROI Color", tgt::vec4(0.f, 1.f, 0.f, 1.f), tgt::vec4(1.f, 0.f, 0.f, 1.f), Processor::VALID, Property::LOD_ADVANCED)
{
    // ports
    addPort(streamlineInport_);
    addPort(streamlineOutport_);
    addPort(geometryOutport_);

    //properties
    addProperty(enabled_);
        //general
    addProperty(clearSelection_);
        ON_CHANGE(clearSelection_, StreamlineSelector, clearSelectionOnChange);
        clearSelection_.setGroupID("select");
    setPropertyGroupGuiName("select","Start Selection");
        //config
    addProperty(inside_);
        inside_.addOption("inside" , "Select inside ROI", true);
        inside_.addOption("outside", "Select outside ROI", false);
        inside_.setGroupID("config");
    addProperty(selectionMode_);
        selectionMode_.addOption("intersect", "Line must intersect ROI", STREAM_INTERSECT);
        selectionMode_.addOption("begin"    , "Line must begin in ROI", STREAM_BEGIN);
        selectionMode_.addOption("end"      , "Line must end in ROI", STREAM_END);
        selectionMode_.setGroupID("config");
    addProperty(roi_);
        roi_.setGroupID("config");
    setPropertyGroupGuiName("config","Selection Configuration");
        //geometry output
    addProperty(color_);
        color_.setGroupID("geometry");
    setPropertyGroupGuiName("geometry","Geometry Configuration");
}

StreamlineSelector::~StreamlineSelector() {
}

StreamlineSelectorComputeInput StreamlineSelector::prepareComputeInput() {
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

    std::unique_ptr<StreamlineListBase> streamlines(streamlineInport_.getData()->clone());

    // Transform roi from voxel into world space.
    tgt::Bounds roi(streamlines->getVoxelToWorldMatrix() * tgt::vec3(roi_.get().getLLF()),
                    streamlines->getVoxelToWorldMatrix() * tgt::vec3(roi_.get().getURB()));

    bool inside = inside_.getValue();

    std::function<bool(const Streamline&)> filter;

    switch (selectionMode_.getValue()) {
    case StreamlineSelector::STREAM_BEGIN:
        filter = [=] (const Streamline& streamline) -> bool {
            return inside != roi.inside(streamline.getFirstElement().position_);
        };
        break;
    case StreamlineSelector::STREAM_INTERSECT:
        filter = [=] (const Streamline& streamline) -> bool {
            bool elementInside = false;
            size_t numElements = streamline.getNumElements();
            for (size_t k = 0; k < numElements; k++) {
                // Is any element inside the roi?
                if (roi.inside(streamline.getElementAt(k).position_)) {
                    elementInside = true;
                    break;
                }
            }
            return elementInside != inside;
        };
        break;
    case StreamlineSelector::STREAM_END:
        filter = [=] (const Streamline& streamline) -> bool {
            return inside != roi.inside(streamline.getLastElement().position_);
        };
        break;
    default:
        tgtAssert(false, "Unimplemented selection mode");
        break;
    }

    return StreamlineSelectorComputeInput {
            std::move(streamlines),
            filter
    };
}

StreamlineSelectorComputeOutput StreamlineSelector::compute(ComputeInput input, ProgressReporter& progressReporter) const {

    std::unique_ptr<StreamlineListBase> streamlines = std::move(input.streamlines);

    size_t numStreamlines = streamlines->getStreamlines().size();
    for (size_t i = 0; i < numStreamlines; i++) {
        size_t index = numStreamlines - i - 1;
        if(input.filter(streamlines->getStreamlines()[index])) {
            streamlines->removeStreamline(index);
        }
    }

    return StreamlineSelectorComputeOutput {
        std::move(streamlines)
    };
}

void StreamlineSelector::processComputeOutput(ComputeOutput output) {
    lastUsedGeometry_ = roi_.get();
    streamlineOutport_.setData(output.streamlines.release());
}

void StreamlineSelector::afterProcess() {
    AsyncComputeProcessor::afterProcess();

    if(streamlineOutport_.hasData()) {
        color_.setUseActiveColor(lastUsedGeometry_.isDefined() ? roi_.get() == lastUsedGeometry_ : false);
    } else {
        color_.setUseActiveColor(false);
    }

    if(enabled_.get()) {
        // Update bounding box geometry.
        PointListGeometryVec3* list = new PointListGeometryVec3();
        list->addPoint(roi_.get().getLLF());
        list->addPoint(roi_.get().getURB() + tgt::ivec3::one);
        list->setTransformationMatrix(streamlineInport_.getData()->getListTransformMatrix() *
                                      streamlineInport_.getData()->getOriginalVoxelToWorldMatrix());
        geometryOutport_.setData(list);
    }
    else {
        geometryOutport_.clear();
    }
}

void StreamlineSelector::adjustPropertiesToInput() {
    lastUsedGeometry_ = tgt::IntBounds();

    const StreamlineListBase* streamlines = streamlineInport_.getData();
    if(streamlines) {
        tgt::ivec3 dim = streamlines->getOriginalDimensions() - tgt::svec3::one;
        roi_.setMaxValue(dim);
        roi_.setMinValue(tgt::ivec3::zero);
        roi_.invalidate();
    }
}

void StreamlineSelector::dataWillChange(const Port* source) {
    streamlineOutport_.clear();
    AsyncComputeProcessor::dataWillChange(source);
}

void StreamlineSelector::clearSelectionOnChange() {
    if(lastUsedGeometry_.isDefined())
        roi_.set(lastUsedGeometry_);
}

}   // namespace
