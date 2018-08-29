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

#include "startuplinkingprocessor.h"

namespace voreen {

StartUpLinkingProcessor::StartUpLinkingProcessor()
    : Processor()
    //property gui output
    , globalSliceAlignment_("globalSliceAlignment", "Slice Alignment", Processor::VALID)
    , globalSliceNumber_("globalSliceNumber", "Slice Number", 0, 0, 1000, Processor::VALID)
    //linked properties
        //compositors
    , compositingMode_Y_X_("compositingMode_Y_X", "To be linked 'Blend Mode' Y_X", Processor::VALID)
    , compositingMode_Z_YX_("compositingMode_Z_YX", "To be linked 'Blend Mode' Z_YX", Processor::VALID)
        //slice pos
            //x slice
    , renderXSlice_("renderXSlice", "To be linked 'Render X Slice'", false, Processor::VALID)
    , xSliceIndexProp_("xSliceIndex", "To be linked 'X Slice Number'", 0, 0, 10000, Processor::VALID)
            //y slice
    , renderYSlice_("renderYSlice", "To be linked 'Render Y Slice'", false, Processor::VALID)
    , ySliceIndexProp_("ySliceIndex", "To be linked 'Y Slice Number'", 0, 0, 10000, Processor::VALID)
            //z slice
    , renderZSlice_("renderZSlice", "To be linked 'Render Z Slice'", true, Processor::VALID)
    , zSliceIndexProp_("zSliceIndex", "To be linked 'Z Slice Number'", 0, 0, 10000, Processor::VALID)
    {
    globalSliceAlignment_.addOption("xy-plane", "XY-Plane (axial)", XY_PLANE);
    globalSliceAlignment_.addOption("xz-plane", "XZ-Plane (coronal)", XZ_PLANE);
    globalSliceAlignment_.addOption("yz-plane", "YZ-Plane (sagittal)", YZ_PLANE);
    globalSliceAlignment_.onChange(MemberFunctionCallback<StartUpLinkingProcessor>(this, &StartUpLinkingProcessor::onGlobalAlignmentChange) );
    globalSliceAlignment_.setGroupID("gui");
    addProperty(globalSliceAlignment_);

    // x, y, and z slice number changes should be set to the global slice number if the global slice alignment is set correctly
    xSliceIndexProp_.onChange(MemberFunctionCallback<StartUpLinkingProcessor>(this, &StartUpLinkingProcessor::onXNumberChange) );
    ySliceIndexProp_.onChange(MemberFunctionCallback<StartUpLinkingProcessor>(this, &StartUpLinkingProcessor::onYNumberChange) );
    zSliceIndexProp_.onChange(MemberFunctionCallback<StartUpLinkingProcessor>(this, &StartUpLinkingProcessor::onZNumberChange) );

    globalSliceNumber_.onChange(MemberFunctionCallback<StartUpLinkingProcessor>(this, &StartUpLinkingProcessor::onGlobalNumberChange) );
    globalSliceNumber_.setGroupID("gui");
    addProperty(globalSliceNumber_);

    compositingMode_Y_X_.addOption("take-first",  "Take First",  "MODE_TAKE_FIRST");
    compositingMode_Y_X_.addOption("take-second", "Take Second", "MODE_TAKE_SECOND");
    compositingMode_Y_X_.setGroupID("link");
    addProperty(compositingMode_Y_X_);
    compositingMode_Z_YX_.addOption("take-first",  "Take First",  "MODE_TAKE_FIRST");
    compositingMode_Z_YX_.addOption("take-second", "Take Second", "MODE_TAKE_SECOND");
    compositingMode_Z_YX_.setGroupID("link");
    addProperty(compositingMode_Z_YX_);

    renderXSlice_.setGroupID("link");
    addProperty(renderXSlice_);
    xSliceIndexProp_.setGroupID("link");
    addProperty(xSliceIndexProp_);
    renderYSlice_.setGroupID("link");
    addProperty(renderYSlice_);
    ySliceIndexProp_.setGroupID("link");
    addProperty(ySliceIndexProp_);
    renderZSlice_.setGroupID("link");
    addProperty(renderZSlice_);
    zSliceIndexProp_.setGroupID("link");
    addProperty(zSliceIndexProp_);

    setPropertyGroupGuiName("gui","visible in app mode");
    setPropertyGroupGuiName("link","to be linked properties");
}

StartUpLinkingProcessor::~StartUpLinkingProcessor() {
}

void StartUpLinkingProcessor::initialize() {
    Processor::initialize();
    //set property range right
    onGlobalAlignmentChange();
}

void StartUpLinkingProcessor::onGlobalAlignmentChange() {
    switch(globalSliceAlignment_.getValue()) {
    case XY_PLANE:
        //compositingMode_Y_X_.select("take-first/second");
        compositingMode_Z_YX_.select("take-first");
        renderXSlice_.set(false);
        renderYSlice_.set(false);
        renderZSlice_.set(true);
        globalSliceNumber_.setMaxValue(zSliceIndexProp_.getMaxValue());
        globalSliceNumber_.setMinValue(zSliceIndexProp_.getMinValue());
        globalSliceNumber_.set(zSliceIndexProp_.get());
        break;
    case XZ_PLANE:
        compositingMode_Y_X_.select("take-first");
        compositingMode_Z_YX_.select("take-second");
        renderXSlice_.set(false);
        renderYSlice_.set(true);
        renderZSlice_.set(false);
        globalSliceNumber_.setMaxValue(ySliceIndexProp_.getMaxValue());
        globalSliceNumber_.setMinValue(ySliceIndexProp_.getMinValue());
        globalSliceNumber_.set(ySliceIndexProp_.get());
        break;
    case YZ_PLANE:
        compositingMode_Y_X_.select("take-second");
        compositingMode_Z_YX_.select("take-second");
        renderXSlice_.set(true);
        renderYSlice_.set(false);
        renderZSlice_.set(false);
        globalSliceNumber_.setMaxValue(xSliceIndexProp_.getMaxValue());
        globalSliceNumber_.setMinValue(xSliceIndexProp_.getMinValue());
        globalSliceNumber_.set(xSliceIndexProp_.get());
        break;
    default:
        tgtAssert(false, "Unknown slice alignment");
    }
}

void StartUpLinkingProcessor::onXNumberChange() {
    if (!globalSliceAlignment_.getValue() == YZ_PLANE)
        return;

    globalSliceNumber_.set(xSliceIndexProp_.get());
}

void StartUpLinkingProcessor::onYNumberChange() {
    if (!globalSliceAlignment_.getValue() == XZ_PLANE)
        return;

    globalSliceNumber_.set(ySliceIndexProp_.get());
}

void StartUpLinkingProcessor::onZNumberChange() {
    if (!globalSliceAlignment_.getValue() == XY_PLANE)
        return;

    globalSliceNumber_.set(zSliceIndexProp_.get());
}

void StartUpLinkingProcessor::onGlobalNumberChange() {
    switch(globalSliceAlignment_.getValue()) {
    case XY_PLANE:
        zSliceIndexProp_.set(globalSliceNumber_.get());
        break;
    case XZ_PLANE:
        ySliceIndexProp_.set(globalSliceNumber_.get());
        break;
    case YZ_PLANE:
        xSliceIndexProp_.set(globalSliceNumber_.get());
        break;
    default:
        tgtAssert(false, "Unknown slice alignment");
    }
}

}   // namespace
