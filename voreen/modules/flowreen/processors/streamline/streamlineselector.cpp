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

#include "streamlineselector.h"

#include "../../datastructures/streamlinelist.h"

#include "voreen/core/datastructures/geometry/pointlistgeometry.h"

namespace voreen {

StreamlineSelector::StreamlineSelector()
    : Processor()
    //ports
    , streamlineInport_(Port::INPORT, "streamlineInport", "Streamline Input")
    , streamlineOutport_(Port::OUTPORT, "streamlineOutport", "Streamline Output")
    , geometryOutport_(Port::OUTPORT, "geometryOutport", "ROI Output")
    //properties
    , enableProp_("enableProp", "Enable", false)
        //generate
    , selectStreamlinesProp_("selectStreamlinesProp","Update Selection")
    , autoGenerateProp_("autoGenerateProp", "Auto Update",false)
    , clearSelectionProp_("clearSelectionProp","Clear Selection")
    , progressProp_("progressProp","Progress")
        //config
    , insideProp_("insideProp", "Inside?")
    , selectionModeProp_("selectionModeProp", "Selection Mode")
    , roiProp_("roiProp","Region of Interest")
        //roi settings
    , showGeometryProp_("showGeometryProp","Show?",true,Processor::INVALID_RESULT,Property::LOD_ADVANCED)
    , colorProp_("usedGeometryColorProp","ROI Color",tgt::vec4(0.f,1.f,0.f,1.f),tgt::vec4(1.f,0.f,0.f,1.f),Processor::INVALID_RESULT,Property::LOD_ADVANCED)
    , resetToLastUsedGeometryProp_("resetToLastUsedGeometryProp","Reset to last Used",Processor::INVALID_RESULT,Property::LOD_ADVANCED)
        //member
    , backgroundThread_(0), streamlineListThreadOutput_(0), lastSelectedList_(0), lastUsedGeometry_(0), reselectStreamlines_(false)
{
    // ports
    addPort(streamlineInport_);
        streamlineInport_.onChange(MemberFunctionCallback<StreamlineSelector>(this,&StreamlineSelector::inportHasChanged));
    addPort(streamlineOutport_);
    addPort(geometryOutport_);

    //properties
    addProperty(enableProp_);
        enableProp_.onChange(MemberFunctionCallback<StreamlineSelector>(this,&StreamlineSelector::enableOnChange));
        //general
    addProperty(selectStreamlinesProp_);
        selectStreamlinesProp_.setGroupID("select");
        selectStreamlinesProp_.onChange(MemberFunctionCallback<StreamlineSelector>(this,&StreamlineSelector::selectStreamlinesOnChange));
    addProperty(autoGenerateProp_);
        autoGenerateProp_.setGroupID("select");
    addProperty(clearSelectionProp_);
        clearSelectionProp_.onChange(MemberFunctionCallback<StreamlineSelector>(this,&StreamlineSelector::clearSelectionOnChange));
        clearSelectionProp_.setGroupID("select");
    addProperty(progressProp_);
        addProgressBar(&progressProp_);
        progressProp_.setGroupID("select");
        progressProp_.onChange(MemberFunctionCallback<StreamlineSelector>(this,&StreamlineSelector::anythingHasBeenChanged));
    setPropertyGroupGuiName("select","Start Selection");
        //config
    addProperty(insideProp_);
        insideProp_.addOption("inside" , "Select inside ROI", true);
        insideProp_.addOption("outside", "Select outside ROI", false);
        insideProp_.onChange(MemberFunctionCallback<StreamlineSelector>(this,&StreamlineSelector::anythingHasBeenChanged));
        insideProp_.setGroupID("config");
    addProperty(selectionModeProp_);
        selectionModeProp_.addOption("intersect", "Line must intersect ROI",STREAM_INTERSECT);
        selectionModeProp_.addOption("begin"    , "Line must begin in ROI",STREAM_BEGIN);
        selectionModeProp_.addOption("end"      , "Line must end in ROI",STREAM_END);
        selectionModeProp_.onChange(MemberFunctionCallback<StreamlineSelector>(this,&StreamlineSelector::anythingHasBeenChanged));
        selectionModeProp_.setGroupID("config");
    addProperty(roiProp_);
        roiProp_.setGroupID("config");
        roiProp_.onChange(MemberFunctionCallback<StreamlineSelector>(this,&StreamlineSelector::anythingHasBeenChanged));
    setPropertyGroupGuiName("config","Selection Configuration");
        //geometry output
    addProperty(showGeometryProp_);
        showGeometryProp_.setGroupID("geometry");
    addProperty(colorProp_);
        colorProp_.setGroupID("geometry");
    addProperty(resetToLastUsedGeometryProp_);
        resetToLastUsedGeometryProp_.setReadOnlyFlag(true);
        resetToLastUsedGeometryProp_.onChange(MemberFunctionCallback<StreamlineSelector>(this,&StreamlineSelector::resetToLastUsedOnChange));
        resetToLastUsedGeometryProp_.setGroupID("geometry");
    setPropertyGroupGuiName("geometry","Geometry Configuration");

    enableOnChange();
}

StreamlineSelector::~StreamlineSelector() {
    stopBackgroundThread();
    delete streamlineListThreadOutput_;
    delete lastSelectedList_;
    delete lastUsedGeometry_;
}

void StreamlineSelector::process() {
    //if disabled, pass through
    if(!enableProp_.get()) {
        stopBackgroundThread();
        streamlineOutport_.setData(streamlineInport_.getData(),false);
        geometryOutport_.setData(0);
        return;
    }

    //check, if background thread is finished
    if(backgroundThread_ && backgroundThread_->isFinished()) {
        if(!backgroundThread_->getProgress().message_.empty()) {
            LWARNING(backgroundThread_->getProgress().message_);
        }
        delete lastUsedGeometry_;
        lastUsedGeometry_ = new tgt::Bounds(backgroundThread_->roi_);
        delete backgroundThread_;
        backgroundThread_ = 0;
        //outport takes ownership
        streamlineOutport_.setData(streamlineListThreadOutput_);
        delete lastSelectedList_;   //clone will always return a list and not listbase
        lastSelectedList_ = dynamic_cast<StreamlineList*>(streamlineListThreadOutput_->clone());
        streamlineListThreadOutput_ = 0;
        // set progress to 100
        setProgress(1.f);
        // activate recalculation again
        selectStreamlinesProp_.setReadOnlyFlag(false);
    }

    //update roi output
    if(showGeometryProp_.get()) {
        PointListGeometryVec3* list = new PointListGeometryVec3();
        list->addPoint(roiProp_.get().getLLF());
        list->addPoint(roiProp_.get().getURB() + tgt::vec3::one); // +1 it is so!!!!
        list->setTransformationMatrix(streamlineInport_.getData()->getListTransformMatrix()*streamlineInport_.getData()->getOriginalVoxelToWorldMatrix());
        geometryOutport_.setData(list);
        if(lastSelectedList_) {
            colorProp_.setUseActiveColor(lastUsedGeometry_ ? roiProp_.get() == *lastUsedGeometry_ : false);
        } else {
            colorProp_.setUseActiveColor(false);
        }
    }

    // return, if nothing has changed
    if(!reselectStreamlines_) {
        if(lastSelectedList_)
            streamlineOutport_.setData(lastSelectedList_,false);
        else
            streamlineOutport_.setData(streamlineInport_.getData(),false);
        return;
    }

    // otherwise clean up
    stopBackgroundThread();
    delete streamlineListThreadOutput_; // it the outport has ownership, the pointer is 0
    reselectStreamlines_ = false;
    selectStreamlinesProp_.setReadOnlyFlag(true);
    // and start a new calculation
    streamlineListThreadOutput_ = static_cast<StreamlineList*>(streamlineInport_.getData()->clone());
    backgroundThread_ = new StreamlineSelectorBackgroundThread(this,streamlineListThreadOutput_,insideProp_.getValue(),
                                                                selectionModeProp_.getValue(), roiProp_.get());
    backgroundThread_->run();
}

//---------------------------------------------------------------------------
//          Observers
//---------------------------------------------------------------------------
void StreamlineSelector::afterConnectionAdded(const Port* source, const Port* connectedPort) {
}

void StreamlineSelector::beforeConnectionRemoved(const Port* source, const Port*) {
    stopBackgroundThread();
}

//---------------------------------------------------------------------------
//          Thread Handling
//---------------------------------------------------------------------------
void StreamlineSelector::stopBackgroundThread() {
    //stop and delete thread
    if (backgroundThread_) {
        backgroundThread_->interruptAndJoin();
        delete backgroundThread_;
        backgroundThread_ = 0;
    }
    //set outport to zero
    streamlineOutport_.setData(lastSelectedList_,false);
    // activate the calculate button
    selectStreamlinesProp_.setReadOnlyFlag(!(enableProp_.get() && streamlineInport_.hasData()));
    // set progress to zero
    setProgress(0.f);
}

//---------------------------------------------------------------------------
//          Callbacks
//---------------------------------------------------------------------------
void StreamlineSelector::inportHasChanged() {
    delete lastSelectedList_;
    lastSelectedList_ = 0;
    delete lastUsedGeometry_;
    lastUsedGeometry_ = 0;

    if(streamlineInport_.hasData()) {
        if(autoGenerateProp_.get()) reselectStreamlines_ = true;
        tgt::vec3 dim = streamlineInport_.getData()->getOriginalDimensions() - tgt::svec3::one;
        roiProp_.setMaxValue(dim);
        roiProp_.setMinValue(tgt::vec3::zero);
        lastSelectedList_ = static_cast<StreamlineList*>(streamlineInport_.getData()->clone());
        //lastUsedGeometry_ = new tgt::Bounds(tgt::vec3::zero,dim);
    }
    stopBackgroundThread();
    setProgress(1.f);
}

void StreamlineSelector::enableOnChange() {
    stopBackgroundThread();
    autoGenerateProp_.setReadOnlyFlag(!enableProp_.get());
    if(autoGenerateProp_.get()) reselectStreamlines_ = true;
}

void StreamlineSelector::selectStreamlinesOnChange() {
    stopBackgroundThread();
    reselectStreamlines_ = true;
}

void StreamlineSelector::anythingHasBeenChanged() {
    stopBackgroundThread();
    if(autoGenerateProp_.get()) reselectStreamlines_ = true;
}

void StreamlineSelector::clearSelectionOnChange() {
    inportHasChanged();
    if(lastUsedGeometry_)
        roiProp_.set(*lastUsedGeometry_);
    setProgress(1.f);
}

void StreamlineSelector::resetToLastUsedOnChange() {
    stopBackgroundThread();
    if(lastUsedGeometry_)
        roiProp_.set(*lastUsedGeometry_);
    setProgress(1.f);
}

//***********************************************************************************************************
//***********************************************************************************************************
//******               Background Thread                                                               ******
//***********************************************************************************************************
//***********************************************************************************************************
StreamlineSelectorBackgroundThread::StreamlineSelectorBackgroundThread(StreamlineSelector* processor,  StreamlineList* output, bool inside,
                                                                       StreamlineSelector::StreamlineSelectionMode selectionMode, tgt::Bounds roi)
    :ProcessorBackgroundThread<StreamlineSelector>(processor)
    , output_(output), inside_(inside), selectionMode_(selectionMode), roi_(roi)
{
    //progress is used to get Error-Messages
    setProgress("",0.f);
}

StreamlineSelectorBackgroundThread::~StreamlineSelectorBackgroundThread() {
}

void StreamlineSelectorBackgroundThread::handleInterruption() {
}

void StreamlineSelectorBackgroundThread::threadMain() {
    //return if empty
    if(output_->getStreamlines().empty()) return;
    switch (selectionMode_)
    {
    case StreamlineSelector::STREAM_BEGIN:
        selectBegin();
        break;
    case StreamlineSelector::STREAM_INTERSECT:
        selectIntersect();
        break;
    case StreamlineSelector::STREAM_END:
        selectEnd();
        break;
    default:
        //used as error massage
        setProgress("Unexpected selection mode!",0.f);
        break;
    }
}

//---------------------------------------------------------------------------
//          Helpers
//---------------------------------------------------------------------------
void StreamlineSelectorBackgroundThread::selectBegin() {

    size_t numStreamlines = output_->getStreamlines().size();
    size_t numStreamlineBundles = output_->getStreamlineBundles().size();
    size_t totalAmount = numStreamlines + numStreamlineBundles;

    size_t updateProcess = std::max<size_t>(1, std::min<size_t>(totalAmount / 9, 100));
    for (size_t i = 0, j = 0; i < numStreamlines + numStreamlineBundles; i++, j++) {
        //handle interruption and progress
        interruptionPoint();
        if(j % updateProcess == 0) {
            processor_->setProgress(0.95f * j / totalAmount);
        }

        //check, if inside for
        if (i < numStreamlines) { // a) streamlines
            bool inside = roi_.inside(output_->getStreamlines()[i].getFirstElement().position_);
            if ((inside_ ? !inside : inside)) {
                numStreamlines = output_->removeStreamline(i).size();
                i--;
                continue;
            }
        } else { // b) bundles
            size_t index = i - numStreamlines;
            bool inside = roi_.inside(output_->getStreamlineBundles()[index].getCentroid().getFirstElement().position_);
            if ((inside_ ? !inside : inside)) {
                numStreamlineBundles = output_->removeStreamlineBundle(index).size();
                i--;
                continue;
            }
        }
    }
}

void StreamlineSelectorBackgroundThread::selectEnd() {

    size_t numStreamlines = output_->getStreamlines().size();
    size_t numStreamlineBundles = output_->getStreamlineBundles().size();
    size_t totalAmount = numStreamlines + numStreamlineBundles;

    size_t updateProcess = std::max<size_t>(1, std::min<size_t>(totalAmount / 9, 100));
    for(size_t i = 0, j = 0; i < numStreamlines + numStreamlineBundles; i++,j++) {
        //handle interruption and progress
        interruptionPoint();
        if(j % updateProcess == 0) {
            processor_->setProgress(0.95f * j / totalAmount);
        }

        //check, if inside for
        if (i < numStreamlines) { // a) streamlines
            bool inside = roi_.inside(output_->getStreamlines()[i].getLastElement().position_);
            if ((inside_ ? !inside : inside)) {
                numStreamlines = output_->removeStreamline(i).size();
                i--;
                continue;
            }
        } else { // b) bundles
            size_t index = i - numStreamlines;
            bool inside = roi_.inside(output_->getStreamlineBundles()[index].getCentroid().getLastElement().position_);
            if ((inside_ ? !inside : inside)) {
                numStreamlineBundles = output_->removeStreamlineBundle(index).size();
                i--;
                continue;
            }
        }
    }
}

void StreamlineSelectorBackgroundThread::selectIntersect() {

    size_t numStreamlines = output_->getStreamlines().size();
    size_t numStreamlineBundles = output_->getStreamlineBundles().size();
    size_t totalAmount = numStreamlines + numStreamlineBundles;

    size_t updateProcess = std::max<size_t>(1, std::min<size_t>(totalAmount / 9, 100));
    for (size_t i = 0, j = 0; i < numStreamlines + numStreamlineBundles; i++, j++) {
        //handle interruption and progress
        interruptionPoint();
        if(j % updateProcess == 0) {
            processor_->setProgress(0.95f * j / totalAmount);
        }

        Streamline streamline = (i < numStreamlines) ? output_->getStreamlines()[i] : output_->getStreamlineBundles()[i - numStreamlines].getCentroid();

        //check, if inside
        bool removeStreamline;
        if(inside_) { //streamline must intersect
            removeStreamline = true;

            for(size_t k = 0; k < streamline.getNumElements(); k++) {
                bool inside = roi_.inside(streamline.getElementAt(k).position_);

                if(inside) {
                    removeStreamline = false;
                    break;
                }
            }
        } else { //streamline must not intersect
            removeStreamline = false;

            for(size_t k = 0; k < streamline.getNumElements(); k++) {
                bool inside = roi_.inside(streamline.getElementAt(k).position_);

                if(inside) {
                    removeStreamline = true;
                    break;
                }
            }
        }
        //remove streamline if needed
        if(removeStreamline) {
            if (i < numStreamlines)
                numStreamlines = output_->removeStreamline(i).size();
            else
                numStreamlineBundles = output_->removeStreamlineBundle(i - numStreamlines).size();
            i--;
        }
    }
}

}   // namespace
