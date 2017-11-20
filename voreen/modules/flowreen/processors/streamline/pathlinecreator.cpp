/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2017 University of Muenster, Germany.                        *
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

#include "pathlinecreator.h"

#include "voreen/core/datastructures/volume/volumeminmaxmagnitude.h"

#include "modules/flowreen/datastructures/streamlinebundle.h"

#include "modules/flowreen/utils/pathlinecreatorbackgroundthread.h"

namespace voreen {

PathlineCreator::PathlineCreator()
    : Processor()
    // ports
    , volInport_(Port::INPORT, "volInport", "Flow Volume Input (vec3)")
    , streamlineOutport_(Port::OUTPORT, "streamlineOutport", "Pathlines Output")
    //general
    , calculatePathlinesProp_("calculatePathlineProp","Create Pathlines")
    , autoGenerateProp_("autoGenerateProp", "Auto Creation",false)
    , waitForThreadFinishedProp_("waitForThreadFinishedProp","Wait for Creation",false,Processor::INVALID_RESULT,Property::LOD_DEBUG)
    , progressProp_("progressProp","Progress:")
    // streamline config
    , maxNumPathlinesProp_("maxNumPathlinesprop", "Maximal Pathlines: ", 5000, 1, 100000)
    , streamlineLengthThresholdProp_("streamlineLengthThresholdProp", "Threshold of Pathline Length: ", tgt::ivec2(10, 100), 2, 1000)
    , absoluteMagnitudeThresholdProp_("absoluteMagnitudeThreshold", "Threshold of Magnitude (absolute)", tgt::vec2(0.0f, 1000.0f), 0.0f, 9999.99f)
    , relativeMagnitudeThresholdProp_("relativeMagnitudeThreshold", "Threshold of Magnitude (relative)", tgt::vec2(0.0f, 100.0f), 0.0f, 100.0f, Processor::VALID)
    , filterModeProp_("filterModeProp","Filtering:",Processor::INVALID_RESULT,false,Property::LOD_DEVELOPMENT)
    // streamlinebundle config
    , detectPathlineBundlesProp_("generateTubesProp", "Detect Pathline Bundles", false)
    , bundleDetectionProgressProp_("detectBundlesProgressProp", "Progress:")
    , maxAverageDistanceThresholdProp_("maxAverageDistanceThreshold", "Max. Average Distance Threshold (mm)", 1.0f, 0.0f, 100.0f)
    , minNumPathlinesPerBundleProp_("minNumPathlinesPerBundle", "Minimal number of Pathlines per Bundle (%)", 1.0f, 0.0f, 100.0f, Processor::INVALID_RESULT, NumericProperty<float>::STATIC, Property::LOD_ADVANCED)
    , resampleSizeProp_("resampleSize", "Pathline Resample Size", 20, 2, 100, Processor::INVALID_RESULT, NumericProperty<int>::STATIC, Property::LOD_ADVANCED)
    // seeds config
    , seedTimeProp_("seedTimeProp", "Current Random Seed", static_cast<int>(time(0)), std::numeric_limits<int>::min(), std::numeric_limits<int>::max(),
                    Processor::INVALID_RESULT, NumericProperty<int>::STATIC, Property::LOD_ADVANCED)
    , resetSeedProp_("resetSeedProp", "Reset Seed")
    //background thread
    , backgroundThread_(0)
    , streamlineListThreadOutput_(0)
    , dirtyFlag_(NONE)
    , backgroundThreadIsPathlineCreator_(false)
{
    //add ports
    //volInport_.onNewData(MemberFunctionCallback<PathlineCreator>(this,&PathlineCreator::volumePortHasChanged));
    static_cast<Observable<PortObserver>*>(&volInport_)->addObserver(static_cast<PortObserver*>(this));
    addPort(volInport_);
    addPort(streamlineOutport_);
    //add properties
        //general
    addProperty(calculatePathlinesProp_);
        calculatePathlinesProp_.onChange(MemberFunctionCallback<PathlineCreator>(this,&PathlineCreator::calculatePathlinesOnChange));
    addProperty(autoGenerateProp_);
    addProperty(waitForThreadFinishedProp_);
    addProperty(progressProp_);
        addProgressBar(&progressProp_);
        //streamlines
    addProperty(maxNumPathlinesProp_);
        maxNumPathlinesProp_.onChange(MemberFunctionCallback<PathlineCreator>(this,&PathlineCreator::streamlineSettingsHaveBeenChanged));
        maxNumPathlinesProp_.setGroupID("streamline");
    addProperty(streamlineLengthThresholdProp_);
        streamlineLengthThresholdProp_.onChange(MemberFunctionCallback<PathlineCreator>(this,&PathlineCreator::streamlineSettingsHaveBeenChanged));
        streamlineLengthThresholdProp_.setGroupID("streamline");
    addProperty(absoluteMagnitudeThresholdProp_);
        absoluteMagnitudeThresholdProp_.onChange(MemberFunctionCallback<PathlineCreator>(this,&PathlineCreator::streamlineSettingsHaveBeenChanged));
        absoluteMagnitudeThresholdProp_.onChange(MemberFunctionCallback<PathlineCreator>(this,&PathlineCreator::adjustRelativeThreshold));
        absoluteMagnitudeThresholdProp_.setGroupID("streamline");
    addProperty(relativeMagnitudeThresholdProp_);
        relativeMagnitudeThresholdProp_.setReadOnlyFlag(true);
        relativeMagnitudeThresholdProp_.setGroupID("streamline");
    addProperty(filterModeProp_);
        filterModeProp_.addOption("linear","Linear",LINEAR);
        filterModeProp_.addOption("nearest","Nearest",NEAREST);
        filterModeProp_.setGroupID("streamline");
        filterModeProp_.onChange(MemberFunctionCallback<PathlineCreator>(this,&PathlineCreator::streamlineSettingsHaveBeenChanged));
    setPropertyGroupGuiName("streamline", "Pathline Settings");
        //streamline bundles
    addProperty(detectPathlineBundlesProp_);
        detectPathlineBundlesProp_.onChange(MemberFunctionCallback<PathlineCreator>(this, &PathlineCreator::generateBundlesHasBeenChanged));
        detectPathlineBundlesProp_.setGroupID("streamlinebundles");
    addProperty(bundleDetectionProgressProp_);
        bundleDetectionProgressProp_.setGroupID("streamlinebundles");
    addProperty(maxAverageDistanceThresholdProp_);
        maxAverageDistanceThresholdProp_.onChange(MemberFunctionCallback<PathlineCreator>(this, &PathlineCreator::streamlineBundleSettingsHaveBeenChanged));
        maxAverageDistanceThresholdProp_.setGroupID("streamlinebundles");
    addProperty(minNumPathlinesPerBundleProp_);
        minNumPathlinesPerBundleProp_.onChange(MemberFunctionCallback<PathlineCreator>(this, &PathlineCreator::streamlineBundleSettingsHaveBeenChanged));
        minNumPathlinesPerBundleProp_.setGroupID("streamlinebundles");
    addProperty(resampleSizeProp_);
        resampleSizeProp_.onChange(MemberFunctionCallback<PathlineCreator>(this, &PathlineCreator::streamlineBundleSettingsHaveBeenChanged));
        resampleSizeProp_.setGroupID("streamlinebundles");
    setPropertyGroupGuiName("streamlinebundles", "Pathline Bundle Settings");
        //seed
    addProperty(seedTimeProp_);
        //seedTimeProp_.setReadOnlyFlag(true);
        seedTimeProp_.onChange(MemberFunctionCallback<PathlineCreator>(this,&PathlineCreator::streamlineSettingsHaveBeenChanged));
        seedTimeProp_.setGroupID("seed");
    addProperty(resetSeedProp_);
        resetSeedProp_.onChange(MemberFunctionCallback<PathlineCreator>(this,&PathlineCreator::resetSeedsOnChange));
        resetSeedProp_.setGroupID("seed");
    setPropertyGroupGuiName("seed", "Seed Configuration");

    // adjust properties
    adjustRelativeThreshold();
    generateBundlesHasBeenChanged();

}

PathlineCreator::~PathlineCreator() {
    stopBackgroundThread();
    delete streamlineListThreadOutput_;
}

void PathlineCreator::process() {
    //check, if background thread is finished
    if (backgroundThread_ && backgroundThread_->isFinished()) {
        // Retrieve error message.
        std::string message = backgroundThread_->getProgress().message_;

        // Delete thread.
        delete backgroundThread_;
        backgroundThread_ = 0;

        // reactivate recalculation
        calculatePathlinesProp_.setReadOnlyFlag(false);

        // Error handling.
        if (!message.empty()) {
            setProgress(0.0f);
            dirtyFlag_ = NONE;
            LWARNING(message);
        }
        else {
            // set progress to 100
            setProgress(1.f);

            // set the update flag for the streamline bundle detection
            if (backgroundThreadIsPathlineCreator_ && detectPathlineBundlesProp_.get())
                dirtyFlag_ = PATHLINEBUNDLES;
            else
                dirtyFlag_ = NONE;

            // send the calculated data to the outport (owning the data)
            streamlineOutport_.setData(streamlineListThreadOutput_);
            streamlineListThreadOutput_ = 0;
        }

        // Reset background thread indicator
        backgroundThreadIsPathlineCreator_ = false;
    }

    // Return, if no data is available or calculation finished.
    if (dirtyFlag_ == NONE || !volInport_.hasData())
        return;

    // Stop Thread.
    unlockMutex();
    stopBackgroundThread();
    lockMutex();

    // Clean up.
    delete streamlineListThreadOutput_; // if the outport has ownership, the pointer is 0
    streamlineListThreadOutput_ = new StreamlineList(volInport_.getData()->at(0));

    switch (dirtyFlag_) {
    case PATHLINES:
    {
        // update state
        calculatePathlinesProp_.setReadOnlyFlag(true);
        backgroundThreadIsPathlineCreator_ = true;

        // create a new thread
        backgroundThread_ = new PathlineCreatorBackgroundThread(
                    this,
                    seedTimeProp_.get(),
                    volInport_.getData(),
                    streamlineListThreadOutput_,
                    maxNumPathlinesProp_.get(),
                    streamlineLengthThresholdProp_.get(),
                    absoluteMagnitudeThresholdProp_.get(),
                    filterModeProp_.getValue()
                    );
        break;
    }
    case PATHLINEBUNDLES:
    {
        break;
    }
    default:
        tgtAssert(false, "unexpected dirty flag");
        return;
    }

    // update flag
    dirtyFlag_ = NONE;

    // Start the background thread
    backgroundThread_->run();

    // wait for background thread to finish computation
    if(waitForThreadFinishedProp_.get()) {
        unlockMutex();
        backgroundThread_->join();
        lockMutex();
    }
}

void PathlineCreator::afterProcess() {
    Processor::setValid();
    if(waitForThreadFinishedProp_.get() && backgroundThread_)
         invalidate();
}

void PathlineCreator::afterConnectionAdded(const Port* source, const Port* connectedPort) {
}

void PathlineCreator::beforeConnectionRemoved(const Port* source, const Port*) {
    stopBackgroundThread();
}

//---------------------------------------------------------------------------
//          Thread Handling
//---------------------------------------------------------------------------

void PathlineCreator::stopBackgroundThread() {
    //stop and delete thread
    if (backgroundThread_) {
        backgroundThread_->interruptAndJoin();
        delete backgroundThread_;
        backgroundThread_ = 0;
    }

    if(backgroundThreadIsPathlineCreator_) {
        // Reset data.
        streamlineOutport_.setData(0);
        // activate the calculate button
        calculatePathlinesProp_.setReadOnlyFlag(false);
        // set progress to zero
        setProgress(0.f);
    }
    bundleDetectionProgressProp_.setProgress(0.f);
}

//---------------------------------------------------------------------------
//          Callbacks
//---------------------------------------------------------------------------

void PathlineCreator::adjustRelativeThreshold() {
    //adjust read only property
    tgt::vec2 range(absoluteMagnitudeThresholdProp_.getMinValue(), absoluteMagnitudeThresholdProp_.getMaxValue());
    tgt::vec2 value = absoluteMagnitudeThresholdProp_.get();

    relativeMagnitudeThresholdProp_.set(  (value - tgt::vec2(range.x)) / tgt::vec2(range.y - range.x) * 100.f);
}

void PathlineCreator::calculatePathlinesOnChange() {
    dirtyFlag_ = PATHLINES;
}

void PathlineCreator::resetSeedsOnChange() {
    seedTimeProp_.set(time(0));
    if(autoGenerateProp_.get())
        dirtyFlag_ = PATHLINES;
}

void PathlineCreator::streamlineSettingsHaveBeenChanged() {
    if(autoGenerateProp_.get())
        dirtyFlag_ = PATHLINES;
}

void PathlineCreator::streamlineBundleSettingsHaveBeenChanged() {

    // Dont't set flag if initialization is not done.
    if (!isInitialized() || !streamlineOutport_.hasData())
        return;

    // Dont't overwrite the streamline dirty flag here
    if (detectPathlineBundlesProp_.get() && dirtyFlag_ == NONE && !backgroundThreadIsPathlineCreator_)
        dirtyFlag_ = PATHLINEBUNDLES;
}

void PathlineCreator::generateBundlesHasBeenChanged() {

    // Update UI.
    bundleDetectionProgressProp_.setVisibleFlag(detectPathlineBundlesProp_.get());
    maxAverageDistanceThresholdProp_.setVisibleFlag(detectPathlineBundlesProp_.get());
    minNumPathlinesPerBundleProp_.setVisibleFlag(detectPathlineBundlesProp_.get());
    resampleSizeProp_.setVisibleFlag(detectPathlineBundlesProp_.get());

    // Trigger a calculation of bundles, if necessary.
    streamlineBundleSettingsHaveBeenChanged();
}

}   // namespace
