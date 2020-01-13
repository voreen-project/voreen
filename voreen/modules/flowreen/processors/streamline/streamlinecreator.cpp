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

#include "streamlinecreator.h"

#include "voreen/core/ports/conditions/portconditionvolumetype.h"
#include "voreen/core/datastructures/volume/volumeminmaxmagnitude.h"

#include "modules/flowreen/datastructures/streamlinebundle.h"

#include "modules/flowreen/utils/streamlinecreatorbackgroundthread.h"
#include "modules/flowreen/utils/streamlinebundledetectorbackgroundthread.h"

namespace voreen {

StreamlineCreator::StreamlineCreator()
    : Processor()
    // ports
    , volInport_(Port::INPORT, "volInport", "Flow Volume Input (vec3)")
    , streamlineOutport_(Port::OUTPORT, "streamlineOutport", "Streamlines Output")
    //general
    , calculateStreamlinesProp_("calculateStreamlineProp","Create Streamlines")
    , autoGenerateProp_("autoGenerateProp", "Auto Creation",false)
    , waitForThreadFinishedProp_("waitForThreadFinishedProp","Wait for Creation",false,Processor::INVALID_RESULT,Property::LOD_DEBUG)
    , progressProp_("progressProp","Progress:")
    // streamline config
    , maxNumStreamlinesProp_("maxNumStreamlinesprop", "Maximal Streamlines: ", 5000, 1, 100000)
    , streamlineLengthThresholdProp_("streamlineLengthThresholdProp", "Threshold of Streamline Length: ", tgt::ivec2(10, 100), 2, 1000)
    , absoluteMagnitudeThresholdProp_("absoluteMagnitudeThreshold", "Threshold of Magnitude (absolute)", tgt::vec2(0.0f, 1000.0f), 0.0f, 9999.99f)
    , fitAbsoluteMagnitudeProp_("fitAbsoluteMagnitude", "Fit absolute Threshold to Input", false)
    , relativeMagnitudeThresholdProp_("relativeMagnitudeThreshold", "Threshold of Magnitude (relative)", tgt::vec2(0.0f, 100.0f), 0.0f, 100.0f, Processor::VALID)
    , filterModeProp_("filterModeProp","Filtering:",Processor::INVALID_RESULT,false,Property::LOD_DEVELOPMENT)
    // streamlinebundle config
    , detectStreamlineBundlesProp_("generateTubesProp", "Detect Streamline Bundles", false)
    , bundleDetectionProgressProp_("detectBundlesProgressProp", "Progress:")
    , maxAverageDistanceThresholdProp_("maxAverageDistanceThreshold", "Max. Average Distance Threshold (mm)", 1.0f, 0.0f, 100.0f)
    , minNumStreamlinesPerBundleProp_("minNumStreamlinesPerBundle", "Minimal number of Streamlines per Bundle (%)", 1.0f, 0.0f, 100.0f, Processor::INVALID_RESULT, NumericProperty<float>::STATIC, Property::LOD_ADVANCED)
    , resampleSizeProp_("resampleSize", "Streamline Resample Size", 20, 2, 100, Processor::INVALID_RESULT, NumericProperty<int>::STATIC, Property::LOD_ADVANCED)
    // seeds config
    , seedTimeProp_("seedTimeProp", "Current Random Seed", static_cast<int>(time(0)), std::numeric_limits<int>::min(), std::numeric_limits<int>::max(),
                    Processor::INVALID_RESULT, NumericProperty<int>::STATIC, Property::LOD_ADVANCED)
    , resetSeedProp_("resetSeedProp", "Reset Seed")
    //background thread
    , backgroundThread_(0)
    , streamlineListThreadOutput_(0)
    , dirtyFlag_(NONE)
    , backgroundThreadIsStreamlineCreator_(false)
{
    //add ports
    volInport_.addCondition(new PortConditionVolumeChannelCount(3));
    volInport_.onNewData(MemberFunctionCallback<StreamlineCreator>(this,&StreamlineCreator::volumePortHasChanged));
    static_cast<Observable<PortObserver>*>(&volInport_)->addObserver(static_cast<PortObserver*>(this));
    addPort(volInport_);
    addPort(streamlineOutport_);
    //add properties
        //general
    addProperty(calculateStreamlinesProp_);
        calculateStreamlinesProp_.onChange(MemberFunctionCallback<StreamlineCreator>(this,&StreamlineCreator::calculateStreamlinesOnChange));
    addProperty(autoGenerateProp_);
    addProperty(waitForThreadFinishedProp_);
    addProperty(progressProp_);
        addProgressBar(&progressProp_);
        //streamlines
    addProperty(maxNumStreamlinesProp_);
        maxNumStreamlinesProp_.onChange(MemberFunctionCallback<StreamlineCreator>(this,&StreamlineCreator::streamlineSettingsHaveBeenChanged));
        maxNumStreamlinesProp_.setGroupID("streamline");
    addProperty(streamlineLengthThresholdProp_);
        streamlineLengthThresholdProp_.onChange(MemberFunctionCallback<StreamlineCreator>(this,&StreamlineCreator::streamlineSettingsHaveBeenChanged));
        streamlineLengthThresholdProp_.setGroupID("streamline");
    addProperty(absoluteMagnitudeThresholdProp_);
        absoluteMagnitudeThresholdProp_.onChange(MemberFunctionCallback<StreamlineCreator>(this,&StreamlineCreator::streamlineSettingsHaveBeenChanged));
        absoluteMagnitudeThresholdProp_.onChange(MemberFunctionCallback<StreamlineCreator>(this,&StreamlineCreator::adjustRelativeThreshold));
        absoluteMagnitudeThresholdProp_.setGroupID("streamline");
    addProperty(fitAbsoluteMagnitudeProp_);
        fitAbsoluteMagnitudeProp_.setGroupID("streamline");
    addProperty(relativeMagnitudeThresholdProp_);
        relativeMagnitudeThresholdProp_.setReadOnlyFlag(true);
        relativeMagnitudeThresholdProp_.setGroupID("streamline");
    addProperty(filterModeProp_);
        filterModeProp_.addOption("linear","Linear",LINEAR);
        filterModeProp_.addOption("nearest","Nearest",NEAREST);
        filterModeProp_.setGroupID("streamline");
        filterModeProp_.onChange(MemberFunctionCallback<StreamlineCreator>(this,&StreamlineCreator::streamlineSettingsHaveBeenChanged));
    setPropertyGroupGuiName("streamline", "Streamline Settings");
        //streamline bundles
    addProperty(detectStreamlineBundlesProp_);
        detectStreamlineBundlesProp_.onChange(MemberFunctionCallback<StreamlineCreator>(this, &StreamlineCreator::generateBundlesHasBeenChanged));
        detectStreamlineBundlesProp_.setGroupID("streamlinebundles");
    addProperty(bundleDetectionProgressProp_);
        bundleDetectionProgressProp_.setGroupID("streamlinebundles");
    addProperty(maxAverageDistanceThresholdProp_);
        maxAverageDistanceThresholdProp_.onChange(MemberFunctionCallback<StreamlineCreator>(this, &StreamlineCreator::streamlineBundleSettingsHaveBeenChanged));
        maxAverageDistanceThresholdProp_.setGroupID("streamlinebundles");
    addProperty(minNumStreamlinesPerBundleProp_);
        minNumStreamlinesPerBundleProp_.onChange(MemberFunctionCallback<StreamlineCreator>(this, &StreamlineCreator::streamlineBundleSettingsHaveBeenChanged));
        minNumStreamlinesPerBundleProp_.setGroupID("streamlinebundles");
    addProperty(resampleSizeProp_);
        resampleSizeProp_.onChange(MemberFunctionCallback<StreamlineCreator>(this, &StreamlineCreator::streamlineBundleSettingsHaveBeenChanged));
        resampleSizeProp_.setGroupID("streamlinebundles");
    setPropertyGroupGuiName("streamlinebundles", "Streamline Bundle Settings");
        //seed
    addProperty(seedTimeProp_);
        //seedTimeProp_.setReadOnlyFlag(true);
        seedTimeProp_.onChange(MemberFunctionCallback<StreamlineCreator>(this,&StreamlineCreator::streamlineSettingsHaveBeenChanged));
        seedTimeProp_.setGroupID("seed");
    addProperty(resetSeedProp_);
        resetSeedProp_.onChange(MemberFunctionCallback<StreamlineCreator>(this,&StreamlineCreator::resetSeedsOnChange));
        resetSeedProp_.setGroupID("seed");
    setPropertyGroupGuiName("seed", "Seed Configuration");

    // adjust properties
    adjustRelativeThreshold();
    generateBundlesHasBeenChanged();

}

StreamlineCreator::~StreamlineCreator() {
    stopBackgroundThread();
    delete streamlineListThreadOutput_;
    if(volInport_.hasData())
        volInport_.getData()->Observable<VolumeObserver>::removeObserver(static_cast<VolumeObserver*>(this));
}

void StreamlineCreator::process() {
    //check, if background thread is finished
    if (backgroundThread_ && backgroundThread_->isFinished()) {
        // Retrieve error message.
        std::string message = backgroundThread_->getProgress().message_;

        // Delete thread.
        delete backgroundThread_;
        backgroundThread_ = 0;

        // reactivate recalculation
        calculateStreamlinesProp_.setReadOnlyFlag(false);

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
            if (backgroundThreadIsStreamlineCreator_ && detectStreamlineBundlesProp_.get())
                dirtyFlag_ = STREAMLINEBUNDLES;
            else
                dirtyFlag_ = NONE;

            // send the calculated data to the outport (owning the data)
            streamlineOutport_.setData(streamlineListThreadOutput_);
            streamlineListThreadOutput_ = 0;
        }

        // Reset background thread indicator
        backgroundThreadIsStreamlineCreator_ = false;
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
    streamlineListThreadOutput_ = new StreamlineList(volInport_.getData());

    switch (dirtyFlag_) {
    case STREAMLINES:
    {
        // update state
        calculateStreamlinesProp_.setReadOnlyFlag(true);
        backgroundThreadIsStreamlineCreator_ = true;

        // create a new thread
        backgroundThread_ = new StreamlineCreatorBackgroundThread(
                    this,
                    seedTimeProp_.get(),
                    volInport_.getThreadSafeData(),
                    streamlineListThreadOutput_,
                    maxNumStreamlinesProp_.get(),
                    streamlineLengthThresholdProp_.get(),
                    absoluteMagnitudeThresholdProp_.get(),
                    filterModeProp_.getValue()
                    );
        break;
    }
    case STREAMLINEBUNDLES:
    {
        // We don't need to calculate the streamlines again, they didn't change.
        // Thus, we add all old streamlines to our new streamlinelist.
        for (const Streamline& streamline : streamlineOutport_.getData()->getStreamlines())
            streamlineListThreadOutput_->addStreamline(streamline);

        // Calculate the noise threshold.
        size_t minNumStreamlines = static_cast<size_t>(minNumStreamlinesPerBundleProp_.get() * streamlineListThreadOutput_->getStreamlines().size() / 100.0f);

        // Create a new thread.
        backgroundThread_ = new StreamlineBundleDetectorBackgroundThread(
                    this,
                    streamlineListThreadOutput_,
                    maxAverageDistanceThresholdProp_.get(),
                    minNumStreamlines,
                    resampleSizeProp_.get()
            );
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

void StreamlineCreator::afterProcess() {
    Processor::setValid();
    if(waitForThreadFinishedProp_.get() && backgroundThread_)
         invalidate();
}

//---------------------------------------------------------------------------
//          Observers
//---------------------------------------------------------------------------

void StreamlineCreator::volumeDelete(const VolumeBase* source) {
    stopBackgroundThread();
    source->Observable<VolumeObserver>::removeObserver(static_cast<VolumeObserver*>(this));
}

void StreamlineCreator::volumeChange(const VolumeBase* source) {
    stopBackgroundThread();
    if(autoGenerateProp_.get()) dirtyFlag_ = STREAMLINES;

    VolumeMinMaxMagnitude* data = volInport_.getData()->getDerivedData<VolumeMinMaxMagnitude>();

    absoluteMagnitudeThresholdProp_.setMinValue(data->getMinMagnitude());
    absoluteMagnitudeThresholdProp_.setMaxValue(data->getMaxMagnitude());
    if(fitAbsoluteMagnitudeProp_.get()) {
        absoluteMagnitudeThresholdProp_.set(tgt::vec2(data->getMinMagnitude(), data->getMaxMagnitude()));
    }

    tgt::vec3 length = volInport_.getData()->getSpacing() * tgt::vec3(volInport_.getData()->getDimensions());
    maxAverageDistanceThresholdProp_.setMaxValue(tgt::length(length));

    // Invalidate processor
    invalidate();
}

void StreamlineCreator::afterConnectionAdded(const Port* source, const Port* connectedPort) {
}

void StreamlineCreator::beforeConnectionRemoved(const Port* source, const Port*) {
    stopBackgroundThread();
    if(volInport_.hasData())
        volInport_.getData()->Observable<VolumeObserver>::removeObserver(static_cast<VolumeObserver*>(this));
}

//---------------------------------------------------------------------------
//          Thread Handling
//---------------------------------------------------------------------------

void StreamlineCreator::stopBackgroundThread() {
    //stop and delete thread
    if (backgroundThread_) {
        backgroundThread_->interruptAndJoin();
        delete backgroundThread_;
        backgroundThread_ = 0;
    }

    if(backgroundThreadIsStreamlineCreator_) {
        // Reset data.
        streamlineOutport_.setData(0);
        // activate the calculate button
        calculateStreamlinesProp_.setReadOnlyFlag(false);
        // set progress to zero
        setProgress(0.f);
    }
    bundleDetectionProgressProp_.setProgress(0.f);
}

//---------------------------------------------------------------------------
//          Callbacks
//---------------------------------------------------------------------------
void StreamlineCreator::volumePortHasChanged() {
    // Add volume observer
    volInport_.getData()->Observable<VolumeObserver>::addObserver(static_cast<VolumeObserver*>(this));
    // handle new input like a change of old input
    volumeChange(0);
}

void StreamlineCreator::adjustRelativeThreshold() {
    //adjust read only property
    tgt::vec2 range(absoluteMagnitudeThresholdProp_.getMinValue(), absoluteMagnitudeThresholdProp_.getMaxValue());
    tgt::vec2 value = absoluteMagnitudeThresholdProp_.get();

    relativeMagnitudeThresholdProp_.set(  (value - tgt::vec2(range.x)) / tgt::vec2(range.y - range.x) * 100.f);
}

void StreamlineCreator::calculateStreamlinesOnChange() {
    dirtyFlag_ = STREAMLINES;
}

void StreamlineCreator::resetSeedsOnChange() {
    seedTimeProp_.set(time(0));
    if(autoGenerateProp_.get())
        dirtyFlag_ = STREAMLINES;
}

void StreamlineCreator::streamlineSettingsHaveBeenChanged() {
    if(autoGenerateProp_.get())
        dirtyFlag_ = STREAMLINES;
}

void StreamlineCreator::streamlineBundleSettingsHaveBeenChanged() {

    // Dont't set flag if initialization is not done.
    if (!isInitialized() || !streamlineOutport_.hasData())
        return;

    // Dont't overwrite the streamline dirty flag here
    if (detectStreamlineBundlesProp_.get() && dirtyFlag_ == NONE && !backgroundThreadIsStreamlineCreator_)
        dirtyFlag_ = STREAMLINEBUNDLES;
}

void StreamlineCreator::generateBundlesHasBeenChanged() {

    // Update UI.
    bundleDetectionProgressProp_.setVisibleFlag(detectStreamlineBundlesProp_.get());
    maxAverageDistanceThresholdProp_.setVisibleFlag(detectStreamlineBundlesProp_.get());
    minNumStreamlinesPerBundleProp_.setVisibleFlag(detectStreamlineBundlesProp_.get());
    resampleSizeProp_.setVisibleFlag(detectStreamlineBundlesProp_.get());

    // Trigger a calculation of bundles, if necessary.
    streamlineBundleSettingsHaveBeenChanged();
}

}   // namespace
