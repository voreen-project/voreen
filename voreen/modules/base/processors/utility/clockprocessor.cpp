/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2024 University of Muenster, Germany,                        *
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

#include "clockprocessor.h"
#include "voreen/core/voreenapplication.h"
#include "voreen/core/datastructures/callback/memberfunctioncallback.h"
#include "tgt/event/timeevent.h"
#include <limits.h>

namespace voreen {

ClockProcessor::ClockProcessor()
    : Processor(),
    counter_(0),
    timerIsActive_(false),
    timer_(0),
    eventHandler_(),
    enableTimerProp_("enableTimerProp", "Enable Timer", false),
    intervalProp_("intervalProp", "Interval (ms)", 1000, 1, 10000),
    counterStyleProp_("counterStyle", "Counter Style"),
    beginProp_("beginProp","Start Count (including)",0,std::numeric_limits<int32_t>::min(),std::numeric_limits<int32_t>::max()),
    endProp_("endProp","End Count (including)",10,std::numeric_limits<int32_t>::min(),std::numeric_limits<int32_t>::max()),
    tickCounterProp_("tickCounterProp", "Counter",0,std::numeric_limits<int32_t>::min(),std::numeric_limits<int32_t>::max()),
    useHighResCounterProp_("useHighResCounterProp", "Use High Resolution Counter", false),
    resolutionProp_("resolutionProp", "Resolution (*counts)", 0.5f, 0.0f, 1.0f),
    highResCounterProp_("highResCounterProp", "High Resolution Counter", 0.0f,-std::numeric_limits<float>::max(),std::numeric_limits<float>::max()),
    resetProp_("resetProp", "Reset Counters")
{
    //set timer
    eventHandler_.addListenerToBack(this);
    timer_ = VoreenApplication::app()->createTimer(&eventHandler_);

    //enable timer
    addProperty(enableTimerProp_);
        enableTimerProp_.onChange(MemberFunctionCallback<ClockProcessor>(this, &ClockProcessor::toggleTimer));
    //settings
    addProperty(intervalProp_);
        intervalProp_.setGroupID("settings");
        intervalProp_.onChange(MemberFunctionCallback<ClockProcessor>(this, &ClockProcessor::onIntervalChange));
    addProperty(counterStyleProp_);
        counterStyleProp_.addOption("linear", "Linear", ClockProcessor::COUNTER_LINEAR);
        counterStyleProp_.addOption("linear_reverse", "Linear (Reverse)", ClockProcessor::COUNTER_LINEAR_REVERSE);
        counterStyleProp_.addOption("cyclic", "Cyclic", ClockProcessor::COUNTER_CYCLIC);
        counterStyleProp_.addOption("cyclic_reverse", "Cyclic (Reverse)", ClockProcessor::COUNTER_CYCLIC_REVERSE);
        counterStyleProp_.set("cyclic");
        counterStyleProp_.onChange(MemberFunctionCallback<ClockProcessor>(this, &ClockProcessor::onCounterStyleChange));
        counterStyleProp_.setGroupID("settings");
    addProperty(beginProp_);
        beginProp_.onChange(MemberFunctionCallback<ClockProcessor>(this, &ClockProcessor::onBeginChange));
        beginProp_.setGroupID("settings");
    addProperty(endProp_);
        endProp_.onChange(MemberFunctionCallback<ClockProcessor>(this, &ClockProcessor::onEndChange));
        endProp_.setGroupID("settings");
    addProperty(tickCounterProp_);
        tickCounterProp_.setGroupID("settings");
    //res counter
    addProperty(useHighResCounterProp_);
        useHighResCounterProp_.setGroupID("res");
        useHighResCounterProp_.onChange(MemberFunctionCallback<ClockProcessor>(this, &ClockProcessor::onUseHighResCounterChange));
    addProperty(resolutionProp_);
        resolutionProp_.setGroupID("res");
        resolutionProp_.setNumDecimals(3);
        resolutionProp_.onChange(MemberFunctionCallback<ClockProcessor>(this, &ClockProcessor::onResolutionChange));
    addProperty(highResCounterProp_);
        highResCounterProp_.setGroupID("res");
    //reset button
    addProperty(resetProp_);
        resetProp_.onChange(MemberFunctionCallback<ClockProcessor>(this, &ClockProcessor::resetCounter));

    //modify properties
    tickCounterProp_.setReadOnlyFlag(true);
    highResCounterProp_.setReadOnlyFlag(true);
    setPropertyGroupGuiName("settings","Settings");
    setPropertyGroupGuiName("res","High Resolution Counter");

    onCounterStyleChange();
    //onPeriodChange();
    onBeginChange();
    onEndChange();
    onUseHighResCounterChange();
}

ClockProcessor::~ClockProcessor() {
    delete timer_;
}

void ClockProcessor::process() {}

void ClockProcessor::initialize() {
    Processor::initialize();

    // Startup the timer if desired
    if (enableTimerProp_.get()) {
        if (timer_) {
            startTimer();
        }
        else {
            LWARNING("Timer object could not be created. Disabling processor.");
            enableTimerProp_.set(false);
            enableTimerProp_.setReadOnlyFlag(true);
            timerIsActive_ = false;
        }
    }
}

void ClockProcessor::timerEvent(tgt::TimeEvent* te) {
    if (counter_ < 0)   // in case of overflow of int, reset counter to zero to prevent strange results
        counter_ = 0;

    if (te)
        te->accept();

    switch (counterStyleProp_.getValue()) {
        case COUNTER_LINEAR:
            if(counter_ == endProp_.get()) {
                enableTimerProp_.set(false);
                return;
            }
            tickCounterProp_.set(++counter_);
            break;
        case COUNTER_LINEAR_REVERSE:
            if(counter_ == beginProp_.get()) {
                enableTimerProp_.set(false);
                return;
            }
            tickCounterProp_.set(--counter_);
            break;
        case COUNTER_CYCLIC:
            counter_++;
            if(counter_ > endProp_.get())
                counter_ = beginProp_.get();
            tickCounterProp_.set(counter_);
            break;
        case COUNTER_CYCLIC_REVERSE:
            counter_--;
            if(counter_ < beginProp_.get())
                counter_ = endProp_.get();
            tickCounterProp_.set(counter_);
            break;
    }

    if (useHighResCounterProp_.get())
        triggerHighResCounter();


}

bool ClockProcessor::isTimerActive() const {
    return timerIsActive_;
}

void ClockProcessor::startTimer() {

    if (!timer_) {
        LWARNING("No timer.");
        return;
    }

    if (!timerIsActive_) {
        timerIsActive_ = true;
        timer_->start(intervalProp_.get());
    }
}

void ClockProcessor::stopTimer() {
    if (timer_) {
        timer_->stop();
    }
    else {
        LWARNING("No timer.");
    }
    timerIsActive_ = false;
}

// protected methods
//

void ClockProcessor::onCounterStyleChange() {
    resetCounter();
}

void ClockProcessor::onIntervalChange() {
    // If the timer has been active and the interval changes, the timer has to be
    // restarted.
    //
    if (timerIsActive_) {
        stopTimer();
        startTimer();
    }
}

void ClockProcessor::onBeginChange() {
    if(timerIsActive_)
        enableTimerProp_.set(false);

    if(beginProp_.get() > endProp_.get())
        endProp_.set(beginProp_.get());

    resetCounter();
}

void ClockProcessor::onEndChange() {
    if(timerIsActive_)
        enableTimerProp_.set(false);

    if(beginProp_.get() > endProp_.get())
        beginProp_.set(endProp_.get());

    resetCounter();
}

void ClockProcessor::resetCounter() {
    if(timerIsActive_)
        enableTimerProp_.set(false);

    if(counterStyleProp_.getValue() == COUNTER_LINEAR ||
       counterStyleProp_.getValue() == COUNTER_CYCLIC )
        counter_ = beginProp_.get();
    else
        counter_ = endProp_.get();
    tickCounterProp_.set(counter_);
    onResolutionChange();
}

void ClockProcessor::onResolutionChange() {
    highResCounterProp_.set(static_cast<float>(counter_)*resolutionProp_.get());
}

void ClockProcessor::onUseHighResCounterChange() {
    if (useHighResCounterProp_.get() == true) {
        resolutionProp_.setVisibleFlag(true);
        highResCounterProp_.setVisibleFlag(true);
        onResolutionChange();
    }
    else {
        resolutionProp_.setVisibleFlag(false);
        highResCounterProp_.setVisibleFlag(false);
        highResCounterProp_.set(0.0f);
    }

}

void ClockProcessor::toggleTimer() {
    if (enableTimerProp_.get())
        startTimer();
    else
        stopTimer();
}

void ClockProcessor::triggerHighResCounter() {
    highResCounterProp_.set(static_cast<float>(counter_)*resolutionProp_.get());
}

}   // namespace
