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

#include "cosmologytimestepanimationprocessor.h"
#include "voreen/core/voreenapplication.h"
#include "voreen/core/datastructures/callback/memberfunctioncallback.h"
#include "tgt/event/timeevent.h"
#include <limits.h>

namespace voreen {

CosmologyTimeStepAnimationProcessor::CosmologyTimeStepAnimationProcessor()
    : Processor()
    , counter_(0)
    , timerIsActive_(false)
    , timer_(0)
    , eventHandler_()
    , enableTimerProp_("enableTimerProp", "Enable Timer", false)
    , intervalProp_("intervalProp", "Interval (ms)", 1000, 1, 10000)
    , counterStyleProp_("counterStyle", "Counter Style")
    , tickCounterProp_("tickCounterProp", "Counter",0,std::numeric_limits<int32_t>::min(),std::numeric_limits<int32_t>::max())
    , stepProp_("stepProp_", "Stepping", 0.0f, 0.0f, 1.0f)
    , animationInterval_("animationInterval_", "Animation Interval", 0.0f, 0.0f, 624.0f)
    , resetProp_("resetProp", "Reset Counters")
{
    //set timer
    eventHandler_.addListenerToBack(this);
    timer_ = VoreenApplication::app()->createTimer(&eventHandler_);

    //enable timer
    addProperty(enableTimerProp_);
        enableTimerProp_.onChange(MemberFunctionCallback<CosmologyTimeStepAnimationProcessor>(this, &CosmologyTimeStepAnimationProcessor::toggleTimer));
    //settings
    addProperty(intervalProp_);
        intervalProp_.setGroupID("settings");
        intervalProp_.onChange(MemberFunctionCallback<CosmologyTimeStepAnimationProcessor>(this, &CosmologyTimeStepAnimationProcessor::onIntervalChange));
    addProperty(counterStyleProp_);
        counterStyleProp_.addOption("linear", "Linear", CosmologyTimeStepAnimationProcessor::COUNTER_LINEAR);
        counterStyleProp_.addOption("linear_reverse", "Linear (Reverse)", CosmologyTimeStepAnimationProcessor::COUNTER_LINEAR_REVERSE);
        counterStyleProp_.addOption("cyclic", "Cyclic", CosmologyTimeStepAnimationProcessor::COUNTER_CYCLIC);
        counterStyleProp_.addOption("cyclic_reverse", "Cyclic (Reverse)", CosmologyTimeStepAnimationProcessor::COUNTER_CYCLIC_REVERSE);
        counterStyleProp_.set("cyclic");
        counterStyleProp_.onChange(MemberFunctionCallback<CosmologyTimeStepAnimationProcessor>(this, &CosmologyTimeStepAnimationProcessor::onCounterStyleChange));
        counterStyleProp_.setGroupID("settings");
    addProperty(animationInterval_);
        ON_CHANGE(animationInterval_, CosmologyTimeStepAnimationProcessor, resetCounter)
        animationInterval_.setGroupID("settings");
    addProperty(stepProp_);
        stepProp_.setGroupID("settings");
    addProperty(tickCounterProp_);
        tickCounterProp_.setGroupID("settings");

    addProperty(resetProp_);
        resetProp_.onChange(MemberFunctionCallback<CosmologyTimeStepAnimationProcessor>(this, &CosmologyTimeStepAnimationProcessor::resetCounter));

    //modify properties
    tickCounterProp_.setReadOnlyFlag(true);
    setPropertyGroupGuiName("settings","Settings");
    setPropertyGroupGuiName("res","High Resolution Counter");

    onCounterStyleChange();
    resetCounter();
}

CosmologyTimeStepAnimationProcessor::~CosmologyTimeStepAnimationProcessor() {
    delete timer_;
}

void CosmologyTimeStepAnimationProcessor::process() {}

void CosmologyTimeStepAnimationProcessor::initialize() {
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

void CosmologyTimeStepAnimationProcessor::timerEvent(tgt::TimeEvent* te) {
    const float min = animationInterval_.get().x;
    const float max = animationInterval_.get().y;
    const float step = stepProp_.get();

    if (counter_ < min)   // in case of overflow of int, reset counter to zero to prevent strange results
        counter_ = min;
    if (counter_ > max)
        counter_ = max;

    if (te)
        te->accept();

    switch (counterStyleProp_.getValue()) {
    case COUNTER_LINEAR:
            counter_ += step;
            if(counter_ >= max) {
                enableTimerProp_.set(false);
                return;
            }
            tickCounterProp_.set(counter_);
            break;
        case COUNTER_LINEAR_REVERSE:
            counter_ -= step;
            if(counter_ <= min) {
                enableTimerProp_.set(false);
                return;
            }
            tickCounterProp_.set(counter_);
            break;
        case COUNTER_CYCLIC:
            counter_ += step;
            if(counter_ >= max)
                counter_ = min;
            tickCounterProp_.set(counter_);
            break;
        case COUNTER_CYCLIC_REVERSE:
            counter_ -= step;
            if(counter_ < min)
                counter_ = max;
            tickCounterProp_.set(counter_);
            break;
    }

}

bool CosmologyTimeStepAnimationProcessor::isTimerActive() const {
    return timerIsActive_;
}

void CosmologyTimeStepAnimationProcessor::startTimer() {

    if (!timer_) {
        LWARNING("No timer.");
        return;
    }

    if (!timerIsActive_) {
        timerIsActive_ = true;
        timer_->start(intervalProp_.get());
    }
}

void CosmologyTimeStepAnimationProcessor::stopTimer() {
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

void CosmologyTimeStepAnimationProcessor::onCounterStyleChange() {
    resetCounter();
}

void CosmologyTimeStepAnimationProcessor::onIntervalChange() {
    // If the timer has been active and the interval changes, the timer has to be
    // restarted.
    //
    if (timerIsActive_) {
        stopTimer();
        startTimer();
    }
}

void CosmologyTimeStepAnimationProcessor::resetCounter() {
    if(timerIsActive_)
        enableTimerProp_.set(false);

    if(counterStyleProp_.getValue() == COUNTER_LINEAR ||
       counterStyleProp_.getValue() == COUNTER_CYCLIC )
       counter_ = animationInterval_.get().x;
    else
        counter_ = animationInterval_.get().y;
    tickCounterProp_.set(counter_);
}


void CosmologyTimeStepAnimationProcessor::toggleTimer() {
    if (enableTimerProp_.get())
        startTimer();
    else
        stopTimer();
}


}   // namespace
