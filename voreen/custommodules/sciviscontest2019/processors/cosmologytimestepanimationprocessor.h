/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2014 University of Muenster, Germany.                        *
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

#ifndef VRN_COSMOLOGYTIMESTEPANIMATIONPROCESSOR_H
#define VRN_COSMOLOGYTIMESTEPANIMATIONPROCESSOR_H

#include "voreen/core/processors/processor.h"

#include "tgt/timer.h"
#include "tgt/event/eventhandler.h"
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/buttonproperty.h"
#include "voreen/core/properties/floatproperty.h"
#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/optionproperty.h"
#include "voreen/core/properties/numeric/intervalproperty.h"

namespace voreen {
class VRN_CORE_API CosmologyTimeStepAnimationProcessor : public Processor {
public:
    CosmologyTimeStepAnimationProcessor();
    virtual ~CosmologyTimeStepAnimationProcessor();
    virtual Processor* create() const { return new CosmologyTimeStepAnimationProcessor(); }

    virtual std::string getClassName() const    { return "CosmologyTimeStepAnimationProcessor"; }
    virtual std::string getCategory() const     { return "Viscontest201x";                      }
    virtual CodeState getCodeState() const      { return CODE_STATE_EXPERIMENTAL;               }
    virtual bool isUtility() const              { return true;                                  }

    virtual void timerEvent(tgt::TimeEvent* te);

    /**
     * Indicates whether the timer is currently running.
     */
    bool isTimerActive() const;

    /**
     * Starts the timer with the intervals contained in intervalProp_ property.
     */
    void startTimer();

    /**
     * Stops the timer.
     */
    void stopTimer();

    void resetCounter();

protected:
    virtual void setDescriptions() {
        setDescription("Performs clocking for time-based manipulation of properties in the network."\
                       "A counter represented by a property is increased on every time event."\
                       "This counter is intended to be linked with other properties in order to propagate the clock tick."\
                       "This processor has neither any ports nor does it require any connections, instead it is used as a"\
                       "'floating' processor in the networks, solely connected by property links.");
    }

    virtual void process();
    virtual void initialize();

    void toggleTimer();
    void onCounterStyleChange();
    void onIntervalChange();

protected:
    enum CounterStyle {
        COUNTER_LINEAR,
        COUNTER_LINEAR_REVERSE,
        COUNTER_CYCLIC,
        COUNTER_CYCLIC_REVERSE
    };
    friend class OptionProperty<CounterStyle>;

    float counter_;
    //int period_;
    bool timerIsActive_;                /** Indicates whether the timer is active. */
    tgt::Timer* timer_;                 /** Timer object. */
    tgt::EventHandler eventHandler_;    /** A local eventhanlde which is added to the timer. */
    //enable
    BoolProperty enableTimerProp_;          /** Determines whether the clocking is active. */
    //general
    IntProperty intervalProp_;              /** Determines the inveral for the clock in ms. */
    OptionProperty<CounterStyle> counterStyleProp_;
    FloatProperty stepProp_;
    ButtonProperty resetProp_;

    FloatProperty tickCounterProp_;
    FloatIntervalProperty animationInterval_;

};

}   // namespace

#endif
