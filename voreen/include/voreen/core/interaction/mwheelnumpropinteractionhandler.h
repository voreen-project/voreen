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

#ifndef VRN_MWHEELNUMPROPINTERACTIONHANDLER_H
#define VRN_MWHEELNUMPROPINTERACTIONHANDLER_H

#include "voreen/core/interaction/interactionhandler.h"
#include "voreen/core/properties/numericproperty.h"
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/utils/voreenqualitymode.h"

#include "tgt/event/mouseevent.h"

namespace voreen {

template<typename T>
class MWheelNumPropInteractionHandler: public InteractionHandler {

public:
    /// Default constructor needed for serialization. Do not call it directly!
    MWheelNumPropInteractionHandler();

    /**
     * Constructor.
     *
     * @param id Identifier that must be unique across all interaction handlers
     *  of a processor. Must not be empty.
     * @param name the string that is to be displayed in the GUI
     * @param numProp numeric property the handler will operate on. Must not be null.
     * @param The modifier that has to be active for allowing the handler to events
     * @param sharing @see InteractionHandler
     * @param enabled @see InteractionHandler
     */
    MWheelNumPropInteractionHandler(const std::string& id, const std::string& name, NumericProperty<T>* numProp,
        tgt::Event::Modifier modifier = tgt::Event::MODIFIER_NONE, bool sharing = false, bool enabled = true);

    virtual ~MWheelNumPropInteractionHandler() {}

    virtual std::string getClassName() const   { return "MWheelNumPropInteractionHandler";     }
    virtual InteractionHandler* create() const { return new MWheelNumPropInteractionHandler(); }

    virtual void setProperty(NumericProperty<T>* prop) {
        numProp_ = prop;
    }

    virtual NumericProperty<T>* setProperty() {
        return numProp_;
    }

private:
    virtual void onEvent(tgt::Event* e);

    NumericProperty<T>* numProp_;
    tgt::ivec2 lastMousePos_;
};

//---------------------------------------------------------------------------


template<typename T>
voreen::MWheelNumPropInteractionHandler<T>::MWheelNumPropInteractionHandler() :
    InteractionHandler("dummy", "dummy"),
    lastMousePos_(-1)
{ }

template<typename T>
MWheelNumPropInteractionHandler<T>::MWheelNumPropInteractionHandler(const std::string& id, const std::string& name,
        NumericProperty<T>* numProp,  tgt::Event::Modifier modifier, bool sharing, bool enabled)
    : InteractionHandler(id, name,
        tgt::MouseEvent::MouseButtons(tgt::MouseEvent::MOUSE_WHEEL | tgt::MouseEvent::MOUSE_BUTTON_MIDDLE),
        static_cast<tgt::MouseEvent::MouseAction>(tgt::MouseEvent::WHEEL | tgt::MouseEvent::MOTION | tgt::MouseEvent::PRESSED | tgt::MouseEvent::RELEASED), modifier,
        sharing, enabled)
{
    tgtAssert(numProp, "No numeric property");
    numProp_ = numProp;
}

template<typename T>
void MWheelNumPropInteractionHandler<T>::onEvent(tgt::Event* eve) {
    tgtAssert(numProp_, "No numeric property");

    const float motionScale = -0.33f;

    if (tgt::MouseEvent* mouseEve = dynamic_cast<tgt::MouseEvent*>(eve)) {

        bool accept = false;

        // wheel up/down => discrete increment/decrement
        if (mouseEve->action() == tgt::MouseEvent::WHEEL) {
            if (mouseEve->button() == tgt::MouseEvent::MOUSE_WHEEL_UP) {
                accept = true;
                if (numProp_->get() != numProp_->getMaxValue()) {
                    if ((numProp_->get() + numProp_->getStepping()) > numProp_->getMaxValue())
                        numProp_->set(numProp_->getMaxValue());
                    else
                        numProp_->increase();
                }
            }
            else if (mouseEve->button() == tgt::MouseEvent::MOUSE_WHEEL_DOWN) {
                if (numProp_->get() != numProp_->getMinValue()) {
                    if ((numProp_->get() - numProp_->getStepping()) < numProp_->getMinValue())
                        numProp_->set(numProp_->getMinValue());
                    else
                        numProp_->decrease();
                }
                accept = true;
            }
        }

        // mouse action with wheel button => continous increment/decrement
        if (mouseEve->button() == tgt::MouseEvent::MOUSE_BUTTON_MIDDLE) {
            if (mouseEve->action() == tgt::MouseEvent::PRESSED) {
                // store initial mouse position, start interaction mode
                lastMousePos_ = mouseEve->coord();
                QualityMode.requestQualityMode(VoreenQualityMode::RQ_INTERACTIVE, this);
                accept = true;
            }
            else if (mouseEve->action() == tgt::MouseEvent::RELEASED) {
                // clear mouse position, leave interaction mode
                lastMousePos_ = tgt::ivec2(-1);
                QualityMode.requestQualityMode(VoreenQualityMode::RQ_DEFAULT, this);
                accept = true;
            }
            else if (mouseEve->action() == tgt::MouseEvent::MOTION) {
                if (lastMousePos_.y != -1) {
                    // take vertical offset for determining increment
                    float offset = (mouseEve->coord().y - lastMousePos_.y) * motionScale;
                    T newVal = static_cast<T>(numProp_->get() + offset*numProp_->getStepping());
                    if (newVal < numProp_->getMinValue()) {
                        newVal = numProp_->getMinValue();
                        lastMousePos_ = mouseEve->coord();
                    }
                    if (newVal > numProp_->getMaxValue()) {
                        newVal = numProp_->getMaxValue();
                        lastMousePos_ = mouseEve->coord();
                    }
                    if (newVal != numProp_->get()) {
                        numProp_->set(newVal);
                        lastMousePos_ = mouseEve->coord();
                    }
                }
               accept = true;
            }

        }

        if (accept)
            eve->accept();

    }
}
} // namespace

#endif // VRN_MWHEELNUMPROPINTERACTIONHANDLER_H
