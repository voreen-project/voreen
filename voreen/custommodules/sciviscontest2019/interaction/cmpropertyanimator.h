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

#ifndef VRN_CMPROPERTYANIMATOR_H
#define VRN_CMPROPERTYANIMATOR_H

#include "voreen/core/properties/numericproperty.h"
#include "voreen/core/voreenapplication.h"
#include "voreen/core/properties/cameraproperty.h"
#include "tgt/event/timeevent.h"
#include "tgt/event/eventlistener.h"
#include "tgt/timer.h"

namespace voreen {

/**
 * Abstract base class for interpolators. Implementing classes should update properties
 * during the interpolate method to make use of CMPropertyAnimator's animation features.
 */
class CMPropertyInterpolator {
public:
    /**
     * Indicates method relation between time spent during animation and the a value, used
     * to interpolate between property values.
     */
    enum InterpolationMode {
        LINEAR, /// progress = a
        SMOOTH  /// Cut sigmoid function, so that f(0) = and f(1) = 1
    };

    /**
     * Implement and set property values relative to a.
     * @param a Progress indicator for animation. a has to be in [0,1].
     *          a = 0 should be the start value, a = 1 the final value of the animation.
     *          Subsequent calls to interpolate from the same animation sequence are
     *          guaranteed give values of a greater or equal to previous values.
     */
    virtual void interpolate(float a) = 0;

protected:
    /**
     * Concrete implementation of interpolation modes.
     * Should be called first when implementing interpolate to get a value for a according
     * to the interpolation mode.
     * @param a has to be in [0,1]
     * @return value in [0,1]
     */
    static float transform(float a, InterpolationMode mode);
};

/**
 * Interpolator for camera properties
 */
class CMCameraPropertyInterpolator : public CMPropertyInterpolator {
public:
    /**
     * Constructor.
     * @param property Camera to be transitioned to new orientation and position
     * @param targetPos Final focus point of camera
     * @param targetDistance Final distance between camera center and targetPos
     */
    CMCameraPropertyInterpolator(CameraProperty* property, tgt::vec3 targetPos, float targetDistance);
    /**
     * Destructor.
     */
    ~CMCameraPropertyInterpolator();

    // See CMPropertyInterpolator.
    virtual void interpolate(float a);

protected:
    /// Camera to be transitioned to new orientation and position
    CameraProperty* property_;
    /// Initial camera (a=0)
    const tgt::Camera initial_;
    /// initial distance between camera and focus
    float initialDistance_;
    /// Final focus point of camera
    tgt::vec3 targetPos_;
    /// Final distance between camera center and targetPos
    float targetDistance_;
    /// Difference between first and final orientation of the camera
    tgt::quat offsetRotation_;
};

/**
 * Generic animator for NumericProperties.
 * Animation yields in a direct transition between the current value and a target value.
 */
template<typename T>
class CMNumericPropertyInterpolator : public CMPropertyInterpolator {
public:
    /**
     * Constructor.
     * @param InterpolationMode Desired mode used for interpolation
     * @param property NumericProperty to be animated
     * @param target Value for property at a=1
     */
    CMNumericPropertyInterpolator(InterpolationMode interpolationMode, NumericProperty<T>* property, T target)
        : interpolationMode_(interpolationMode)
        , property_(property)
        , initial_(property->get())
        , target_(target)
    {
    }

    /**
     * Destructor.
     */
    ~CMNumericPropertyInterpolator() {}


protected:
    // See CMPropertyInterpolator.
    virtual void interpolate(float a) {
        a = transform(a, interpolationMode_);
        property_->set(initial_*(1-a)+a*target_);
    }

    /// Interpolation mode used for animation
    InterpolationMode interpolationMode_;
    /// NumericProperty to be animated
    NumericProperty<T>* property_;
    /// Value for property at a=0
    T initial_;
    /// Value for property at a=1
    T target_;
};

/**
 * Generic animator for NumericProperties.
 * Animation yields in a step wise transition between the current value,  a target value, and an intermediate value.
 */
template<typename T>
class CMTwoStepNumericPropertyInterpolator : public CMPropertyInterpolator {
public:
    /**
     * Constructor.
     * @param InterpolationMode Desired mode used for interpolation
     * @param property NumericProperty to be animated
     * @param interstation Value for property at a=0.5
     * @param target Value for property at a=1
     */
    CMTwoStepNumericPropertyInterpolator(InterpolationMode interpolationMode, NumericProperty<T>* property, T interstation, T target)
        : interpolationMode_(interpolationMode)
        , property_(property)
        , initial_(property_->get())
        , interstation_(interstation)
        , target_(target)
    {
    }

    /**
     * Destructor.
     */
    ~CMTwoStepNumericPropertyInterpolator() {}


protected:
    // See CMPropertyInterpolator.
    virtual void interpolate(float a) {
        a = transform(a, interpolationMode_);
        if(a<=0.5f) {
            a=a*2;
            property_->set(initial_*(1-a)+a*interstation_);
        } else {
            a=a*2-1;
            property_->set(interstation_*(1-a)+a*target_);
        }
    }

    /// Interpolation mode used for animation
    InterpolationMode interpolationMode_;
    /// NumericProperty to be animated
    NumericProperty<T>* property_;
    /// Value for property at a=0
    T initial_;
    /// Value for property at a=0.5
    T interstation_;
    /// Value for property at a=1
    T target_;
};



/**
 * Class used to animate arbitrary properties.
 */
class CMPropertyAnimator : public tgt::EventListener {
public:
    /**
     * Constructor.
     * @param maxUPS Maximum number of updates per second during animation.
     */
    CMPropertyAnimator(unsigned int maxUPS);
    /**
     * Destructor.
     */
    ~CMPropertyAnimator();

    /**
     * Animate the properties using the interpolators for the given time
     * @param lengthInMS Animation duration
     * @param interpolators Interpolators of the properties subject to animation.
     */
    void animate(unsigned int lengthInMS, std::vector<CMPropertyInterpolator*> interpolators);

protected:
    virtual void onEvent(tgt::Event* e);

private:
    /// Used to generate events
    tgt::Timer* timer_;
    /// Handles timer's events
    tgt::EventHandler eventHandler_;
    /// Maximum number of updates per second during animation.
    unsigned int maxUPS_;
    /// Interpolators of the properties subject to animation.
    std::vector<CMPropertyInterpolator*> interpolators_;
};
} //namespace
#endif // VRN_CMPROPERTYANIMATOR_H
