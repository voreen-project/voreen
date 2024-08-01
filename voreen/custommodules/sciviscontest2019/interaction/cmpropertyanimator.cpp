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

#include "cmpropertyanimator.h"
namespace tgt {
template<class T>
Quaternion<T> gqft(const Vector3<T> src, const Vector3<T> dest) {
    // This is taken from Game Programming Gems 1 and Real Time Rendering

    Quaternion<T> q;

    Vector3<T> v1 = normalize(src);
    Vector3<T> v2 = normalize(dest);

    Vector3<T> crs = cross(v1, v2);

    T v1v2dot = dot(v1, v2);

    if (std::abs(v1v2dot) >= 1) // the vectors are identical
        return Quaternion<T>(0, 0, 0, 1); // ... so we return a rotation that does nothing

    T root = std::sqrt((1 + v1v2dot) * 2);

    if (root < T(1e-6)) { // do this for numerical stability
        Vector3<T> axis = cross(Vector3<T>(1, 0, 0), src);
        if (length(axis) == 0) {
            axis = cross(Vector3<T>(0, 1, 0), src);
        }
        axis = normalize(axis);
        q = Quaternion<T>::createQuat(PIf, axis);
    }
    else {
        q.x = crs.x;
        q.y = crs.y;
        q.z = crs.z;
        q.w = T(1) + v1v2dot;
    }

    return normalize(q);
}
}

namespace voreen {

float CMPropertyInterpolator::transform(float a, InterpolationMode mode) {
    switch(mode) {
        case LINEAR:
            return a;
            break;
        case SMOOTH:
            if(a==0) {
                return 0;
            } else if(a==1) {
                return 1;
            } else {
                static const float lambda = 10;
                static const float d = 1.0f/(1+exp(lambda*(0.5)));
                return (1.0f/(1+exp(lambda*(0.5f-a)))-d)/(1-2*d);
            }
    }
}



CMCameraPropertyInterpolator::CMCameraPropertyInterpolator(CameraProperty* property, tgt::vec3 targetPos, float targetDistance)
    : property_(property)
    , initial_(property->get())
    , initialDistance_(tgt::distance(initial_.getPosition(), initial_.getFocus()))
    , targetPos_(targetPos)
    , targetDistance_(targetDistance)
    //, offsetRotation_(tgt::gqft(initial_.getPosition() - initial_.getFocus(), property->get().getFocus() - targetPos_))
    , offsetRotation_(tgt::gqft(initial_.getPosition() - initial_.getFocus(), property->get().getPosition() - targetPos_))
{
}
CMCameraPropertyInterpolator::~CMCameraPropertyInterpolator() {
}
void CMCameraPropertyInterpolator::interpolate(float a) {
    float aCam = transform(a, SMOOTH);
    float aOr = transform(a, SMOOTH);
    float aDist = transform(a, SMOOTH);
    tgt::vec3 pos = initial_.getFocus()*(1-aCam)+aCam*targetPos_;
    tgt::quat orientationRotation = tgt::slerpQuat(tgt::quat(0,0,0,1), offsetRotation_, aOr);
    float dist = initialDistance_*(1-aDist)+aDist*targetDistance_;
    if(aOr == 1.0f) {
        orientationRotation = offsetRotation_;
    }
    tgt::Camera newCam = initial_;
    tgt::vec3 camOffset = tgt::normalize(tgt::quat::rotate(initial_.getPosition()-initial_.getFocus(), orientationRotation))*dist;
    newCam.setFocus(pos);
    newCam.setUpVector(tgt::quat::rotate(initial_.getUpVector(), orientationRotation));
    newCam.setPosition(pos+camOffset);
    property_->set(newCam);
}


CMPropertyAnimator::CMPropertyAnimator(unsigned int maxUPS)
    : timer_(nullptr)
    , maxUPS_(maxUPS)
    , eventHandler_()
{
    eventHandler_.addListenerToBack(this);
}
CMPropertyAnimator::~CMPropertyAnimator() {
    if (timer_)
        //Apparently deleted somewhere else?
        //delete timer_;
        timer_ = nullptr;
    for(auto i : interpolators_) {
        delete i;
    }
}

void CMPropertyAnimator::animate(unsigned int lengthInMS, std::vector<CMPropertyInterpolator*> interpolators) {
#if 0
    for (auto& interpolator : interpolators) {
        interpolator->interpolate(1.0f);
    }
    return;
#endif
    timer_ = VoreenApplication::app()->createTimer(&eventHandler_);
    //Get rid of old interpolators
    for(auto i : interpolators_) {
        delete i;
    }
    if(timer_) {
        float updateInterval = 1000.0f/maxUPS_;
        interpolators_ = interpolators;
        //if limit is 0, we still want at least one update:
        float limit = std::max(1.0f, lengthInMS/updateInterval);
        timer_->start(updateInterval, limit);
    }
}
void CMPropertyAnimator::onEvent(tgt::Event* e) {
    tgt::TimeEvent* te = dynamic_cast<tgt::TimeEvent*>(e);
    if (te){
        float progress = (float)(timer_->getCount())/timer_->getLimit();
        for(auto& interpolator : interpolators_) {
            interpolator->interpolate(progress);
        }
    }
}
}
