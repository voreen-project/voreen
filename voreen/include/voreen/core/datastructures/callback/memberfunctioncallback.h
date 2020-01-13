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

#ifndef VRN_MEMBERFUNCTIONCALLBACK_H
#define VRN_MEMBERFUNCTIONCALLBACK_H

#include "voreen/core/datastructures/callback/callback.h"

namespace voreen {

#define ON_PROPERTY_CHANGE(PROPERTY,PROPERTYOWNER,FUNCTION) PROPERTY.onChange(MemberFunctionCallback<PROPERTYOWNER>(this, &PROPERTYOWNER::FUNCTION));
#define ON_CHANGE(CHANGEABLE, RECEIVER, FUNCTION) CHANGEABLE.onChange(MemberFunctionCallback<RECEIVER>(this, &RECEIVER::FUNCTION));
/**
 * This Callback can call a member function of an Object with signature void mem()
 *
 * Useage for Property::onChange e.g.:
 * myProperty.onChange(MemberFunctionCallback<myProcessor>(this, &myProcessor::myFunction));
 */
template<class T>
class MemberFunctionCallback : public Callback {
public:

    MemberFunctionCallback(T* target, void (T::*fpt)())
        : target_(target), fpt_(fpt)
    {
    }

    virtual ~MemberFunctionCallback() {}
    virtual MemberFunctionCallback* clone() const { return new MemberFunctionCallback(*this); }

    virtual void exec() {
        if ((target_ != 0) && (fpt_ != 0))
            (target_->*fpt_)();
    }

private:
    T*    target_;
    void (T::*fpt_)();
};

// ============================================================================

/**
 * This Callback can call a member function of an Object with signature void mem(P param)
 */
template<class T, class P>
class MemberFunctionOneParameterCallback : public Callback {

public:
    MemberFunctionOneParameterCallback(T* target, void (T::*fpt)(P), P param)
        : target_(target), fpt_(fpt), param_(param)
    {
    }

    virtual ~MemberFunctionOneParameterCallback() {}
    virtual MemberFunctionOneParameterCallback* clone() const { return new MemberFunctionOneParameterCallback(*this); }

    virtual void exec() {
        if ((target_ != 0) && (fpt_ != 0))
            (target_->*fpt_)(param_);
    }

private:
    T* target_;
    void (T::*fpt_)(P);
    P param_;
};

// ============================================================================

}   // namespace

#endif
