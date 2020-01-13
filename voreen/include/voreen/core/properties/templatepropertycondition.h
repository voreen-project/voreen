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

#ifndef VRN_TEMPLATEPROPERTYCONDITION_H
#define VRN_TEMPLATEPROPERTYCONDITION_H

#include "voreen/core/properties/condition.h"

namespace voreen {

template<typename T> class TemplateProperty;

/**
 * Base for Conditions that observe a TemplateProperty and exec actions depending on that
 * Property's value
 */
template<class T, class S>
class TemplatePropertyCondition : public Condition {
public:
    TemplatePropertyCondition(TemplateProperty<T>* observed, const S& value)
        : observed_(observed), value_(value)
    {}

    virtual ~TemplatePropertyCondition() {}

    virtual bool met() const throw() = 0;

    void set(const S& value) { value_ = value; }
    const S& get() const { return value_; }

protected:
    TemplateProperty<T>* observed_;
    S value_;
};

// TODO find any use for this or throw it away...
// This is an attempt to use <functional> in Conditions
// BinaryOperatorCondition<T, T>(..., equal_to<T>, ...)
// would be the same as
// EqualCondition<T>(...)

template<class T, class S>
class BinaryOperatorCondition : public TemplatePropertyCondition<T, S> {
protected:
    using TemplatePropertyCondition<T, S>::observed_;
    using TemplatePropertyCondition<T, S>::value_;

public:
    BinaryOperatorCondition(TemplateProperty<T>* observed, const S& value,
                            std::binary_function<T, S, bool> function, const Callback& action,
                            const Callback& elseaction)
        : TemplatePropertyCondition<T, S>(observed, value, action, elseaction), function_(function)
    {}

    virtual ~BinaryOperatorCondition() {}
    virtual BinaryOperatorCondition* clone() const { return new BinaryOperatorCondition(*this); }

    virtual bool met() const throw() { return function_(observed_->get(), value_); }

private:
    std::binary_function<T, T, bool> function_;
};

/*template<class T>
class EqualCondition : public TemplatePropertyCondition<T, T> {
protected:
    using TemplatePropertyCondition<T, T>::observed_;
    using TemplatePropertyCondition<T, T>::value_;

public:
    EqualCondition(TemplateProperty<T>* observed, T value,
                   const Action& action, const Action& elseaction)
        : TemplatePropertyCondition<T, T>(observed, value, action, elseaction)
    {}

    virtual ~EqualCondition() {}
    virtual EqualCondition* clone() const { return new EqualCondition(*this); }

    virtual bool met() const throw() {
        if (observed_ == 0)
            return false;
        else
            return (value_ == observed_->get());
    }
};*/

template<class T>
class LessCondition : public TemplatePropertyCondition<T, T> {
protected:
    using TemplatePropertyCondition<T, T>::observed_;
    using TemplatePropertyCondition<T, T>::value_;

public:
    LessCondition(TemplateProperty<T>* observed, T value,
                  const Callback& action, const Callback& elseaction)
        : TemplatePropertyCondition<T, T>(observed, value, action, elseaction)
    {}

    virtual ~LessCondition() {}
    virtual LessCondition* clone() const { return new LessCondition(*this); }
    virtual bool met() const throw() {
        if (this->observed_ == 0)
            return false;
        else
            return (observed_->get() < value_);
    }
};

template<class T>
class LessEqualCondition : public TemplatePropertyCondition<T, T> {
protected:
    using TemplatePropertyCondition<T, T>::observed_;
    using TemplatePropertyCondition<T, T>::value_;

public:
    LessEqualCondition(TemplateProperty<T>* observed, T value,
                       const Callback& action, const Callback& elseaction)
        : TemplatePropertyCondition<T, T>(observed, value, action, elseaction)
    {}

    virtual ~LessEqualCondition() { }
    virtual LessEqualCondition* clone() const { return new LessEqualCondition(*this); }

    virtual bool met() const throw() {
        if (observed_ == 0)
            return false;
        else
            return (observed_->get() <= value_);
    }
};

template<class T>
class GreaterCondition : public TemplatePropertyCondition<T, T> {
protected:
    using TemplatePropertyCondition<T, T>::observed_;
    using TemplatePropertyCondition<T, T>::value_;

public:
    GreaterCondition(TemplateProperty<T>* observed, T value,
                     const Callback& action, const Callback& elseaction)
        : TemplatePropertyCondition<T, T>(observed, value, action, elseaction)
    {}

    virtual ~GreaterCondition() {}
    virtual GreaterCondition* clone() const { return new GreaterCondition(*this); }

    virtual bool met() const throw() {
        if (observed_ == 0)
            return false;
        else
            return (observed_->get() > value_);
    }
};

template<class T>
class GreaterEqualCondition : public TemplatePropertyCondition<T, T> {
protected:
    using TemplatePropertyCondition<T, T>::observed_;
    using TemplatePropertyCondition<T, T>::value_;

public:
    GreaterEqualCondition(TemplateProperty<T>* observed, T value,
                          const Callback& action, const Callback& elseaction)
        : TemplatePropertyCondition<T, T>(observed, value, action, elseaction)
    {}

    virtual ~GreaterEqualCondition() {}
    virtual GreaterEqualCondition* clone() const { return new GreaterEqualCondition(*this); }

    virtual bool met() const throw() {
        if (observed_ == 0)
            return false;
        else
            return (observed_->get() >= value_);
    }
};

template<class T>
class InCondition : public TemplatePropertyCondition<T, std::set<T> > {
protected:
    using TemplatePropertyCondition<T, std::set<T> >::observed_;
    using TemplatePropertyCondition<T, std::set<T> >::value_;

public:
    InCondition(TemplateProperty<T>* observed, std::set<T> value, const Callback& action,
                const Callback& elseaction)
        : TemplatePropertyCondition<T, std::set<T> >(observed, value, action, elseaction)
    {}

    virtual ~InCondition() {}
    virtual InCondition* clone() const { return new InCondition(*this); }

    virtual bool met() const throw() {
        return (value_.find(observed_->get()) != value_.end());
    }
};

}   // namespace

#endif
