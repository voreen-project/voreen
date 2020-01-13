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

#ifndef VRN_LINKEVALUATORIDGENERIC_H
#define VRN_LINKEVALUATORIDGENERIC_H

#include "voreen/core/properties/link/linkevaluatorbase.h"
#include "voreen/core/properties/templateproperty.h"

namespace voreen {

template<class T>
class LinkEvaluatorIdGeneric : public LinkEvaluatorBase {
public:
    virtual void eval(Property* src, Property* dst);

    virtual std::string getGuiName() const { return "Identity (value)"; }

    virtual bool arePropertiesLinkable(const Property* p1, const Property* p2) const;
};

template<class T>
void LinkEvaluatorIdGeneric<T>::eval(Property* src, Property* dst) {
    static_cast<TemplateProperty<T>*>(dst)->set(static_cast<TemplateProperty<T>*>(src)->get());
}

template<class T>
bool LinkEvaluatorIdGeneric<T>::arePropertiesLinkable(const Property* p1, const Property* p2) const {
    tgtAssert(p1, "null pointer");
    tgtAssert(p2, "null pointer");

    return (dynamic_cast<const TemplateProperty<T>*>(p1) && dynamic_cast<const TemplateProperty<T>*>(p2));
}

// ----------------------------------------------------------------------------

template<class T>
class LinkEvaluatorIdNumeric : public LinkEvaluatorBase {
public:
    virtual void eval(Property* src, Property* dst);

    virtual std::string getGuiName() const { return "Identity (with bounds)"; }

    virtual bool arePropertiesLinkable(const Property* p1, const Property* p2) const;
};

} // namespace

//HACK: resolving circular dependencie
#include "voreen/core/properties/numericproperty.h"

namespace voreen {

template<class T>
void LinkEvaluatorIdNumeric<T>::eval(Property* src, Property* dst) {
    static_cast<NumericProperty<T>*>(dst)->setMaxValue(static_cast<NumericProperty<T>*>(src)->getMaxValue());
    static_cast<NumericProperty<T>*>(dst)->setMinValue(static_cast<NumericProperty<T>*>(src)->getMinValue());
    static_cast<NumericProperty<T>*>(dst)->set(static_cast<NumericProperty<T>*>(src)->get());
}

template<class T>
bool LinkEvaluatorIdNumeric<T>::arePropertiesLinkable(const Property* p1, const Property* p2) const {
    tgtAssert(p1, "null pointer");
    tgtAssert(p2, "null pointer");

    return (dynamic_cast<const NumericProperty<T>*>(p1) && dynamic_cast<const NumericProperty<T>*>(p2));
}

// ----------------------------------------------------------------------------


template<class T, class R>
class LinkEvaluatorIdGenericConversion : public LinkEvaluatorBase {
public:
    virtual void eval(Property* src, Property* dst);

    virtual std::string getGuiName() const { return "Identity (type conversion)"; }

    virtual bool arePropertiesLinkable(const Property* p1, const Property* p2) const;
};

template<class T, class R>
void LinkEvaluatorIdGenericConversion<T, R>::eval(Property* src, Property* dst) {
    //Find out direction:
    TemplateProperty<R>* srcR = dynamic_cast<TemplateProperty<R>*>(src);
    TemplateProperty<T>* dstT = dynamic_cast<TemplateProperty<T>*>(dst);

    TemplateProperty<T>* srcT = dynamic_cast<TemplateProperty<T>*>(src);
    TemplateProperty<R>* dstR = dynamic_cast<TemplateProperty<R>*>(dst);

    if(srcR && dstT) {
        dstT->set(static_cast<T>(srcR->get()));
    }
    else if(srcT && dstR) {
        dstR->set(static_cast<R>(srcT->get()));
    }
    else {
        tgtAssert(false, "Should not get here!");
    }
}

template<class T, class R>
bool LinkEvaluatorIdGenericConversion<T, R>::arePropertiesLinkable(const Property* p1, const Property* p2) const {
    tgtAssert(p1, "null pointer");
    tgtAssert(p2, "null pointer");

    return ( (dynamic_cast<const TemplateProperty<T>*>(p1) && dynamic_cast<const TemplateProperty<R>*>(p2))
          || (dynamic_cast<const TemplateProperty<R>*>(p1) && dynamic_cast<const TemplateProperty<T>*>(p2)) );
}

// ----------------------------------------------------------------------------

template<class T>
class LinkEvaluatorIdNormalizedGeneric : public LinkEvaluatorBase {
public:
    virtual void eval(Property* src, Property* dst);

    virtual std::string getGuiName() const { return "Normalization"; }

    virtual bool arePropertiesLinkable(const Property* p1, const Property* p2) const;
};

template<class T>
void LinkEvaluatorIdNormalizedGeneric<T>::eval(Property* src, Property* dst) {
    NumericProperty<T>* srcConv = static_cast<NumericProperty<T>*>(src);
    NumericProperty<T>* dstConv = static_cast<NumericProperty<T>*>(dst);
    T diff_s = srcConv->getMaxValue() - srcConv->getMinValue();
    T diff_d = dstConv->getMaxValue() - dstConv->getMinValue();
    T result = dstConv->getMinValue() + (srcConv->get() - srcConv->getMinValue()) * diff_d / diff_s;
    dstConv->set(result);
}

template<class T>
bool LinkEvaluatorIdNormalizedGeneric<T>::arePropertiesLinkable(const Property* p1, const Property* p2) const {
    tgtAssert(p1, "null pointer");
    tgtAssert(p2, "null pointer");

    return (dynamic_cast<const NumericProperty<T>*>(p1) && dynamic_cast<const NumericProperty<T>*>(p2));
}

// ----------------------------------------------------------------------------

template<class T, class R>
class LinkEvaluatorIdNormalizedGenericConversion : public LinkEvaluatorBase {
public:
    virtual void eval(Property* src, Property* dst);

    virtual std::string getGuiName() const { return "Normalization (type conversion)"; }

    virtual bool arePropertiesLinkable(const Property* p1, const Property* p2) const;
};

template<class T, class R>
void LinkEvaluatorIdNormalizedGenericConversion<T, R>::eval(Property* src, Property* dst) {
    //Find out direction:
    NumericProperty<R>* srcR = dynamic_cast<NumericProperty<R>*>(src);
    NumericProperty<T>* dstT = dynamic_cast<NumericProperty<T>*>(dst);

    NumericProperty<T>* srcT = dynamic_cast<NumericProperty<T>*>(src);
    NumericProperty<R>* dstR = dynamic_cast<NumericProperty<R>*>(dst);

    if (srcR && dstT) {
        R diff_s = srcR->getMaxValue() - srcR->getMinValue();
        T diff_d = dstT->getMaxValue() - dstT->getMinValue();
        T result = dstT->getMinValue() + static_cast<T>(srcR->get() - srcR->getMinValue()) * diff_d / static_cast<T>(diff_s);
        dstT->set(result);
    }
    else if (srcT && dstR) {
        T diff_s = srcT->getMaxValue() - srcT->getMinValue();
        R diff_d = dstR->getMaxValue() - dstR->getMinValue();
        R result = dstR->getMinValue() + static_cast<R>(srcT->get() - srcT->getMinValue()) * diff_d / static_cast<R>(diff_s);
        dstR->set(result);
    }
    else {
        tgtAssert(false, "Should not get here!");
    }
}

template<class T, class R>
bool LinkEvaluatorIdNormalizedGenericConversion<T, R>::arePropertiesLinkable(const Property* p1, const Property* p2) const {
    tgtAssert(p1, "null pointer");
    tgtAssert(p2, "null pointer");

    return ( (dynamic_cast<const NumericProperty<T>*>(p1) && dynamic_cast<const NumericProperty<R>*>(p2))
        || (dynamic_cast<const NumericProperty<R>*>(p1) && dynamic_cast<const NumericProperty<T>*>(p2)) );
}

// ----------------------------------------------------------------------------

template<class T>
class LinkEvaluatorGenericStringId : public LinkEvaluatorBase {

public:
    virtual void eval(Property* src, Property* dst);

    virtual std::string getGuiName() const { return "Value"; }

    virtual bool arePropertiesLinkable(const Property* p1, const Property* p2) const;
};

template<class T>
void LinkEvaluatorGenericStringId<T>::eval(Property* src, Property* dst) {
    tgtAssert(src, "null pointer");
    tgtAssert(dst, "null pointer");
    dynamic_cast<TemplateProperty<std::string>*>(dst)->set(std::to_string(dynamic_cast<TemplateProperty<T>*>(src)->get()));
}

template<class T>
bool LinkEvaluatorGenericStringId<T>::arePropertiesLinkable(const Property* src, const Property* dst) const {
    tgtAssert(src, "null pointer");
    tgtAssert(dst, "null pointer");
    return (dynamic_cast<const TemplateProperty<T>*>(src) && dynamic_cast<const TemplateProperty<std::string>*>(dst));
}

// ----------------------------------------------------------------------------

template<class Vec, class Scalar, size_t dim>
class LinkEvaluatorIdGenericVecComponentConversion: public LinkEvaluatorBase {

    // ensure that dim is not larger that vector size
    static_assert(dim < Vec::size, "dim is larger than vector");

public:
    virtual void eval(Property* src, Property* dst);

    virtual std::string getGuiName() const { return "Component[" + std::to_string(dim) + "] (Value)"; }

    virtual bool arePropertiesLinkable(const Property* p1, const Property* p2) const;
};

template<class Vec, class Scalar, size_t dim>
void LinkEvaluatorIdGenericVecComponentConversion<Vec,Scalar,dim>::eval(Property* src, Property* dst) {
    //Find out direction:
    TemplateProperty<Vec>* srcVec = dynamic_cast<TemplateProperty<Vec>*>(src);
    TemplateProperty<Scalar>* dstScalar = dynamic_cast<TemplateProperty<Scalar>*>(dst);

    TemplateProperty<Scalar>* srcScalar = dynamic_cast<TemplateProperty<Scalar>*>(src);
    TemplateProperty<Vec>* dstVec = dynamic_cast<TemplateProperty<Vec>*>(dst);

    if(srcVec && dstScalar) {
        dstScalar->set(static_cast<Scalar>(srcVec->get()[dim]));
    }
    else if(srcScalar && dstVec) {
        Vec v = dstVec->get();
        v[dim] = static_cast<typename Vec::ElemType>(srcScalar->get());
        dstVec->set(v);
    }
    else {
        tgtAssert(false, "Should not get here!");
    }
}

template<class Vec, class Scalar, size_t dim>
bool LinkEvaluatorIdGenericVecComponentConversion<Vec,Scalar,dim>::arePropertiesLinkable(const Property* src, const Property* dst) const {
    tgtAssert(src, "null pointer");
    tgtAssert(dst, "null pointer");

    return ( (dynamic_cast<const TemplateProperty<Vec>*>(src) && dynamic_cast<const TemplateProperty<Scalar>*>(dst))
          || (dynamic_cast<const TemplateProperty<Scalar>*>(src) && dynamic_cast<const TemplateProperty<Vec>*>(dst)) );
}

// ----------------------------------------------------------------------------

template<class Vec, class Scalar, size_t dim>
class LinkEvaluatorIdNormalizedGenericVecComponentConversion: public LinkEvaluatorBase {

    // ensure that dim is not larger that vector size
    static_assert(dim < Vec::size, "dim is larger than vector");

public:
    virtual void eval(Property* src, Property* dst);

    virtual std::string getGuiName() const { return "Component[" + std::to_string(dim) + "] (Normalized)"; }

    virtual bool arePropertiesLinkable(const Property* p1, const Property* p2) const;
};

template<class Vec, class Scalar, size_t dim>
void LinkEvaluatorIdNormalizedGenericVecComponentConversion<Vec,Scalar,dim>::eval(Property* src, Property* dst) {
    //Find out direction:
    NumericProperty<Vec>* srcVec = dynamic_cast<NumericProperty<Vec>*>(src);
    NumericProperty<Scalar>* dstScalar = dynamic_cast<NumericProperty<Scalar>*>(dst);

    NumericProperty<Scalar>* srcScalar = dynamic_cast<NumericProperty<Scalar>*>(src);
    NumericProperty<Vec>* dstVec = dynamic_cast<NumericProperty<Vec>*>(dst);

    if(srcVec && dstScalar) {
        Scalar diff_src = static_cast<Scalar>(srcVec->getMaxValue().elem[dim] - srcVec->getMinValue().elem[dim]);
        Scalar diff_dst = dstScalar->getMaxValue() - dstScalar->getMinValue();
        Scalar result = dstScalar->getMinValue() + static_cast<Scalar>(srcVec->get().elem[dim] - srcVec->getMinValue().elem[dim]) * diff_dst / diff_src;
        dstScalar->set(result);
    }
    else if(srcScalar && dstVec) {
        Vec v = dstVec->get();
        typename Vec::ElemType diff_src = static_cast<typename Vec::ElemType>(srcScalar->getMaxValue() - srcScalar->getMinValue());
        typename Vec::ElemType diff_dst = dstVec->getMaxValue().elem[dim] - dstVec->getMinValue().elem[dim];
        v[dim] = dstVec->getMinValue().elem[dim] + static_cast<Scalar>(srcScalar->get() - srcScalar->getMinValue()) * diff_dst / diff_src;
        dstVec->set(v);
    }
    else {
        tgtAssert(false, "Should not get here!");
    }
}

template<class Vec, class Scalar, size_t dim>
bool LinkEvaluatorIdNormalizedGenericVecComponentConversion<Vec,Scalar,dim>::arePropertiesLinkable(const Property* src, const Property* dst) const {
    tgtAssert(src, "null pointer");
    tgtAssert(dst, "null pointer");

    return ( (dynamic_cast<const TemplateProperty<Vec>*>(src) && dynamic_cast<const TemplateProperty<Scalar>*>(dst))
          || (dynamic_cast<const TemplateProperty<Scalar>*>(src) && dynamic_cast<const TemplateProperty<Vec>*>(dst)) );
}

} // namespace

#endif
