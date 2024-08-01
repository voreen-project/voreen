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

#ifndef VRN_LINKEVALUATORMATRIXINVERT_H
#define VRN_LINKEVALUATORMATRIXINVERT_H

#include "voreen/core/properties/link/linkevaluatorbase.h"
#include "voreen/core/properties/templateproperty.h"

namespace voreen {

template<class T>
class LinkEvaluatorInvertGeneric : public LinkEvaluatorBase {
public:
    virtual void eval(Property* src, Property* dst);

    virtual std::string getGuiName() const { return "Inversion (value)"; }

    virtual bool arePropertiesLinkable(const Property* p1, const Property* p2) const;
};

template<class T>
void LinkEvaluatorInvertGeneric<T>::eval(Property* src, Property* dst) {
    auto val = static_cast<TemplateProperty<T>*>(src)->get();
    T linkValue;
    if(val.invert(linkValue)) {
        static_cast<TemplateProperty<T>*>(dst)->set(linkValue);
    } else {
        LWARNINGC("voreen.linkevaluatorinvert", "Inversion of linked property failed.");
    }
}

template<class T>
bool LinkEvaluatorInvertGeneric<T>::arePropertiesLinkable(const Property* p1, const Property* p2) const {
    tgtAssert(p1, "null pointer");
    tgtAssert(p2, "null pointer");

    return (dynamic_cast<const TemplateProperty<T>*>(p1) && dynamic_cast<const TemplateProperty<T>*>(p2));
}

/*
// tgt::mat2 currently does not have an invert-function.
class VRN_CORE_API LinkEvaluatorMat2Invert : public LinkEvaluatorInvertGeneric<tgt::mat2> {
public:
    virtual std::string getClassName() const { return "LinkEvaluatorMat2Invert"; }
    virtual LinkEvaluatorBase* create() const { return new LinkEvaluatorMat2Invert(); }
};
*/

class VRN_CORE_API LinkEvaluatorMat3Invert : public LinkEvaluatorInvertGeneric<tgt::mat3> {
public:
    virtual std::string getClassName() const { return "LinkEvaluatorMat3Invert"; }
    virtual LinkEvaluatorBase* create() const { return new LinkEvaluatorMat3Invert(); }
};

class VRN_CORE_API LinkEvaluatorMat4Invert : public LinkEvaluatorInvertGeneric<tgt::mat4> {
public:
    virtual std::string getClassName() const { return "LinkEvaluatorMat4Invert"; }
    virtual LinkEvaluatorBase* create() const { return new LinkEvaluatorMat4Invert(); }
};

} // namespace voreen

#endif // VRN_LINKEVALUATORMATRIXINVERT_H
