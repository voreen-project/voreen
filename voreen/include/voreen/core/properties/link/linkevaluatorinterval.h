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

#ifndef VRN_LINKEVALUATORINTERVAL_H
#define VRN_LINKEVALUATORINTERVAL_H

#include "voreen/core/properties/link/linkevaluatorbase.h"
#include "voreen/core/properties/numeric/intervalproperty.h"
#include "voreen/core/properties/floatproperty.h"
#include "voreen/core/properties/intproperty.h"

namespace voreen {


template<typename TYPE, int COMP>
class TemplateLinkEvaluatorIntervalProjection : public LinkEvaluatorBase{
public:
    virtual void eval(Property* src, Property* dst);
    virtual bool arePropertiesLinkable(const Property* p1, const Property* p2) const;
};


class LinkEvaluatorFloatIntervalMinProjection
    :public TemplateLinkEvaluatorIntervalProjection<float, 0> {
public:
    virtual std::string getClassName()  const { return "LinkEvaluatorFloatIntervalMinProjection";     }
    virtual std::string getGuiName()    const { return "Minimum component";                           }
    virtual LinkEvaluatorBase* create() const { return new LinkEvaluatorFloatIntervalMinProjection(); }
};

class LinkEvaluatorFloatIntervalMaxProjection
    :public TemplateLinkEvaluatorIntervalProjection<float, 1> {
public:
    virtual std::string getClassName()  const { return "LinkEvaluatorFloatIntervalMaxProjection";     }
    virtual std::string getGuiName()    const { return "Maximum component";                           }
    virtual LinkEvaluatorBase* create() const { return new LinkEvaluatorFloatIntervalMaxProjection(); }
};

class LinkEvaluatorIntIntervalMinProjection
    :public TemplateLinkEvaluatorIntervalProjection<int, 0> {
public:
    virtual std::string getClassName()  const { return "LinkEvaluatorIntIntervalMinProjection";     }
    virtual std::string getGuiName()    const { return "Minimum component";                         }
    virtual LinkEvaluatorBase* create() const { return new LinkEvaluatorIntIntervalMinProjection(); }
};

class LinkEvaluatorIntIntervalMaxProjection
    :public TemplateLinkEvaluatorIntervalProjection<int, 1> {
public:
    virtual std::string getClassName()  const { return "LinkEvaluatorIntIntervalMaxProjection";     }
    virtual std::string getGuiName()    const { return "Maximum component";                         }
    virtual LinkEvaluatorBase* create() const { return new LinkEvaluatorIntIntervalMaxProjection(); }
};

template<typename TYPE>
class TemplateLinkEvaluatorIntervalIdentityWithBounds : public LinkEvaluatorBase{
public:
    virtual void eval(Property* src, Property* dst);
    virtual bool arePropertiesLinkable(const Property* p1, const Property* p2) const;
};

class TemplateLinkEvaluatorIntIntervalIdentityWithBounds
    :public TemplateLinkEvaluatorIntervalIdentityWithBounds<int> {
public:
    virtual std::string getClassName()  const { return "TemplateLinkEvaluatorIntIntervalIdentityWithBounds";     }
    virtual std::string getGuiName()    const { return "Identity (With Bounds)";                                 }
    virtual LinkEvaluatorBase* create() const { return new TemplateLinkEvaluatorIntIntervalIdentityWithBounds(); }
};

class TemplateLinkEvaluatorFloatIntervalIdentityWithBounds
    :public TemplateLinkEvaluatorIntervalIdentityWithBounds<float> {
public:
    virtual std::string getClassName()  const { return "TemplateLinkEvaluatorFloatIntervalIdentityWithBounds";     }
    virtual std::string getGuiName()    const { return "Identity (With Bounds)";                                   }
    virtual LinkEvaluatorBase* create() const { return new TemplateLinkEvaluatorFloatIntervalIdentityWithBounds(); }
};


//////////////////////////////////////////////////////////////////////////
//                                                                      //
//                       Template Implementation                        //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////
//                 TemplateLinkEvaluatorIntervalProjection              //
//////////////////////////////////////////////////////////////////////////
template<typename TYPE, int COMP>
void TemplateLinkEvaluatorIntervalProjection<TYPE, COMP>::eval(Property* src, Property* dst) {
    IntervalProperty<TYPE>* intervalProperty = dynamic_cast<IntervalProperty<TYPE>*>(src);
    if (intervalProperty){
        // Interval -> Float
        TemplateProperty<TYPE> * scalarProperty = static_cast<TemplateProperty<TYPE>*>(dst);
        scalarProperty->set(intervalProperty->get()[COMP]);
    }else{
        // Float -> Interval
        intervalProperty = dynamic_cast<IntervalProperty<TYPE>*>(dst);
        TemplateProperty<TYPE>* scalarProperty = static_cast<TemplateProperty<TYPE>*>(src);
        typename tgt::Vector2<TYPE> val = intervalProperty->get();
        val[COMP] = scalarProperty->get();
        intervalProperty->set(val);
    }

}

template<typename TYPE, int COMP>
bool TemplateLinkEvaluatorIntervalProjection<TYPE, COMP>::arePropertiesLinkable(const Property* p1, const Property* p2) const {
    tgtAssert(p1, "null pointer");
    tgtAssert(p2, "null pointer");

    if (dynamic_cast<const IntervalProperty<TYPE>*>(p1) && dynamic_cast<const TemplateProperty<TYPE>*>(p2))
        return true;
    if (dynamic_cast<const IntervalProperty<TYPE>*>(p2) && dynamic_cast<const TemplateProperty<TYPE>*>(p1))
        return true;
    return false;
}


//////////////////////////////////////////////////////////////////////////
//            TemplateLinkEvaluatorIntervalIdentityWithBounds           //
//////////////////////////////////////////////////////////////////////////
template<typename TYPE>
void TemplateLinkEvaluatorIntervalIdentityWithBounds<TYPE>::eval(Property* src, Property* dst) {
    IntervalProperty<TYPE>* srcProp = dynamic_cast<IntervalProperty<TYPE>*>(src);
    IntervalProperty<TYPE>* dstProp = dynamic_cast<IntervalProperty<TYPE>*>(dst);

    dstProp->set(srcProp->get());
    dstProp->setMaxValue(srcProp->getMaxValue());
    dstProp->setMinValue(srcProp->getMinValue());

}

template<typename TYPE>
bool TemplateLinkEvaluatorIntervalIdentityWithBounds<TYPE>::arePropertiesLinkable(const Property* p1, const Property* p2) const {
    tgtAssert(p1, "null pointer");
    tgtAssert(p2, "null pointer");

    return dynamic_cast<const IntervalProperty<TYPE>*>(p1) && dynamic_cast<const IntervalProperty<TYPE>*>(p2);
}

} // namespace

#endif // VRN_LINKEVALUATORBOOLINVERT_H
