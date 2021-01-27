/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2021 University of Muenster, Germany,                        *
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

#ifndef VRN_LINKEVALUATORBOUNDINGBOX_H
#define VRN_LINKEVALUATORBOUNDINGBOX_H

#include "voreen/core/properties/link/linkevaluatorbase.h"
#include "voreen/core/properties/boundingboxproperty.h"
#include "voreen/core/properties/numeric/intervalproperty.h"
#include "voreen/core/properties/floatproperty.h"
#include "voreen/core/properties/intproperty.h"

namespace voreen {

//////////////////////////////////////////////////////////////////////////
//             TemplateLinkEvaluatorBoundingBoxComponent                //
//////////////////////////////////////////////////////////////////////////
template<typename T, int SIDE, int COMP>
class TemplateLinkEvaluatorBoundingBoxComponent : public LinkEvaluatorBase{
public:
    virtual void eval(Property* src, Property* dst);
    virtual bool arePropertiesLinkable(const Property* p1, const Property* p2) const;
};

class LinkEvaluatorIntBoundingBoxMinXComponent
    :public TemplateLinkEvaluatorBoundingBoxComponent<int, 0, 0> {
public:
    virtual std::string getClassName()  const { return "LinkEvaluatorIntBoundingBoxMinXComponent";     }
    virtual std::string getGuiName()    const { return "X component minimum";                       }
    virtual LinkEvaluatorBase* create() const { return new LinkEvaluatorIntBoundingBoxMinXComponent(); }
};

class LinkEvaluatorIntBoundingBoxMaxXComponent
    :public TemplateLinkEvaluatorBoundingBoxComponent<int, 1, 0> {
public:
    virtual std::string getClassName()  const { return "LinkEvaluatorIntBoundingBoxMaxXComponent";     }
    virtual std::string getGuiName()    const { return "X component maximum";                       }
    virtual LinkEvaluatorBase* create() const { return new LinkEvaluatorIntBoundingBoxMaxXComponent(); }
};

class LinkEvaluatorIntBoundingBoxMinYComponent
    :public TemplateLinkEvaluatorBoundingBoxComponent<int, 0, 1> {
public:
    virtual std::string getClassName()  const { return "LinkEvaluatorIntBoundingBoxMinYComponent";     }
    virtual std::string getGuiName()    const { return "Y component minimum";                       }
    virtual LinkEvaluatorBase* create() const { return new LinkEvaluatorIntBoundingBoxMinYComponent(); }
};

class LinkEvaluatorIntBoundingBoxMaxYComponent
    :public TemplateLinkEvaluatorBoundingBoxComponent<int, 1, 1> {
public:
    virtual std::string getClassName()  const { return "LinkEvaluatorIntBoundingBoxMaxYComponent";     }
    virtual std::string getGuiName()    const { return "Y component maximum";                       }
    virtual LinkEvaluatorBase* create() const { return new LinkEvaluatorIntBoundingBoxMaxYComponent(); }
};

class LinkEvaluatorIntBoundingBoxMinZComponent
    :public TemplateLinkEvaluatorBoundingBoxComponent<int, 0, 2> {
public:
    virtual std::string getClassName()  const { return "LinkEvaluatorIntBoundingBoxMinZComponent";     }
    virtual std::string getGuiName()    const { return "Z component minimum";                       }
    virtual LinkEvaluatorBase* create() const { return new LinkEvaluatorIntBoundingBoxMinZComponent(); }
};

class LinkEvaluatorIntBoundingBoxMaxZComponent
    :public TemplateLinkEvaluatorBoundingBoxComponent<int, 1, 2> {
public:
    virtual std::string getClassName()  const { return "LinkEvaluatorIntBoundingBoxMaxZComponent";     }
    virtual std::string getGuiName()    const { return "Z component maximum";                       }
    virtual LinkEvaluatorBase* create() const { return new LinkEvaluatorIntBoundingBoxMaxZComponent(); }
};

class LinkEvaluatorFloatBoundingBoxMinXComponent
    :public TemplateLinkEvaluatorBoundingBoxComponent<float, 0, 0> {
public:
    virtual std::string getClassName()  const { return "LinkEvaluatorFloatBoundingBoxMinXComponent";     }
    virtual std::string getGuiName()    const { return "X component minimum";                       }
    virtual LinkEvaluatorBase* create() const { return new LinkEvaluatorFloatBoundingBoxMinXComponent(); }
};

class LinkEvaluatorFloatBoundingBoxMaxXComponent
    :public TemplateLinkEvaluatorBoundingBoxComponent<float, 1, 0> {
public:
    virtual std::string getClassName()  const { return "LinkEvaluatorFloatBoundingBoxMaxXComponent";     }
    virtual std::string getGuiName()    const { return "X component maximum";                       }
    virtual LinkEvaluatorBase* create() const { return new LinkEvaluatorFloatBoundingBoxMaxXComponent(); }
};

class LinkEvaluatorFloatBoundingBoxMinYComponent
    :public TemplateLinkEvaluatorBoundingBoxComponent<float, 0, 1> {
public:
    virtual std::string getClassName()  const { return "LinkEvaluatorFloatBoundingBoxMinYComponent";     }
    virtual std::string getGuiName()    const { return "Y component minimum";                       }
    virtual LinkEvaluatorBase* create() const { return new LinkEvaluatorFloatBoundingBoxMinYComponent(); }
};

class LinkEvaluatorFloatBoundingBoxMaxYComponent
    :public TemplateLinkEvaluatorBoundingBoxComponent<float, 1, 1> {
public:
    virtual std::string getClassName()  const { return "LinkEvaluatorFloatBoundingBoxMaxYComponent";     }
    virtual std::string getGuiName()    const { return "Y component maximum";                       }
    virtual LinkEvaluatorBase* create() const { return new LinkEvaluatorFloatBoundingBoxMaxYComponent(); }
};

class LinkEvaluatorFloatBoundingBoxMinZComponent
    :public TemplateLinkEvaluatorBoundingBoxComponent<float, 0, 2> {
public:
    virtual std::string getClassName()  const { return "LinkEvaluatorFloatBoundingBoxMinZComponent";     }
    virtual std::string getGuiName()    const { return "Z component minimum";                       }
    virtual LinkEvaluatorBase* create() const { return new LinkEvaluatorFloatBoundingBoxMinZComponent(); }
};

class LinkEvaluatorFloatBoundingBoxMaxZComponent
    :public TemplateLinkEvaluatorBoundingBoxComponent<float, 1, 2> {
public:
    virtual std::string getClassName()  const { return "LinkEvaluatorFloatBoundingBoxMaxZComponent";     }
    virtual std::string getGuiName()    const { return "Z component maximum";                       }
    virtual LinkEvaluatorBase* create() const { return new LinkEvaluatorFloatBoundingBoxMaxZComponent(); }
};

//////////////////////////////////////////////////////////////////////////
//           TemplateLinkEvaluatorBoundingBoxComponentInterval          //
//////////////////////////////////////////////////////////////////////////
template<typename T, int COMP>
class TemplateLinkEvaluatorBoundingBoxComponentInterval : public LinkEvaluatorBase{
public:
    virtual void eval(Property* src, Property* dst);
    virtual bool arePropertiesLinkable(const Property* p1, const Property* p2) const;
};

class LinkEvaluatorFloatBoundingBoxXComponentInterval
    :public TemplateLinkEvaluatorBoundingBoxComponentInterval<float, 0> {
public:
    virtual std::string getClassName()  const { return "LinkEvaluatorFloatBoundingBoxXComponentInterval";     }
    virtual std::string getGuiName()    const { return "X component";                                    }
    virtual LinkEvaluatorBase* create() const { return new LinkEvaluatorFloatBoundingBoxXComponentInterval(); }
};

class LinkEvaluatorFloatBoundingBoxYComponentInterval
    :public TemplateLinkEvaluatorBoundingBoxComponentInterval<float, 1> {
public:
    virtual std::string getClassName()  const { return "LinkEvaluatorFloatBoundingBoxYComponentInterval";     }
    virtual std::string getGuiName()    const { return "Y component";                                    }
    virtual LinkEvaluatorBase* create() const { return new LinkEvaluatorFloatBoundingBoxYComponentInterval(); }
};

class LinkEvaluatorFloatBoundingBoxZComponentInterval
    :public TemplateLinkEvaluatorBoundingBoxComponentInterval<float, 2> {
public:
    virtual std::string getClassName()  const { return "LinkEvaluatorFloatBoundingBoxZComponentInterval";     }
    virtual std::string getGuiName()    const { return "Z component";                                    }
    virtual LinkEvaluatorBase* create() const { return new LinkEvaluatorFloatBoundingBoxZComponentInterval(); }
};

class LinkEvaluatorIntBoundingBoxXComponentInterval
    :public TemplateLinkEvaluatorBoundingBoxComponentInterval<int, 0> {
public:
    virtual std::string getClassName()  const { return "LinkEvaluatorIntBoundingBoxXComponentInterval";     }
    virtual std::string getGuiName()    const { return "X component";                                    }
    virtual LinkEvaluatorBase* create() const { return new LinkEvaluatorIntBoundingBoxXComponentInterval(); }
};

class LinkEvaluatorIntBoundingBoxYComponentInterval
    :public TemplateLinkEvaluatorBoundingBoxComponentInterval<int, 1> {
public:
    virtual std::string getClassName()  const { return "LinkEvaluatorIntBoundingBoxYComponentInterval";     }
    virtual std::string getGuiName()    const { return "Y component";                                    }
    virtual LinkEvaluatorBase* create() const { return new LinkEvaluatorIntBoundingBoxYComponentInterval(); }
};

class LinkEvaluatorIntBoundingBoxZComponentInterval
    :public TemplateLinkEvaluatorBoundingBoxComponentInterval<int, 2> {
public:
    virtual std::string getClassName()  const { return "LinkEvaluatorIntBoundingBoxZComponentInterval";     }
    virtual std::string getGuiName()    const { return "Z component";                                    }
    virtual LinkEvaluatorBase* create() const { return new LinkEvaluatorIntBoundingBoxZComponentInterval(); }
};

//////////////////////////////////////////////////////////////////////////
//     TemplateLinkEvaluatorBoundingBoxComponentIntervalWithBounds      //
//////////////////////////////////////////////////////////////////////////
template<typename T, int COMP>
class TemplateLinkEvaluatorBoundingBoxComponentIntervalWithBounds : public LinkEvaluatorBase{
public:
    virtual void eval(Property* src, Property* dst);
    virtual bool arePropertiesLinkable(const Property* p1, const Property* p2) const;
};

class LinkEvaluatorFloatBoundingBoxXComponentIntervalWithBounds
    :public TemplateLinkEvaluatorBoundingBoxComponentIntervalWithBounds<float, 0> {
public:
    virtual std::string getClassName()  const { return "LinkEvaluatorFloatBoundingBoxXComponentIntervalWithBounds";     }
    virtual std::string getGuiName()    const { return "X component with bounds";                                    }
    virtual LinkEvaluatorBase* create() const { return new LinkEvaluatorFloatBoundingBoxXComponentIntervalWithBounds(); }
};

class LinkEvaluatorFloatBoundingBoxYComponentIntervalWithBounds
    :public TemplateLinkEvaluatorBoundingBoxComponentIntervalWithBounds<float, 1> {
public:
    virtual std::string getClassName()  const { return "LinkEvaluatorFloatBoundingBoxYComponentIntervalWithBounds";     }
    virtual std::string getGuiName()    const { return "Y component with bounds";                                    }
    virtual LinkEvaluatorBase* create() const { return new LinkEvaluatorFloatBoundingBoxYComponentIntervalWithBounds(); }
};

class LinkEvaluatorFloatBoundingBoxZComponentIntervalWithBounds
    :public TemplateLinkEvaluatorBoundingBoxComponentIntervalWithBounds<float, 2> {
public:
    virtual std::string getClassName()  const { return "LinkEvaluatorFloatBoundingBoxZComponentIntervalWithBounds";     }
    virtual std::string getGuiName()    const { return "Z component with bounds";                                    }
    virtual LinkEvaluatorBase* create() const { return new LinkEvaluatorFloatBoundingBoxZComponentIntervalWithBounds(); }
};

class LinkEvaluatorIntBoundingBoxXComponentIntervalWithBounds
    :public TemplateLinkEvaluatorBoundingBoxComponentIntervalWithBounds<int, 0> {
public:
    virtual std::string getClassName()  const { return "LinkEvaluatorIntBoundingBoxXComponentIntervalWithBounds";     }
    virtual std::string getGuiName()    const { return "X component with bounds";                                    }
    virtual LinkEvaluatorBase* create() const { return new LinkEvaluatorIntBoundingBoxXComponentIntervalWithBounds(); }
};

class LinkEvaluatorIntBoundingBoxYComponentIntervalWithBounds
    :public TemplateLinkEvaluatorBoundingBoxComponentIntervalWithBounds<int, 1> {
public:
    virtual std::string getClassName()  const { return "LinkEvaluatorIntBoundingBoxYComponentIntervalWithBounds";     }
    virtual std::string getGuiName()    const { return "Y component with bounds";                                    }
    virtual LinkEvaluatorBase* create() const { return new LinkEvaluatorIntBoundingBoxYComponentIntervalWithBounds(); }
};

class LinkEvaluatorIntBoundingBoxZComponentIntervalWithBounds
    :public TemplateLinkEvaluatorBoundingBoxComponentIntervalWithBounds<int, 2> {
public:
    virtual std::string getClassName()  const { return "LinkEvaluatorIntBoundingBoxZComponentIntervalWithBounds";     }
    virtual std::string getGuiName()    const { return "Z component with bounds";                                    }
    virtual LinkEvaluatorBase* create() const { return new LinkEvaluatorIntBoundingBoxZComponentIntervalWithBounds(); }
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
//                       Template Implementation                        //
//                                                                      //
//////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////
//             TemplateLinkEvaluatorBoundingBoxComponent                //
//////////////////////////////////////////////////////////////////////////
template<typename T, int SIDE, int COMP>
void TemplateLinkEvaluatorBoundingBoxComponent<T, SIDE, COMP>::eval(Property* src, Property* dst) {
    TemplateBoundingBoxProperty<T>* boundingBoxProperty;
    NumericProperty<T>*       floatProperty;
    bool boundingBoxSource;

    boundingBoxProperty = dynamic_cast<TemplateBoundingBoxProperty<T>*>(src);
    if (boundingBoxProperty){
        floatProperty = static_cast<NumericProperty<T>*>(dst);
        boundingBoxSource = true;
    }else{
        boundingBoxProperty = dynamic_cast<TemplateBoundingBoxProperty<T>*>(dst);
        floatProperty = static_cast<NumericProperty<T>*>(src);
        boundingBoxSource = false;
    }

    typename tgt::Vector3<T> bounds[2];
    bounds[0] = boundingBoxProperty->get().getLLF();
    bounds[1] = boundingBoxProperty->get().getURB();

    if (boundingBoxSource){
        // BoundingBox -> Float
        floatProperty->set(bounds[SIDE][COMP]);
    }else{
        // Float -> BoundingBox
        bounds[SIDE][COMP] = floatProperty->get();
        tgt::TemplateBounds<T> boundsTgt(bounds[0], bounds[1]);
        boundingBoxProperty->set(boundsTgt);
    }
}

template<typename T, int SIDE, int COMP>
bool TemplateLinkEvaluatorBoundingBoxComponent<T, SIDE, COMP>::arePropertiesLinkable(const Property* p1, const Property* p2) const {
    tgtAssert(p1, "null pointer");
    tgtAssert(p2, "null pointer");

    if (dynamic_cast<const TemplateBoundingBoxProperty<T>*>(p1) && dynamic_cast<const NumericProperty<T>*>(p2))
        return true;
    if (dynamic_cast<const TemplateBoundingBoxProperty<T>*>(p2) && dynamic_cast<const NumericProperty<T>*>(p1))
        return true;
    return false;
}

//////////////////////////////////////////////////////////////////////////
//           TemplateLinkEvaluatorBoundingBoxComponentInterval          //
//////////////////////////////////////////////////////////////////////////
template<typename T, int COMP>
void TemplateLinkEvaluatorBoundingBoxComponentInterval<T, COMP>::eval(Property* src, Property* dst) {
    TemplateBoundingBoxProperty<T>*   boundingBoxProperty;
    IntervalProperty<T>* floatProperty;
    bool boundingBoxSource;

    boundingBoxProperty = dynamic_cast<TemplateBoundingBoxProperty<T>*>(src);
    if (boundingBoxProperty){
        floatProperty = static_cast<IntervalProperty<T>*>(dst);
        boundingBoxSource = true;
    }else{
        boundingBoxProperty = dynamic_cast<TemplateBoundingBoxProperty<T>*>(dst);
        floatProperty = static_cast<IntervalProperty<T>*>(src);
        boundingBoxSource = false;
    }

    typename tgt::Vector3<T> bounds[2];
    bounds[0] = boundingBoxProperty->get().getLLF();
    bounds[1] = boundingBoxProperty->get().getURB();

    if (boundingBoxSource){
        // BoundingBox -> Float
        floatProperty->set(tgt::Vector2<T>(bounds[0][COMP], bounds[1][COMP]));
    }else{
        // Float -> BoundingBox
        bounds[0][COMP] = floatProperty->get().x;
        bounds[1][COMP] = floatProperty->get().y;
        tgt::TemplateBounds<T> boundsTgt(bounds[0], bounds[1]);
        boundingBoxProperty->set(boundsTgt);
    }
}

template<typename T, int COMP>
bool TemplateLinkEvaluatorBoundingBoxComponentInterval<T, COMP>::arePropertiesLinkable(const Property* p1, const Property* p2) const {
    tgtAssert(p1, "null pointer");
    tgtAssert(p2, "null pointer");

    if (dynamic_cast<const TemplateBoundingBoxProperty<T>*>(p1) && dynamic_cast<const IntervalProperty<T>*>(p2))
        return true;
    if (dynamic_cast<const TemplateBoundingBoxProperty<T>*>(p2) && dynamic_cast<const IntervalProperty<T>*>(p1))
        return true;
    return false;
}

//////////////////////////////////////////////////////////////////////////
//     TemplateLinkEvaluatorBoundingBoxComponentIntervalWithBounds      //
//////////////////////////////////////////////////////////////////////////
template<typename T, int COMP>
void TemplateLinkEvaluatorBoundingBoxComponentIntervalWithBounds<T, COMP>::eval(Property* src, Property* dst) {
    TemplateBoundingBoxProperty<T>*   boundingBoxProperty;
    IntervalProperty<T>* floatProperty;
    bool boundingBoxSource;

    boundingBoxProperty = dynamic_cast<TemplateBoundingBoxProperty<T>*>(src);
    if (boundingBoxProperty){
        floatProperty = static_cast<IntervalProperty<T>*>(dst);
        boundingBoxSource = true;
    }else{
        boundingBoxProperty = dynamic_cast<TemplateBoundingBoxProperty<T>*>(dst);
        floatProperty = static_cast<IntervalProperty<T>*>(src);
        boundingBoxSource = false;
    }

    typename tgt::Vector3<T> bounds[2];
    bounds[0] = boundingBoxProperty->get().getLLF();
    bounds[1] = boundingBoxProperty->get().getURB();

    if (boundingBoxSource){
        // BoundingBox -> Float
        floatProperty->set(tgt::Vector2<T>(bounds[0][COMP], bounds[1][COMP]));
        floatProperty->setMaxValue(boundingBoxProperty->getMaxValue()[COMP]);
        floatProperty->setMinValue(boundingBoxProperty->getMinValue()[COMP]);
    }else{
        // Float -> BoundingBox
        bounds[0][COMP] = floatProperty->get().x;
        bounds[1][COMP] = floatProperty->get().y;
        tgt::TemplateBounds<T> boundsTgt(bounds[0], bounds[1]);
        typename tgt::Vector3<T> minValue = boundingBoxProperty->getMinValue();
        typename tgt::Vector3<T> maxValue = boundingBoxProperty->getMaxValue();
        minValue[COMP] = floatProperty->getMinValue();
        maxValue[COMP] = floatProperty->getMaxValue();

        boundingBoxProperty->setMaxValue(maxValue);
        boundingBoxProperty->setMinValue(minValue);
        boundingBoxProperty->set(boundsTgt);
    }
}

template<typename T, int COMP>
bool TemplateLinkEvaluatorBoundingBoxComponentIntervalWithBounds<T, COMP>::arePropertiesLinkable(const Property* p1, const Property* p2) const {
    tgtAssert(p1, "null pointer");
    tgtAssert(p2, "null pointer");

    if (dynamic_cast<const TemplateBoundingBoxProperty<T>*>(p1) && dynamic_cast<const IntervalProperty<T>*>(p2))
        return true;
    if (dynamic_cast<const TemplateBoundingBoxProperty<T>*>(p2) && dynamic_cast<const IntervalProperty<T>*>(p1))
        return true;
    return false;
}

} // namespace

#endif // VRN_LINKEVALUATORBOUNDINGBOX_H
