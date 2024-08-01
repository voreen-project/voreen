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

#ifndef VRN_LINKEVALUATORVECTORELEMENT_H
#define VRN_LINKEVALUATORVECTORELEMENT_H

#include "voreen/core/properties/link/linkevaluatorbase.h"
#include "voreen/core/properties/numericproperty.h"

namespace voreen {

template<typename V, typename E=V>
class LinkEvaluatorVectorElementFront : public LinkEvaluatorBase {

    using ElementType = NumericProperty<E>;
    using VectorType  = TemplateProperty<std::vector<V>>;

public:
    virtual ~LinkEvaluatorVectorElementFront() {}

    virtual std::string getGuiName() const { return "VectorElement (Front)"; }

    virtual void eval(Property* src, Property* dst) {
        //Find out direction:
        auto* srcR = dynamic_cast<VectorType*>(src);
        auto* dstT = dynamic_cast<ElementType*>(dst);

        auto* srcT = dynamic_cast<ElementType*>(src);
        auto* dstR = dynamic_cast<VectorType*>(dst);

        if (srcR && dstT) {
            if(!srcR->get().empty()) {
                dstT->set(srcR->get().front());
            }
        }
        else if (srcT && dstR) {
            dstR->set(std::vector<V>{srcT->get()});
        }
        else {
            tgtAssert(false, "Should not get here!");
        }
    };

    //Returns true if the LinkEvaluator can link the two properties.
    virtual bool arePropertiesLinkable(const Property* src, const Property* dst) const {
        return (dynamic_cast<const VectorType*>(src) && dynamic_cast<const ElementType*>(dst)) ||
            (dynamic_cast<const VectorType*>(dst) && dynamic_cast<const ElementType*>(src));
    }

};

template<typename V, typename E=V>
class LinkEvaluatorVectorElementBack : public LinkEvaluatorBase {

    using ElementType = NumericProperty<E>;
    using VectorType  = TemplateProperty<std::vector<V>>;

public:
    virtual ~LinkEvaluatorVectorElementBack() {}

    virtual std::string getGuiName() const { return "VectorElement (Back)"; }

    virtual void eval(Property* src, Property* dst) {
        //Find out direction:
        auto* srcR = dynamic_cast<VectorType*>(src);
        auto* dstT = dynamic_cast<ElementType*>(dst);

        auto* srcT = dynamic_cast<ElementType*>(src);
        auto* dstR = dynamic_cast<VectorType*>(dst);

        if (srcR && dstT) {
            if(!srcR->get().empty()) {
                dstT->set(srcR->get().back());
            }
        }
        else if (srcT && dstR) {
            dstR->set(std::vector<V>{srcT->get()});
        }
        else {
            tgtAssert(false, "Should not get here!");
        }
    };

    //Returns true if the LinkEvaluator can link the two properties.
    virtual bool arePropertiesLinkable(const Property* src, const Property* dst) const {
        return (dynamic_cast<const VectorType*>(src) && dynamic_cast<const ElementType*>(dst)) ||
               (dynamic_cast<const VectorType*>(dst) && dynamic_cast<const ElementType*>(src));
    }

};

class LinkEvaluatorIntVectorElementFront : public LinkEvaluatorVectorElementFront<int> {
public:
    virtual std::string getClassName() const { return "LinkEvaluatorIntVectorElementFront"; }
    virtual LinkEvaluatorBase* create() const { return new LinkEvaluatorIntVectorElementFront(); }
};

class LinkEvaluatorIntVectorElementBack : public LinkEvaluatorVectorElementBack<int> {
public:
    virtual std::string getClassName() const { return "LinkEvaluatorIntVectorElementBack"; }
    virtual LinkEvaluatorBase* create() const { return new LinkEvaluatorIntVectorElementBack(); }
};

} // namespace

#endif // VRN_LINKEVALUATORVECTORELEMENT_H
