/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany,                        *
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

#include "voreen/core/properties/link/propertylink.h"

#include "voreen/core/voreenapplication.h"
#include "voreen/core/voreenmodule.h"
#include "voreen/core/properties/link/linkevaluatorbase.h"
#include "voreen/core/properties/link/linkevaluatorhelper.h"
#include "voreen/core/properties/link/linkevaluatorid.h"
#include "voreen/core/properties/property.h"
#include "voreen/core/properties/link/linkevaluatorhelper.h"
#include <vector>
#include <map>
#include <typeinfo>

namespace voreen {

const std::string PropertyLink::loggerCat_("voreen.PropertyLink");

std::vector<Property*> PropertyLink::visitedProperties_;

PropertyLink::PropertyLink(Property* src, Property* dest, LinkEvaluatorBase* linkEvaluator)
    : src_(src)
    , dest_(dest)
    , evaluator_(linkEvaluator)
{
    tgtAssert(src_, "No source property");
    tgtAssert(dest_, "No destination property");

    if (linkEvaluator)
        evaluator_ = linkEvaluator;
    else {
        std::vector<std::pair<std::string, std::string> > availableFunctions = LinkEvaluatorHelper::getCompatibleLinkEvaluators(src, dest);
        std::string evalType;
        for(std::vector<std::pair<std::string, std::string> >::iterator i=availableFunctions.begin(); i!=availableFunctions.end(); i++) {
            if(evalType == "")
                evalType = i->first;
            else {
                if(i->second == "id")
                    evalType = i->first;
            }
        }
        evaluator_ = LinkEvaluatorHelper::createEvaluator(evalType);
    }

    evaluator_->propertiesChanged(src_, dest_);
}

PropertyLink::PropertyLink()
    : src_(0)
    , dest_(0)
    , evaluator_(0)
{}

PropertyLink::~PropertyLink() {
    if (src_)
        src_->removeLink(this);

    if (!dest_->getLink(src_) || dest_->getLink(src_)->getLinkEvaluator() != evaluator_)
        delete evaluator_;
    else {
        // should never happen, but if it does, prevent crash
        LWARNING("destructor(): two links with the same evaluator instance: "
            << "src=" << src_->getFullyQualifiedID() << ", dest=" << dest_->getFullyQualifiedID());
    }
}

/*
 * Call to test if the link will work to find out about failures in an early state.
 */
bool PropertyLink::testPropertyLink() {
    tgtAssert(src_, "No source property");
    tgtAssert(dest_, "No destination property");
    tgtAssert(visitedProperties_.empty(), "Test can only be called, if visitedProperties_ are empty.");

    //push back the first property
    visitedProperties_.push_back(src_);

    try {
        executeLink();
    }
    catch (...) {
        //clean-up vector
        visitedProperties_.clear();
        return false;
    }

    //clean-up vector
    visitedProperties_.clear();

    return true;
}


void PropertyLink::executeLink() {

    tgtAssert(src_, "No source property");
    tgtAssert(dest_, "No destination property");

    // links without evaluator must not exist, but do a graceful check
    if (!evaluator_) {
        //tgtAssert(false, "PropertyLink::onChange(): no link evaluator");
        LERROR("onChange(): no link evaluator" << " (src='" + src_->getFullyQualifiedID()
            << "', dest='" << dest_->getFullyQualifiedID() << "')");
        return;
    }

    // check, if destintion has been visited
    if (std::find(visitedProperties_.begin(), visitedProperties_.end(), dest_) != visitedProperties_.end()) {
        return;
    }

    //always push destination
    visitedProperties_.push_back(dest_);

    try {
        evaluator_->eval(src_, dest_);
    }
    catch (const VoreenException& e) {
        throw VoreenException(std::string(e.what()) +
            " (src='" + src_->getFullyQualifiedID() +
            "', dest='" + dest_->getFullyQualifiedID() + "')");
    }
}

void PropertyLink::setLinkEvaluator(LinkEvaluatorBase* evaluator) {
    tgtAssert(evaluator, "Null pointer passed");

    if(evaluator_ && (evaluator_ != evaluator))
        delete evaluator_;

    evaluator_ = evaluator;

    if(evaluator_) {
        evaluator_->propertiesChanged(src_, dest_);
    }
}

LinkEvaluatorBase* PropertyLink::getLinkEvaluator() const {
    return evaluator_;
}

Property* PropertyLink::getSourceProperty() const {
    return src_;
}

Property* PropertyLink::getDestinationProperty() const {
    return dest_;
}

void PropertyLink::serialize(Serializer& s) const {
    // Serialize source property reference...
    s.serialize("SourceProperty", src_);

    // Serialize destination property reference...
    s.serialize("DestinationProperty", dest_);

    // Serialize link evaluator...
    s.serialize("Evaluator", evaluator_);
}

void PropertyLink::deserialize(Deserializer& s) {
    // Deserialize source property reference...
    s.deserialize("SourceProperty", src_);

    // Deserialize destination property reference...
    s.deserialize("DestinationProperty", dest_);

    // Was either the source or the destination property not deserialized?
    if (!src_ || !dest_) {
        std::string addOn;
        if (!src_ && !dest_) {
            addOn = "No source and destination.";
        }
        else if (src_) {
            addOn = "Link source: '";
            if (src_->getOwner())
                addOn += src_->getOwner()->getID() + "::";
            addOn += src_->getGuiName() + "'";
        }
        else if (dest_) {
            addOn = "Link dest: '";
            if (dest_->getOwner())
                addOn += dest_->getOwner()->getID() + "::";
            addOn += dest_->getGuiName() + "'";
        }
        s.raise(SerializationMemoryAllocationException("Property link could not be established. " + addOn));
    }

    src_->registerLink(this);

    // Deserialize link evaluator...
    s.deserialize("Evaluator", evaluator_);

    if (dest_->getLink(src_) && dest_->getLink(src_)->getLinkEvaluator() == evaluator_) {
        // this should never happen, but if it does, replace evaluator with fresh instance
        LWARNING("deserialize(): link has been assigned the same evaluator as its reverse link: "
            << "src=" << src_->getFullyQualifiedID() << ", dest=" << dest_->getFullyQualifiedID());
        evaluator_ = dynamic_cast<LinkEvaluatorBase*>(evaluator_->create());
        if (evaluator_) {
            LERROR(evaluator_->getClassName() << "::create() " << " did not return a LinkEvaluatorBase");
            delete evaluator_;
            evaluator_ = 0;
        }
    }

    if (evaluator_) {
        // auto-convert old LinkEvaluatorId:
        if (evaluator_->getClassName() == "LinkEvaluatorId") {
            std::vector<std::pair<std::string, std::string> > availableFunctions = LinkEvaluatorHelper::getCompatibleLinkEvaluators(src_, dest_);
            std::string evalType = "";
            for(std::vector<std::pair<std::string, std::string> >::iterator i=availableFunctions.begin(); i!=availableFunctions.end(); i++) {
                if(i->second == "id")
                    evalType = i->first;
            }
            if(!evalType.empty()) {
                //delete evaluator_;
                evaluator_ = LinkEvaluatorHelper::createEvaluator(evalType);
                LINFO("Replaced deprecated link evaluator with " << evaluator_->getClassName());
            }
            else {
                LERROR("Could not find and alternative for old LinkEvaluatorId between " << src_->getTypeDescription() << " and " << dest_->getTypeDescription());
            }
        }
        // --------------------------------
        // auto-convert old LinkEvaluatorIdNormalized:
        if (evaluator_->getClassName() == "LinkEvaluatorIdNormalized") {
            std::vector<std::pair<std::string, std::string> > availableFunctions = LinkEvaluatorHelper::getCompatibleLinkEvaluators(src_, dest_);
            std::string evalType = "";
            for(std::vector<std::pair<std::string, std::string> >::iterator i=availableFunctions.begin(); i!=availableFunctions.end(); i++) {
                if(i->second == "id normalized")
                    evalType = i->first;
            }
            if(!evalType.empty()) {
                //delete evaluator_;
                evaluator_ = LinkEvaluatorHelper::createEvaluator(evalType);
                LINFO("Replaced deprecated link evaluator with " << evaluator_->getClassName());
            }
            else {
                LERROR("Could not find and alternative for old LinkEvaluatorIdNormalized between " << src_->getTypeDescription() << " and " << dest_->getTypeDescription());
            }
        }

        evaluator_->propertiesChanged(src_, dest_);
    }
}


} // namespace
