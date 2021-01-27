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

#include "plotpredicateproperty.h"

#include "../datastructures/plotpredicate.h"
#include <sstream>


namespace voreen {

const std::string PlotPredicateProperty::loggerCat_("voreen.Properties.PlotPredicateProperty");

PlotPredicateProperty::PlotPredicateProperty(const std::string& id, const std::string& guiText, Processor::InvalidationLevel invalidationLevel)
    : TemplateProperty<std::vector<std::pair<int, PlotPredicate*> > >(id, guiText, std::vector<std::pair<int, PlotPredicate*> >(), invalidationLevel)
{
    value_ = std::vector<std::pair<int, PlotPredicate*> >();
}

PlotPredicateProperty::PlotPredicateProperty() {
}

PlotPredicateProperty::~PlotPredicateProperty() {
    clear();
}

Property* PlotPredicateProperty::create() const {
    return new PlotPredicateProperty();
}

void PlotPredicateProperty::clear() {
    for (size_t i=0; i<value_.size(); ++i)
        delete value_[i].second;
    value_.clear();
    notifyAll();;
}

const std::vector<std::pair<int, PlotPredicate*> >& PlotPredicateProperty::get() const {
    return value_;
}

std::vector<std::pair<int,PlotPredicate*> > PlotPredicateProperty::getCloned() const {
    std::vector<std::pair<int,PlotPredicate*> > toReturn;
    for (size_t i = 0; i < value_.size(); ++i) {
        toReturn.push_back(std::make_pair(value_.at(i).first, value_.at(i).second->clone()));
    }
    return toReturn;
}

void PlotPredicateProperty::set(const std::vector<std::pair<int, PlotPredicate*> >& value) {
    clear();
    for (size_t i = 0; i < value.size(); ++i) {
        value_.push_back(std::make_pair(value.at(i).first, value.at(i).second->clone()));
    }
    notifyAll();
}

int PlotPredicateProperty::size() const {
    return static_cast<int>(value_.size());
}

void PlotPredicateProperty::select(const std::vector<int>& columns){
    std::vector<std::pair<int, PlotPredicate*> >::iterator predit;
    bool found;
    int pos;
    predit = value_.begin();
    while (predit < value_.end()) {
        found = false;
        pos = -1;
        for (size_t i = 0; i < columns.size(); ++i) {
            if (predit->first == columns[i]) {
                found = true;
                pos = static_cast<int>(i);
                break;
            }
        }
        if (found) {
            predit->first = pos;
            ++predit;
        }
        else {
            delete predit->second;
            predit = value_.erase(predit);
        }
    }
}

void PlotPredicateProperty::notifyAll() {
    std::vector<std::pair<int, PlotPredicate*> > cpy = value_;

    onChangeCallbacks_.execute();

    updateWidgets();

    // invalidate owner:
    invalidateOwner();
}


void PlotPredicateProperty::serialize(Serializer& s) const {
    Property::serialize(s);

    //std::vector<SerializablePair<int, PlotPredicate*> > serializableVector;
    //for (size_t i = 0; i < value_.size(); ++i) {
    //    serializableVector.push_back(SerializablePair<int,PlotPredicate*>(value_.at(i).first,value_.at(i).second));
    //}
    //s.serialize("PredicateList",serializableVector);

    /*    // serialize the properties of the processor
        typedef std::map<std::string, Property*> PropertyMapType;

        PropertyMapType propertyMap;
        for (Properties::const_iterator it = properties_.begin(); it != properties_.end(); ++it)
            propertyMap[(*it)->getID()] = *it;

        const bool usePointerContentSerialization = s.getUsePointerContentSerialization();
        s.setUsePointerContentSerialization(true);
        s.serialize("ElementProperties", propertyMap, "Property", "name");
        s.setUsePointerContentSerialization(usePointerContentSerialization);*/
}

void PlotPredicateProperty::deserialize(Deserializer& s) {
    Property::deserialize(s);

    /*    // deserialize the properties of the processor
        typedef std::map<std::string, Property*> PropertyMapType;

        PropertyMapType propertyMap;
        for (Properties::const_iterator it = properties_.begin(); it != properties_.end(); ++it)
            propertyMap[(*it)->getID()] = *it;

        const bool usePointerContentSerialization = s.getUsePointerContentSerialization();
        s.setUsePointerContentSerialization(true);
        s.deserialize("ElementProperties", propertyMap, "Property", "name");
        s.setUsePointerContentSerialization(usePointerContentSerialization);*/
}

}   // namespace
