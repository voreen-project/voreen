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

#include "voreen/core/properties/stringexpressionproperty.h"
#include "voreen/core/utils/stringutils.h"

namespace voreen {

StringExpressionProperty::StringExpressionProperty(const std::string& id, const std::string& guiText,
                       const std::string& value, int invalidationLevel, Property::LevelOfDetail lod)
    : StringProperty(id, guiText, value, invalidationLevel, lod)
{
    placeHolderMarkers_.first = "{";
    placeHolderMarkers_.second = "}";
}

StringExpressionProperty::StringExpressionProperty()
{
    placeHolderMarkers_.first = "{";
    placeHolderMarkers_.second = "}";
}

StringExpressionProperty::~StringExpressionProperty() {

}

Property* StringExpressionProperty::create() const {
    return new StringExpressionProperty();
}

const std::vector<std::string> StringExpressionProperty::getGuiStrings() const {
    std::vector<std::string> result;
    std::map<std::string, std::string>::const_iterator iter;
    for (iter = items_.begin(); iter != items_.end(); ++iter) {
        result.push_back(iter->first);
    }
    return result;
}

std::string StringExpressionProperty::getPlaceHolder(const std::string& guiName) const {
    std::string result = "";
    std::map<std::string, std::string>::const_iterator iter;
    iter = items_.find(guiName);
    if (iter != items_.end())
        result = iter->second;

    return result;
}

void StringExpressionProperty::addPlaceHolder(const std::string& guiName, const std::string& placeholder) {

    if (guiName.empty() || placeholder.empty())
        return;
    else {
        items_.insert(std::make_pair(guiName, placeholder));
        //update the Widgets
        updateWidgets();
    }
}

bool StringExpressionProperty::setPlaceHolderMarkers(const std::string& beginMarker, const std::string& endMarker) {
    if (beginMarker.empty() || endMarker.empty() || beginMarker == endMarker)
        return false;
    else {
        placeHolderMarkers_.first = beginMarker;
        placeHolderMarkers_.second = endMarker;
    }

    return true;
}

std::string StringExpressionProperty::makePlaceHolder(const std::string& text) const {
    std::stringstream s;
    s << placeHolderMarkers_.first << text << placeHolderMarkers_.second;
    return s.str();
}

std::string StringExpressionProperty::getPlaceHolderText(const std::string& placeholder) const {
    if (!((placeholder.find(placeHolderMarkers_.first) == 0) && (placeholder.find(placeHolderMarkers_.second) == placeholder.length() - (placeHolderMarkers_.second.length()))))
        return placeholder;
    else
        return placeholder.substr(placeHolderMarkers_.first.length(),placeholder.length()-(placeHolderMarkers_.first.length()+placeHolderMarkers_.second.length()));
}

void StringExpressionProperty::eraseItems() {
    items_.clear();
    updateWidgets();
}

void StringExpressionProperty::itemWasSelected(const std::string& guiName) {
    invalidateOwner();

    set(get() + guiName + ": " + getPlaceHolder(guiName) + "\n");
}

const std::pair<std::string, std::string> StringExpressionProperty::getPlaceHolderMarkers() const {
    return placeHolderMarkers_;
}

std::set<std::string> StringExpressionProperty::getPlaceholdersInText() const {

    std::set<std::string> result;

    //get the expression
    std::string expression = get();

    size_t pos = 0; //starting position

    while (true) {
        //find placeholder begin marker
        size_t beginPos = expression.find(placeHolderMarkers_.first, pos);
        if (beginPos == std::string::npos)
            break;

        //find placeholder end marker
        size_t endPos = expression.find(placeHolderMarkers_.second, beginPos);
        if (endPos == std::string::npos)
            break;

        //check if there is no begin marker in between these two
        //if there are begin markers in between, the last begin marker is used
        size_t betweenPos = beginPos;
        while (true) {
            size_t bPos = expression.find(placeHolderMarkers_.first, betweenPos + placeHolderMarkers_.first.length());

            if ((bPos == std::string::npos) || (bPos >= endPos))
                break;

            if (bPos < endPos)
                betweenPos = bPos;
        }

        beginPos = betweenPos;

        //add placeholder to list
        std::string key = getPlaceHolderText(expression.substr(beginPos, (endPos+placeHolderMarkers_.second.length())-beginPos));
        result.insert(key);

        //change starting position for parsing
        pos = endPos + placeHolderMarkers_.second.length();

    }
    return result;
}

std::string StringExpressionProperty::replacePlaceHoldersInText(const std::map<std::string, std::string>& replacements) const {
    std::stringstream output;

    //get the expression
    std::string expression = get();

    size_t pos = 0; //starting position

    while (true) {
        //find placeholder begin marker
        size_t beginPos = expression.find(placeHolderMarkers_.first, pos);
        if (beginPos == std::string::npos) {
            output << expression.substr(pos);
            break;
        }
        //find placeholder end marker
        size_t endPos = expression.find(placeHolderMarkers_.second, beginPos);
        if (endPos == std::string::npos) {
            output << expression.substr(pos);
            break;
        }

        //check if there is no begin marker in between these two
        //if there are begin markers in between, the last begin marker is used
        size_t betweenPos = beginPos;
        while (true) {
            size_t bPos = expression.find(placeHolderMarkers_.first, betweenPos + placeHolderMarkers_.first.length());

            if ((bPos == std::string::npos) || (bPos >= endPos))
                break;

            if (bPos < endPos)
                betweenPos = bPos;
        }

        beginPos = betweenPos;

        //text before the new placeholder
        output << expression.substr(pos, beginPos-pos);

        //process placeholder
        std::string placeholder = expression.substr(beginPos, (endPos+placeHolderMarkers_.second.length())-beginPos);
        std::string key = getPlaceHolderText(placeholder);

        //get replacement value for this placeholder
        std::string rString = "";
        std::map<std::string, std::string>::const_iterator iter;
        iter = replacements.find(key);
        if (iter != replacements.end())
            rString = iter->second;

        output << rString;

        //change starting position for parsing
        pos = endPos + placeHolderMarkers_.second.length();

    }

    return output.str();
}

}   // namespace
