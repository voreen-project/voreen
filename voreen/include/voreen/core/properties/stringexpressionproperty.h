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

#ifndef VRN_STRINGEXPRESSIONPROPERTY_H
#define VRN_STRINGEXPRESSIONPROPERTY_H

#include "voreen/core/properties/stringproperty.h"
#include <map>

namespace voreen {

class VRN_CORE_API StringExpressionProperty : public StringProperty {
public:
    StringExpressionProperty(const std::string& id, const std::string& guiText, const std::string& value = "",
        int invalidationLevel=Processor::INVALID_RESULT, Property::LevelOfDetail lod = Property::LOD_DEFAULT);
    StringExpressionProperty();
    virtual ~StringExpressionProperty();

    virtual Property* create() const;

    virtual std::string getClassName() const       { return "StringExpressionProperty"; }
    virtual std::string getTypeDescription() const { return "StringExpression"; }

    /**
     * Returns a list of the item names to be shown to the user.
     */
    virtual const std::vector<std::string> getGuiStrings() const;

    /**
     * Returns the placeholder assigned to the given item name.
     * If no item with this name ist found, an empty string is returned.
     */
    virtual std::string getPlaceHolder(const std::string& guiName) const;

    /**
     * Adds a new item (a name/placeholder combination) to the list.
     * If an item with this name already exists, it will be overwritten.
     * If at least one of the two strings is empty, nothing will be done.
     */
    virtual void addPlaceHolder(const std::string& guiName, const std::string& placeholder);

    /**
     * Sets the strings to mark the beginning and the end of a placeholder string.
     * The standard setting is "{" and "}".
     * The beginSign and the endSign have to be different and should be chosen well, since they
     * will always be interpreted as marking placeholders and may therefore not be used in actual text.
     * If the given markers could not be set (e.g. empty strings or the same string for beginMarker and endMarker), false will be returned.
     * Otherwise true is returned and the new markers are set.
     */
    virtual bool setPlaceHolderMarkers(const std::string& beginMarker, const std::string& endMarker);

    /**
     * Get the strings that mark the beginning and end of a placeholder string.
     */
    const std::pair<std::string, std::string> getPlaceHolderMarkers() const;

    /**
     * Constructs a placeholder string out of the given text by concatenating begin marker, the given text and end marker.
     */
    virtual std::string makePlaceHolder(const std::string& text) const;

    /**
     * Returns the text of a placeholder without the beginning and end marker
     * If not both markers are present in the string and at begin and end of the string, it is returned without change
     */
    virtual std::string getPlaceHolderText(const std::string& placeholder) const;

    /**
     * Deletes all guiName/placeholder combinations from the list
     */
    virtual void eraseItems();

    /**
     * To be called by Widget when Button has been clicked -> add GuiName and Placeholder to text field
     */
    void itemWasSelected(const std::string& guiName);

    /**
     * Returns a list of the placeholders that have been found in the text field
     */
    std::set<std::string> getPlaceholdersInText() const;

    /**
     * Gets a map of Placeholder-Texts and the corresponding text to replace them with
     * and returns the text, in which all placeholders have been replaced.
     */
    std::string replacePlaceHoldersInText(const std::map<std::string, std::string>&) const;


private:
    std::map<std::string, std::string> items_; ///< a map with guiNames and the corresponding placeHolders
    std::pair<std::string, std::string> placeHolderMarkers_; ///< the begin- and end-markers for placeholders
};

}

#endif // VRN_STRINGEXPRESSIONPROPERTY_H
