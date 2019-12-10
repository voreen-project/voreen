/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2019 University of Muenster, Germany,                        *
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

#ifndef VRN_STRINGLISTPROPERTY_H
#define VRN_STRINGLISTPROPERTY_H

#include "voreen/core/properties/templateproperty.h"

#include <vector>
#include <string>

namespace voreen {

#ifdef DLL_TEMPLATE_INST
template class VRN_CORE_API TemplateProperty<std::vector<int>>;
#endif

class VRN_CORE_API StringListProperty : public TemplateProperty<std::vector<int>> {

    friend class StringListPropertyWidget;

public:
    StringListProperty(const std::string& id, const std::string& guiText,
                       int invalidationLevel=Processor::INVALID_RESULT, Property::LevelOfDetail lod = Property::LOD_DEFAULT);
    StringListProperty();
    virtual ~StringListProperty() {}

    virtual Property* create() const {return new StringListProperty();}
    virtual std::string getClassName() const       { return "StringListProperty"; }
    virtual std::string getTypeDescription() const { return "String"; }

    /** @see Property::serialize */
    virtual void serialize(Serializer& s) const;
    /** @see Property::deserialize */
    virtual void deserialize(Deserializer& s);

    /** To set the TableUpdateFlag. */
    virtual void addWidget(PropertyWidget* widget);
    /** To reset the TableUpdateFlag. */
    virtual void updateWidgets();

    /** Enum used to update only needed parts of the widget */
    enum TableUpdateFlag {
        UPDATE_NONE      = 0,
        UPDATE_ROWS      = 1,
        UPDATE_SELECTION = 2,

        UPDATE_ALL = UPDATE_SELECTION | UPDATE_ROWS
    };

public:
    /** Returns the number of rows. */
    int getNumRows() const;
    /** Clears the rows. */
    virtual void reset();
    /** Adds a new row to the end of the table. */
    void addRow(const std::string& row, const tgt::vec3& color = tgt::vec3::zero);
    /** Get the selected row index or -1 if nothing is selected */
    const std::vector<int>& getSelectedRowIndices() const;
    inline const std::vector<int>& get() const { return getSelectedRowIndices(); }
    /** Selects a new row. */
    void setSelectedRowIndices(const std::vector<int>& rowIndices);
    inline void set(const std::vector<int>& value) { setSelectedRowIndices(value); }

private:
    //------------------
    //  Member
    //------------------
    std::vector<int> selectedRows_;                         ///< member containing the selected row
    std::vector<std::pair<std::string, tgt::vec3>> values_; ///< The stored value (table content)
    TableUpdateFlag neededTableUpdates_;                    ///< Enum containing th current table changes
};

    inline StringListProperty::TableUpdateFlag operator|(StringListProperty::TableUpdateFlag a, StringListProperty::TableUpdateFlag b)
    {return static_cast<StringListProperty::TableUpdateFlag>(static_cast<int>(a) | static_cast<int>(b));}
    inline StringListProperty::TableUpdateFlag operator&(StringListProperty::TableUpdateFlag a, StringListProperty::TableUpdateFlag b)
    {return static_cast<StringListProperty::TableUpdateFlag>(static_cast<int>(a) & static_cast<int>(b));}
    inline StringListProperty::TableUpdateFlag& operator|=(StringListProperty::TableUpdateFlag& a, StringListProperty::TableUpdateFlag b)
    { a = a | b; return a;}

} // namespace voreen

#endif //VRN_STRINGLISTPROPERTY_H
