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

#ifndef VRN_STRINGTABLEPROPERTY_H
#define VRN_STRINGTABLEPROPERTY_H

#include "voreen/core/properties/property.h"

#include <vector>
#include <string>

namespace voreen {

    /**
     * Property used as user information for list data.
     * For each element of a list, the user can add a new row to the property.
     * The number of columns is defined by the constructor.
     * An example of usage looks like this:
     *
     * StringTableProperty property("id","gui",2);
     *
     * property.changeColumnLabel(0,"First");
     * property.changeColumnLabel(1,"Second");
     *
     * std::vector<std::string> vec(2);
     * vec[0] = "top left"; vec[0] = "top right";
     * property.addRow(vec);
     * vec[0] = "bottom left"; vec[0] = "bottom right";
     * property.addRow(vec);
     *
     * TODO: Add dynamic behavior. Add and remove columns etc..
     */
class VRN_CORE_API StringTableProperty : public Property {

    friend class StringTablePropertyWidget;

public:
    StringTableProperty(const std::string& id, const std::string& guiText, int columnCount,
        int invalidationLevel=Processor::INVALID_RESULT, Property::LevelOfDetail lod = Property::LOD_DEFAULT);
    StringTableProperty();
    virtual ~StringTableProperty() {}

    virtual Property* create() const {return new StringTableProperty();}
    virtual std::string getClassName() const       { return "StringtableProperty"; }
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
        UPDATE_LABELS    = 2,
        UPDATE_SELECTION = 4,

        UPDATE_ALL = UPDATE_NONE | UPDATE_SELECTION | UPDATE_ROWS | UPDATE_LABELS
    };

    //------------------
    //  Column Functions
    //------------------
public:
    /** Returns the number of columns created in the constructor. */
    int getNumColumns() const;
    /** Change the label of a column. Starts at 0. */
    void setColumnLabel(size_t number, const std::string& label);
protected:
    /** Returns all column labels. */
    const std::vector<std::string>& getColumnLabels() const;

    //------------------
    //  Row    Functions
    //------------------
public:
    /** Returns the number of rows. */
    int getNumRows() const;
    /** Clears the rows. */
    virtual void reset();
    /** Adds a new row to the end of the table. */
    void addRow(const std::vector<std::string>& row);
    /** Removes the specified row. */
    void removeRow(int rowIndex);
    /** Get the selected row index or -1 if nothing is selected */
    int getSelectedRowIndex() const;
    /** Selects a new row. */
    void setSelectedRowIndex(int rowIndex);


private:
    //------------------
    //  Member
    //------------------
    int columnCount_;                      ///< number of columns
    int selectedRow_;                               ///< member containing the selected row
    std::vector<std::string> columnLabels_;         ///< labels of each column
    std::vector<std::vector<std::string> > values_; ///< The stored value (table content)
    TableUpdateFlag neededTableUpdates_;            ///< Enum containing th current table changes
};

    inline StringTableProperty::TableUpdateFlag operator|(StringTableProperty::TableUpdateFlag a, StringTableProperty::TableUpdateFlag b)
    {return static_cast<StringTableProperty::TableUpdateFlag>(static_cast<int>(a) | static_cast<int>(b));}
    inline StringTableProperty::TableUpdateFlag operator&(StringTableProperty::TableUpdateFlag a, StringTableProperty::TableUpdateFlag b)
    {return static_cast<StringTableProperty::TableUpdateFlag>(static_cast<int>(a) & static_cast<int>(b));}
    inline StringTableProperty::TableUpdateFlag& operator|=(StringTableProperty::TableUpdateFlag& a, StringTableProperty::TableUpdateFlag b)
    { a = a | b; return a;}

} // namespace voreen

#endif //VRN_STRINGTABLEPROPERTY_H
