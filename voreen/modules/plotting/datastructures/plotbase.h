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

#ifndef VRN_PLOTBASE_H
#define VRN_PLOTBASE_H

#include <string>
#include <vector>

#include "voreen/core/voreencoreapi.h"
#include "voreen/core/utils/exception.h"
#include "tgt/vector.h"

namespace voreen {

typedef double plot_t;

class VRN_CORE_API PlotBase {
public:
    /// type for the columns of a plotdata
    enum ColumnType {
        EMPTY = 0,
        STRING = 1,
        NUMBER = 2,
    };

    /**
     * \brief initialises keyColumnCount, dataColumnCount by given values and columnLabels with empty strings.
     *
     * \param keyColumnCount   count of key columns for this PlotBase
     * \param dataColumnCount  count of data columns for this PlotBase
     * \note    Be sure the parameter bounds are (>= 0) or this function might crash!
     **/
    PlotBase(int keyColumnCount, int dataColumnCount);


    /// copy constructor
    PlotBase(const PlotBase& rhs);

    /// destructor
    virtual ~PlotBase();

    /// returns the count of key columns
    int getKeyColumnCount() const;

    /// returns the count of data columns
    int getDataColumnCount() const;

    /// returns the total column count
    int getColumnCount() const;

    /**
     * Returns the label of column with index \a column, ordering: first key columns then data columns.
     *
     * \note    Be sure \a column is within bounds (>= 0 && < getColumnCount()) or this function might crash!
     **/
    const std::string& getColumnLabel(int column) const;

    /**
     * Sets the label of column with index \a column to \a label, ordering: first key columns then data columns.
     *
     * \note    Be sure \a column is within bounds (>= 0 && < getColumnCount()) or this function might crash!
     **/
    void setColumnLabel(int column, std::string label);

    void setColumnColorHint(int column, tgt::Color colorHint);
    void removeColumnColorHint(int column);
    bool hasColumnColorHint(int column) const;
    tgt::Color getColumnColorHint(int column) const;

    /**
     * Returns the Type of column with index \a column, ordering: first key columns then data columns.
     *
     * \note    Be sure \a column is within bounds (>= 0 && < getColumnCount()) or this function might crash!
     **/
    PlotBase::ColumnType getColumnType(int column) const;

    /// assignment operator which correctly handles the pointers
    PlotBase& operator=(const PlotBase& rhs);

    /**
     * \brief   Compares the structure with \a other.
     *
     * \param   other   PlotBase to compare with.
     **/
    bool compareStructure(const PlotBase* other) const;

protected:
    /// clears all data and resets key- and data column count
    void reset(int keyColumnCount, int dataColumnCount);

    /**
     * Sets the type of column with index \a column to \a type, ordering: first key columns then data columns.
     *
     * \note    Be sure \a column is within bounds (>= 0 && < getColumnCount()) or this function might crash!
     **/
    void setColumnType(int column, ColumnType type);

    int keyColumnCount_;          ///< count of key columns
    int dataColumnCount_;         ///< count of data columns

    struct Column {
        std::string label_;
        ColumnType type_;
        bool hasColorHint_;
        tgt::Color colorHint_;

        Column(std::string label, ColumnType type) : label_(label), type_(type), hasColorHint_(false) {}
        Column(std::string label, ColumnType type, tgt::Color colorHint) : label_(label), type_(type), hasColorHint_(true), colorHint_(colorHint) {}
        bool operator==(const Column& c) const
        {
            if((label_ == c.label_) && (type_ == c.type_) && (hasColorHint_ == c.hasColorHint_)) {
                if(hasColorHint_)
                    return colorHint_ == c.colorHint_;
                else
                    return true;
            }
            else
                return false;
        }
    };

    std::vector<Column> columns_;
};

} // namespace voreen

#endif // VRN_PLOTBASE_H
