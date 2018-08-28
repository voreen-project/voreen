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

#ifndef VRN_PLOTCELL_H
#define VRN_PLOTCELL_H

#include "plotbase.h"
#include "voreen/core/io/serialization/serializable.h"

namespace voreen {

class AggregationFunction;

/**
 * \brief Single data cell storing either a plot_t or a string pointer.
 *
 * A PlotCellValue stores one piece of data in a union which is either a plot_t or a string pointer.
 **/
class VRN_CORE_API PlotCellValue : public Serializable {
public:
    /// standard constructor: value and tag will be 0, null-flag will be set
    PlotCellValue();

    /// copy constructor
    PlotCellValue(const PlotCellValue& value);
    /// creates a PlotCellValue from a plot_t
    PlotCellValue(plot_t value);
    /// creates a PlotCellValue from a string
    PlotCellValue(const std::string& tag);

    ~PlotCellValue();

    /// assignment operator which correctly handles the pointers
    PlotCellValue& operator=(PlotCellValue rhs);
    /// comparison = operator
    bool operator==(const PlotCellValue& rhs) const;

    /**
     * Comparison operator. Values are smaller than strings.
     */
    bool operator<(const PlotCellValue& rhs) const;
    /// comparison > operator
    bool operator>(const PlotCellValue& rhs) const;
    /// stream operator
    friend std::ostream& operator<<(std::ostream& ostr, const PlotCellValue& rhs);

    /// checks whether a String is stored in this PlotCell
    bool isTag() const {
        return isTag_;
    };

    /// checks whether a plot_t value type is stored in this PlotCell
    bool isValue() const {
        return isValue_;
    };

    /// checks if there is no data in this PlotCell
    bool isNull() const {
        return !(isTag_ || isValue_);
    };

    /// checks if this PlotCell is highlighted (selected)
    bool isHighlighted() const {
        return isHighlighted_;
    };

    /// sets the highlighted state
    void setHighlighted(bool highlightedState);

    /// gets the value stored in this PlotCell
    plot_t getValue() const;

    /// gets the pointer to the tag stored in this PlotCell
    std::string getTag() const;

    /// stores the value \a value in this cell
    void setValue(plot_t value);

    /// stores the string pointer \a tag in this cell
    void setTag(const std::string& tag);

    /// clears all data stored in this cell and sets isNull flag
    void clear();

    /// @see Serializable
    virtual void serialize(Serializer& s) const;

    /// @see Serializable
    virtual void deserialize(Deserializer& s);

private:
    /// union storing a plot_t value and a string pointer in one memory block
    union {
        plot_t value_;
        std::string* tag_;
    };

    bool isValue_;          ///< flag whether cell contains a plot_t value
    bool isTag_;            ///< flag whether cell contains a string
    bool isHighlighted_;    ///< flag whether cell is highlighted
};

/// forward declaration for PlotData class
class PlotData;

/**
* \brief Single data cell containing implicit data generated by AggregationFunctions.
**/
class PlotCellImplicit {
public:
    PlotCellImplicit(AggregationFunction* func, const PlotData* data, int column);

    ~PlotCellImplicit();

    /// gets the AggregationFunction stored in this PlotCell
    AggregationFunction* getAggregationFunction() const;

    /// returns a clone of the stored AggregationFunction in this PlotCell
    AggregationFunction* getAggregationFunctionClone() const;

    /// gets the value stored in this PlotCell
    plot_t getValue() const;

    /// checks if there is no data in this PlotCell
    bool isNull() const;

    /// clears all data stored in this cell and sets isNull flag
    void clear();

private:

    AggregationFunction* aggregationFunction_;      ///< AggregationFunction which specifies the value of this cell

    const PlotData* plotData_;                      ///< Pointer to const PlotData which contains the data to aggregate
    int column_;                                    ///< column of the PlotData which is to aggregate
};

}

#endif // VRN_PLOTCELL_H
