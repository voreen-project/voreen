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

#ifndef VRN_PLOTSELECTION_H
#define VRN_PLOTSELECTION_H

#include "plotpredicate.h"
#include "voreen/core/io/serialization/serializable.h"

#include "tgt/vector.h"

#include <vector>
#include <utility>
#include <string>

namespace voreen {

/**
 * \brief   Encapsulates either a table position or a bunch of PlotPredicates.
 *
 * \note    The contained PlotPredicates may be 0!
 **/
class VRN_CORE_API PlotSelection : public Serializable {
public:
    /**
     * Constructs a new empty PlotSelection
     **/
    PlotSelection();

    /**
     * Constructs a new PlotSelection pointing to the table position \a tablePosition.
     **/
    PlotSelection(tgt::ivec2 tablePosition);

    /**
     * Copy constructor
     **/
    PlotSelection(const PlotSelection& rhs);

    /**
     * Destructor
     **/
    ~PlotSelection();

    /**
     * Assignment operator
     **/
    PlotSelection& operator=(PlotSelection rhs);

    /**
    * @see Property::serialize
    */
    virtual void serialize(Serializer& s) const;

    /**
    * @see Property::deserialize
    */
    virtual void deserialize(Deserializer& s);

    /**
     * \brief   Adds predicate \a pred associated with column \a column to this property.
     *
     * \param   column  column index to which \a pred shall later be applied
     * \param   pred    PlotPredicate to add, may be 0
     *
     * \note    Makes a deep copy of \a pred, previously contained predicate will be properly deleted.
     **/
    void addPredicate(int column, const PlotPredicate* pred);

    /**
     * \brief   Removes PlotPredicate in row number \a index.
     *
     * \param   index   row to remove
     **/
    void removeRow(size_t index);

    /**
     * Checks if this selection contains only a simple table position.
     **/
    bool isTablePosition() const;

    /**
     * Returns the selected table position.
     **/
    const tgt::ivec2& getTablePosition() const;

    /**
     * \brief   Returns all held selections via PlotPredicates, each associated with a column.
     *
     * \note    The pointers to the PlotPredicates may be 0!
     **/
    const std::vector<std::pair<int, PlotPredicate*> >& getSelection() const;

    /**
     * \brief   Sets the column index in row \a index to \a column.
     *
     * \param   index   row index of predicate to change
     * \param   column  new column index
     **/
    void setColumn(size_t index, int column);

    /**
     * \brief   Sets the PlotPredicate in row \a index to \a pred.
     *
     * \param   index   row index of predicate to change
     * \param   pred    PlotPredicate to set, may be 0
     *
     * \note    Makes a deep copy of \a pred, previously contained predicate will be properly deleted.
     **/
    void setPredicate(size_t index, const PlotPredicate* pred);

    /**
     * Returns a std::string representation of this selection using the column labels of \a pData.
     *
     * \param   pData   PlotData which will be used for getting the column labels.
     **/
    std::string toString(const PlotData& pData) const;

    /**
     * Equals predicate
     **/
    bool operator==(const PlotSelection& rhs) const;

private:

    bool isTablePositionFlag_;
    tgt::ivec2 tablePosition_;
    std::vector<std::pair<int, PlotPredicate*> > selection_;    ///< Note: held PlotPredicates may be 0!

};

/**
 * Small struct encapsulating a PlotSelection together with flags declaring how to apply it.
 **/
class VRN_CORE_API PlotSelectionEntry : public Serializable {
public:
    PlotSelection selection_; ///< the PlotSelection
    bool highlight_;          ///< flag whether this selection shall be highlighted
    bool renderLabel_;        ///< flag whether there shall be rendered labels for this selection
    bool zoomTo_;             ///< flag whether it shall be zoomed to this selection

    PlotSelectionEntry();
    PlotSelectionEntry(PlotSelection selection, bool highlight, bool renderLabel, bool zoomTo);

    bool operator==(const PlotSelectionEntry& rhs) const {
        return (   selection_ == rhs.selection_
                && renderLabel_ == rhs.renderLabel_
                && zoomTo_ == rhs.zoomTo_
                && highlight_ == rhs.highlight_);
    }

    virtual void serialize(Serializer& s) const;

    virtual void deserialize(Deserializer& s);
};

} // namespace voreen

#endif // VRN_PLOTSELECTION_H
