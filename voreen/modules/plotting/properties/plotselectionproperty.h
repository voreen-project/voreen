/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany.                        *
 * Visualization and Computer Graphics Group <http://viscg.uni-muenster.de>        *
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

#ifndef VRN_PLOTSELECTIONPROPERTY_H
#define VRN_PLOTSELECTIONPROPERTY_H


#include "plotentitiesproperty.h"

#include "voreen/core/properties/templateproperty.h"
#include "../datastructures/plotselection.h"
#include "../datastructures/plotzoomstate.h"

#include "tgt/vector.h"
#include "tgt/logmanager.h"

#include <algorithm>
#include <iterator>
#include <string>
#include <utility>
#include <vector>

namespace voreen {

// forward declaration
class PlotData;

#ifdef DLL_TEMPLATE_INST
template class VRN_CORE_API TemplateProperty<std::vector<PlotSelectionEntry> >;
#endif

/**
 * This property holds a vector of plotzoomstates.
 */
class VRN_CORE_API PlotSelectionProperty : public TemplateProperty< std::vector< PlotSelectionEntry > > {
public:
    /**
     * Iterator adapter for the held PlotSelectionEntries iterating just over entries with renderLabel flag set.
     **/
    class LabelSelectionIterator : public std::iterator<std::bidirectional_iterator_tag, PlotSelection> {
    public:
        /**
         * Constructor for a new LabelSelectionIterator.
         * Current position will be the first element with renderLabel flag set in \a selectionEntries
         * or selectionEntries.end() if no such element exists.
         *
         * \param   selectionEntries    vector of PlotSelectionEntries to iterate over
         **/
        LabelSelectionIterator(const std::vector<PlotSelectionEntry>& selectionEntries);

        /**
         * Constructor for a new LabelSelectionIterator.
         *
         * \param   selectionEntries    vector of PlotSelectionEntries to iterate over
         * \param   position            start position of the iterator in \a selectionEntries
         **/
        LabelSelectionIterator(const std::vector<PlotSelectionEntry>& selectionEntries, std::vector<PlotSelectionEntry>::const_iterator position);

        bool operator==(LabelSelectionIterator const& rhs) const;
        bool operator!=(LabelSelectionIterator const& rhs) const;

        LabelSelectionIterator& operator++();
        LabelSelectionIterator operator++(int);

        LabelSelectionIterator& operator--();
        LabelSelectionIterator operator--(int);

        const PlotSelection& operator* () const;
        const PlotSelection* operator-> () const;

    private:
        const std::vector<PlotSelectionEntry>& selectionEntries_;
        std::vector<PlotSelectionEntry>::const_iterator position_;
    };


    /**
     * Constructor for a PlotSelectionProperty
     *
     * A PlotSelectionProperty can only exist in context with a PlotData object and some
     * information telling how to interprete this data. This is implemented by having references
     * to a PlotEntitiesProperty property defining the axis- and data columns and and to a bool
     * flag for the plot dimensionality.
     *
     * \param   pData               PlotData object containing the data to plot
     * \param   entitiesProp        PlotEntitiesProperty specifying how to interprete the data in \a pData
     * \param   threeDimensional    flag whether it is a 2D- or 3D-plot
     **/
    PlotSelectionProperty(const std::string& id, const std::string& guiText,
                          PlotData* pData, const PlotEntitiesProperty& entitiesProp, bool threeDimensional,
                          Processor::InvalidationLevel invalidationLevel = Processor::INVALID_RESULT);
    PlotSelectionProperty();

    virtual Property* create() const;

    virtual std::string getClassName() const       { return "PlotSelectionProperty"; }
    virtual std::string getTypeDescription() const { return "PlotSelection"; }

    /**
     * Adds the PlotSelection \a value to the property.
     * PlotSelection will be appended to the end of current PlotSelection vector.
     *
     * \param   value       PlotSelection to add
     * \param   additive    flag whether selection shall be added to (true) or replace (false) current selection
     **/
    void add(const PlotSelectionEntry& value, bool additive = false);

    /**
     * Removes PlotSelection with index \a index from property.
     *
     * \param   index   index in vector of PlotSelection to remove
     **/
    void remove(size_t index);

    /**
     * Removes unused PlotSelections from property.
     **/
    void removeUnused();

    /**
     * Removes all PlotSelection from property.
     **/
    void clear();

    void clearFlags(bool highlight = true, bool label = true, bool zoomTo = true);    ///< clears all flags
    void clearLastFlag();           ///< clears flag with highest index
    void clearZoomFlags();          ///< clears all zoom flags
    void clearLastZoomFlag();       ///< clears zoom flag with highest index
    void clearHighlightFlags();     ///< clears all highlight flags
    void clearLastHighlightFlag();  ///< clears highlight flag with highest index
    void clearLabelFlags();         ///< clears all label flags
    void clearLastLabelFlag();      ///< clears label flag with highest index

    /**
     * Returns const vector of current PlotSelections.
     **/
    const std::vector<PlotSelectionEntry>& get() const;

    /**
     * Returns the selection with index \a index from property.
     **/
    const PlotSelectionEntry& getSelectionAt(size_t index);

    /**
     * Returns a pointer to held PlotData
     **/
    const PlotData* getPlotData() const;

    /**
     * \brief   Sets to PlotData to \a plotData.
     *
     * \note    Be aware that \a plotData might be changed in this class (e.g. highlight flags).
     **/
    void setPlotData(PlotData* plotData);

    /**
     * Returns current PlotEntitiesProperty.
     **/
    const PlotEntitiesProperty& getEntitiesProperty() const;

    /**
     * Returns threeDimensional flag.
     **/
    bool getThreeDimensional() const;

    /**
     * Sets vector of PlotSelections.
     *
     * \param   value   new vector of PlotSelections of this property
     */
    void set(const std::vector<PlotSelectionEntry>& value);

    /**
     * Sets PlotSelectionsEntry.
     *
     * \param   entry   PlotSelectionEntry
     * \param   index   position in vector
     */
    void set(PlotSelectionEntry entry, int index);

    /**
     * Returns the number of PlotSelections in property.
     **/
    int size() const;

    /// Sets the dimension of the plot.
    void setDimension(bool threeDimensional);

    /// Sets the base zoom state.
    void setBaseZoomState(PlotZoomState state);

    /// Sets highlight flag at index \a index to \a value.
    void setHighlight(size_t index, bool value);

    /// Sets renderLabel flag at index \a index to \a value.
    void setRenderLabel(size_t index, bool value);

    /// Sets zoomTo flag at index \a index to \a value.
    void setZoomTo(size_t index, bool value);

    /// Returns whether to enable highlight flags and show them in the widget.
    bool getEnableHighlightFlags() const;

    /// Sets whether to enable highlight flags and show them in the widget.
    void setEnableHighlightFlags(bool value);

    /// Returns whether to enable label flags and show them in the widget.
    bool getEnableLabelFlags() const;

    /// Sets whether to enable label flags and show them in the widget.
    void setEnableLabelFlags(bool value);

    /// Returns whether to enable ZoomTo flags and show them in the widget.
    bool getEnableZoomToFlags() const;

    /// Sets whether to enable ZoomTo flags and show them in the widget.
    void setEnableZoomToFlags(bool value);

    /**
     * Returns the union of all current PlotSelections with zoomTo flag as PlotZoomState.
     **/
    const PlotZoomState& getZoom() const;

    /**
     * Sets highlight flag for all cells in currently held  PlotData object matching any PlotSelection with highlight flag.
     * Beware that for all other cells the highlight flag will be set to false.
     **/
    void applyHighlights();

    /**
     * Returns a LabelSelectionIterator to the begin of the labels.
     **/
    LabelSelectionIterator getLabelsBegin() const;

    /**
     * Returns a LabelSelectionIterator to the end of the labels.
     **/
    LabelSelectionIterator getLabelsEnd() const;

    /**
     * @see Property::serialize
     */
    virtual void serialize(Serializer& s) const;

    /**
     * @see Property::deserialize
     */
    virtual void deserialize(Deserializer& s);

    static const int X_AXIS_COLUMN = -1;        ///< sometimes you want to define the x axis instead of a column in PlotData - if so use this define instead
    static const int Y_AXIS_COLUMN = -2;        ///< sometimes you want to define the y axis instead of a column in PlotData - if so use this define instead
    static const int Z_AXIS_COLUMN = -3;        ///< sometimes you want to define the z axis instead of a column in PlotData - if so use this define instead

private:
    void updateZoomState();                     ///< recalculates currentZoomState_

    PlotEntitiesProperty dummyEntitiesProp_;

    PlotData* plotData_;                        ///< pointer to PlotData object of which this property holds selections - handle with care as it is non-const!
    const PlotEntitiesProperty& entitiesProp_;  ///< const reference to according plot settings property (e.g. needed for the meaning of the different PlotData columns)
    bool threeDimensional_;                     ///< flag specifying whether this are selections of a 2D- or 3D-plot

    PlotZoomState currentZoomState_;            ///< PlotZoomState built out of current PlotSelectionEntries with zoomTo flag set
    PlotZoomState baseZoomState_;               ///< automaticly calculated PlotZoomState, used if no zoomTo flags are set

    bool enableHighlightFlags_;                 ///< enable highlight flags and show them in widget
    bool enableLabelFlags_;                     ///< enable label flags and show them in widget
    bool enableZoomToFlags_;                    ///< enable zoom to flags and show them in widget

    static const std::string loggerCat_;

};

}   // namespace

#endif //VRN_PLOTSELECTIONPROPERTY_H
