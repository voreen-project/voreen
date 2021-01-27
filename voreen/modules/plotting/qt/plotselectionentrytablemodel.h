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

#ifndef VRN_PLOTSELECTIONENTRYTABLEMODEL_H
#define VRN_PLOTSELECTIONENTRYTABLEMODEL_H

#include <QItemDelegate>
#include <QModelIndex>
#include <QAbstractTableModel>
#include <QObject>
#include <QVBoxLayout>
#include <QHeaderView>

#include "../properties/plotentitiesproperty.h"

namespace voreen {

// forward declaration
class PlotSelectionEntry;
class PlotData;

/**
 * Table model for PlotSelectionEntry
 **/
class PlotSelectionEntryTableModel : public QAbstractTableModel {
    Q_OBJECT

public:
    /**
     * Constructs a new PlotSelectionEntryTableModel
     *
     * \param   prop        PlotSelectionEntry whose data shall be visualized
     * \param   plotData    pointer to PlotData (needed for column labels)
     * \param   parent      parent widget
     **/
    PlotSelectionEntryTableModel(PlotSelectionEntry& entry, const PlotData* plotData, QObject *parent=0);

    /**
     * Returns number of rows.
     **/
    int rowCount(const QModelIndex &parent) const;

    /**
     * Returns number of columns.
     **/
    int columnCount(const QModelIndex &parent) const;

    /**
     * Returns the data at index \a index for role \a role
     **/
    QVariant data(const QModelIndex &index, int role) const;

    /**
     * Sets the data at index \a index to \a value.
     *
     * \param   index   position where to set data
     * \param   value   new data
     * \param   role    edit role
     **/
    bool setData(const QModelIndex& index, const QVariant& value, int role);

    /**
     * Returns the header data.
     **/
    QVariant headerData(int section, Qt::Orientation orientation, int role) const;

    /**
     * Returns the item flags at index \a index
     **/
    Qt::ItemFlags flags(const QModelIndex& index) const;

    /**
     * Force reinitialization of all data in table model
     **/
    void invalidate();

private:
    PlotSelectionEntry& entry_;   ///< PlotSelectionEntry to edit
    const PlotData* plotData_;          ///< according PlotData with metadata
};


/**
 * ItemDelegate for PlotSelectionEntry
 **/
class PlotSelectionEntryItemDelegate : public QItemDelegate {
    Q_OBJECT

public:
    /**
     * Constructs a new PlotSelectionEntryItemDelegate
     *
     * \param   prop                PlotSelectionEntry whose data shall be visualized
     * \param   plotData            pointer to PlotData (needed for column labels)
     * \param   entitiesProp        according PlotEntitiesProperty (needed for determining x,y and z column indices)
     * \param   threeDimensional    flag if plot is threedimensional or not
     * \param   parent              parent widget
     **/
    PlotSelectionEntryItemDelegate(PlotSelectionEntry& entry, const PlotData* plotData, const PlotEntitiesProperty& entitiesProp, bool threeDimensional, QObject* parent = 0);

    /**
     * Creates and returns an editor widget for the cell at \a index.
     *
     * \param   parent  parent widget
     * \param   option  style options
     * \param   index   table position
     **/
    QWidget* createEditor(QWidget* parent, const QStyleOptionViewItem& option, const QModelIndex& index) const;

    /**
     * Fills the editor widget \a editor with appropriate data for table cell \a index.
     *
     * \param   editor  editor widget to fill
     * \param   index   table position
     **/
    void setEditorData(QWidget* editor, const QModelIndex& index) const;

    /**
     * Retrieves the edited data from \a editor and stores it in the model (usually a PlotSelectionEntryTableModel).
     *
     * \param   editor  editor widget to get the data drom
     * \param   model   model to store the data in
     * \param   index   table position
     **/
    void setModelData(QWidget* editor, QAbstractItemModel* model, const QModelIndex& index) const;

    /**
     * Recalculates the dimensions and geometry of the editor widget.
     **/
    void updateEditorGeometry(QWidget* editor, const QStyleOptionViewItem& option, const QModelIndex& index) const;

private:
    PlotSelectionEntry& entry_;                 ///< PlotSelectionEntry to edit
    const PlotData* plotData_;                  ///< according PlotData with metadata
    const PlotEntitiesProperty& entitiesProp_;  ///< according PlotEntitiesProperty (needed for determining x,y and z column indices)
    bool threeDimensional_;                     ///< flag if plot is threedimensional or not
};

}

#endif // VRN_PLOTSELECTIONENTRYTABLEMODEL_H
