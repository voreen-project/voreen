/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2024 University of Muenster, Germany,                        *
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

#ifndef VRN_PLOTSELECTIONDIALOG_H
#define VRN_PLOTSELECTIONDIALOG_H

#include "plotselectionentrytablemodel.h"
#include "../datastructures/plotselection.h"

#include <QDialog>
#include <QWidget>
#include <QTableView>
#include <QPushButton>

namespace voreen {

// forward declarations
class PlotSelectionProperty;
class PlotData;

/**
 * Modal dialog for editing a PlotSelectionEntry
 **/
class PlotSelectionDialog : public QDialog {
    Q_OBJECT
public:
    /**
     * \brief   Constructor for a new Dialog offering editing opportunities for a PlotSelectionEntry.
     *
     * \param   entry               PlotSelectionEntry to edit (initial value), dialog will edit a copy
     * \param   plotData            PlotData reference for determining the column names
     * \param   entitiesProp        according PlotEntitiesProperty (needed for determining x,y and z column indices)
     * \param   threeDimensional    flag if plot is threedimensional or not
     * \param   parent              parent widget
     *
     * \note    This dialog will edit a copy of \a entry.
     **/
    PlotSelectionDialog(const PlotSelectionEntry& entry, const PlotData* plotData, const PlotEntitiesProperty& entitiesProp, bool threeDimensional, QWidget* parent = 0);

    /**
     * Returns edited PlotSelectionEntry.
     **/
    const PlotSelectionEntry& getEntry() const;

private slots:
    void clickedAddPredicate();
    void clickedRemovePredicate();

private:
    /// creates and layouts all widgets
    void createWidgets();

    // private members
    PlotSelectionEntry entry_;                      ///< PlotSelectionEntry to edit
    const PlotData* plotData_;                      ///< reference to current PlotData (used for column names)
    const PlotEntitiesProperty& entitiesProp_;      ///< according PlotEntitiesProperty (needed for determining x,y and z column indices)
    bool threeDimensional_;                         ///< flag if plot is threedimensional or not

    // Qt widgets:
    PlotSelectionEntryTableModel* proxyModel_;      ///< proxy table model for table_
    PlotSelectionEntryItemDelegate* itemDelegate_;  ///< item delegate for table_
    QTableView* table_;                             ///< table showing the predicates and their columns
    QPushButton* btnAddPredicate_;                  ///< button for adding a row
    QPushButton* btnRemovePredicate_;               ///< button for removing a row
    QPushButton* btnOK_;                            ///< button to accept changes
    QPushButton* btnCancel_;                        ///< button to reject changes
};


} // namespace

#endif // VRN_PLOTSELECTIONDIALOG_H
