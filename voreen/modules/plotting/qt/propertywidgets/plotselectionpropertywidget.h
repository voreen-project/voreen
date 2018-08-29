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

#ifndef VRN_PLOTSELECTIONPROPERTYWIDGET_H
#define VRN_PLOTSELECTIONPROPERTYWIDGET_H

#include "voreen/qt/widgets/property/qpropertywidget.h"

#include "../plotselectiontablemodel.h"

class QDoubleSpinBox;
class QLabel;
class QPushButton;
class QTableView;
class QModelIndex;

namespace voreen {

// forward declaration
class PlotSelectionProperty;

/**
 * Widget showing the selections held by PlotSelectionProperty.
 *
 * \note The selection vector is shown in reversed order.
 **/
class PlotSelectionPropertyWidget : public QPropertyWidget {
Q_OBJECT
public:
    /**
     * Constructs a new PlotSelectionPropertyWidget showing the selections of \a prop.
     *
     * \param   prop    parent property
     * \param   parent  parent widget
     **/
    PlotSelectionPropertyWidget(PlotSelectionProperty* prop, QWidget* parent = 0);

protected slots:
    virtual void updateFromPropertySlot();

private slots:
    void btnAddClicked();
    void btnRemoveClicked();
    void btnRemoveUnusedClicked();
    void tableViewDoubleClicked(const QModelIndex& index);

private:
    PlotSelectionProperty* property_;       ///< parent property

    QPushButton* btnAdd_;
    QPushButton* btnRemove_;
    QPushButton* btnRemoveUnused_;

    PlotSelectionTableModel* proxyModel_;   ///< proxy model of the QTableView
    QTableView* table_;                     ///< extende QTableView widget which shows the data in data_

};

} // namespace

#endif // VRN_PLOTSELECTIONPROPERTYWIDGET_H

