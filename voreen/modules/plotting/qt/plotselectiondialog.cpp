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

#include "plotselectiondialog.h"

#include "plotpredicatedialog.h"
#include "../properties/plotselectionproperty.h"

#include <QLabel>
#include <QGridLayout>
#include <QPushButton>
#include <QDoubleSpinBox>

namespace voreen {

PlotSelectionDialog::PlotSelectionDialog(const PlotSelectionEntry& entry,
                                         const PlotData* plotData,
                                         const PlotEntitiesProperty& entitiesProp,
                                         bool threeDimensional,
                                         QWidget* parent)
    : QDialog(parent)
    , entry_(entry)
    , plotData_(plotData)
    , entitiesProp_(entitiesProp)
    , threeDimensional_(threeDimensional)
    , proxyModel_(0)
    , itemDelegate_(0)
    , table_(0)
    , btnAddPredicate_(0)
    , btnRemovePredicate_(0)
    , btnOK_(0)
    , btnCancel_(0)
{
    tgtAssert(plotData_, "PlotSelectionDialog(): plotData is NULL");
    setObjectName("PlotSelectionDialog");
    createWidgets();
}

void PlotSelectionDialog::createWidgets() {
    QGridLayout* layout = new QGridLayout();

    proxyModel_ = new PlotSelectionEntryTableModel(entry_, plotData_);
    itemDelegate_ = new PlotSelectionEntryItemDelegate(entry_, plotData_, entitiesProp_, threeDimensional_);

    table_ = new QTableView();
    table_->setModel(proxyModel_);
    table_->setItemDelegate(itemDelegate_);
    table_->setSortingEnabled(false);
    table_->verticalHeader()->hide();
    table_->resizeColumnsToContents();
    layout->addWidget(table_, 0, 0, 1, 3);

    btnAddPredicate_ = new QPushButton(tr("Add Row"));
    layout->addWidget(btnAddPredicate_, 1, 1);
    btnRemovePredicate_ = new QPushButton(tr("Remove Row"));
    layout->addWidget(btnRemovePredicate_, 1, 2);

    btnOK_ = new QPushButton(tr("OK"));
    layout->addWidget(btnOK_, 3, 1);
    btnCancel_ = new QPushButton(tr("Cancel"));
    layout->addWidget(btnCancel_, 3, 2);

    connect(btnAddPredicate_, SIGNAL(clicked()), this, SLOT(clickedAddPredicate()));
    connect(btnRemovePredicate_, SIGNAL(clicked()), this, SLOT(clickedRemovePredicate()));
    connect(btnOK_, SIGNAL(clicked()), this, SLOT(accept()));
    connect(btnCancel_, SIGNAL(clicked()), this, SLOT(reject()));

    setLayout(layout);
}

const PlotSelectionEntry& PlotSelectionDialog::getEntry() const {
    return entry_;
}

void PlotSelectionDialog::clickedAddPredicate() {
    entry_.selection_.addPredicate(0, 0);
    proxyModel_->invalidate();
}

void PlotSelectionDialog::clickedRemovePredicate() {
    QModelIndexList mil = table_->selectionModel()->selectedIndexes();
    if (! mil.empty())
        entry_.selection_.removeRow(mil.first().row());
    proxyModel_->invalidate();
}

} // namespace
