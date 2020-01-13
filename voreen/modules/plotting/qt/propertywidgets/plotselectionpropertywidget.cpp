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

#include "plotselectionpropertywidget.h"

#include "../plotselectiondialog.h"
#include "../../datastructures/plotpredicate.h"
#include "../../properties/plotselectionproperty.h"

#include "voreen/qt/widgets/voreentoolwindow.h"

#include <QDoubleSpinBox>
#include <QPushButton>
#include <QVBoxLayout>
#include <QGridLayout>
#include <QTableView>

#include <sstream>
#include <limits>

namespace voreen {

PlotSelectionPropertyWidget::PlotSelectionPropertyWidget(PlotSelectionProperty* prop, QWidget* parent)
    : QPropertyWidget(prop, parent)
    , property_(prop)
    , btnAdd_(new QPushButton(tr("Add")))
    , btnRemove_(new QPushButton(tr("Remove")))
    , btnRemoveUnused_(new QPushButton(tr("Remove unused")))
    , proxyModel_(0)
    , table_(new QTableView())
    {
    QWidget* mainWidget = new QWidget();
    QGridLayout* gl = new QGridLayout();

    table_->setMinimumHeight(160);
    table_->setSelectionBehavior(QAbstractItemView::SelectRows);
    gl->addWidget(table_, 0, 0, 1, 4);

    gl->addWidget(btnAdd_, 1, 1);
    gl->addWidget(btnRemove_, 1, 2);
    gl->addWidget(btnRemoveUnused_, 1, 3);

    connect(btnAdd_, SIGNAL(clicked()), this, SLOT(btnAddClicked()));
    connect(btnRemove_, SIGNAL(clicked()), this, SLOT(btnRemoveClicked()));
    connect(btnRemoveUnused_, SIGNAL(clicked()), this, SLOT(btnRemoveUnusedClicked()));
    connect(table_, SIGNAL(doubleClicked(const QModelIndex&)), this, SLOT(tableViewDoubleClicked(const QModelIndex&)));
    updateFromPropertySlot();

    mainWidget->setLayout(gl);
    addWidget(mainWidget);
}

void PlotSelectionPropertyWidget::updateFromPropertySlot() {
    delete proxyModel_;

    proxyModel_ = new PlotSelectionTableModel(property_, this);
    table_->setModel(proxyModel_);
    table_->verticalHeader()->hide();
    table_->resizeColumnsToContents();
    table_->resizeRowsToContents();
}

void PlotSelectionPropertyWidget::btnAddClicked() {
    // create and add new PlotSelectionEntry
    if (property_->getEntitiesProperty().dataValid()) {
        PlotSelectionEntry tmp(PlotSelection(), false, false, false);
        property_->add(tmp);
        // open dialog to edit created PlotSelectionEntry
        PlotSelectionDialog dialog(property_->getSelectionAt(property_->size()-1), property_->getPlotData(), property_->getEntitiesProperty(), property_->getThreeDimensional(), this);
        if (dialog.exec() == QDialog::Accepted) {
            property_->set(dialog.getEntry(), property_->size()-1);
        }
        proxyModel_->invalidate();
    }
}

void PlotSelectionPropertyWidget::btnRemoveClicked() {
    QModelIndexList mil = table_->selectionModel()->selectedRows();
    if (! mil.empty())
        property_->remove(property_->size()-1-mil.first().row());
    proxyModel_->invalidate();
}

void PlotSelectionPropertyWidget::btnRemoveUnusedClicked() {
    property_->removeUnused();
}

void PlotSelectionPropertyWidget::tableViewDoubleClicked(const QModelIndex& index) {
    int row = property_->size() - 1 - index.row();
    if (!property_->getSelectionAt(row).selection_.isTablePosition() && property_->getEntitiesProperty().dataValid()) {
        PlotSelectionDialog dialog(property_->getSelectionAt(row), property_->getPlotData(), property_->getEntitiesProperty(), property_->getThreeDimensional(), this);
        if (dialog.exec() == QDialog::Accepted) {
            property_->set(dialog.getEntry(), row);
        }
    }
}

} // namespace
