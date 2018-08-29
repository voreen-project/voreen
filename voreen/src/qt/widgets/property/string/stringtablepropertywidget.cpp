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

#include "voreen/qt/widgets/property/string/stringtablepropertywidget.h"

#include "voreen/qt/widgets/customlabel.h"
#include "voreen/core/utils/stringutils.h"

#include <QTableWidget>
#include <QHeaderView>

namespace voreen {

StringTablePropertyWidget::StringTablePropertyWidget(StringTableProperty* prop, QWidget* parent)
    : QPropertyWidget(prop, parent)
    , table_(0)
    , property_(prop)
{
    tgtAssert(property_,"null pointer");

    table_ = new QTableWidget(0,property_->getNumColumns());
    addWidget(table_);
    table_->setMaximumHeight(150); table_->setMinimumHeight(150);   // set to nice size
    table_->setEditTriggers(QAbstractItemView::NoEditTriggers);     // removes edit
    table_->verticalHeader()->hide();                               //removes vertical headers
    table_->setSelectionBehavior(QAbstractItemView::SelectRows);    //select entire row
    table_->setSelectionMode( QAbstractItemView::SingleSelection ); //and only one at a time
    table_->horizontalHeader()->setHighlightSections(false);        //headers are not selected
    table_->horizontalHeader()->setSectionsClickable(false);                //headers are not clickable
    table_->horizontalHeader()->setSectionResizeMode(QHeaderView::ResizeToContents); // adjust headers to content size
    table_->horizontalHeader()->setStretchLastSection(true);        // except the last one...

    connect(table_,SIGNAL(itemSelectionChanged()),this,SLOT(selectionOnChange()));

    //update
    updateTable();
}

StringTablePropertyWidget::~StringTablePropertyWidget() {
}


//----------------------------------------------------------
//  Overrides
//----------------------------------------------------------
CustomLabel* StringTablePropertyWidget::getOrCreateNameLabel() const {
    // prevents name label from being shown left of the property widget
    return 0;
}

void StringTablePropertyWidget::updateFromPropertySlot() {
    updateTable();
}

//----------------------------------------------------------
//  Helpers
//----------------------------------------------------------
void StringTablePropertyWidget::updateTable() {
    tgtAssert(property_,"no property");
    table_->blockSignals(true);

    if(property_->neededTableUpdates_ & StringTableProperty::UPDATE_LABELS) {
        QStringList labels;
        for(size_t i = 0; i < property_->getColumnLabels().size(); i++)
            labels << property_->getColumnLabels()[i].c_str();
        table_->setHorizontalHeaderLabels(labels);
    }

    if(property_->neededTableUpdates_ & StringTableProperty::UPDATE_ROWS) {
        table_->setRowCount(0); //clears rows
        for(size_t r = 0; r < property_->values_.size(); r++) {
            table_->insertRow(r);
            for(int c = 0; c < (int) property_->columnCount_; c++) {
                table_->setItem(r,c,new QTableWidgetItem(property_->values_[r][c].c_str()));
            }
        }
    }

    if(property_->neededTableUpdates_ & StringTableProperty::UPDATE_SELECTION) {
        tgtAssert(table_->rowCount() > property_->getSelectedRowIndex(), "not enough rows");
        table_->selectRow(property_->getSelectedRowIndex());
    }

    table_->blockSignals(false);
    //property_->neededTableUpdates_ = StringTableProperty::UPDATE_NONE;// done in updateWidgets
}


void StringTablePropertyWidget::selectionOnChange() {
    if(!disconnected_) {
        int index = -1;
        if(!table_->selectedItems().empty()) {
            index = table_->selectedItems().at(0)->row();
        }
        property_->setSelectedRowIndex(index);
    }
}

} // namespace
