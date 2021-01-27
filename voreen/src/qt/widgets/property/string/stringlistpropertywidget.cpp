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

#include "include/voreen/qt/widgets/property/string/stringlistpropertywidget.h"

#include "voreen/qt/widgets/customlabel.h"

#include <QListWidget>

namespace voreen {

StringListPropertyWidget::StringListPropertyWidget(StringListProperty* prop, QWidget* parent)
    : QPropertyWidget(prop, parent)
    , table_(0)
    , property_(prop)
{
    tgtAssert(property_,"null pointer");

    table_ = new QListWidget();
    addWidget(table_);
    table_->setMaximumHeight(150); table_->setMinimumHeight(150);   // set to nice size
    table_->setEditTriggers(QAbstractItemView::NoEditTriggers);     // removes edit
    table_->setSelectionBehavior(QAbstractItemView::SelectRows);    //select entire row
    table_->setSelectionMode( QAbstractItemView::ExtendedSelection );  //allow multiple selection

    connect(table_,SIGNAL(itemSelectionChanged()),this,SLOT(selectionOnChange()));

    //update
    updateTable();
}

StringListPropertyWidget::~StringListPropertyWidget() {
}


//----------------------------------------------------------
//  Overrides
//----------------------------------------------------------
void StringListPropertyWidget::updateFromPropertySlot() {
    updateTable();
}

//----------------------------------------------------------
//  Helpers
//----------------------------------------------------------
void StringListPropertyWidget::updateTable() {
    tgtAssert(property_,"no property");
    table_->blockSignals(true);

    if(property_->neededTableUpdates_ & StringListProperty::UPDATE_ROWS) {
        table_->clear(); //clears rows
        for (const auto& value : property_->values_) {
            QListWidgetItem* item = new QListWidgetItem(value.first.c_str());
            tgt::ivec4 color(value.second.r * 255, value.second.g * 255, value.second.b * 255, 255);
            item->setForeground(QBrush(QColor(color.r, color.g, color.b, color.a)));
            table_->addItem(item);
        }
    }

    if(property_->neededTableUpdates_ & StringListProperty::UPDATE_SELECTION) {
        table_->clearSelection();
        for(int row : property_->getSelectedRowIndices())
            table_->item(row)->setSelected(true);
    }

    table_->blockSignals(false);
    //property_->neededTableUpdates_ = StringListProperty::UPDATE_NONE;// done in updateWidgets
}

void StringListPropertyWidget::selectionOnChange() {
    if(!disconnected_) {
        std::vector<int> indices;
        for(int i = 0; i < table_->count(); i++) {
            if(table_->item(i)->isSelected())
                indices.push_back(i);
        }
        property_->setSelectedRowIndices(indices);
    }
}

} // namespace
