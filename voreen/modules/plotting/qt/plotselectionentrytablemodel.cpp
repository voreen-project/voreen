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

#include "plotselectionentrytablemodel.h"

#include "plotpredicatedialog.h"
#include "../properties/plotselectionproperty.h"
#include "../datastructures/plotselection.h"
#include "../datastructures/plotdata.h"

#include "tgt/logmanager.h"

#include <string>
#include <QComboBox>

namespace voreen {

PlotSelectionEntryTableModel::PlotSelectionEntryTableModel(PlotSelectionEntry& entry, const PlotData* plotData, QObject* parent)
    : QAbstractTableModel(parent)
    , entry_(entry)
    , plotData_(plotData)
{
    tgtAssert(plotData_, "PlotSelectionEntryTableModel(): plotData is NULL");
}

int PlotSelectionEntryTableModel::rowCount(const QModelIndex& /*parent*/) const {
    if (entry_.selection_.isTablePosition())
        return 0;
    else
        return static_cast<int>(entry_.selection_.getSelection().size());
}

int PlotSelectionEntryTableModel::columnCount(const QModelIndex& /*parent*/) const {
    return 2;
}

QVariant PlotSelectionEntryTableModel::data(const QModelIndex &index, int role) const {
    if (!index.isValid())
        return QVariant();

    if (index.row() >= rowCount(index) || index.row() < 0)
        return QVariant();

    if (role == Qt::DisplayRole) {
        const std::pair<int, const PlotPredicate*>& pair = entry_.selection_.getSelection()[index.row()];
        if (index.column() == 0)
            switch (pair.first) {
                case PlotSelectionProperty::X_AXIS_COLUMN:
                    return QVariant(tr("X Axis"));
                    break;
                case PlotSelectionProperty::Y_AXIS_COLUMN:
                    return QVariant(tr("Y Axis"));
                    break;
                case PlotSelectionProperty::Z_AXIS_COLUMN:
                    return QVariant(tr("Z Axis"));
                    break;
                default:
                    if (pair.first < plotData_->getColumnCount())
                        return QVariant(QString::fromStdString(plotData_->getColumnLabel(pair.first)));
                    else
                        return QVariant(pair.first);
                    break;
            }
        else if (index.column() == 1) {
            if (pair.second != 0)
                return QVariant(QString::fromStdString(pair.second->toString()));
            else
                return QVariant(tr("No Predicate"));
        }
    }

    return QVariant();
}

bool PlotSelectionEntryTableModel::setData(const QModelIndex& index, const QVariant& value, int role) {
    if (role == Qt::EditRole && index.column() == 0) {
        entry_.selection_.setColumn(index.row(), value.toInt());
    }
    return false;
}

QVariant PlotSelectionEntryTableModel::headerData(int section, Qt::Orientation /*orientation*/, int role) const {
    if (role == Qt::TextAlignmentRole && section > 0 && section <= 3) {
        return QVariant(Qt::AlignHCenter);
    }

    if (role != Qt::DisplayRole && role != Qt::ToolTipRole)
        return QVariant();

    switch (section) {
        case 0:
            return QVariant(tr("Column"));
            break;
        case 1:
            return QVariant(tr("Predicate"));
            break;
        default:
            return QVariant();
    }
}

Qt::ItemFlags PlotSelectionEntryTableModel::flags(const QModelIndex &index) const {
    if (!index.isValid())
        return Qt::ItemIsEnabled;

    return QAbstractTableModel::flags(index) | Qt::ItemIsSelectable | Qt::ItemIsEnabled | Qt::ItemIsEditable;
}

void PlotSelectionEntryTableModel::invalidate() {
    layoutChanged();
}

// ----------------------------------------------------------------------------

PlotSelectionEntryItemDelegate::PlotSelectionEntryItemDelegate(PlotSelectionEntry& entry,
                                                               const PlotData* plotData,
                                                               const PlotEntitiesProperty& entitiesProp,
                                                               bool threeDimensional,
                                                               QObject* parent)
    : QItemDelegate(parent)
    , entry_(entry)
    , plotData_(plotData)
    , entitiesProp_(entitiesProp)
    , threeDimensional_(threeDimensional)
{
    tgtAssert(plotData_, "PlotSelectionEntryItemDelegate(): plotData is NULL");
}

QWidget *PlotSelectionEntryItemDelegate::createEditor(QWidget* parent, const QStyleOptionViewItem& /*option*/, const QModelIndex& index) const {
    if (index.column() == 0) {
        QComboBox* editor = new QComboBox(parent);
        editor->addItem(tr("X Axis"), QVariant(PlotSelectionProperty::X_AXIS_COLUMN));
        editor->addItem(tr("Y Axis"), QVariant(PlotSelectionProperty::Y_AXIS_COLUMN));
        if (threeDimensional_)
            editor->addItem(tr("Z Axis"), QVariant(PlotSelectionProperty::Z_AXIS_COLUMN));

        for (int i = 0; i < plotData_->getColumnCount(); ++i) {
            const std::string& label = plotData_->getColumnLabel(i);
            if (!label.empty())
                editor->addItem(QString::fromStdString(label), i);
            else
                editor->addItem(QString(i), i);
        }

        return editor;
    }
    if (index.column() == 1) {
        // we do not offer a real editor widget here
        // all we do is starting a PlotPredicateDialog and on accept setting the new predicate

        const std::pair<int, const PlotPredicate*>& pair = entry_.selection_.getSelection()[index.row()];
        int realColumnIndex = pair.first;
        if (realColumnIndex < 0) {
            if (realColumnIndex == PlotSelectionProperty::X_AXIS_COLUMN)
                realColumnIndex = entitiesProp_.getXColumnIndex();
            else if (realColumnIndex == PlotSelectionProperty::Y_AXIS_COLUMN && threeDimensional_)
                realColumnIndex = entitiesProp_.getYColumnIndex();
            else // must be data axis
                realColumnIndex = (entitiesProp_.getDataColumnIndices().empty() ? 0 : entitiesProp_.getDataColumnIndices().front());
        }
        PlotPredicateDialog dialog(pair.second, (plotData_->getColumnType(realColumnIndex) == PlotData::STRING));
        if (dialog.exec() == QDialog::Accepted) {
            const PlotPredicate* p = 0;
            if (dialog.getPlotPredicate() != 0)
                p = dialog.getPlotPredicate();
            entry_.selection_.setPredicate(index.row(), p);
        }

        // as mentioned above no widget here - we're already finished
        return NULL;
    }
    return NULL;
}

void PlotSelectionEntryItemDelegate::setEditorData(QWidget *editor, const QModelIndex &index) const {
    const std::pair<int, const PlotPredicate*>& pair = entry_.selection_.getSelection()[index.row()];
    if (index.column() == 0) {
        QComboBox *cb = static_cast<QComboBox*>(editor);
        switch (pair.first) {
            case PlotSelectionProperty::X_AXIS_COLUMN:
                cb->setCurrentIndex(0);
                break;
            case PlotSelectionProperty::Y_AXIS_COLUMN:
                cb->setCurrentIndex(1);
                break;
            case PlotSelectionProperty::Z_AXIS_COLUMN:
                cb->setCurrentIndex(2);
                break;
            default:
                int item = pair.first + 3;
                if (cb->count() > item)
                    cb->setCurrentIndex(item);
                break;
        }
    }
}

void PlotSelectionEntryItemDelegate::setModelData(QWidget *editor, QAbstractItemModel *model, const QModelIndex &index) const {
    if (index.column() == 0) {
        QComboBox *cb = static_cast<QComboBox*>(editor);
        switch (cb->currentIndex()) {
            case 0:
                model->setData(index, PlotSelectionProperty::X_AXIS_COLUMN);
                break;
            case 1:
                model->setData(index, PlotSelectionProperty::Y_AXIS_COLUMN);
                break;
            case 2:
                if (threeDimensional_)
                    model->setData(index, PlotSelectionProperty::Z_AXIS_COLUMN);
                else
                    model->setData(index, cb->currentIndex()-2);
                break;
            default:
                model->setData(index, cb->currentIndex()-(threeDimensional_ ? 3 : 2));
                break;
        }
    }
}

void PlotSelectionEntryItemDelegate::updateEditorGeometry(QWidget* editor, const QStyleOptionViewItem& option, const QModelIndex& /* index */) const {
    editor->setGeometry(option.rect);
}

}
