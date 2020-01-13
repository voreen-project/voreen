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

#include "plotdatasimpletablemodel.h"

#include "../datastructures/plotdata.h"

#include "tgt/logmanager.h"

namespace voreen {

PlotDataSimpleTableModel::PlotDataSimpleTableModel(const voreen::PlotData* data, QObject* parent,
        const QColor& keyColumnColor, const QColor& dataColumnColor, const QColor& highlightedColor,
        const QColor& selectColumnColor, int selectedColumn)
    : QAbstractTableModel(parent)
    , pData_(data)
    , keyColumnColor_(keyColumnColor)
    , dataColumnColor_(dataColumnColor)
    , highlightedColor_(highlightedColor)
    , selectColumnColor_(selectColumnColor)
    , selectedColumn_(selectedColumn)
{
}

int PlotDataSimpleTableModel::rowCount(const QModelIndex& /*parent*/) const {
    if (! pData_)
        return 0;

    return pData_->getRowsCount();
}

int PlotDataSimpleTableModel::columnCount(const QModelIndex& /*parent*/) const {
    if (! pData_)
        return 0;

    return pData_->getColumnCount();
}

QVariant PlotDataSimpleTableModel::data(const QModelIndex &index, int role) const {
    if (! pData_)
        return QVariant();


    if (!index.isValid())
        return QVariant();

    if (index.row() >= rowCount(index) || index.row() < 0)
        return QVariant();

    if (role == Qt::TextAlignmentRole) {
        if (pData_->getColumnType(index.column()) == PlotBase::NUMBER)
            return QVariant(Qt::AlignRight | Qt::AlignVCenter);
        if (pData_->getColumnType(index.column()) == PlotBase::STRING)
            return QVariant(Qt::AlignLeft | Qt::AlignVCenter);
        return QVariant(Qt::AlignCenter);
    }

    if (role == Qt::DisplayRole)
        return getCellAt(index.row(), index.column());

    if (role == Qt::BackgroundColorRole && getPlotCellAt(index.row(),index.column()).isHighlighted()) {
        return highlightedColor_;
    }
    else if (role == Qt::BackgroundColorRole && selectedColumn_ == index.column()) {
        return selectColumnColor_;
    }
    else if (role == Qt::BackgroundColorRole && pData_->getKeyColumnCount() > index.column()) {
        return keyColumnColor_;
    }
    else if (role == Qt::BackgroundColorRole && pData_->getKeyColumnCount() <= index.column()) {
        return dataColumnColor_;
    }

    return QVariant();
}

QVariant PlotDataSimpleTableModel::getCellAt(int row, int column) const {
    if (!pData_ || pData_->getRowsCount() == 0)
        return QVariant();
    std::vector<PlotRowValue>::const_iterator it = pData_->getRowsBegin();
    it += row;

    PlotCellValue cell = it->getCellAt(column);
    return (cell.isValue() ? QVariant(cell.getValue()) : QVariant(QString::fromStdString(cell.getTag())));
}

PlotCellValue PlotDataSimpleTableModel::getPlotCellAt(int row, int column) const {
    if (!pData_ || pData_->getRowsCount() == 0)
        return PlotCellValue();
    std::vector<PlotRowValue>::const_iterator it = pData_->getRowsBegin();
    it += row;

    return it->getCellAt(column);
}


QVariant PlotDataSimpleTableModel::headerData(int section, Qt::Orientation /*orientation*/, int role) const {
    if (! pData_)
        return QVariant();

    if (role != Qt::DisplayRole)
        return QVariant();

    if (section < 0 || section >= pData_->getColumnCount())
        return QVariant();

    return QVariant(QString::fromStdString(pData_->getColumnLabel(section)));
}

Qt::ItemFlags PlotDataSimpleTableModel::flags(const QModelIndex &index) const {
    if (!index.isValid())
        return Qt::ItemIsEnabled;

    return QAbstractTableModel::flags(index) | Qt::ItemIsSelectable | Qt::ItemIsEnabled ;
}

void PlotDataSimpleTableModel::setKeyColumnColor(const QColor &keyColumnColor) {
    keyColumnColor_ = keyColumnColor;
}

void PlotDataSimpleTableModel::setDataColumnColor(const QColor &dataColumnColor) {
    dataColumnColor_ = dataColumnColor;
}

}
