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

#include "plotdataextendedtablemodel.h"
#include "tgt/logmanager.h"

namespace voreen {

    PlotDataExtendedTableModel::PlotDataExtendedTableModel(PlotData* data,
        const std::vector<std::pair<int, std::string> >& stringVector, QObject* parent,
        const QColor& keyColumnColor, const QColor& dataColumnColor,
        const QColor& selectColumnColor, int selectedColumn)
    : PlotDataSimpleTableModel(data,parent,keyColumnColor,dataColumnColor,QColor(200,240,240,255),selectColumnColor,selectedColumn)
    , data_(data)
    , stringVector_(stringVector)
    , addLines_(0)
    , type_(FunctionLibrary::NONE)
{
    showColumnType_ = false;
}

void PlotDataExtendedTableModel::setFunctionType(FunctionLibrary::ProcessorFunctionalityType type, bool showColumnType) {
    switch (type) {
        case FunctionLibrary::NONE : addLines_ = 0;
                              break;
        case FunctionLibrary::SELECT : addLines_ = 1;
                                break;
        case FunctionLibrary::COLUMNORDER : addLines_ = 0;
                                     break;
        case FunctionLibrary::FILTER : addLines_ = 0;
                                break;
        case FunctionLibrary::GROUPBY : addLines_ = 1;
                                 break;
        default : addLines_ = 0;
    }
    if (showColumnType) {
        addLines_ += 1;
    }
    showColumnType_ = showColumnType;
    type_ = type;
}

int PlotDataExtendedTableModel::getAddLines() const {
    return addLines_;
}

int PlotDataExtendedTableModel::rowCount(const QModelIndex& /*parent*/) const {
    if (! data_)
        return 0;

    return data_->getRowsCount()+addLines_;
}

QVariant PlotDataExtendedTableModel::getCellAt(int row, int column) const {
    return PlotDataSimpleTableModel::getCellAt(row-addLines_,column);
}

PlotCellValue PlotDataExtendedTableModel::getPlotCellAt(int row, int column) const {
    return PlotDataSimpleTableModel::getPlotCellAt(row-addLines_,column);
}

QVariant PlotDataExtendedTableModel::data(const QModelIndex &index, int role) const {
    if (! data_)
        return QVariant();

    if (!index.isValid())
        return QVariant();

    if (index.row() >= rowCount(index) || index.row() < 0)
        return QVariant();

    if (index.row() == addLines_-1 && role == Qt::DisplayRole && addLines_ > 0 && showColumnType_) {
        int columntype = static_cast<int>(data_->getColumnType(index.column()));
        std::string ausgabe;
        switch (columntype) {
            case 1 : ausgabe = "STRING"; break;
            case 2 : ausgabe = "NUMBER"; break;
            default: ausgabe = "EMPTY"; break;
        }
        return QVariant(QString::fromStdString(ausgabe));
    }
    if (index.row() == 0 && role == Qt::DisplayRole && (type_ == FunctionLibrary::SELECT || type_ == FunctionLibrary::GROUPBY)) {
        for (size_t i = 0; i < stringVector_.size(); ++i) {
            if (stringVector_.at(i).first == index.column()) {
                return QVariant(QString::fromStdString(stringVector_.at(i).second));
                break;
            }
        }
        return QVariant("");
    }

    if (role == Qt::TextAlignmentRole && index.row() < addLines_) {
        return QVariant(Qt::AlignRight | Qt::AlignVCenter);
    }

    if (index.row() < addLines_ && (role == Qt::DisplayRole || role == Qt::BackgroundColorRole))
        return QVariant();

    if ((index.row() >=  addLines_) && (role == Qt::DisplayRole))
        return getCellAt(index.row(), index.column());

    return PlotDataSimpleTableModel::data(index,role);
}

Qt::ItemFlags PlotDataExtendedTableModel::flags(const QModelIndex &index) const {
    if (!index.isValid())
        return Qt::ItemIsEnabled;

    return PlotDataSimpleTableModel::flags(index);
}

}
