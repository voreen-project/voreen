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

#include "plotselectiontablemodel.h"
#include "../properties/plotselectionproperty.h"

#include "tgt/logmanager.h"

namespace voreen {

PlotSelectionTableModel::PlotSelectionTableModel(PlotSelectionProperty* prop, QObject* parent)
    : QAbstractTableModel(parent)
    , prop_(prop)
{
    tgtAssert(prop_, "no PlotData");
}

int PlotSelectionTableModel::rowCount(const QModelIndex& /*parent*/) const {
    return prop_->size();
}

int PlotSelectionTableModel::columnCount(const QModelIndex& /*parent*/) const {
    return 4;
}

QVariant PlotSelectionTableModel::data(const QModelIndex &index, int role) const {
    if (!index.isValid())
        return QVariant();

    if (index.row() >= rowCount(index) || index.row() < 0)
        return QVariant();

    if (role == Qt::DisplayRole) {
        if (index.column() == 0)
            return QVariant(QString::fromStdString(prop_->getSelectionAt(prop_->size()-1-index.row()).selection_.toString(*(prop_->getPlotData()))));
    }

    if (role == Qt::CheckStateRole) {
        switch (index.column()) {
            case 1:
                if (!prop_->getEnableHighlightFlags())
                    return QVariant();
                else
                    return (prop_->getSelectionAt(prop_->size()-1-index.row()).highlight_) ? QVariant(Qt::Checked) : QVariant(Qt::Unchecked);
                break;
            case 2:
                if (!prop_->getEnableLabelFlags() || !prop_->getSelectionAt(prop_->size()-1-index.row()).selection_.isTablePosition())
                    return QVariant();
                else
                    return (prop_->getSelectionAt(prop_->size()-1-index.row()).renderLabel_) ? QVariant(Qt::Checked) : QVariant(Qt::Unchecked);
                break;
            case 3:
                if (!prop_->getEnableZoomToFlags())
                    return QVariant();
                else
                    return (prop_->getSelectionAt(prop_->size()-1-index.row()).zoomTo_) ? QVariant(Qt::Checked) : QVariant(Qt::Unchecked);
                break;
            default:
                return QVariant();
        }
    }
    return QVariant();
}

bool PlotSelectionTableModel::setData(const QModelIndex& index, const QVariant& value, int role) {
    if (role == Qt::CheckStateRole) {
        switch (index.column()) {
            case 1:
                prop_->setHighlight(prop_->size()-1-index.row(), (value == Qt::Checked ? true : false));
                return true;
                break;
            case 2:
                prop_->setRenderLabel(prop_->size()-1-index.row(), (value == Qt::Checked ? true : false));
                return true;
                break;
            case 3:
                prop_->setZoomTo(prop_->size()-1-index.row(), (value == Qt::Checked ? true : false));
                return true;
                break;
            default:
                return false;
        }
    }
    return false;
}

QVariant PlotSelectionTableModel::headerData(int section, Qt::Orientation /*orientation*/, int role) const {
    if (role == Qt::TextAlignmentRole && section > 0 && section <= 3) {
        return QVariant(Qt::AlignHCenter);
    }

    if (role != Qt::DisplayRole && role != Qt::ToolTipRole)
        return QVariant();

    switch (section) {
        case 0:
            return QVariant(tr("Selection"));
            break;
        case 1:
            if (role == Qt::DisplayRole) {
                if (!prop_->getEnableHighlightFlags())
                    return QVariant();
                else
                    return QVariant(tr("H"));
            }
            else {
                return QVariant(tr("set Highlight"));
            }
            break;
        case 2:
            if (role == Qt::DisplayRole) {
                if (!prop_->getEnableLabelFlags())
                    return QVariant();
                else
                    return QVariant(tr("L"));
            }
            else {
                return QVariant(tr("set Label"));
            }
            break;
        case 3:
            if (role == Qt::DisplayRole) {
                if (!prop_->getEnableZoomToFlags())
                    return QVariant();
                else
                    return QVariant(tr("Z"));
            }
            else {
                return QVariant(tr("Zoom"));
            }
            break;
        default:
            return QVariant();
    }
}

Qt::ItemFlags PlotSelectionTableModel::flags(const QModelIndex &index) const {
    if (!index.isValid())
        return Qt::ItemIsEnabled;

    if (index.column() == 0)
        return QAbstractTableModel::flags(index) | Qt::ItemIsSelectable | Qt::ItemIsEnabled ;
    else if (index.column() == 1 && prop_->getEnableHighlightFlags())
        return QAbstractTableModel::flags(index) | Qt::ItemIsSelectable | Qt::ItemIsEnabled | Qt::ItemIsUserCheckable;
    else if (index.column() == 2 && prop_->getEnableLabelFlags())
        return QAbstractTableModel::flags(index) | Qt::ItemIsSelectable | Qt::ItemIsEnabled | Qt::ItemIsUserCheckable;
    else if (index.column() == 3 && prop_->getEnableZoomToFlags())
        return QAbstractTableModel::flags(index) | Qt::ItemIsSelectable | Qt::ItemIsEnabled | Qt::ItemIsUserCheckable;

    return Qt::ItemIsEnabled;
}

void PlotSelectionTableModel::invalidate() {
    layoutChanged();
}

}
