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

#ifndef VRN_PLOTSELECTIONTABLEMODEL_H
#define VRN_PLOTSELECTIONTABLEMODEL_H

#include <QAbstractTableModel>
#include <QObject>
#include <QVBoxLayout>
#include <QHeaderView>

namespace voreen {

// forward declaration
class PlotSelectionProperty;

/**
 * Table model for PlotSelectionProperty
 **/
class PlotSelectionTableModel : public QAbstractTableModel {
    Q_OBJECT

public:
    /**
     * Constructor
     *
     * \param   prop    PlotSelectionProperty whose data shall be visualized
     **/
    PlotSelectionTableModel(PlotSelectionProperty* prop, QObject *parent=0);

    /// Returns number of rows.
    int rowCount(const QModelIndex &parent) const;
    /// Returns number of columns.
    int columnCount(const QModelIndex &parent) const;

    /**
     * Returns the data at index \a index for role \a role
     **/
    QVariant data(const QModelIndex &index, int role) const;

    /**
     * Sets the data at index \a index to \a value.
     *
     * \param   index   position where to set data
     * \param   value   new data
     * \param   role    edit role
     **/
    bool setData(const QModelIndex& index, const QVariant& value, int role);

    /**
     * Returns the header data.
     **/
    QVariant headerData(int section, Qt::Orientation orientation, int role) const;

    /**
     * Returns the item flags at index \a index
     **/
    Qt::ItemFlags flags(const QModelIndex& index) const;

    /**
     * Force reinitialization of all data in table model
     **/
    void invalidate();

private:
    PlotSelectionProperty* prop_;
};

}

#endif // VRN_PLOTSELECTIONTABLEMODEL_H
