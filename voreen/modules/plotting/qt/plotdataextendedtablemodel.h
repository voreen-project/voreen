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

#ifndef VRN_PLOTDATAEXTENDEDTABLEMODEL_H
#define VRN_PLOTDATAEXTENDEDTABLEMODEL_H

#include "../datastructures/plotdata.h"
#include "../utils/functionlibrary.h"

#include "plotdatasimpletablemodel.h"

#include <QModelIndex>
#include <QObject>

namespace voreen {

class PlotDataExtendedTableModel : public PlotDataSimpleTableModel {
    Q_OBJECT
public:
    PlotDataExtendedTableModel(PlotData* data,
        const std::vector<std::pair<int, std::string> >& StringVector,  QObject* parent=0,
         const QColor& keyColumnColor = QColor(200,200,200,255),
         const QColor& dataColumnColor = QColor(240,240,240,255),
         const QColor& selectColumnColor = QColor(255,180,142,255),
         int selectedColumn = -1);
    /// set what is able to show in the widget
    void setFunctionType(FunctionLibrary::ProcessorFunctionalityType type, bool showColumnType = false);
    /// how many additional lines will be shown above the datalines
    int getAddLines() const;
    /// @see PlotDataSimpleTableModel
    int rowCount(const QModelIndex &parent) const;
    /// @see PlotDataSimpleTableModel
    QVariant data(const QModelIndex &index, int role) const;
    /// @see PlotDataSimpleTableModel
    Qt::ItemFlags flags(const QModelIndex &index) const;
    /// @see PlotDataSimpleTableModel
    QVariant getCellAt(int row, int column) const;
    /// @see PlotDataSimpleTableModel
    PlotCellValue getPlotCellAt(int row, int column) const;

private:

    PlotData* data_;
    std::vector<std::pair<int, std::string > > stringVector_;
    int addLines_;
    FunctionLibrary::ProcessorFunctionalityType type_;
    bool showColumnType_;
};

}
#endif // VRN_PLOTSELSORTFILTERTABLEMODEL_H
