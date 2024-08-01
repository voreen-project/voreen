/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2024 University of Muenster, Germany,                        *
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

#ifndef VRN_PLOTDATAFITFUNCTIONWIDGET_H
#define VRN_PLOTDATAFITFUNCTIONWIDGET_H

#include "../../processors/plotdatafitfunction.h"
#include "voreen/qt/widgets/processor/qprocessorwidget.h"

class QMenu;

namespace voreen {

class ExtendedTable;
class PlotDataExtendedTableModel;

class PlotDataFitFunctionWidget : public QProcessorWidget {
    Q_OBJECT
public:
    PlotDataFitFunctionWidget(QWidget* parent, PlotDataFitFunction* PlotDataFitFunction);
    virtual ~PlotDataFitFunctionWidget();

    void initialize();

    virtual void updateFromProcessor();

private:
    PlotDataExtendedTableModel* proxyModel_;
    ExtendedTable* table_;
    QMenu* contextMenuTable_;
    QMenu* fittingFunctionMenu_;
    QMenu* multiRegessionMenu_;
    QMenu* splineFunctionMenu_;
    std::vector<QWidget*> widgetVector_;

private slots:

    void tableContextMenu(const QPoint& pos);
    void addFunction();
    void ignoreFalseValues();

};

} // namespace voreen

#endif // VRN_PLOTDATAFITFUNCTIONWIDGET_H

