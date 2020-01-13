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

#ifndef VRN_PLOTDATASELECTWIDGET_H
#define VRN_PLOTDATASELECTWIDGET_H

#include "../../processors/plotdataselect.h"
#include "voreen/qt/widgets/processor/qprocessorwidget.h"
//#include "voreen/core/plotting/functionlibrary.h"

class QMenu;

namespace voreen {


class ExtendedTable;
class PlotDataExtendedTableModel;


class PlotDataSelectWidget : public QProcessorWidget {
    Q_OBJECT
public:
    PlotDataSelectWidget(QWidget* parent, PlotDataSelect* plotDataSelect);
    virtual ~PlotDataSelectWidget();

    void initialize();

    virtual void updateFromProcessor();

public slots:
//    void checkboxchange(int state);

protected:

    //void createSnapshotTool();

private:
    PlotDataExtendedTableModel* proxyModel_;
    ExtendedTable* table_;
    QMenu* contextMenuTable_;
    QMenu* beforeMenu_;
    QMenu* afterMenu_;
    QMenu* switchMenu_;
    QMenu* copyColumnMenu_;
    FunctionLibrary::ProcessorFunctionalityType function_;
    std::vector<QWidget*> widgetVector_;
    PlotPredicateProperty* pPVProperty_;

    static const std::string loggerCat_;

private slots:
    void functionSelect();
    void tableContextMenu(const QPoint& pos);
    void selectChange();
    void selectColumnChange();
    void selectColumnRowChange();
    void selectRowChange();
    void selectAreaChange();
    void selectResetAll();
    void selectResetLast();
    void sortChangeBefore();
    void sortChangeAfter();
    void sortChangeSwitch();
    void predicateChange();
    void filterEmptyCells();
    void filterAlphaNumeric();
    void copyColumn();
    void renameColumn();
};

} // namespace voreen

#endif // VRN_PLOTDATASELECTWIDGET_H

