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

#include "plotdatamergewidget.h"

#include "../plotdataextendedtablemodel.h"
#include "../extendedtable.h"
#include "../../utils/functionlibrary.h"

#include "voreen/qt/voreenapplicationqt.h"
#include "voreen/core/properties/boolproperty.h"

#include <QGridLayout>
#include <QMainWindow>
#include <QLabel>
#include <QMenu>


namespace voreen {

PlotDataMergeWidget::PlotDataMergeWidget(QWidget* parent, PlotDataMerge* plotDataMerge)
    : QProcessorWidget(plotDataMerge, parent)
    , tableFirst_()
    , tableSecond_()
    , tableOut_()
{
    tgtAssert(plotDataMerge, "No PlotDataMerge Processor");
    proxyModelFirst_ = new PlotDataExtendedTableModel(NULL, std::vector<std::pair<int,std::string> >(0), this);
    proxyModelSecond_ = new PlotDataExtendedTableModel(NULL, std::vector<std::pair<int,std::string> >(0), this);
    proxyModelOut_ = new PlotDataExtendedTableModel(NULL, std::vector<std::pair<int,std::string> >(0), this);
    setWindowTitle(QString::fromStdString(plotDataMerge->getGuiName()));
    contextMenu_ = new QMenu(this);
    resize(512, 512);
}

PlotDataMergeWidget::~PlotDataMergeWidget() {
    delete proxyModelFirst_;
    delete proxyModelSecond_;
    delete proxyModelOut_;
}

void PlotDataMergeWidget::initialize() {

    QProcessorWidget::initialize();
    QGridLayout* mainLayout = new QGridLayout(this);
    QLabel* label = new QLabel(this);
    QFont font = label->font();
    font.setBold(true);
    font.setPointSize(10);

    label->setText(QString(tr("First PlotData",0)));
    label->setFont(font);
    mainLayout->addWidget(label,0,0);
    label = new QLabel(this);
    label->setText(QString(tr("Second PlotData",0)));
    label->setFont(font);
    mainLayout->addWidget(label,0,2);
    label = new QLabel(this);
    label->setText(QString(tr("Merged PlotData",0)));
    label->setFont(font);
    mainLayout->addWidget(label,2,0);

    delete proxyModelFirst_;
    delete proxyModelSecond_;
    delete proxyModelOut_;
    proxyModelFirst_ = new PlotDataExtendedTableModel(NULL, std::vector<std::pair<int,std::string> >(0), this);
    proxyModelSecond_ = new PlotDataExtendedTableModel(NULL, std::vector<std::pair<int,std::string> >(0), this);
    proxyModelOut_ = new PlotDataExtendedTableModel(NULL, std::vector<std::pair<int,std::string> >(0), this);

    tableFirst_ = new ExtendedTable(this);
    tableFirst_->setModel(proxyModelFirst_);
    tableFirst_->setSortingEnabled(false);
    tableFirst_->verticalHeader()->hide();
    tableFirst_->resizeColumnsToContents();
    tableFirst_->setContextMenuPolicy(Qt::CustomContextMenu);
    tableFirst_->setMinimumHeight(150);
    tableFirst_->setMinimumWidth(150);

    tableSecond_ = new ExtendedTable(this);
    tableSecond_->setModel(proxyModelSecond_);
    tableSecond_->setSortingEnabled(false);
    tableSecond_->verticalHeader()->hide();
    tableSecond_->resizeColumnsToContents();
    tableSecond_->setContextMenuPolicy(Qt::CustomContextMenu);
    tableSecond_->setMinimumHeight(150);
    tableSecond_->setMinimumWidth(150);

    tableOut_ = new ExtendedTable(this);
    tableOut_->setModel(proxyModelOut_);
    tableOut_->setSortingEnabled(false);
    tableOut_->verticalHeader()->hide();
    tableOut_->resizeColumnsToContents();
    tableOut_->setContextMenuPolicy(Qt::CustomContextMenu);
    tableOut_->setMinimumHeight(150);
    tableOut_->setMinimumWidth(150);

    mainLayout->addWidget(tableFirst_, 1, 0, 1, 1);
    mainLayout->addWidget(tableSecond_, 1, 2, 1, 1);
    mainLayout->addWidget(tableOut_, 3, 0, 1, 3);

    setContextMenuPolicy(Qt::CustomContextMenu);
    QObject::connect(this,SIGNAL(customContextMenuRequested(const QPoint&)),this,SLOT(tableContextMenu(const QPoint&)));
    QObject::connect(tableFirst_,SIGNAL(customContextMenuRequested(const QPoint&)),this,SLOT(tableContextMenu(const QPoint&)));
    QObject::connect(tableSecond_,SIGNAL(customContextMenuRequested(const QPoint&)),this,SLOT(tableContextMenu(const QPoint&)));
    QObject::connect(tableOut_,SIGNAL(customContextMenuRequested(const QPoint&)),this,SLOT(tableContextMenu(const QPoint&)));

    setLayout(mainLayout);

    //show();

    initialized_ = true;
}

void PlotDataMergeWidget::updateFromProcessor(){
    PlotDataMerge* plotDataMerge = dynamic_cast<PlotDataMerge*>(processor_);
    delete proxyModelFirst_;
    delete proxyModelSecond_;
    delete proxyModelOut_;
    PlotData* dataFirst = const_cast<PlotData*>(plotDataMerge->getPlotDataFirst());
    PlotData* dataSecond = const_cast<PlotData*>(plotDataMerge->getPlotDataSecond());
    PlotData* dataOut = const_cast<PlotData*>(plotDataMerge->getPlotDataOut());
    const std::vector<Property*>& properties = plotDataMerge->getProperties();

    tgt::vec4 colorValue;
    std::vector<QColor> qColorValues;
    int maxColumnSize = 0;
    for (size_t i = 0; i < properties.size(); ++i) {
        if (dynamic_cast<FloatVec4Property*>(properties[i])) {
            colorValue = static_cast<FloatVec4Property*>(properties[i])->get();
            qColorValues.push_back(QColor(int(colorValue.r*255),
                int(colorValue.g*255),int(colorValue.b*255),
                int(colorValue.a*255)));
        }
        else if (dynamic_cast<IntProperty*>(properties[i])) {
            maxColumnSize = static_cast<IntProperty*>(properties[i])->get();
        }
    }

    std::vector<std::pair<int,std::string> > stringvector;

    proxyModelFirst_ = new PlotDataExtendedTableModel(dataFirst, stringvector, this,qColorValues[0],qColorValues[1]);
    proxyModelSecond_ = new PlotDataExtendedTableModel(dataSecond, stringvector, this,qColorValues[0],qColorValues[1]);
    proxyModelOut_ = new PlotDataExtendedTableModel(dataOut, stringvector, this,qColorValues[0],qColorValues[1]);

    tableFirst_->setModel(proxyModelFirst_);
    tableFirst_->setSortingEnabled(false);
    tableFirst_->verticalHeader()->hide();
    tableFirst_->resizeColumnsToContents();

    tableSecond_->setModel(proxyModelSecond_);
    tableSecond_->setSortingEnabled(false);
    tableSecond_->verticalHeader()->hide();
    tableSecond_->resizeColumnsToContents();

    tableOut_->setModel(proxyModelOut_);
    tableOut_->setSortingEnabled(false);
    tableOut_->verticalHeader()->hide();
    tableOut_->resizeColumnsToContents();
    if (maxColumnSize > 0) {
        for (int i = 0; i < dataFirst->getColumnCount(); ++i) {
            if (tableFirst_->columnWidth(i) > maxColumnSize)
                tableFirst_->setColumnWidth(i,maxColumnSize);
        }
        for (int i = 0; i < dataSecond->getColumnCount(); ++i) {
            if (tableSecond_->columnWidth(i) > maxColumnSize)
                tableSecond_->setColumnWidth(i,maxColumnSize);
        }
        for (int i = 0; i < dataOut->getColumnCount(); ++i) {
            if (tableOut_->columnWidth(i) > maxColumnSize)
                tableOut_->setColumnWidth(i,maxColumnSize);
        }
    }
}

void PlotDataMergeWidget::tableContextMenu(const QPoint& pos) {
    PlotDataMerge* plotDataMerge = dynamic_cast<PlotDataMerge*>(processor_);
    const std::vector<Property*>& properties = plotDataMerge->getProperties();
    const PlotData* dataFirst = plotDataMerge->getPlotDataFirst();
    const PlotData* dataSecond = plotDataMerge->getPlotDataSecond();
    //PlotData* dataOut = const_cast<PlotData*>(plotDataMerge->getPlotDataOut());

    QAction* newAct;
    std::string s;
    contextMenu_->clear();

    s = "Swap PlotData";
    newAct = new QAction(QString::fromStdString(s),this);
    newAct->setCheckable(true);
    newAct->setChecked(dynamic_cast<BoolProperty*>(properties[0])->get());
    contextMenu_->addAction(newAct);
    QObject::connect(newAct,SIGNAL(triggered()),this,SLOT(swapPlotData()));

    if (dataFirst && dataSecond && dataFirst->getKeyColumnCount() == dataSecond->getKeyColumnCount()) {
        s = "Ignore Column Labels";
        newAct = new QAction(QString::fromStdString(s),this);
        newAct->setCheckable(true);
        newAct->setChecked(dynamic_cast<BoolProperty*>(properties[1])->get());
        contextMenu_->addAction(newAct);
        QObject::connect(newAct,SIGNAL(triggered()),this,SLOT(ignoreColumnLabels()));
    }

    contextMenu_->popup(static_cast<QWidget*>(sender())->mapToGlobal(pos));
}

void PlotDataMergeWidget::swapPlotData() {
    PlotDataMerge* plotDataMerge = dynamic_cast<PlotDataMerge*>(processor_);
    plotDataMerge->swapPlotData();
}

void PlotDataMergeWidget::ignoreColumnLabels() {
    PlotDataMerge* plotDataMerge = dynamic_cast<PlotDataMerge*>(processor_);
    plotDataMerge->ignoreColumnLabels();
}

} //namespace voreen

