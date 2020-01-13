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

#include "plotdatagroupwidget.h"

#include "../plotpredicatedialog.h"
#include "../plotdataextendedtablemodel.h"
#include "../extendedtable.h"

#include "voreen/qt/voreenapplicationqt.h"

#include <QGridLayout>
#include <QMainWindow>
#include <QLabel>
#include <QMenu>


namespace voreen {

const std::string PlotDataGroupWidget::loggerCat_("voreen.PlotDataGroupWidget");

PlotDataGroupWidget::PlotDataGroupWidget(QWidget* parent, PlotDataGroup* plotDataGroup)
    : QProcessorWidget(plotDataGroup, parent)
    , table_()
    , widgetVector_(0)
{
    tgtAssert(plotDataGroup, "No PlotDataGroup processor");
    proxyModel_ = new PlotDataExtendedTableModel(NULL, std::vector<std::pair<int,std::string> >(0), this);
    setWindowTitle(QString::fromStdString(plotDataGroup->getGuiName()));
    contextMenuTable_ = new QMenu(this);
    addFunctionMenu_ = new QMenu(tr("Add Function",0),this);
    averageFunctionMenu_ = new QMenu(tr("Average/Mean",0),this);
    histogramFunctionMenu_ = new QMenu(tr("Histogram-Function",0),this);
    resize(256, 256);
}

PlotDataGroupWidget::~PlotDataGroupWidget() {
    delete proxyModel_;
    widgetVector_.clear();
}

void PlotDataGroupWidget::initialize() {

    QProcessorWidget::initialize();
    QGridLayout* mainLayout = new QGridLayout(this);
    QLabel* label = new QLabel(this);
    label->setText(QString(tr("Group by one Column",0)));
    label->setAlignment(Qt::AlignHCenter);
    QFont font = label->font();
    font.setBold(true);
    font.setPointSize(14);
    label->setFont(font);
    mainLayout->addWidget(label,0,0,1,3);
    widgetVector_.push_back(label);
    delete proxyModel_;
    proxyModel_ = new PlotDataExtendedTableModel(NULL, std::vector<std::pair<int,std::string> >(0), this);
    table_ = new ExtendedTable(this);

    table_->setModel(proxyModel_);
    table_->setSortingEnabled(false);
    table_->verticalHeader()->hide();
    table_->resizeColumnsToContents();
    mainLayout->addWidget(table_, 3, 0, 1, 3);
    table_->setContextMenuPolicy(Qt::CustomContextMenu);
    table_->horizontalHeader()->setContextMenuPolicy(Qt::CustomContextMenu);

    QObject::connect(table_,SIGNAL(customContextMenuRequested(const QPoint&)),this,SLOT(tableContextMenu(const QPoint&)));
    QObject::connect(table_->horizontalHeader(),SIGNAL(customContextMenuRequested(const QPoint&)),this,SLOT(tableContextMenu(const QPoint&)));

    setLayout(mainLayout);

    //show();

    initialized_ = true;
}

void PlotDataGroupWidget::updateFromProcessor(){
    PlotDataGroup* plotDataGroup = dynamic_cast<PlotDataGroup*>(processor_);
    delete proxyModel_;
    PlotData* data = const_cast<PlotData*>(plotDataGroup->getPlotData());
    std::vector<std::pair<int,std::string> > stringvector;
    for (size_t i = 0; i < aggregatedcolumns.size(); ++i) {
        if (aggregatedcolumns[i].second)
            stringvector.push_back(std::pair<int,std::string>(aggregatedcolumns.at(i).first,aggregatedcolumns.at(i).second->toString()));
    }
    for (int i = 0; i < data->getColumnCount(); ++i){
        stringvector.push_back(std::pair<int,std::string>(i,"No Function"));
    }

    const std::vector<Property*>& properties = plotDataGroup->getProperties();

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
    proxyModel_ = new PlotDataExtendedTableModel(data,stringvector,this,qColorValues[0],qColorValues[1]);
    proxyModel_ ->setFunctionType(FunctionLibrary::GROUPBY);
    table_->setModel(proxyModel_);
    table_->setSortingEnabled(false);
    table_->verticalHeader()->hide();
    table_->setSelectionMode(QAbstractItemView::SingleSelection);
    table_->setSelectionBehavior(QAbstractItemView::SelectColumns);
    table_->resizeColumnsToContents();
    table_->setMinimumHeight(150);
    table_->setMinimumWidth(150);
    if (maxColumnSize > 0) {
        for (int i = 0; i < data->getColumnCount(); ++i) {
            if (table_->columnWidth(i) > maxColumnSize)
                table_->setColumnWidth(i,maxColumnSize);
        }
    }
}

void PlotDataGroupWidget::tableContextMenu(const QPoint& pos) {
    PlotDataGroup* plotDataGroup = dynamic_cast<PlotDataGroup*>(processor_);
    QModelIndex index = table_->indexAt(pos);
    contextMenuTable_->clear();
    addFunctionMenu_->clear();
    averageFunctionMenu_->clear();
    histogramFunctionMenu_->clear();
    PlotData* data = const_cast<PlotData*>(plotDataGroup->getPlotData());
    QAction* newAct;
    QList<QVariant>* qlist = new QList<QVariant>();
    std::string s;

    s = "Reset All";
    newAct = new QAction(QString::fromStdString(s),this);
    contextMenuTable_->addAction(newAct);
    QObject::connect(newAct,SIGNAL(triggered()),this,SLOT(selectResetAll()));

    s = "Reset Last";
    newAct = new QAction(QString::fromStdString(s),this);
    contextMenuTable_->addAction(newAct);
    QObject::connect(newAct,SIGNAL(triggered()),this,SLOT(selectResetLast()));

    contextMenuTable_->addSeparator();
    bool found = false;

    s = "Delete Function";
    newAct = new QAction(QString::fromStdString(s),this);
    qlist->push_back(0);
    qlist->push_back(index.column());
    newAct->setData(QVariant(*qlist));
    addFunctionMenu_->addAction(newAct);
    QObject::connect(newAct,SIGNAL(triggered()),this,SLOT(addFunction()));
    qlist->clear();
    for(size_t i=0; i< aggregatedcolumns.size(); ++i){
        if(index.column() == aggregatedcolumns.at(i).first){
            found = true;
            break;
        }
    }
    if (!found || aggregatedcolumns.size() == 0 || data->getColumnType(index.column()) == PlotBase::EMPTY) {
        newAct->setEnabled(false);
    }
    else {
        newAct->setEnabled(true);
    }
    qlist->clear();

    addFunctionMenu_->addSeparator();

    s = "Arithmetic";
    newAct = new QAction(QString::fromStdString(s),this);
    qlist->push_back(static_cast<int>(AggregationFunction::AVERAGE));
    qlist->push_back(index.column());
    newAct->setData(QVariant(*qlist));
    averageFunctionMenu_->addAction(newAct);
    QObject::connect(newAct,SIGNAL(triggered()),this,SLOT(addFunction()));
    qlist->clear();

    s = "Geometric";
    newAct = new QAction(QString::fromStdString(s),this);
    qlist->push_back(static_cast<int>(AggregationFunction::GEOMETRICMEAN));
    qlist->push_back(index.column());
    newAct->setData(QVariant(*qlist));
    averageFunctionMenu_->addAction(newAct);
    QObject::connect(newAct,SIGNAL(triggered()),this,SLOT(addFunction()));
    qlist->clear();

    s = "Harmonic";
    newAct = new QAction(QString::fromStdString(s),this);
    qlist->push_back(static_cast<int>(AggregationFunction::HARMONICMEAN));
    qlist->push_back(index.column());
    newAct->setData(QVariant(*qlist));
    averageFunctionMenu_->addAction(newAct);
    QObject::connect(newAct,SIGNAL(triggered()),this,SLOT(addFunction()));
    qlist->clear();

    addFunctionMenu_->addMenu(averageFunctionMenu_);

    s = "Count";
    newAct = new QAction(QString::fromStdString(s),this);
    qlist->push_back(static_cast<int>(AggregationFunction::COUNT));
    qlist->push_back(index.column());
    newAct->setData(QVariant(*qlist));
    addFunctionMenu_->addAction(newAct);
    QObject::connect(newAct,SIGNAL(triggered()),this,SLOT(addFunction()));
    qlist->clear();

    s = "Max";
    newAct = new QAction(QString::fromStdString(s),this);
    qlist->push_back(static_cast<int>(AggregationFunction::MAX));
    qlist->push_back(index.column());
    newAct->setData(QVariant(*qlist));
    addFunctionMenu_->addAction(newAct);
    QObject::connect(newAct,SIGNAL(triggered()),this,SLOT(addFunction()));
    qlist->clear();

    s = "Median";
    newAct = new QAction(QString::fromStdString(s),this);
    qlist->push_back(static_cast<int>(AggregationFunction::MEDIAN));
    qlist->push_back(index.column());
    newAct->setData(QVariant(*qlist));
    addFunctionMenu_->addAction(newAct);
    QObject::connect(newAct,SIGNAL(triggered()),this,SLOT(addFunction()));
    qlist->clear();

    s = "Min";
    newAct = new QAction(QString::fromStdString(s),this);
    qlist->push_back(static_cast<int>(AggregationFunction::MIN));
    qlist->push_back(index.column());
    newAct->setData(QVariant(*qlist));
    addFunctionMenu_->addAction(newAct);
    QObject::connect(newAct,SIGNAL(triggered()),this,SLOT(addFunction()));
    qlist->clear();

    s = "Mode";
    newAct = new QAction(QString::fromStdString(s),this);
    qlist->push_back(static_cast<int>(AggregationFunction::MODE));
    qlist->push_back(index.column());
    newAct->setData(QVariant(*qlist));
    addFunctionMenu_->addAction(newAct);
    QObject::connect(newAct,SIGNAL(triggered()),this,SLOT(addFunction()));
    qlist->clear();

    s = "Standard Deviation";
    newAct = new QAction(QString::fromStdString(s),this);
    qlist->push_back(static_cast<int>(AggregationFunction::STANDARDDEVIATION));
    qlist->push_back(index.column());
    newAct->setData(QVariant(*qlist));
    addFunctionMenu_->addAction(newAct);
    QObject::connect(newAct,SIGNAL(triggered()),this,SLOT(addFunction()));
    qlist->clear();

    s = "Sum";
    newAct = new QAction(QString::fromStdString(s),this);
    qlist->push_back(static_cast<int>(AggregationFunction::SUM));
    qlist->push_back(index.column());
    newAct->setData(QVariant(*qlist));
    addFunctionMenu_->addAction(newAct);
    QObject::connect(newAct,SIGNAL(triggered()),this,SLOT(addFunction()));
    qlist->clear();

    s = "Variance";
    newAct = new QAction(QString::fromStdString(s),this);
    qlist->push_back(static_cast<int>(AggregationFunction::VARIANCE));
    qlist->push_back(index.column());
    newAct->setData(QVariant(*qlist));
    addFunctionMenu_->addAction(newAct);
    QObject::connect(newAct,SIGNAL(triggered()),this,SLOT(addFunction()));
    qlist->clear();

    addFunctionMenu_->addSeparator();

    //Histogramm-Menu
    s = "Count";
    newAct = new QAction(QString::fromStdString(s),this);
    qlist->push_back(static_cast<int>(AggregationFunction::COUNTHISTOGRAM));
    qlist->push_back(index.column());
    newAct->setData(QVariant(*qlist));
    histogramFunctionMenu_->addAction(newAct);
    QObject::connect(newAct,SIGNAL(triggered()),this,SLOT(addFunction()));
    qlist->clear();

    s = "Min";
    newAct = new QAction(QString::fromStdString(s),this);
    qlist->push_back(static_cast<int>(AggregationFunction::MINHISTOGRAM));
    qlist->push_back(index.column());
    newAct->setData(QVariant(*qlist));
    histogramFunctionMenu_->addAction(newAct);
    QObject::connect(newAct,SIGNAL(triggered()),this,SLOT(addFunction()));
    qlist->clear();

    s = "Max";
    newAct = new QAction(QString::fromStdString(s),this);
    qlist->push_back(static_cast<int>(AggregationFunction::MAXHISTOGRAM));
    qlist->push_back(index.column());
    newAct->setData(QVariant(*qlist));
    histogramFunctionMenu_->addAction(newAct);
    QObject::connect(newAct,SIGNAL(triggered()),this,SLOT(addFunction()));
    qlist->clear();

    s = "Median";
    newAct = new QAction(QString::fromStdString(s),this);
    qlist->push_back(static_cast<int>(AggregationFunction::MEDIANHISTOGRAM));
    qlist->push_back(index.column());
    newAct->setData(QVariant(*qlist));
    histogramFunctionMenu_->addAction(newAct);
    QObject::connect(newAct,SIGNAL(triggered()),this,SLOT(addFunction()));
    qlist->clear();

    s = "Sum";
    newAct = new QAction(QString::fromStdString(s),this);
    qlist->push_back(static_cast<int>(AggregationFunction::SUMHISTOGRAM));
    qlist->push_back(index.column());
    newAct->setData(QVariant(*qlist));
    histogramFunctionMenu_->addAction(newAct);
    QObject::connect(newAct,SIGNAL(triggered()),this,SLOT(addFunction()));
    qlist->clear();

    addFunctionMenu_->addMenu(histogramFunctionMenu_);

    if (data->getColumnType(index.column()) != PlotBase::NUMBER) {
        addFunctionMenu_->setEnabled(false);
    }
    else {
        addFunctionMenu_->setEnabled(true);
    }

    contextMenuTable_->addMenu(addFunctionMenu_);

    s = "Group By";
    newAct = new QAction(QString::fromStdString(s),this);
    newAct->setData(QVariant(index.column()));
    contextMenuTable_->addAction(newAct);
    QObject::connect(newAct,SIGNAL(triggered()),this,SLOT(groupby()));
    for(size_t i=0; i< aggregatedcolumns.size(); ++i){
        if(index.column() == aggregatedcolumns.at(i).first){
            found = true;
            break;
        }
    }
    if (found || aggregatedcolumns.size() == 0 || data->getColumnType(index.column()) == PlotBase::EMPTY) {
        newAct->setEnabled(false);
    }
    else {
        newAct->setEnabled(true);
    }
    qlist->clear();


    contextMenuTable_->popup(table_->mapToGlobal(pos));
    delete qlist;

}

void PlotDataGroupWidget::addFunction() {
    QVariant x = dynamic_cast<QAction*>(sender())->data();
    AggregationFunction::AggregationFunctionEnum function = static_cast<AggregationFunction::AggregationFunctionEnum>(x.toList().at(0).toInt());
    AggregationFunction* aggf;
    switch (function) {
        case AggregationFunction::AVERAGE : aggf = new AggregationFunctionAverage(); break;
        case AggregationFunction::COUNT : aggf = new AggregationFunctionCount(); break;
        case AggregationFunction::COUNTHISTOGRAM : aggf = new AggregationFunctionCountHistogram(); break;
        case AggregationFunction::GEOMETRICMEAN : aggf = new AggregationFunctionGeometricAverage(); break;
        case AggregationFunction::HARMONICMEAN : aggf = new AggregationFunctionHarmonicAverage(); break;
        case AggregationFunction::MAX : aggf = new AggregationFunctionMax(); break;
        case AggregationFunction::MAXHISTOGRAM : aggf = new AggregationFunctionMaxHistogram(); break;
        case AggregationFunction::MEDIAN : aggf = new AggregationFunctionMedian(); break;
        case AggregationFunction::MEDIANHISTOGRAM : aggf = new AggregationFunctionMedianHistogram(); break;
        case AggregationFunction::MIN : aggf = new AggregationFunctionMin(); break;
        case AggregationFunction::MINHISTOGRAM : aggf = new AggregationFunctionMinHistogram(); break;
        case AggregationFunction::MODE : aggf = new AggregationFunctionMode(); break;
        case AggregationFunction::SUM : aggf = new AggregationFunctionSum(); break;
        case AggregationFunction::SUMHISTOGRAM : aggf = new AggregationFunctionSumHistogram(); break;
        case AggregationFunction::STANDARDDEVIATION : aggf = new AggregationFunctionStandardDeviation(); break;
        case AggregationFunction::VARIANCE : aggf = new AggregationFunctionVariance(); break;
        default : aggf = NULL;
    }
    int column = x.toList().at(1).toInt();
    std::vector<std::pair<int, AggregationFunction*> >::iterator it;
    bool found = false;

    for (it = aggregatedcolumns.begin(); it < aggregatedcolumns.end(); ++it){
        if (it->first == column) {
            if (function != 0) {
                it->second = aggf;
            }
            else {
                delete it->second;
                aggregatedcolumns.erase(it);
            }
            found = true;
            break;
        }
    }
    if (!found && aggf != NULL) {
        aggregatedcolumns.push_back(std::pair<int,AggregationFunction*>(column,aggf));
    }
    updateFromProcessor();
}

void PlotDataGroupWidget:: groupby(){
    PlotDataGroup* plotDataGroup = dynamic_cast<PlotDataGroup*>(processor_);
    QVariant x = dynamic_cast<QAction*>(sender())->data();
    std::vector<std::pair<int, AggregationFunction*> > aggregatedcolumns2;
    for (size_t i = 0; i < aggregatedcolumns.size(); ++i) {
        aggregatedcolumns2.push_back(std::pair<int, AggregationFunction*>(aggregatedcolumns.at(i).first,aggregatedcolumns.at(i).second));
    }
    aggregatedcolumns.clear();
    plotDataGroup->groupbyWidget(x.toInt(),aggregatedcolumns2);
}

void PlotDataGroupWidget::selectResetAll() {
    PlotDataGroup* plotDataGroup = dynamic_cast<PlotDataGroup*>(processor_);
    aggregatedcolumns.clear();
    plotDataGroup->resetAllData();
}

void PlotDataGroupWidget::selectResetLast() {
    PlotDataGroup* plotDataGroup = dynamic_cast<PlotDataGroup*>(processor_);
    aggregatedcolumns.clear();
    plotDataGroup->resetLastData();
}


} //namespace voreen

