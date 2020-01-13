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

#include "plotdataselectwidget.h"

#include "../plotpredicatedialog.h"
#include "../plotdataextendedtablemodel.h"
#include "../extendedtable.h"
#include "../../datastructures/plotdata.h"
#include "../../datastructures/plotpredicate.h"

#include "voreen/qt/voreenapplicationqt.h"

#include <utility>
#include <QGridLayout>
#include <QMainWindow>
#include <QLabel>
#include <QCheckBox>
#include <QLineEdit>
#include <QInputDialog>
#include <QMenu>


namespace voreen {

const std::string PlotDataSelectWidget::loggerCat_("voreen.PlotDataSelectWidget");

PlotDataSelectWidget::PlotDataSelectWidget(QWidget* parent, PlotDataSelect* plotDataSelect)
    : QProcessorWidget(plotDataSelect, parent)
    , table_()
    , widgetVector_(0)
{
    function_ = FunctionLibrary::NONE;
    tgtAssert(plotDataSelect, "No PlotDataSelect processor");
    proxyModel_ = new PlotDataExtendedTableModel(NULL, std::vector<std::pair<int,std::string> >(0), this);
    setWindowTitle(QString::fromStdString(plotDataSelect->getGuiName()));
    contextMenuTable_ = new QMenu(this);
    beforeMenu_ = new QMenu(tr("Insert Before",0),this);
    afterMenu_ = new QMenu(tr("Insert After",0),this);
    switchMenu_ = new QMenu(tr("Switch With",0),this);
    copyColumnMenu_ = new QMenu(tr("Copy Column",0),this);
    resize(256, 256);
    pPVProperty_ = NULL;
}

PlotDataSelectWidget::~PlotDataSelectWidget() {
    delete proxyModel_;
    widgetVector_.clear();
}

void PlotDataSelectWidget::initialize() {
    QProcessorWidget::initialize();
    QGridLayout* mainLayout = new QGridLayout(this);
    QLabel* label = new QLabel(this);
    label->setText(QString(tr("Function Type",0)));
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
    table_->setMinimumHeight(150);
    table_->setMinimumWidth(150);
    mainLayout->addWidget(table_, 3, 0, 1, 3);
    table_->setContextMenuPolicy(Qt::CustomContextMenu);
    table_->horizontalHeader()->setContextMenuPolicy(Qt::CustomContextMenu);

    QObject::connect(table_,SIGNAL(customContextMenuRequested(const QPoint&)),this,SLOT(tableContextMenu(const QPoint&)));
    QObject::connect(table_->horizontalHeader(),SIGNAL(customContextMenuRequested(const QPoint&)),this,SLOT(tableContextMenu(const QPoint&)));

    for (size_t i = 0; i < widgetVector_.size(); ++i) {
        widgetVector_.at(i)->setVisible(false);
    }

    setLayout(mainLayout);

    //show();

    initialized_ = true;
}

void PlotDataSelectWidget::updateFromProcessor(){
    PlotDataSelect* plotDataSelect = dynamic_cast<PlotDataSelect*>(processor_);
    const std::vector<Property*>& properties = plotDataSelect->getProperties();
    function_ = FunctionLibrary::NONE;
    pPVProperty_ = NULL;
    tgt::vec4 colorValue;
    std::vector<QColor> qColorValues;
    int maxColumnSize = -1;
    for (size_t i = 0; i < properties.size(); ++i) {
        if (dynamic_cast<BoolProperty*>(properties[i])) {
            BoolProperty* bProperty = dynamic_cast<BoolProperty*>(properties[i]);
            if (bProperty->get()) {
                if (bProperty->getID() == "select") {
                    function_ = FunctionLibrary::SELECT;
                    dynamic_cast<QLabel*>(widgetVector_.at(0))->setText(tr("Select",0));
                    for (size_t i = 0; i < widgetVector_.size(); ++i) {
                        widgetVector_.at(i)->setVisible(true);
                    }
                }
                else if (bProperty->getID() == "columnorder") {
                    function_ = FunctionLibrary::COLUMNORDER;
                    dynamic_cast<QLabel*>(widgetVector_.at(0))->setText(tr("Column Order",0));
                    widgetVector_.at(0)->setVisible(true);
                    for (size_t i = 1; i < widgetVector_.size(); ++i) {
                        widgetVector_.at(i)->setVisible(false);
                    }
                }
                else if (bProperty->getID() == "filter") {
                    function_ = FunctionLibrary::FILTER;
                    dynamic_cast<QLabel*>(widgetVector_.at(0))->setText(tr("Filter",0));
                    widgetVector_.at(0)->setVisible(true);
                    for (size_t i = 1; i < widgetVector_.size(); ++i) {
                        widgetVector_.at(i)->setVisible(false);
                    }
                }
            }
        }
        else if (dynamic_cast<PlotPredicateProperty*>(properties[i])) {
            pPVProperty_ = dynamic_cast<PlotPredicateProperty*>(properties[i]);
        }
        else if (dynamic_cast<FloatVec4Property*>(properties[i])) {
            colorValue = static_cast<FloatVec4Property*>(properties[i])->get();
            qColorValues.push_back(QColor(int(colorValue.r*255),
                int(colorValue.g*255),int(colorValue.b*255),
                int(colorValue.a*255)));
        }
        else if (dynamic_cast<IntProperty*>(properties[i])) {
            maxColumnSize = static_cast<IntProperty*>(properties[i])->get();
        }
    }
    if (function_ == 0) {
        for (size_t i = 0; i < widgetVector_.size(); ++i) {
            widgetVector_.at(i)->setVisible(false);
        }
    }
    delete proxyModel_;
    PlotData* data = const_cast<PlotData*>(plotDataSelect->getPlotData());

    std::vector<std::pair<int,std::string> > stringvector;
    for (size_t i = 0; i < pPVProperty_->get().size(); ++i) {
        stringvector.push_back(std::pair<int,std::string>(pPVProperty_->get().at(i).first,pPVProperty_->get().at(i).second->toString()));
    }
    for (int i = 0; i < data->getColumnCount(); ++i) {
        stringvector.push_back(std::pair<int,std::string>(i,"No Predicate"));
    }

    proxyModel_ = new PlotDataExtendedTableModel(data,stringvector, this,qColorValues[0],qColorValues[1]);
    proxyModel_ ->setFunctionType(function_);
    table_->setModel(proxyModel_);
    table_->setSortingEnabled(false);
    table_->verticalHeader()->hide();
    if (function_ == FunctionLibrary::NONE) {
        table_->setSelectionMode(QAbstractItemView::NoSelection);
        table_->setSelectionBehavior(QAbstractItemView::SelectItems);
    }
    else if (function_ == FunctionLibrary::SELECT) {
        table_->setSelectionMode(QAbstractItemView::ExtendedSelection);
        table_->setSelectionBehavior(QAbstractItemView::SelectItems);
    }
    else if (function_ == FunctionLibrary::COLUMNORDER || function_ == FunctionLibrary::FILTER) {
        table_->setSelectionMode(QAbstractItemView::SingleSelection);
        table_->setSelectionBehavior(QAbstractItemView::SelectColumns);
    }
    table_->resizeColumnsToContents();
    if (maxColumnSize > 0)
        for (int i = 0; i < data->getColumnCount(); ++i) {
            if (table_->columnWidth(i) > maxColumnSize)
                table_->setColumnWidth(i,maxColumnSize);
        }
}

void PlotDataSelectWidget::functionSelect() {
   PlotDataSelect* plotDataSelect = dynamic_cast<PlotDataSelect*>(processor_);
   QVariant x = dynamic_cast<QAction*>(sender())->data();
   plotDataSelect->selectFunctionalityType(static_cast<FunctionLibrary::ProcessorFunctionalityType>(x.toInt()));
}

void PlotDataSelectWidget::tableContextMenu(const QPoint& pos) {
    PlotDataSelect* plotDataSelect = dynamic_cast<PlotDataSelect*>(processor_);
    QModelIndex index = table_->indexAt(pos);
    std::vector<int> selectedColumns = table_->selectedColumns();
    std::vector<int> selectedRows = table_->selectedRows();
    std::pair<QPoint,QPoint> selectedArea = table_->selectedArea();
    contextMenuTable_->clear();
    PlotData* data = const_cast<PlotData*>(plotDataSelect->getPlotData());
    QAction* newAct;
    QList<QVariant>* qlist = new QList<QVariant>();
    std::stringstream str;
    std::string menuString;
    // create popupmenu
    if (function_ == FunctionLibrary::NONE) {
        str.str("");
        str.clear();
        str << "Select";
        newAct = new QAction(QString::fromStdString(str.str()),this);
        newAct->setData(QVariant(static_cast<int>(FunctionLibrary::SELECT)));
        contextMenuTable_->addAction(newAct);
        QObject::connect(newAct,SIGNAL(triggered()),this,SLOT(functionSelect()));

        menuString = "Column Order";
        newAct = new QAction(QString::fromStdString(menuString),this);
        newAct->setData(QVariant(static_cast<int>(FunctionLibrary::COLUMNORDER)));
        contextMenuTable_->addAction(newAct);
        QObject::connect(newAct,SIGNAL(triggered()),this,SLOT(functionSelect()));

        menuString = "Filter";
        newAct = new QAction(QString::fromStdString(menuString),this);
        newAct->setData(QVariant(static_cast<int>(FunctionLibrary::FILTER)));
        contextMenuTable_->addAction(newAct);
        QObject::connect(newAct,SIGNAL(triggered()),this,SLOT(functionSelect()));

        contextMenuTable_->popup(table_->mapToGlobal(pos));
    }
    else if (function_ != FunctionLibrary::NONE) {

        menuString = "Deactivate Function";
        newAct = new QAction(QString::fromStdString(menuString),this);
        newAct->setData(QVariant(0));
        contextMenuTable_->addAction(newAct);
        QObject::connect(newAct,SIGNAL(triggered()),this,SLOT(functionSelect()));

        if (function_ == FunctionLibrary::SELECT || function_ == FunctionLibrary::COLUMNORDER) {
            contextMenuTable_->addSeparator();

            str.str("");
            str.clear();
            if (plotDataSelect->getResetListSize() > 0) {
                str << plotDataSelect->getResetListSize();
                menuString = "Reset All (" + str.str() + " Actions)";
            }
            else
                menuString = "Reset All";
            newAct = new QAction(QString::fromStdString(menuString),this);
            if (plotDataSelect->getResetListSize() == 0)
                newAct->setEnabled(false);
            contextMenuTable_->addAction(newAct);
            QObject::connect(newAct,SIGNAL(triggered()),this,SLOT(selectResetAll()));

            menuString = "Reset Last";
            newAct = new QAction(QString::fromStdString(menuString),this);
            if (plotDataSelect->getResetListSize() == 0)
                newAct->setEnabled(false);
            contextMenuTable_->addAction(newAct);
            QObject::connect(newAct,SIGNAL(triggered()),this,SLOT(selectResetLast()));

            contextMenuTable_->addSeparator();

            if (function_ == FunctionLibrary::SELECT) {
                if (index.column() >= 0) {
                    menuString = "Disable \""+data->getColumnLabel(index.column())+"\" Column";
                }
                else {
                    menuString = "Selecting not Well-Defined";
                }
                newAct = new QAction(QString::fromStdString(menuString),this);
                if (index.column() >= 0) {
                    newAct->setData(QVariant(index.column()));
                }
                else {
                    newAct->setEnabled(false);
                }
                contextMenuTable_->addAction(newAct);
                QObject::connect(newAct,SIGNAL(triggered()),this,SLOT(selectChange()));

                // only show this part of the Menu, if something is selected
                if (selectedColumns.size() > 0 || selectedRows.size() > 0) {
                    if (selectedColumns.size() > 1) {
                        qlist->clear();
                        str.str("");
                        str.clear();
                        str << selectedColumns.size();
                        menuString = "Disable Selected Columns ("+str.str()+")";
                        newAct = new QAction(QString::fromStdString(menuString),this);
                        for (size_t j = 0; j < selectedColumns.size(); ++j) {
                            qlist->push_back(QVariant(selectedColumns.at(j)));
                        }
                        newAct->setData(QVariant(*qlist));
                        contextMenuTable_->addAction(newAct);
                        QObject::connect(newAct,SIGNAL(triggered()),this,SLOT(selectChange()));
                    }
                    qlist->clear();
                    str.str("");
                    str.clear();
                    str << selectedColumns.size();
                    menuString = "Select Columns ("+str.str()+")";
                    newAct = new QAction(QString::fromStdString(menuString),this);
                    for (size_t j = 0; j < selectedColumns.size(); ++j) {
                        qlist->push_back(QVariant(selectedColumns.at(j)));
                    }
                    newAct->setData(QVariant(*qlist));
                    contextMenuTable_->addAction(newAct);
                    QObject::connect(newAct,SIGNAL(triggered()),this,SLOT(selectColumnChange()));

                    qlist->clear();
                    str.str("");
                    str.clear();
                    str << selectedRows.size();
                    menuString = "Select Rows ("+str.str()+")";
                    newAct = new QAction(QString::fromStdString(menuString),this);
                    for (size_t j = 0; j < selectedRows.size(); ++j) {
                        qlist->push_back(QVariant(selectedRows.at(j)-proxyModel_->getAddLines()));
                    }
                    newAct->setData(QVariant(*qlist));
                    contextMenuTable_->addAction(newAct);
                    QObject::connect(newAct,SIGNAL(triggered()),this,SLOT(selectRowChange()));

                    qlist->clear();
                    str.str("");
                    str.clear();
                    str << selectedRows.size();
                    menuString = "Select Columns && Rows ("+str.str()+")";
                    newAct = new QAction(QString::fromStdString(menuString),this);
                    QList<QVariant>* qlist1 = new QList<QVariant>();
                    QList<QVariant>* qlist2 = new QList<QVariant>();
                    for (size_t j = 0; j < selectedColumns.size(); ++j) {
                        qlist1->push_back(QVariant(selectedColumns.at(j)));
                    }
                    for (size_t j = 0; j < selectedRows.size(); ++j) {
                        qlist2->push_back(QVariant(selectedRows.at(j)-proxyModel_->getAddLines()));
                    }
                    qlist->push_back(QVariant(*qlist1));
                    qlist->push_back(QVariant(*qlist2));
                    newAct->setData(QVariant(*qlist));
                    contextMenuTable_->addAction(newAct);
                    QObject::connect(newAct,SIGNAL(triggered()),this,SLOT(selectColumnRowChange()));

                    qlist->clear();
                    str.str("");
                    str.clear();
                    str << selectedArea.second.x()- selectedArea.first.x() + 1  << "x" << selectedArea.second.y() - selectedArea.first.y() + 1;
                    menuString = "Select Area ("+str.str()+")";
                    newAct = new QAction(QString::fromStdString(menuString),this);
                    qlist->push_back(QVariant(QPoint(selectedArea.first.x(),selectedArea.first.y()-proxyModel_->getAddLines())));
                    qlist->push_back(QVariant(QPoint(selectedArea.second.x(),selectedArea.second.y()-proxyModel_->getAddLines())));
                    newAct->setData(QVariant(*qlist));
                    contextMenuTable_->addAction(newAct);
                    QObject::connect(newAct,SIGNAL(triggered()),this,SLOT(selectAreaChange()));

                }

                contextMenuTable_->addSeparator();

                if (index.column() >= 0) {
                    menuString = "Rename \""+data->getColumnLabel(index.column())+"\"";
                    newAct = new QAction(QString::fromStdString(menuString),this);
                    newAct->setData(QVariant(index.column()));
                    contextMenuTable_->addAction(newAct);
                    QObject::connect(newAct,SIGNAL(triggered()),this,SLOT(renameColumn()));
                }

                const std::vector<Property*>& properties = plotDataSelect->getProperties();
                PlotPredicateProperty* pPVProperty = NULL;
                for (size_t i = 0; i < properties.size(); ++i) {
                    if (dynamic_cast<PlotPredicateProperty*>(properties[i])) {
                        pPVProperty = dynamic_cast<PlotPredicateProperty*>(properties[i]);
                        break;
                    }
                }
                menuString = "Set ";
                if (index.column() >= 0) {
                    for (size_t i = 0; i < pPVProperty->get().size(); ++i) {
                        if (index.column() == pPVProperty->get().at(i).first) {
                            menuString = "Change ";
                        }
                    }
                }
                menuString += "Predicate";
                newAct = new QAction(QString::fromStdString(menuString),this);
                if (index.column() >= 0) {
                    newAct->setData(QVariant(index.column()));
                }
                else {
                    newAct->setEnabled(false);
                }
                contextMenuTable_->addAction(newAct);
                QObject::connect(newAct,SIGNAL(triggered()),this,SLOT(predicateChange()));
            }
            beforeMenu_->clear();
            for (int i = 0; i < data->getColumnCount(); ++i) {
                if (i != index.column()) {
                    qlist->clear();
                    menuString = "\""+data->getColumnLabel(i)+"\" Column";
                    newAct = new QAction(QString::fromStdString(menuString),this);
                    qlist->push_back(QVariant(index.column()));
                    qlist->push_back(QVariant(i));
                    newAct->setData(QVariant(*qlist));
                    beforeMenu_->addAction(newAct);
                    QObject::connect(newAct,SIGNAL(triggered()),this,SLOT(sortChangeBefore()));
                }
            }
            if (index.column() >= 0) {
                beforeMenu_->setEnabled(true);
            }
            else {
                beforeMenu_->setEnabled(false);
            }
            contextMenuTable_->addMenu(beforeMenu_);
            afterMenu_->clear();
            for (int i = 0; i < data->getColumnCount(); ++i) {
                if (i != index.column()) {
                    qlist->clear();
                    menuString = "\""+data->getColumnLabel(i)+"\" Column";
                    newAct = new QAction(QString::fromStdString(menuString),this);
                    qlist->push_back(QVariant(index.column()));
                    qlist->push_back(QVariant(i));
                    newAct->setData(QVariant(*qlist));
                    afterMenu_->addAction(newAct);
                    QObject::connect(newAct,SIGNAL(triggered()),this,SLOT(sortChangeAfter()));
                }
            }
            if (index.column() >= 0) {
                afterMenu_->setEnabled(true);
            }
            else {
                afterMenu_->setEnabled(false);
            }
            contextMenuTable_->addMenu(afterMenu_);
            switchMenu_->clear();
            for (int i = 0; i < data->getColumnCount(); ++i) {
                if (i != index.column()) {
                    qlist->clear();
                    menuString = "\""+data->getColumnLabel(i)+"\" Column";
                    newAct = new QAction(QString::fromStdString(menuString),this);
                    qlist->push_back(QVariant(index.column()));
                    qlist->push_back(QVariant(i));
                    newAct->setData(QVariant(*qlist));
                    switchMenu_->addAction(newAct);
                    QObject::connect(newAct,SIGNAL(triggered()),this,SLOT(sortChangeSwitch()));
                }
            }
            if (index.column() >= 0) {
                switchMenu_->setEnabled(true);
            }
            else {
                switchMenu_->setEnabled(false);
            }
            contextMenuTable_->addMenu(switchMenu_);

            if (function_ == FunctionLibrary::SELECT) {
                copyColumnMenu_->clear();
                for (int i = 0; i < data->getColumnCount(); ++i) {
                    qlist->clear();
                    menuString = "\""+data->getColumnLabel(i)+"\" Column";
                    newAct = new QAction(QString::fromStdString(menuString),this);
                    qlist->push_back(QVariant(index.column()));
                    qlist->push_back(QVariant(i));
                    newAct->setData(QVariant(*qlist));
                    copyColumnMenu_->addAction(newAct);
                    QObject::connect(newAct,SIGNAL(triggered()),this,SLOT(copyColumn()));
                }

                contextMenuTable_->addMenu(copyColumnMenu_);
            }
        }
        else if (function_ == FunctionLibrary::FILTER) {

            contextMenuTable_->addSeparator();

            str.str("");
            str.clear();
            if (plotDataSelect->getResetListSize() > 0) {
                str << plotDataSelect->getResetListSize();
                menuString = "Reset All (" + str.str() + " Actions)";
            }
            else
                menuString = "Reset All";
            newAct = new QAction(QString::fromStdString(menuString),this);
            if (plotDataSelect->getResetListSize() == 0)
                newAct->setEnabled(false);
            contextMenuTable_->addAction(newAct);
            QObject::connect(newAct,SIGNAL(triggered()),this,SLOT(selectResetAll()));

            menuString = "Reset Last";
            newAct = new QAction(QString::fromStdString(menuString),this);
            if (plotDataSelect->getResetListSize() == 0)
                newAct->setEnabled(false);
            contextMenuTable_->addAction(newAct);
            QObject::connect(newAct,SIGNAL(triggered()),this,SLOT(selectResetLast()));


            contextMenuTable_->addSeparator();

            menuString = "Delete All Rows With EMPTY Cells";
            newAct = new QAction(QString::fromStdString(menuString),this);
            contextMenuTable_->addAction(newAct);
            QObject::connect(newAct,SIGNAL(triggered()),this,SLOT(filterEmptyCells()));

            menuString = "Delete All Rows With Non Alpha-Numeric Cells";
            newAct = new QAction(QString::fromStdString(menuString),this);
            contextMenuTable_->addAction(newAct);
            QObject::connect(newAct,SIGNAL(triggered()),this,SLOT(filterAlphaNumeric()));
        }
        contextMenuTable_->popup(table_->mapToGlobal(pos));
    }
    delete qlist;
}

void PlotDataSelectWidget::selectChange() {
    PlotDataSelect* plotDataSelect = dynamic_cast<PlotDataSelect*>(processor_);
    QVariant x = dynamic_cast<QAction*>(sender())->data();
    if (x.canConvert<QList<QVariant> >()) {
        QList<QVariant> qlist = x.toList();
        std::vector<int> cVector;
        for(int i = 0; i < qlist.size(); ++i) {
            cVector.push_back(qlist.at(i).toInt());
        }
        plotDataSelect->disSelectColumns(cVector);
        cVector.clear();
    }
    else {
        plotDataSelect->disSelectColumn(x.toInt());
    }
}

void PlotDataSelectWidget::selectColumnChange() {
    PlotDataSelect* plotDataSelect = dynamic_cast<PlotDataSelect*>(processor_);
    QVariant x = dynamic_cast<QAction*>(sender())->data();
    if (x.canConvert<QList<QVariant> >()) {
        QList<QVariant> qlist = x.toList();
        std::vector<int> cVector;
        for(int i = 0; i < qlist.size(); ++i) {
            cVector.push_back(qlist.at(i).toInt());
        }
        plotDataSelect->selectColumns(cVector);
        cVector.clear();
    }
}

void PlotDataSelectWidget::selectRowChange() {
    PlotDataSelect* plotDataSelect = dynamic_cast<PlotDataSelect*>(processor_);
    PlotData* data = const_cast<PlotData*>(plotDataSelect->getPlotData());
    QVariant x = dynamic_cast<QAction*>(sender())->data();
    if (x.canConvert<QList<QVariant> >()) {
        QList<QVariant> qlist = x.toList();
        std::vector<int> cVector;
        for(int i = 0; i < data->getColumnCount(); ++i) {
            cVector.push_back(i);
        }
        std::vector<int> rVector;
        for(int i = 0; i < qlist.size(); ++i) {
            rVector.push_back(qlist.at(i).toInt());
        }
        plotDataSelect->selectColumnsRows(cVector,rVector);
        cVector.clear();
        rVector.clear();
    }
}

void PlotDataSelectWidget::selectColumnRowChange() {
    PlotDataSelect* plotDataSelect = dynamic_cast<PlotDataSelect*>(processor_);
    QVariant x = dynamic_cast<QAction*>(sender())->data();
    if (x.canConvert<QList<QVariant> >()) {
        QList<QVariant> qlist = x.toList();
        std::vector<int> cVector;
        for(int i = 0; i < qlist.at(0).toList().size(); ++i) {
            cVector.push_back(qlist.at(0).toList().at(i).toInt());
        }
        std::vector<int> rVector;
        for(int i = 0; i < qlist.at(1).toList().size(); ++i) {
            rVector.push_back(qlist.at(1).toList().at(i).toInt());
        }
        plotDataSelect->selectColumnsRows(cVector,rVector);
        cVector.clear();
        rVector.clear();
    }
}

void PlotDataSelectWidget::selectAreaChange() {
    PlotDataSelect* plotDataSelect = dynamic_cast<PlotDataSelect*>(processor_);
    QVariant x = dynamic_cast<QAction*>(sender())->data();
    if (x.canConvert<QList<QVariant> >()) {
        QList<QVariant> qlist = x.toList();
        std::pair<int,int> point1;
        std::pair<int,int> point2;
        QPoint qpoint = qlist.at(0).toPoint();
        point1.first = qpoint.x();
        point1.second = qpoint.y();
        qpoint = qlist.at(1).toPoint();
        point2.first = qpoint.x();
        point2.second = qpoint.y();
        plotDataSelect->selectIncludeArea(point1,point2);
    }
}

void PlotDataSelectWidget::selectResetAll() {
    PlotDataSelect* plotDataSelect = dynamic_cast<PlotDataSelect*>(processor_);
    plotDataSelect->resetAllData();
}

void PlotDataSelectWidget::selectResetLast() {
    PlotDataSelect* plotDataSelect = dynamic_cast<PlotDataSelect*>(processor_);
    plotDataSelect->resetLastData();
}

void PlotDataSelectWidget::sortChangeBefore() {
    PlotDataSelect* plotDataSelect = dynamic_cast<PlotDataSelect*>(processor_);
    QVariant x = dynamic_cast<QAction*>(sender())->data();
    int oldcolumn = x.toList().at(0).toInt();
    int newcolumn = x.toList().at(1).toInt();
    plotDataSelect->slideColumnBefore(oldcolumn, newcolumn);
}

void PlotDataSelectWidget::sortChangeAfter() {
    PlotDataSelect* plotDataSelect = dynamic_cast<PlotDataSelect*>(processor_);
    QVariant x = dynamic_cast<QAction*>(sender())->data();
    int oldcolumn = x.toList().at(0).toInt();
    int newcolumn = x.toList().at(1).toInt();
    plotDataSelect->slideColumnAfter(oldcolumn, newcolumn);
}

void PlotDataSelectWidget::sortChangeSwitch() {
    PlotDataSelect* plotDataSelect = dynamic_cast<PlotDataSelect*>(processor_);
    QVariant x = dynamic_cast<QAction*>(sender())->data();
    int oldcolumn = x.toList().at(0).toInt();
    int newcolumn = x.toList().at(1).toInt();
    plotDataSelect->slideColumnSwitch(oldcolumn, newcolumn);
}

void PlotDataSelectWidget::predicateChange() {
    PlotDataSelect* plotDataSelect = dynamic_cast<PlotDataSelect*>(processor_);
    int x = dynamic_cast<QAction*>(sender())->data().toInt();

    // we edit only a copy of the according predicate, but we have to find it first
    PlotPredicate* pred = 0;
    std::vector< std::pair<int, PlotPredicate*> > predsCopy = pPVProperty_->getCloned();
    std::vector< std::pair<int, PlotPredicate*> >::iterator it = predsCopy.begin();
    for ( ; it != predsCopy.end(); ++it) {
        if (it->first == x) {
            pred = it->second;
            break;
        }
    }

    // instantiate a PlotPredicateDialog
    PlotPredicateDialog preddialog(pred, plotDataSelect->getPlotData()->getColumnType(x) == PlotData::STRING, false, this);
    if (preddialog.exec() == QDialog::Accepted) {
        // replace predicate in predicate vector with deep copy of dialog predicate if that exists
        if (preddialog.getPlotPredicate() != 0) {
            if (it !=predsCopy.end()) {
                delete it->second;
                it->second = preddialog.getPlotPredicate()->clone();
            }
            else {
                predsCopy.push_back(std::make_pair(x, preddialog.getPlotPredicate()->clone()));
            }
        }
        // delete old predicate if existed
        else {
            if (it !=predsCopy.end()) {
                delete it->second;
                predsCopy.erase(it);
            }
        }
        // update predicate property
        pPVProperty_->set(predsCopy);
    }
    // delete all predicates in copy of predicate vector
    for (it = predsCopy.begin() ; it != predsCopy.end(); ++it) {
        delete it->second;
    }

    plotDataSelect->selectPredicate();
}

void PlotDataSelectWidget::filterEmptyCells() {
    PlotDataSelect* plotDataSelect = dynamic_cast<PlotDataSelect*>(processor_);
    plotDataSelect->filterEmptyCells();
}

void PlotDataSelectWidget::filterAlphaNumeric() {
    PlotDataSelect* plotDataSelect = dynamic_cast<PlotDataSelect*>(processor_);
    plotDataSelect->filterAlphaNumeric();
}

void PlotDataSelectWidget::copyColumn() {
    PlotDataSelect* plotDataSelect = dynamic_cast<PlotDataSelect*>(processor_);
    QVariant x = dynamic_cast<QAction*>(sender())->data();
    int newcolumn = x.toList().at(0).toInt();
    int copycolumn = x.toList().at(1).toInt();
    plotDataSelect->copyColumn(copycolumn,newcolumn);
}

void PlotDataSelectWidget::renameColumn() {
    PlotDataSelect* plotDataSelect = dynamic_cast<PlotDataSelect*>(processor_);
    QVariant x = dynamic_cast<QAction*>(sender())->data();
    PlotData* data = const_cast<PlotData*>(plotDataSelect->getPlotData());
    std::string oldtext = data->getColumnLabel(x.toInt());
    bool ok;
    QString text = QInputDialog::getText(this,"Change ColumnLabel", "Enter the new Label:", QLineEdit::Normal,
        QString::fromStdString(oldtext), &ok);
    if (ok && !text.isEmpty())
        plotDataSelect->renameColumn(x.toInt(),text.toStdString());
}

} //namespace voreen

