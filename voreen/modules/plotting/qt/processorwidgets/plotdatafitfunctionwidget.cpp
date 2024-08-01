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

#include "plotdatafitfunctionwidget.h"

#include "../plotpredicatedialog.h"
#include "../plotdataextendedtablemodel.h"
#include "../extendedtable.h"

#include "../../utils/functionlibrary.h"
#include "../../datastructures/plotfunction.h"

#include "voreen/qt/voreenapplicationqt.h"

#include <QGridLayout>
#include <QMainWindow>
#include <QLabel>
#include <QMenu>
#include <QTextEdit>


namespace voreen {

PlotDataFitFunctionWidget::PlotDataFitFunctionWidget(QWidget* parent, PlotDataFitFunction* PlotDataFitFunction)
    : QProcessorWidget(PlotDataFitFunction, parent)
    , table_()
    , widgetVector_(0)
{
    tgtAssert(PlotDataFitFunction, "No PlotDataFitFunction processor");
    proxyModel_ = new PlotDataExtendedTableModel(NULL, std::vector<std::pair<int,std::string> >(0), this);
    setWindowTitle(QString::fromStdString(PlotDataFitFunction->getGuiName()));
    contextMenuTable_ = new QMenu(this);
    fittingFunctionMenu_ = new QMenu(tr("Fitting Functions",0),this);
    multiRegessionMenu_ = new QMenu(tr("Multidimensional Regession",0),this);
    splineFunctionMenu_ = new QMenu(tr("Interpolation Functions",0),this);
    resize(256, 256);
}

PlotDataFitFunctionWidget::~PlotDataFitFunctionWidget() {
    delete proxyModel_;
    widgetVector_.clear();
    delete fittingFunctionMenu_;
    delete multiRegessionMenu_;
    delete splineFunctionMenu_;

}

void PlotDataFitFunctionWidget::initialize() {

    QProcessorWidget::initialize();
    QGridLayout* mainLayout = new QGridLayout(this);
    QLabel* label = new QLabel(this);
    label->setText(QString(tr("Fitting Function",0)));
    label->setAlignment(Qt::AlignHCenter);
    QFont font = label->font();
    font.setBold(true);
    font.setPointSize(14);
    label->setFont(font);
    mainLayout->addWidget(label,0,0,1,3);
    widgetVector_.push_back(label);
    label = new QLabel(this);
    label->setText(QString(tr("Function:",0)));
    label->setFixedWidth(100);
    mainLayout->addWidget(label,1,0,1,1);
    widgetVector_.push_back(label);
    QTextEdit* textedit = new QTextEdit(QString(tr("No FitFunction Selected",0)),this);
    textedit->setHorizontalScrollBarPolicy(Qt::ScrollBarAsNeeded);
    textedit->setReadOnly(true);
    textedit->setMaximumHeight(200);
    textedit->setMinimumHeight(20);
    mainLayout->addWidget(textedit,1,1,1,2);
    widgetVector_.push_back(textedit);
    label = new QLabel(this);
    label->setText(QString(tr("Fittingcolumn:",0)));
    label->setFixedWidth(100);
    mainLayout->addWidget(label,2,0,1,1);
    widgetVector_.push_back(label);
    label = new QLabel(this);
    label->setText(QString(tr("Column-Value",0)));
    mainLayout->addWidget(label,2,1,1,1);
    widgetVector_.push_back(label);
    label = new QLabel(this);
    label->setText(QString(tr("MSE:",0)));
    label->setFixedWidth(100);
    mainLayout->addWidget(label,3,0,1,1);
    widgetVector_.push_back(label);
    label = new QLabel(this);
    label->setText(QString(tr("MSE-Value",0)));
    mainLayout->addWidget(label,3,1,1,1);
    widgetVector_.push_back(label);

    delete proxyModel_;
    proxyModel_ = new PlotDataExtendedTableModel(NULL, std::vector<std::pair<int,std::string> >(0), this);
    table_ = new ExtendedTable(this);

    table_->setModel(proxyModel_);
    table_->setSortingEnabled(false);
    table_->verticalHeader()->hide();
    table_->resizeColumnsToContents();
    table_->setMinimumHeight(150);
    mainLayout->addWidget(table_, 4, 0, 1, 3);
    table_->setContextMenuPolicy(Qt::CustomContextMenu);
    table_->horizontalHeader()->setContextMenuPolicy(Qt::CustomContextMenu);

    QObject::connect(table_,SIGNAL(customContextMenuRequested(const QPoint&)),this,SLOT(tableContextMenu(const QPoint&)));
    QObject::connect(table_->horizontalHeader(),SIGNAL(customContextMenuRequested(const QPoint&)),this,SLOT(tableContextMenu(const QPoint&)));

    setLayout(mainLayout);

    for (size_t i = 1; i < widgetVector_.size(); ++i) {
        widgetVector_[i]->setVisible(false);
    }

    //show();

    initialized_ = true;
}


void PlotDataFitFunctionWidget::updateFromProcessor(){
    PlotDataFitFunction* plotDataFitFunction = dynamic_cast<PlotDataFitFunction*>(processor_);
    delete proxyModel_;
    PlotData* data = const_cast<PlotData*>(plotDataFitFunction->getPlotData());
    PlotFunction* function = const_cast<PlotFunction*>(plotDataFitFunction->getPlotFunction());
    PlotDataFitFunction::FittingValues fittingValues = plotDataFitFunction->getFittingValues();
    const std::vector<Property*> properties = plotDataFitFunction->getProperties();
    std::vector<std::pair<int,std::string> > stringvector;

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
    proxyModel_ = new PlotDataExtendedTableModel(data,stringvector, this,qColorValues[0],qColorValues[1],qColorValues[2],fittingValues.column);
    proxyModel_->setFunctionType(FunctionLibrary::NONE);
    table_->setModel(proxyModel_);
    table_->setSortingEnabled(false);
    table_->verticalHeader()->hide();
    table_->setSelectionMode(QAbstractItemView::SingleSelection);
    table_->setSelectionBehavior(QAbstractItemView::SelectColumns);
    table_->resizeColumnsToContents();
    for (size_t i = 1; i < widgetVector_.size(); ++i) {
        widgetVector_[i]->setVisible(false);
    }
    if (fittingValues.mse > -0.5) {
        for (size_t i = 1; i < widgetVector_.size(); ++i) {
            widgetVector_[i]->setVisible(true);
        }
        static_cast<QTextEdit*>(widgetVector_[2])->setText(QString::fromStdString(function->getPlotExpression().toString()));
        std::stringstream str("");
        str << fittingValues.column;
        static_cast<QLabel*>(widgetVector_[4])->setText(QString::fromStdString(str.str()));
        str.str("");
        str.clear();
        str << fittingValues.mse;
        static_cast<QLabel*>(widgetVector_[6])->setText(QString::fromStdString(str.str()));
    }
    else if (fittingValues.mse > -1.5) {
        widgetVector_[1]->setVisible(true);
        widgetVector_[2]->setVisible(true);
        widgetVector_[3]->setVisible(true);
        widgetVector_[4]->setVisible(true);
        static_cast<QTextEdit*>(widgetVector_[2])->setText(QString::fromStdString(function->getPlotExpression().toString()));
        std::stringstream str("");
        str << fittingValues.column;
        static_cast<QLabel*>(widgetVector_[4])->setText(QString::fromStdString(str.str()));
    }
    else {
        widgetVector_[1]->setVisible(true);
        widgetVector_[2]->setVisible(true);
        widgetVector_[3]->setVisible(true);
        widgetVector_[4]->setVisible(true);
        static_cast<QTextEdit*>(widgetVector_[2])->setText(QString(tr("No FittingFunction calculable",0)));
        std::stringstream str("");
        if (fittingValues.column >= 0) {
            str << fittingValues.column;
            static_cast<QLabel*>(widgetVector_[4])->setText(QString::fromStdString(str.str()));
        }
        else {
            widgetVector_[3]->setVisible(false);
            widgetVector_[4]->setVisible(false);
        }
    }
    if (maxColumnSize > 0) {
        for (int i = 0; i < data->getColumnCount(); ++i) {
            if (table_->columnWidth(i) > maxColumnSize)
                table_->setColumnWidth(i,maxColumnSize);
        }
    }
}


void PlotDataFitFunctionWidget::tableContextMenu(const QPoint& pos) {

    PlotDataFitFunction* plotDataFitFunction = dynamic_cast<PlotDataFitFunction*>(processor_);
    const std::vector<Property*> properties = plotDataFitFunction->getProperties();
    std::vector<bool> boolValues;
    for (size_t i = 0; i < properties.size(); ++i) {
        if (dynamic_cast<BoolProperty*>(properties[i])) {
            boolValues.push_back(static_cast<BoolProperty*>(properties[i])->get());
        }
    }
    QModelIndex index = table_->indexAt(pos);
    contextMenuTable_->clear();
    fittingFunctionMenu_->clear();
    multiRegessionMenu_->clear();
    splineFunctionMenu_->clear();
    PlotData* data = const_cast<PlotData*>(plotDataFitFunction->getPlotData());
    QAction* newAct;
    QList<QVariant>* qlist = new QList<QVariant>();
    std::stringstream s;

    bool numberColumn = (data->getColumnType(index.column()) == PlotBase::NUMBER);

    plot_t minX = 1;
    //plot_t maxX = 1;
    plot_t minY = 1;
    plot_t maxY = 1;
    if (numberColumn) {
        AggregationFunctionMin* aggmin = new AggregationFunctionMin();
        AggregationFunctionMax* aggmax = new AggregationFunctionMax();
        minY = data->aggregate(index.column(),aggmin);
        maxY = data->aggregate(index.column(),aggmax);
        minX = data->aggregate(0,aggmin);
        //maxX = data->aggregate(0,aggmax);
        delete aggmin;
        delete aggmax;
    }

    s.str("");
    s.clear();
    s << "Ignore False Values";
    newAct = new QAction(QString::fromStdString(s.str()),this);
    newAct->setCheckable(true);
    newAct->setChecked(boolValues[0]);
    contextMenuTable_->addAction(newAct);
    QObject::connect(newAct,SIGNAL(triggered()),this,SLOT(ignoreFalseValues()));
    qlist->clear();

    s.str("");
    s.clear();
    s << "Linear Regression";
    newAct = new QAction(QString::fromStdString(s.str()),this);
    qlist->push_back(FunctionLibrary::LINEAR);
    qlist->push_back(index.column());
    newAct->setData(QVariant(*qlist));
    fittingFunctionMenu_->addAction(newAct);
    QObject::connect(newAct,SIGNAL(triggered()),this,SLOT(addFunction()));
    qlist->clear();

    s.str("");
    s.clear();
    s << "Square Regression";
    newAct = new QAction(QString::fromStdString(s.str()),this);
    qlist->push_back(FunctionLibrary::SQUARE);
    qlist->push_back(index.column());
    newAct->setData(QVariant(*qlist));
    fittingFunctionMenu_->addAction(newAct);
    QObject::connect(newAct,SIGNAL(triggered()),this,SLOT(addFunction()));
    qlist->clear();

    s.str("");
    s.clear();
    s << "Cubic Regression";
    newAct = new QAction(QString::fromStdString(s.str()),this);
    qlist->push_back(FunctionLibrary::CUBIC);
    qlist->push_back(index.column());
    newAct->setData(QVariant(*qlist));
    fittingFunctionMenu_->addAction(newAct);
    QObject::connect(newAct,SIGNAL(triggered()),this,SLOT(addFunction()));
    qlist->clear();

    for (int i = 4; i < 11; ++i) {
        s.str("");
        s.clear();
        s << i << "-Dimension Regression";
        newAct = new QAction(QString::fromStdString(s.str()),this);
        qlist->push_back(FunctionLibrary::POLYNOMIAL);
        qlist->push_back(index.column());
        qlist->push_back(i);
        newAct->setData(QVariant(*qlist));
        multiRegessionMenu_->addAction(newAct);
        QObject::connect(newAct,SIGNAL(triggered()),this,SLOT(addFunction()));
        qlist->clear();
    }

    fittingFunctionMenu_->addMenu(multiRegessionMenu_);

    s.str("");
    s.clear();
    s << "Constant-Spline Regression";
    newAct = new QAction(QString::fromStdString(s.str()),this);
    qlist->push_back(FunctionLibrary::CONSTANTSPLINE);
    qlist->push_back(index.column());
    newAct->setData(QVariant(*qlist));
    splineFunctionMenu_->addAction(newAct);
    QObject::connect(newAct,SIGNAL(triggered()),this,SLOT(addFunction()));
    qlist->clear();

    s.str("");
    s.clear();
    s << "Simple-Spline Regression";
    newAct = new QAction(QString::fromStdString(s.str()),this);
    qlist->push_back(FunctionLibrary::SIMPLESPLINE);
    qlist->push_back(index.column());
    newAct->setData(QVariant(*qlist));
    splineFunctionMenu_->addAction(newAct);
    QObject::connect(newAct,SIGNAL(triggered()),this,SLOT(addFunction()));
    qlist->clear();

    s.str("");
    s.clear();
    s << "Square-Spline Regression";
    newAct = new QAction(QString::fromStdString(s.str()),this);
    qlist->push_back(FunctionLibrary::SQUARESPLINE);
    qlist->push_back(index.column());
    newAct->setData(QVariant(*qlist));
    splineFunctionMenu_->addAction(newAct);
    QObject::connect(newAct,SIGNAL(triggered()),this,SLOT(addFunction()));
    qlist->clear();

    s.str("");
    s.clear();
    s << "Cubic-Spline Regression";
    newAct = new QAction(QString::fromStdString(s.str()),this);
    qlist->push_back(FunctionLibrary::CUBICSPLINE);
    qlist->push_back(index.column());
    newAct->setData(QVariant(*qlist));
    splineFunctionMenu_->addAction(newAct);
    QObject::connect(newAct,SIGNAL(triggered()),this,SLOT(addFunction()));
    qlist->clear();

    splineFunctionMenu_->addSeparator();

    s.str("");
    s.clear();
    s << "Power Regression";
    newAct = new QAction(QString::fromStdString(s.str()),this);
    qlist->push_back(FunctionLibrary::POWER);
    qlist->push_back(index.column());
    newAct->setData(QVariant(*qlist));
    if ((minX > 0 && minY > 0) || boolValues[0])
        fittingFunctionMenu_->addAction(newAct);
    QObject::connect(newAct,SIGNAL(triggered()),this,SLOT(addFunction()));
    qlist->clear();

    fittingFunctionMenu_->addSeparator();

    s.str("");
    s.clear();
    s << "SQRT Regression";
    newAct = new QAction(QString::fromStdString(s.str()),this);
    qlist->push_back(FunctionLibrary::SQRT);
    qlist->push_back(index.column());
    newAct->setData(QVariant(*qlist));
    if ((minX >= 0) || boolValues[0])
        fittingFunctionMenu_->addAction(newAct);
    QObject::connect(newAct,SIGNAL(triggered()),this,SLOT(addFunction()));
    qlist->clear();

    fittingFunctionMenu_->addSeparator();

    s.str("");
    s.clear();
    s << "Logarithmic Regression";
    newAct = new QAction(QString::fromStdString(s.str()),this);
    qlist->push_back(FunctionLibrary::LOGARITHMIC);
    qlist->push_back(index.column());
    newAct->setData(QVariant(*qlist));
    if ((minX > 0) || boolValues[0])
        fittingFunctionMenu_->addAction(newAct);
    QObject::connect(newAct,SIGNAL(triggered()),this,SLOT(addFunction()));
    qlist->clear();

    s.str("");
    s.clear();
    s << "Exponential Regression";
    newAct = new QAction(QString::fromStdString(s.str()),this);
    qlist->push_back(FunctionLibrary::EXPONENTIAL);
    qlist->push_back(index.column());
    newAct->setData(QVariant(*qlist));
    if ((minY > 0) || boolValues[0])
        fittingFunctionMenu_->addAction(newAct);
    QObject::connect(newAct,SIGNAL(triggered()),this,SLOT(addFunction()));
    qlist->clear();

    s.str("");
    s.clear();
    s << "Sin Regression";
    newAct = new QAction(QString::fromStdString(s.str()),this);
    qlist->push_back(FunctionLibrary::SIN);
    qlist->push_back(index.column());
    newAct->setData(QVariant(*qlist));
    if ((minY > -1 && maxY < 1) || boolValues[0])
        fittingFunctionMenu_->addAction(newAct);
    QObject::connect(newAct,SIGNAL(triggered()),this,SLOT(addFunction()));
    qlist->clear();

    s.str("");
    s.clear();
    s << "Cos Regression";
    newAct = new QAction(QString::fromStdString(s.str()),this);
    qlist->push_back(FunctionLibrary::COS);
    qlist->push_back(index.column());
    newAct->setData(QVariant(*qlist));
    if ((minY > -1 && maxY < 1) || boolValues[0])
        fittingFunctionMenu_->addAction(newAct);
    QObject::connect(newAct,SIGNAL(triggered()),this,SLOT(addFunction()));
    qlist->clear();

    s.str("");
    s.clear();
    s << "Lagrange-Interpolation";
    newAct = new QAction(QString::fromStdString(s.str()),this);
    qlist->push_back(FunctionLibrary::INTERPOLATION);
    qlist->push_back(index.column());
    newAct->setData(QVariant(*qlist));
    splineFunctionMenu_->addAction(newAct);
    QObject::connect(newAct,SIGNAL(triggered()),this,SLOT(addFunction()));
    qlist->clear();

    if (!numberColumn) {
        fittingFunctionMenu_->setEnabled(false);
        splineFunctionMenu_->setEnabled(false);
        multiRegessionMenu_->setEnabled(false);
    }
    else {
        fittingFunctionMenu_->setEnabled(true);
        splineFunctionMenu_->setEnabled(true);
        multiRegessionMenu_->setEnabled(true);
    }
    contextMenuTable_->addMenu(fittingFunctionMenu_);
    contextMenuTable_->addMenu(splineFunctionMenu_);
    contextMenuTable_->popup(table_->mapToGlobal(pos));
    delete qlist;
}

void PlotDataFitFunctionWidget::addFunction() {
    PlotDataFitFunction* plotDataFitFunction = dynamic_cast<PlotDataFitFunction*>(processor_);
    QVariant x = dynamic_cast<QAction*>(sender())->data();
    FunctionLibrary::RegressionType function = static_cast<FunctionLibrary::RegressionType>(x.toList().at(0).toInt());
    int column = x.toList().at(1).toInt();
    int dimension = -1;
    if (x.toList().size() > 2)
        dimension = x.toList().at(2).toInt();
    if (function == FunctionLibrary::LINEAR)
        plotDataFitFunction->linearRegression(column);
    else if (function == FunctionLibrary::SQRT)
        plotDataFitFunction->sqrtRegression(column);
    else if (function == FunctionLibrary::SIMPLESPLINE)
        plotDataFitFunction->simpleSplineRegression(column);
    else if (function == FunctionLibrary::CONSTANTSPLINE)
        plotDataFitFunction->constantSpline(column);
    else if (function == FunctionLibrary::SQUARESPLINE)
        plotDataFitFunction->squareSplineRegression(column);
    //else if (function == FunctionLibrary::BSPLINE)
    //    plotDataFitFunction->bSplineRegression(column);
    else if (function == FunctionLibrary::CUBICSPLINE)
        plotDataFitFunction->cubicSplineRegression(column);
    else if (function == FunctionLibrary::SQUARE)
        plotDataFitFunction->multiRegression(2,column);
    else if (function == FunctionLibrary::CUBIC)
        plotDataFitFunction->multiRegression(3,column);
    else if (function == FunctionLibrary::POLYNOMIAL)
        plotDataFitFunction->multiRegression(dimension,column);
    //else if (function == FunctionLibrary::SPLINE)
    //    plotDataFitFunction->splineRegression(column);
    else if (function == FunctionLibrary::LOGARITHMIC)
        plotDataFitFunction->logarithmicRegression(column);
    else if (function == FunctionLibrary::POWER)
        plotDataFitFunction->powerRegression(column);
    else if (function == FunctionLibrary::EXPONENTIAL)
        plotDataFitFunction->exponentialRegression(column);
    else if (function == FunctionLibrary::SIN)
       plotDataFitFunction->sinRegression(column);
    else if (function == FunctionLibrary::COS)
        plotDataFitFunction->cosRegression(column);
    else if (function == FunctionLibrary::INTERPOLATION)
        plotDataFitFunction->interpolRegression(column);
}

void PlotDataFitFunctionWidget::ignoreFalseValues() {
    PlotDataFitFunction* plotDataFitFunction = dynamic_cast<PlotDataFitFunction*>(processor_);
    plotDataFitFunction->ignoreFalseValues();
}

} //namespace voreen

