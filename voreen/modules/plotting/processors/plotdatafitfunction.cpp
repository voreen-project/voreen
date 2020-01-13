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

#include "plotdatafitfunction.h"

#include "voreen/core/processors/processorwidgetfactory.h"
#include "voreen/core/processors/processorwidget.h"

#include "../datastructures/plotfunction.h"
#include "../datastructures/plotbase.h"
#include "../datastructures/plotdata.h"
#include "../datastructures/plotrow.h"
#include "../datastructures/interval.h"

namespace voreen {

const std::string PlotDataFitFunction::loggerCat_("voreen.PlotDataFitFunction");


PlotDataFitFunction::PlotDataFitFunction()
    :Processor()
    , outPort_(Port::OUTPORT,"PlotData.OutPort")
    , inPort_(Port::INPORT,"PlotData.InPort")
    , pData_(0)
    , pDataOut_(0)
    , ignoreFalseValues_("Ignore False Values","Ignore False Values",true)
    , expressionNameInput_("expressionNameInput_","Expression-Name","f",Processor::VALID)
    , expressionText_("Expression Format","Expression Format",Processor::VALID)
    , selfDescription_("selfDescription","Custom Text","",Processor::VALID)
    , maxLength_("maxLength","Max Expression Length",200,1,1000,Processor::VALID)
    , recalculate_("recalculate","recalculate",Processor::VALID)
    , maxTableColumnWidth_("Max Column Width","Max Column Width",0,0,500,Processor::VALID)
    , keyColumnColor_("Key-Column Color","Key-Column Color",tgt::vec4(200.0/255.0,200.0/255.0,200.0/255.0,1.0),
        Processor::VALID)
    , dataColumnColor_("Data-Column Color","Data-Column Color",tgt::vec4(240.0/255.0,240.0/255.0,240.0/255.0,1.0),
        Processor::VALID)
    , fittingColumnColor_("Fitting-Column Color","Fitting-Column Color",tgt::vec4(1.0,180.0/255.0,142.0/255.0,1.0),
        Processor::VALID)
{
    keyColumnColor_.onChange(MemberFunctionCallback<PlotDataFitFunction>(this, &PlotDataFitFunction::updateView));
    dataColumnColor_.onChange(MemberFunctionCallback<PlotDataFitFunction>(this, &PlotDataFitFunction::updateView));
    fittingColumnColor_.onChange(MemberFunctionCallback<PlotDataFitFunction>(this, &PlotDataFitFunction::updateView));
    maxTableColumnWidth_.onChange(MemberFunctionCallback<PlotDataFitFunction>(this, &PlotDataFitFunction::updateView));
    maxTableColumnWidth_.onChange(MemberFunctionCallback<PlotDataFitFunction>(this, &PlotDataFitFunction::updateView));
    recalculate_.onChange(MemberFunctionCallback<PlotDataFitFunction>(this, &PlotDataFitFunction::recalculate));

    fittingValues_.column = -1;
    fittingValues_.regressionType = FunctionLibrary::NOREGRESSION;
    fittingValues_.dimension = 2;
    fittingValues_.mse = -2;

    addPort(inPort_);
    addPort(outPort_);

    addProperty(ignoreFalseValues_);

    addProperty(expressionNameInput_);
    addProperty(expressionText_);
    addProperty(selfDescription_);
    addProperty(maxLength_);
    addProperty(recalculate_);

    addProperty(maxTableColumnWidth_);

    addProperty(keyColumnColor_);
    addProperty(dataColumnColor_);
    addProperty(fittingColumnColor_);

    ignoreFalseValues_.setGroupID("Fitting-Options");
    keyColumnColor_.setGroupID("Widget-Configuration");
    dataColumnColor_.setGroupID("Widget-Configuration");
    fittingColumnColor_.setGroupID("Widget-Configuration");
    maxTableColumnWidth_.setGroupID("Widget-Configuration");
    expressionNameInput_.setGroupID("Expression Representation");
    expressionText_.setGroupID("Expression Representation");
    selfDescription_.setGroupID("Expression Representation");
    maxLength_.setGroupID("Expression Representation");

    expressionText_.addOption("full","Full Expression Text",static_cast<int>(PlotFunction::FULL));
    expressionText_.addOption("length","Defined Length",static_cast<int>(PlotFunction::MAXLENGTH));
    expressionText_.addOption("self","Self Defined Expression Name",static_cast<int>(PlotFunction::CUSTOM));
    expressionText_.addOption("name","Expression Name + (..)",static_cast<int>(PlotFunction::ONLYNAME));
    expressionText_.selectByValue(0);
    expressionText_.onChange(MemberFunctionCallback<PlotDataFitFunction>(this, &PlotDataFitFunction::changeText));
    selfDescription_.setVisibleFlag(false);
    maxLength_.setVisibleFlag(false);

    setPropertyGroupGuiName("Fitting-Options","Fitting-Options");
    setPropertyGroupGuiName("Widget-Configuration","Widget-Configuration");
    setPropertyGroupGuiName("Expression Representation","Expression Representation");
}

PlotDataFitFunction::~PlotDataFitFunction() {
}

Processor* PlotDataFitFunction::create() const {
    return new PlotDataFitFunction();
}

void PlotDataFitFunction::initialize() {
    Processor::initialize();

    pDataOut_ = new PlotFunction(PlotExpression(""));

    if (inPort_.hasData()) {
        readData();
        calculate();
    }
    setOutPortData();
}

void PlotDataFitFunction::deinitialize() {
    outPort_.setData(0, false);
    pDataOut_->getPlotExpression().deletePlotExpressionNodes();
    delete pDataOut_;
    pData_ = 0;
    pDataOut_ = 0;

    Processor::deinitialize();
}

bool PlotDataFitFunction::isReady() const {
    return true;
}

bool PlotDataFitFunction::isEndProcessor() const {
    return (!outPort_.isConnected());
}

void PlotDataFitFunction::serialize(Serializer &s) const {
    Processor::serialize(s);
    s.serialize("Fittingvalues",fittingValues_);
}

void PlotDataFitFunction::deserialize(Deserializer &s) {
    Processor::deserialize(s);
    try {
        s.deserialize("Fittingvalues",fittingValues_);
    }
    catch (SerializationNoSuchDataException&) {
        // superfluous or missing XML data item, just ignore
        s.removeLastError();
    }

}

void PlotDataFitFunction::setOutPortData(){
    if (outPort_.isConnected()){
        if ((pDataOut_)){
            outPort_.setData(pDataOut_);
        }
        else
            outPort_.setData(0);
    }
}
void PlotDataFitFunction::readData() {
    if (dynamic_cast<const PlotData*>(inPort_.getData())) {
        pData_ = dynamic_cast<const PlotData*>(inPort_.getData());
    }
    else {
        pData_ = 0;
        LWARNING("PlotDataFitFunction can only handle with PlotData-Objects.");
    }
}

const PlotData* PlotDataFitFunction::getPlotData() const {
    return pData_;
}

const PlotFunction* PlotDataFitFunction::getPlotFunction() const {
    return pDataOut_;
}

const PlotDataFitFunction::FittingValues PlotDataFitFunction::getFittingValues() const {
    return fittingValues_;
}

void PlotDataFitFunction::process() {
    tgtAssert(isInitialized(), "not initialized");
    if (inPort_.hasData() && inPort_.hasChanged()) {
        readData();
        calculate();
    }
    else if (!inPort_.hasData()){
        pData_ = 0;
        PlotFunction* oldPlotData;
        oldPlotData = pDataOut_;
        pDataOut_ = new PlotFunction();
        pDataOut_->getPlotExpression().setExpressionName(expressionNameInput_.get());
        pDataOut_->setExpressionLength(static_cast<PlotFunction::ExpressionDescriptionLengthType>(expressionText_.getValue()),
            maxLength_.get(),selfDescription_.get());
        setOutPortData();
        if (getProcessorWidget()){
            getProcessorWidget()->updateFromProcessor();
        }
        oldPlotData->getPlotExpression().deletePlotExpressionNodes();
        delete oldPlotData;
    }
    else if (inPort_.hasData() && !inPort_.hasChanged()) {
        calculate();
    }
    else {
        if (getProcessorWidget()){
            getProcessorWidget()->updateFromProcessor();
        }
    }
}

void PlotDataFitFunction::updateView() {
    if (getProcessorWidget()){
        getProcessorWidget()->updateFromProcessor();
    }
}

void PlotDataFitFunction::changeText() {
    int value = expressionText_.getValue();
    if (value == static_cast<int>(PlotFunction::MAXLENGTH)) {
        selfDescription_.setVisibleFlag(false);
        maxLength_.setVisibleFlag(true);
    }
    else if (value == static_cast<int>(PlotFunction::CUSTOM)) {
        selfDescription_.setVisibleFlag(true);
        maxLength_.setVisibleFlag(false);
    }
    else {
        selfDescription_.setVisibleFlag(false);
        maxLength_.setVisibleFlag(false);
    }
}

void PlotDataFitFunction::calculate() {
    if (!pData_)
        return;
    if (fittingValues_.column >= 0) {
        fitData();
    }
    else {
        PlotFunction* newPlotData;
        newPlotData = new PlotFunction();
        newPlotData->getPlotExpression().setExpressionName(expressionNameInput_.get());
        newPlotData->setExpressionLength(static_cast<PlotFunction::ExpressionDescriptionLengthType>(expressionText_.getValue()),
            maxLength_.get(),selfDescription_.get());
        PlotFunction* oldData;
        oldData = pDataOut_;
        pDataOut_ = newPlotData;
        setOutPortData();
        if (getProcessorWidget()){
            getProcessorWidget()->updateFromProcessor();
        }
        oldData->getPlotExpression().deletePlotExpressionNodes();
        delete oldData;
    }
}

void PlotDataFitFunction::recalculate() {
    tgtAssert(isInitialized(), "not initialized");
    if (!pData_)
        return;
    if (fittingValues_.column >= 0) {
        PlotFunction* newPlotData;
        newPlotData = new PlotFunction(*pDataOut_);
        newPlotData->getPlotExpression().setExpressionName(expressionNameInput_.get());
        newPlotData->setExpressionLength(static_cast<PlotFunction::ExpressionDescriptionLengthType>(expressionText_.getValue()),
            maxLength_.get(),selfDescription_.get());
        PlotFunction* oldData;
        oldData = pDataOut_;
        pDataOut_ = newPlotData;
        setOutPortData();
        if (getProcessorWidget()){
            getProcessorWidget()->updateFromProcessor();
        }
        oldData->getPlotExpression().deletePlotExpressionNodes();
        delete oldData;
    }
}


void PlotDataFitFunction::fitData() {
    if (fittingValues_.regressionType == FunctionLibrary::LINEAR)
        linearRegression(fittingValues_.column);
    else if (fittingValues_.regressionType == FunctionLibrary::SQRT)
        sqrtRegression(fittingValues_.column);
    else if (fittingValues_.regressionType == FunctionLibrary::SIMPLESPLINE)
        simpleSplineRegression(fittingValues_.column);
    else if (fittingValues_.regressionType == FunctionLibrary::CONSTANTSPLINE)
        constantSpline(fittingValues_.column);
    else if (fittingValues_.regressionType == FunctionLibrary::BSPLINE)
        bSplineRegression(fittingValues_.column);
    else if (fittingValues_.regressionType == FunctionLibrary::CUBICSPLINE)
        cubicSplineRegression(fittingValues_.column);
    else if (fittingValues_.regressionType == FunctionLibrary::SQUARE)
        squareRegression(fittingValues_.column);
    else if (fittingValues_.regressionType == FunctionLibrary::CUBIC)
        cubicRegression(fittingValues_.column);
    else if (fittingValues_.regressionType == FunctionLibrary::POLYNOMIAL)
        multiRegression(fittingValues_.dimension,fittingValues_.column);
    else if (fittingValues_.regressionType == FunctionLibrary::LOGARITHMIC)
        logarithmicRegression(fittingValues_.column);
    else if (fittingValues_.regressionType == FunctionLibrary::POWER)
        powerRegression(fittingValues_.column);
    else if (fittingValues_.regressionType == FunctionLibrary::EXPONENTIAL)
        exponentialRegression(fittingValues_.column);
    else if (fittingValues_.regressionType == FunctionLibrary::SIN)
        sinRegression(fittingValues_.column);
    else if (fittingValues_.regressionType == FunctionLibrary::COS)
        cosRegression(fittingValues_.column);
    else if (fittingValues_.regressionType == FunctionLibrary::INTERPOLATION)
        interpolRegression(fittingValues_.column);
    else {
        if (getProcessorWidget()){
            getProcessorWidget()->updateFromProcessor();
        }
    }
}

void PlotDataFitFunction::linearRegression(int column) {
    std::vector<std::pair<plot_t,plot_t> > points;
    if (pData_->getColumnCount() > column && pData_->getColumnType(column) == PlotBase::NUMBER&& pData_->getColumnType(0) == PlotBase::NUMBER) {
        PlotFunction* oldfunction = pDataOut_;
        bool ignoreFalseValues = ignoreFalseValues_.get();
        plot_t x;
        plot_t y;
        for (int i = 0; i < pData_->getRowsCount(); ++i) {
            PlotRowValue row = pData_->getRow(i);
            x = row.getCellAt(0).getValue();
            y = row.getCellAt(column).getValue();
            if (!ignoreFalseValues || (x == x && y == y))
                points.push_back(std::pair<plot_t,plot_t>(x,y));
        }
        std::pair<PlotExpression,plot_t> result = FunctionLibrary::fittingLinear(points);
        pDataOut_ = new PlotFunction(result.first);
        pDataOut_->getPlotExpression().setExpressionName(expressionNameInput_.get());
        pDataOut_->setExpressionLength(static_cast<PlotFunction::ExpressionDescriptionLengthType>(expressionText_.getValue()),
            maxLength_.get(),selfDescription_.get());
        fittingValues_.mse = result.second;
        fittingValues_.column = column;
        fittingValues_.regressionType = FunctionLibrary::LINEAR;
        setOutPortData();
        oldfunction->getPlotExpression().deletePlotExpressionNodes();
        delete oldfunction;
    }
    if (getProcessorWidget()){
        getProcessorWidget()->updateFromProcessor();
    }
}

void PlotDataFitFunction::sqrtRegression(int column) {
    std::vector<std::pair<plot_t,plot_t> > points;
    if (pData_->getColumnCount() > column && pData_->getColumnType(column) == PlotBase::NUMBER&& pData_->getColumnType(0) == PlotBase::NUMBER) {
        PlotFunction* oldfunction = pDataOut_;
        bool ignoreFalseValues = ignoreFalseValues_.get();
        plot_t x;
        plot_t y;
        for (int i = 0; i < pData_->getRowsCount(); ++i) {
            PlotRowValue row = pData_->getRow(i);
            x = row.getCellAt(0).getValue();
            y = row.getCellAt(column).getValue();
            if (!ignoreFalseValues || (x == x && y == y && x >= 0))
                points.push_back(std::pair<plot_t,plot_t>(x,y));
        }
        std::pair<PlotExpression,plot_t> result = FunctionLibrary::fittingSqrt(points);
        pDataOut_ = new PlotFunction(result.first);
        pDataOut_->getPlotExpression().setExpressionName(expressionNameInput_.get());
        pDataOut_->setExpressionLength(static_cast<PlotFunction::ExpressionDescriptionLengthType>(expressionText_.getValue()),
            maxLength_.get(),selfDescription_.get());
        fittingValues_.mse = result.second;
        fittingValues_.column = column;
        fittingValues_.regressionType = FunctionLibrary::SQRT;
        setOutPortData();
        oldfunction->getPlotExpression().deletePlotExpressionNodes();
        delete oldfunction;
    }
    if (getProcessorWidget()){
        getProcessorWidget()->updateFromProcessor();
    }
}

void PlotDataFitFunction::squareRegression(int column) {
    std::vector<std::pair<plot_t,plot_t> > points;
    if (pData_->getColumnCount() > column && pData_->getColumnType(column) == PlotBase::NUMBER&& pData_->getColumnType(0) == PlotBase::NUMBER) {
        PlotFunction* oldfunction = pDataOut_;
        bool ignoreFalseValues = ignoreFalseValues_.get();
        plot_t x;
        plot_t y;
        for (int i = 0; i < pData_->getRowsCount(); ++i) {
            PlotRowValue row = pData_->getRow(i);
            x = row.getCellAt(0).getValue();
            y = row.getCellAt(column).getValue();
            if (!ignoreFalseValues || (x == x && y == y))
                points.push_back(std::pair<plot_t,plot_t>(x,y));
        }
        std::pair<PlotExpression,plot_t> result = FunctionLibrary::fittingSquare(points);
        pDataOut_ = new PlotFunction(result.first);
        pDataOut_->getPlotExpression().setExpressionName(expressionNameInput_.get());
        pDataOut_->setExpressionLength(static_cast<PlotFunction::ExpressionDescriptionLengthType>(expressionText_.getValue()),
            maxLength_.get(),selfDescription_.get());
        fittingValues_.mse = result.second;
        fittingValues_.column = column;
        fittingValues_.regressionType = FunctionLibrary::SQUARE;
        setOutPortData();
        oldfunction->getPlotExpression().deletePlotExpressionNodes();
        delete oldfunction;
    }
    if (getProcessorWidget()){
        getProcessorWidget()->updateFromProcessor();
    }
}

void PlotDataFitFunction::cubicRegression(int column) {
    std::vector<std::pair<plot_t,plot_t> > points;
    if (pData_->getColumnCount() > column && pData_->getColumnType(column) == PlotBase::NUMBER&& pData_->getColumnType(0) == PlotBase::NUMBER) {
        PlotFunction* oldfunction = pDataOut_;
        bool ignoreFalseValues = ignoreFalseValues_.get();
        plot_t x;
        plot_t y;
        for (int i = 0; i < pData_->getRowsCount(); ++i) {
            PlotRowValue row = pData_->getRow(i);
            x = row.getCellAt(0).getValue();
            y = row.getCellAt(column).getValue();
            if (!ignoreFalseValues || (x == x && y == y))
                points.push_back(std::pair<plot_t,plot_t>(x,y));
        }
        std::pair<PlotExpression,plot_t> result = FunctionLibrary::fittingCubic(points);
        pDataOut_ = new PlotFunction(result.first);
        pDataOut_->getPlotExpression().setExpressionName(expressionNameInput_.get());
        pDataOut_->setExpressionLength(static_cast<PlotFunction::ExpressionDescriptionLengthType>(expressionText_.getValue()),
            maxLength_.get(),selfDescription_.get());
        fittingValues_.mse = result.second;
        fittingValues_.column = column;
        fittingValues_.regressionType = FunctionLibrary::CUBIC;
        setOutPortData();
        oldfunction->getPlotExpression().deletePlotExpressionNodes();
        delete oldfunction;
    }
    if (getProcessorWidget())
        getProcessorWidget()->updateFromProcessor();
}

void PlotDataFitFunction::multiRegression(int dimension,int column) {
    std::vector<std::pair<plot_t,plot_t> > points;
    if (pData_->getColumnCount() > column && pData_->getColumnType(column) == PlotBase::NUMBER&& pData_->getColumnType(0) == PlotBase::NUMBER) {
        PlotFunction* oldfunction = pDataOut_;
        bool ignoreFalseValues = ignoreFalseValues_.get();
        plot_t x;
        plot_t y;
        for (int i = 0; i < pData_->getRowsCount(); ++i) {
            PlotRowValue row = pData_->getRow(i);
            x = row.getCellAt(0).getValue();
            y = row.getCellAt(column).getValue();
            if (!ignoreFalseValues || (x == x && y == y))
                points.push_back(std::pair<plot_t,plot_t>(x,y));
        }
        std::pair<PlotExpression,plot_t> result = FunctionLibrary::fittingMulti(dimension,points);
        pDataOut_ = new PlotFunction(result.first);
        pDataOut_->getPlotExpression().setExpressionName(expressionNameInput_.get());
        pDataOut_->setExpressionLength(static_cast<PlotFunction::ExpressionDescriptionLengthType>(expressionText_.getValue()),
            maxLength_.get(),selfDescription_.get());
        fittingValues_.mse = result.second;
        fittingValues_.column = column;
        fittingValues_.dimension = dimension;
        if (dimension == 2)
            fittingValues_.regressionType = FunctionLibrary::SQUARE;
        else if (dimension == 2)
            fittingValues_.regressionType = FunctionLibrary::CUBIC;
        else
            fittingValues_.regressionType = FunctionLibrary::POLYNOMIAL;
        setOutPortData();
        oldfunction->getPlotExpression().deletePlotExpressionNodes();
        delete oldfunction;
    }
    if (getProcessorWidget())
        getProcessorWidget()->updateFromProcessor();
}

void PlotDataFitFunction::constantSpline(int column) {
    std::vector<std::pair<plot_t,plot_t> > points;
    if (pData_->getColumnCount() > column && pData_->getColumnType(column) == PlotBase::NUMBER&& pData_->getColumnType(0) == PlotBase::NUMBER) {
        PlotFunction* oldfunction = pDataOut_;
        bool ignoreFalseValues = ignoreFalseValues_.get();
        plot_t x;
        plot_t y;
        for (int i = 0; i < pData_->getRowsCount(); ++i) {
            PlotRowValue row = pData_->getRow(i);
            x = row.getCellAt(0).getValue();
            y = row.getCellAt(column).getValue();
            if (!ignoreFalseValues || (x == x && y == y))
                points.push_back(std::pair<plot_t,plot_t>(x,y));
        }
        std::pair<PlotExpression,plot_t> result = FunctionLibrary::constantSpline(points);
        pDataOut_ = new PlotFunction(result.first);
        pDataOut_->getPlotExpression().setExpressionName(expressionNameInput_.get());
        pDataOut_->setExpressionLength(static_cast<PlotFunction::ExpressionDescriptionLengthType>(expressionText_.getValue()),
            maxLength_.get(),selfDescription_.get());
        fittingValues_.mse = result.second;
        fittingValues_.column = column;
        fittingValues_.regressionType = FunctionLibrary::CONSTANTSPLINE;
        setOutPortData();
        oldfunction->getPlotExpression().deletePlotExpressionNodes();
        delete oldfunction;
    }
    if (getProcessorWidget())
       getProcessorWidget()->updateFromProcessor();
}

void PlotDataFitFunction::simpleSplineRegression(int column) {
    std::vector<std::pair<plot_t,plot_t> > points;
    if (pData_->getColumnCount() > column && pData_->getColumnType(column) == PlotBase::NUMBER&& pData_->getColumnType(0) == PlotBase::NUMBER) {
        PlotFunction* oldfunction = pDataOut_;
        bool ignoreFalseValues = ignoreFalseValues_.get();
        plot_t x;
        plot_t y;
        for (int i = 0; i < pData_->getRowsCount(); ++i) {
            PlotRowValue row = pData_->getRow(i);
            x = row.getCellAt(0).getValue();
            y = row.getCellAt(column).getValue();
            if (!ignoreFalseValues || (x == x && y == y))
                points.push_back(std::pair<plot_t,plot_t>(x,y));
        }
        std::pair<PlotExpression,plot_t> result = FunctionLibrary::simpleSpline(points);
        pDataOut_ = new PlotFunction(result.first);
        pDataOut_->getPlotExpression().setExpressionName(expressionNameInput_.get());
        pDataOut_->setExpressionLength(static_cast<PlotFunction::ExpressionDescriptionLengthType>(expressionText_.getValue()),
            maxLength_.get(),selfDescription_.get());
        fittingValues_.mse = result.second;
        fittingValues_.column = column;
        fittingValues_.regressionType = FunctionLibrary::SIMPLESPLINE;
        setOutPortData();
        oldfunction->getPlotExpression().deletePlotExpressionNodes();
        delete oldfunction;
    }
    if (getProcessorWidget())
       getProcessorWidget()->updateFromProcessor();
}

void PlotDataFitFunction::bSplineRegression(int column) {
    std::vector<std::pair<plot_t,plot_t> > points;
    if (pData_->getColumnCount() > column && pData_->getColumnType(column) == PlotBase::NUMBER&& pData_->getColumnType(0) == PlotBase::NUMBER) {
        PlotFunction* oldfunction = pDataOut_;
        bool ignoreFalseValues = ignoreFalseValues_.get();
        plot_t x;
        plot_t y;
        for (int i = 0; i < pData_->getRowsCount(); ++i) {
            PlotRowValue row = pData_->getRow(i);
            x = row.getCellAt(0).getValue();
            y = row.getCellAt(column).getValue();
            if (!ignoreFalseValues || (x == x && y == y))
                points.push_back(std::pair<plot_t,plot_t>(x,y));
        }
        std::pair<PlotExpression,plot_t> result = FunctionLibrary::bSpline(points);
        pDataOut_ = new PlotFunction(result.first);
        pDataOut_->getPlotExpression().setExpressionName(expressionNameInput_.get());
        pDataOut_->setExpressionLength(static_cast<PlotFunction::ExpressionDescriptionLengthType>(expressionText_.getValue()),
            maxLength_.get(),selfDescription_.get());
        fittingValues_.mse = result.second;
        fittingValues_.column = column;
        fittingValues_.regressionType = FunctionLibrary::BSPLINE;
        setOutPortData();
        oldfunction->getPlotExpression().deletePlotExpressionNodes();
        delete oldfunction;
    }
    if (getProcessorWidget())
        getProcessorWidget()->updateFromProcessor();
}

void PlotDataFitFunction::squareSplineRegression(int column) {
    std::vector<std::pair<plot_t,plot_t> > points;
    if (pData_->getColumnCount() > column && pData_->getColumnType(column) == PlotBase::NUMBER&& pData_->getColumnType(0) == PlotBase::NUMBER) {
        PlotFunction* oldfunction = pDataOut_;
        bool ignoreFalseValues = ignoreFalseValues_.get();
        plot_t x;
        plot_t y;
        for (int i = 0; i < pData_->getRowsCount(); ++i) {
            PlotRowValue row = pData_->getRow(i);
            x = row.getCellAt(0).getValue();
            y = row.getCellAt(column).getValue();
            if (!ignoreFalseValues || (x == x && y == y))
                points.push_back(std::pair<plot_t,plot_t>(x,y));
        }
        std::pair<PlotExpression,plot_t> result = FunctionLibrary::squareSpline(points);
        pDataOut_ = new PlotFunction(result.first);
        pDataOut_->getPlotExpression().setExpressionName(expressionNameInput_.get());
        pDataOut_->setExpressionLength(static_cast<PlotFunction::ExpressionDescriptionLengthType>(expressionText_.getValue()),
            maxLength_.get(),selfDescription_.get());
        fittingValues_.mse = result.second;
        fittingValues_.column = column;
        fittingValues_.regressionType = FunctionLibrary::SQUARESPLINE;
        setOutPortData();
        oldfunction->getPlotExpression().deletePlotExpressionNodes();
        delete oldfunction;
    }
    if (getProcessorWidget())
        getProcessorWidget()->updateFromProcessor();
}

void PlotDataFitFunction::cubicSplineRegression(int column) {
    std::vector<std::pair<plot_t,plot_t> > points;
    if (pData_->getColumnCount() > column && pData_->getColumnType(column) == PlotBase::NUMBER&& pData_->getColumnType(0) == PlotBase::NUMBER) {
        PlotFunction* oldfunction = pDataOut_;
        bool ignoreFalseValues = ignoreFalseValues_.get();
        plot_t x;
        plot_t y;
        for (int i = 0; i < pData_->getRowsCount(); ++i) {
            PlotRowValue row = pData_->getRow(i);
            x = row.getCellAt(0).getValue();
            y = row.getCellAt(column).getValue();
            if (!ignoreFalseValues || (x == x && y == y))
                points.push_back(std::pair<plot_t,plot_t>(x,y));
        }
        std::pair<PlotExpression,plot_t> result = FunctionLibrary::cubicSpline(points);
        pDataOut_ = new PlotFunction(result.first);
        pDataOut_->getPlotExpression().setExpressionName(expressionNameInput_.get());
        pDataOut_->setExpressionLength(static_cast<PlotFunction::ExpressionDescriptionLengthType>(expressionText_.getValue()),
            maxLength_.get(),selfDescription_.get());
        fittingValues_.mse = result.second;
        fittingValues_.column = column;
        fittingValues_.regressionType = FunctionLibrary::CUBICSPLINE;
        setOutPortData();
        oldfunction->getPlotExpression().deletePlotExpressionNodes();
        delete oldfunction;
    }
    if (getProcessorWidget())
        getProcessorWidget()->updateFromProcessor();
}

void PlotDataFitFunction::logarithmicRegression(int column) {
    std::vector<std::pair<plot_t,plot_t> > points;
    if (pData_->getColumnCount() > column && pData_->getColumnType(column) == PlotBase::NUMBER&& pData_->getColumnType(0) == PlotBase::NUMBER) {
        PlotFunction* oldfunction = pDataOut_;
        bool ignoreFalseValues = ignoreFalseValues_.get();
        plot_t x;
        plot_t y;
        for (int i = 0; i < pData_->getRowsCount(); ++i) {
            PlotRowValue row = pData_->getRow(i);
            x = row.getCellAt(0).getValue();
            y = row.getCellAt(column).getValue();
            if (!ignoreFalseValues || (x == x && y == y && x > 0))
                points.push_back(std::pair<plot_t,plot_t>(x,y));
        }
        std::pair<PlotExpression,plot_t> result = FunctionLibrary::fittingLogarithmus(points);
        pDataOut_ = new PlotFunction(result.first);
        pDataOut_->getPlotExpression().setExpressionName(expressionNameInput_.get());
        pDataOut_->setExpressionLength(static_cast<PlotFunction::ExpressionDescriptionLengthType>(expressionText_.getValue()),
            maxLength_.get(),selfDescription_.get());
        fittingValues_.mse = result.second;
        fittingValues_.column = column;
        fittingValues_.regressionType = FunctionLibrary::LOGARITHMIC;
        setOutPortData();
        oldfunction->getPlotExpression().deletePlotExpressionNodes();
        delete oldfunction;
    }
    if (getProcessorWidget())
        getProcessorWidget()->updateFromProcessor();
}

void PlotDataFitFunction::powerRegression(int column) {
    std::vector<std::pair<plot_t,plot_t> > points;
    if (pData_->getColumnCount() > column && pData_->getColumnType(column) == PlotBase::NUMBER&& pData_->getColumnType(0) == PlotBase::NUMBER) {
        PlotFunction* oldfunction = pDataOut_;
        bool ignoreFalseValues = ignoreFalseValues_.get();
        plot_t x;
        plot_t y;
        for (int i = 0; i < pData_->getRowsCount(); ++i) {
            PlotRowValue row = pData_->getRow(i);
            x = row.getCellAt(0).getValue();
            y = row.getCellAt(column).getValue();
            if (!ignoreFalseValues || (x == x && y == y && x > 0 && y > 0))
                points.push_back(std::pair<plot_t,plot_t>(x,y));
        }
        std::pair<PlotExpression,plot_t> result = FunctionLibrary::fittingPower(points);
        pDataOut_ = new PlotFunction(result.first);
        pDataOut_->getPlotExpression().setExpressionName(expressionNameInput_.get());
        pDataOut_->setExpressionLength(static_cast<PlotFunction::ExpressionDescriptionLengthType>(expressionText_.getValue()),
            maxLength_.get(),selfDescription_.get());
        fittingValues_.mse = result.second;
        fittingValues_.column = column;
        fittingValues_.regressionType = FunctionLibrary::POWER;
        setOutPortData();
        oldfunction->getPlotExpression().deletePlotExpressionNodes();
        delete oldfunction;
    }
    if (getProcessorWidget())
        getProcessorWidget()->updateFromProcessor();
}

void PlotDataFitFunction::exponentialRegression(int column) {
    std::vector<std::pair<plot_t,plot_t> > points;
    if (pData_->getColumnCount() > column && pData_->getColumnType(column) == PlotBase::NUMBER&& pData_->getColumnType(0) == PlotBase::NUMBER) {
        PlotFunction* oldfunction = pDataOut_;
        bool ignoreFalseValues = ignoreFalseValues_.get();
        plot_t x;
        plot_t y;
        for (int i = 0; i < pData_->getRowsCount(); ++i) {
            PlotRowValue row = pData_->getRow(i);
            x = row.getCellAt(0).getValue();
            y = row.getCellAt(column).getValue();
            if (!ignoreFalseValues || (x == x && y == y && y > 0))
                points.push_back(std::pair<plot_t,plot_t>(x,y));
        }
        std::pair<PlotExpression,plot_t> result = FunctionLibrary::fittingExponential(points);
        pDataOut_ = new PlotFunction(result.first);
        pDataOut_->getPlotExpression().setExpressionName(expressionNameInput_.get());
        pDataOut_->setExpressionLength(static_cast<PlotFunction::ExpressionDescriptionLengthType>(expressionText_.getValue()),
            maxLength_.get(),selfDescription_.get());
        fittingValues_.mse = result.second;
        fittingValues_.column = column;
        fittingValues_.regressionType = FunctionLibrary::EXPONENTIAL;
        setOutPortData();
        oldfunction->getPlotExpression().deletePlotExpressionNodes();
        delete oldfunction;
    }
    if (getProcessorWidget())
        getProcessorWidget()->updateFromProcessor();
}

void PlotDataFitFunction::sinRegression(int column) {
    std::vector<std::pair<plot_t,plot_t> > points;
    if (pData_->getColumnCount() > column && pData_->getColumnType(column) == PlotBase::NUMBER&& pData_->getColumnType(0) == PlotBase::NUMBER) {
        PlotFunction* oldfunction = pDataOut_;
        bool ignoreFalseValues = ignoreFalseValues_.get();
        plot_t x;
        plot_t y;
        for (int i = 0; i < pData_->getRowsCount(); ++i) {
            PlotRowValue row = pData_->getRow(i);
            x = row.getCellAt(0).getValue();
            y = row.getCellAt(column).getValue();
            if (!ignoreFalseValues || (x == x && y == y))
                points.push_back(std::pair<plot_t,plot_t>(x,y));
        }
        std::pair<PlotExpression,plot_t> result = FunctionLibrary::fittingSin(points);
        pDataOut_ = new PlotFunction(result.first);
        pDataOut_->getPlotExpression().setExpressionName(expressionNameInput_.get());
        pDataOut_->setExpressionLength(static_cast<PlotFunction::ExpressionDescriptionLengthType>(expressionText_.getValue()),
            maxLength_.get(),selfDescription_.get());
        fittingValues_.mse = result.second;
        fittingValues_.column = column;
        fittingValues_.regressionType = FunctionLibrary::SIN;
        setOutPortData();
        oldfunction->getPlotExpression().deletePlotExpressionNodes();
        delete oldfunction;
    }
    if (getProcessorWidget())
        getProcessorWidget()->updateFromProcessor();
}

void PlotDataFitFunction::cosRegression(int column) {
    std::vector<std::pair<plot_t,plot_t> > points;
    if (pData_->getColumnCount() > column && pData_->getColumnType(column) == PlotBase::NUMBER&& pData_->getColumnType(0) == PlotBase::NUMBER) {
        PlotFunction* oldfunction = pDataOut_;
        bool ignoreFalseValues = ignoreFalseValues_.get();
        plot_t x;
        plot_t y;
        for (int i = 0; i < pData_->getRowsCount(); ++i) {
            PlotRowValue row = pData_->getRow(i);
            x = row.getCellAt(0).getValue();
            y = row.getCellAt(column).getValue();
            if (!ignoreFalseValues || (x == x && y == y))
                points.push_back(std::pair<plot_t,plot_t>(x,y));
        }
        std::pair<PlotExpression,plot_t> result = FunctionLibrary::fittingCos(points);
        pDataOut_ = new PlotFunction(result.first);
        pDataOut_->getPlotExpression().setExpressionName(expressionNameInput_.get());
        pDataOut_->setExpressionLength(static_cast<PlotFunction::ExpressionDescriptionLengthType>(expressionText_.getValue()),
            maxLength_.get(),selfDescription_.get());
        fittingValues_.mse = result.second;
        fittingValues_.column = column;
        fittingValues_.regressionType = FunctionLibrary::COS;
        setOutPortData();
        oldfunction->getPlotExpression().deletePlotExpressionNodes();
        delete oldfunction;
    }
    if (getProcessorWidget())
        getProcessorWidget()->updateFromProcessor();
 }

void PlotDataFitFunction::interpolRegression(int column) {
    std::vector<std::pair<plot_t,plot_t> > points;
    if (pData_->getColumnCount() > column && pData_->getColumnType(column) == PlotBase::NUMBER&& pData_->getColumnType(0) == PlotBase::NUMBER) {
        PlotFunction* oldfunction = pDataOut_;
        bool ignoreFalseValues = ignoreFalseValues_.get();
        plot_t x;
        plot_t y;
        for (int i = 0; i < pData_->getRowsCount(); ++i) {
            PlotRowValue row = pData_->getRow(i);
            x = row.getCellAt(0).getValue();
            y = row.getCellAt(column).getValue();
            if (!ignoreFalseValues || (x == x && y == y))
                points.push_back(std::pair<plot_t,plot_t>(x,y));
        }
        std::pair<PlotExpression,plot_t> result = FunctionLibrary::fittingInterpol(points);
        pDataOut_ = new PlotFunction(result.first);
        pDataOut_->getPlotExpression().setExpressionName(expressionNameInput_.get());
        pDataOut_->setExpressionLength(static_cast<PlotFunction::ExpressionDescriptionLengthType>(expressionText_.getValue()),
            maxLength_.get(),selfDescription_.get());
        fittingValues_.mse = result.second;
        fittingValues_.column = column;
        fittingValues_.regressionType = FunctionLibrary::INTERPOLATION;
        setOutPortData();
        oldfunction->getPlotExpression().deletePlotExpressionNodes();
        delete oldfunction;
    }
    if (getProcessorWidget())
        getProcessorWidget()->updateFromProcessor();
}

void PlotDataFitFunction::ignoreFalseValues() {
    if (ignoreFalseValues_.get())
        ignoreFalseValues_.set(false);
    else
        ignoreFalseValues_.set(true);
}

void PlotDataFitFunction::FittingValues::serialize(Serializer& s) const {
    s.serialize("regressionType",regressionType);
    s.serialize("column",column);
    s.serialize("Dimension",dimension);
}

void PlotDataFitFunction::FittingValues::deserialize(Deserializer& s){
    int value;
    s.deserialize("regressionType",value);
    regressionType = static_cast<FunctionLibrary::RegressionType>(value);
    s.deserialize("column",column);
    mse = -1;
    s.deserialize("Dimension",dimension);
}



//plot_t PlotDataFitFunction::linear(const std::vector<plot_t>& values) {
//    if (values.size() <= 0)
//        return 0;
//
//    plot_t sum = 0;
//    for (size_t i = 0; i < values.size(); ++i) {
//        if (values[i] == values[i])
//            sum += values[i];
//    }
//    int xMid = sum/values.size();
//    return xMid;
        //unsigned int i = 1;
        //while ((m_x[i] < input) && (i < count))
        //    i++;

        //double dif1 = m_x[i] - m_x[i-1];  //!=0 per definition
        //double dif2 = input - m_x[i-1];

        //return (m_y[i-1] + ((m_y[i] - m_y[i-1])*dif2 / dif1));
//}




//FittingFunctionLinear* FittingFunctionLinear::clone() const {
//    return new FittingFunctionLinear();
//}
//
//std::string FittingFunctionLinear::toString() const {
//    return "Linear";
//}
//plot_t PlotDataFitFunction::square(const std::vector<plot_t>& values) {
//    if (values.size() <= 0)
//        return 0;
//
//    plot_t sum = 0;
//    for (size_t i = 0; i < values.size(); ++i) {
//        if (values[i] == values[i])
//            sum += values[i];
//    }
//    return sum/values.size();
        //double a, b, diff;
        //unsigned int j = 1;

        //while ((m_x[j] < input) && (j < count))
        //    j++;

        //diff = m_x[j] - m_x[j-1];

        //a = (m_x[j] - input) / diff;    //div should not be 0
        //b = (input - m_x[j-1]) / diff;

        //return (a*m_y[j-1] + b*m_y[j] + ((a*a*a - a)*m_der[j - 1] +
        //       (b*b*b - b)*m_der[j])*(diff*diff) / 6);


//}
//
//FittingFunctionSquare* FittingFunctionSquare::clone() const {
//    return new FittingFunctionSquare();
//}
//
//std::string FittingFunctionSquare::toString() const {
//    return "Square";
//}
//plot_t FittingFunctionCubic::evaluate(std::vector<plot_t>& values) const {
//float cubic_interpolate( float y0, float y1, float y2, float y3, float mu ) {
//
//   float a0, a1, a2, a3, mu2;
//
//   mu2 = mu*mu;
//   a0 = y3 - y2 - y0 + y1; //p
//   a1 = y0 - y1 - a0;
//   a2 = y2 - y0;
//   a3 = y1;
//
//   return ( a0*mu*mu2 + a1*mu2 + a2*mu + a3 );
//}
//    if (values.size() <= 0)
//        return 0;
//
//    plot_t sum = 0;
//    for (size_t i = 0; i < values.size(); ++i) {
//        if (values[i] == values[i])
//            sum += values[i];
//    }
//    return sum/values.size();
//}
//
//FittingFunctionCubic* FittingFunctionCubic::clone() const {
//    return new FittingFunctionCubic();
//}
//
//std::string FittingFunctionLinear::toString() const {
//    return "Cubic";
//}



/*
plot_t PlotFunction::evaluateAt(plot_t* value) {
    #define PI  3.1415
    #define N   12

int main()
{
    // Declare and initialize two arrays to hold the coordinates of the initial data points
    double x[N], y[N];

    // Generate the points
    double xx = PI, step = 4 * PI / (N - 1);
    for (int i = 0; i < N; ++i, xx += step) {
        x[i] = xx;
        y[i] = sin(2 * xx) / xx;
    }

    // Initialize the linear interpolation routine with known data points
    Maths::Interpolation::Linear A(N, x, y);

    // Interrogate linear fitting curve to find interpolated values
    int N_out = 20;
    xx = PI, step = (3 * PI) / (N_out - 1);
    for (int i = 0; i < N_out; ++i, xx += step) {
        cout << "x = " << setw(7) << xx << "  y = ";
        cout << setw(13) << A.getValue(xx) << endl;
  }
    return 0;

    return expr_.evaluateAt(value);
}


PlotData* PlotFunction::select(Interval<plot_t>* interval, plot_t* step, int count) {

    //check that count equals the number of variables in expr_
    if (count != expr_.numberOfVariables())
        throw "Wrong number of Variables.";

    //number of key columns equals count, one data column
    PlotData* data = new PlotData(count, 1);

    //the function is of type R->R
    if (count == 1) {
        for (plot_t x = (interval[0].getLeftOpen() ? interval[0].getLeft() + step[0] : interval[0].getLeft());
            interval[0].contains(x) ; x+= step[0]) {
            plot_t value[] = {x};
            plot_t row[] = {x, expr_.evaluateAt((value))};
            data->insert(row, 2);
        }
    }
    //the function is of typre RxR->R
    else if (count == 2) {
        for (plot_t x = (interval[0].getLeftOpen() ? interval[0].getLeft() + step[0] : interval[0].getLeft());
            (x < interval[0].getRight()) || (x == interval[0].getRight() && interval[0].getRightOpen()) ; x+= step[0]) {
            for (plot_t y = (interval[1].getLeftOpen() ? interval[1].getLeft() + step[1] : interval[1].getLeft());
                (y < interval[1].getRight()) || (y == interval[1].getRight() && interval[1].getRightOpen()) ; y+= step[1]) {
            plot_t value[] = {x, y};
            plot_t row[] = {x, y, expr_.evaluateAt((value))};
            data->insert(row, 3);
            }
        }
    }
    return data;

}
*/
}

