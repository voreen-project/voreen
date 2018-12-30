/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany,                        *
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

#include "plotfunctionsource.h"

#include "../datastructures/plotfunction.h"
#include "../datastructures/plotexpression.h"


namespace voreen {

const std::string PlotFunctionSource::loggerCat_("voreen.PlotFunctionSource");

PlotFunctionSource::PlotFunctionSource():
    Processor()
    , outPort_(Port::OUTPORT,"PlotFunctionOutPort")
    , functionInput_("functionInput_","Function Input","x",Processor::VALID)
    , expressionNameInput_("expressionNameInput_","Expression-Name","f",Processor::VALID)
    , expressionText_("Expression Format","Expression Format",Processor::VALID)
    , selfDescription_("selfDescription","Custom Text","",Processor::VALID)
    , maxLength_("maxLength","Max Expression Length",200,1,1000,Processor::VALID)
    , pData_(0)
    , recalculate_("Calculate","Calculate")
{
    recalculate_.onChange(MemberFunctionCallback<PlotFunctionSource>(this, &PlotFunctionSource::recalculate));
    addProperty(&functionInput_);
    addProperty(&expressionNameInput_);
    addProperty(&expressionText_);
    addProperty(&selfDescription_);
    addProperty(&maxLength_);
    addProperty(&recalculate_);
    addPort(&outPort_);

    expressionText_.addOption("full","Full Expression Text",static_cast<int>(PlotFunction::FULL));
    expressionText_.addOption("length","Defined Length",static_cast<int>(PlotFunction::MAXLENGTH));
    expressionText_.addOption("self","Self Defined Expression Name",static_cast<int>(PlotFunction::CUSTOM));
    expressionText_.addOption("name","Expression Name + (..)",static_cast<int>(PlotFunction::ONLYNAME));
    expressionText_.selectByValue(0);
    expressionText_.onChange(MemberFunctionCallback<PlotFunctionSource>(this, &PlotFunctionSource::changeText));
    selfDescription_.setVisibleFlag(false);
    maxLength_.setVisibleFlag(false);

    expressionNameInput_.setGroupID("Expression Representation");
    expressionText_.setGroupID("Expression Representation");
    selfDescription_.setGroupID("Expression Representation");
    maxLength_.setGroupID("Expression Representation");
    setPropertyGroupGuiName("Expression Representation","Expression Representation");
}

PlotFunctionSource::~PlotFunctionSource() {
    // nothing to do here
}

Processor* PlotFunctionSource::create() const {
    return new PlotFunctionSource();
}

void PlotFunctionSource::process() {
    recalculate();
}

void PlotFunctionSource::initialize() {
    Processor::initialize();

    pData_ = new PlotFunction();
}

void PlotFunctionSource::deinitialize() {
    outPort_.setData(0, false);
    pData_->getPlotExpression().deletePlotExpressionNodes();
    delete pData_;
    pData_ = 0;

    Processor::deinitialize();
}

void PlotFunctionSource::recalculate() {
    if (pData_->getPlotExpression().getExpressionString() != functionInput_.get()) {
        PlotFunction* newData = readInputString();
        if (newData->getPlotExpression().getExpressionString().size() > 0)
            functionInput_.set(newData->getPlotExpression().getExpressionString());
        newData->setExpressionLength(static_cast<PlotFunction::ExpressionDescriptionLengthType>(expressionText_.getValue()),
            maxLength_.get(),selfDescription_.get());
        PlotFunction* oldData = pData_;
        pData_ = newData;
        outPort_.setData(pData_);
        oldData->getPlotExpression().deletePlotExpressionNodes();
        delete oldData;
    }
    else if (pData_->getPlotExpression().getExpressionString() != expressionNameInput_.get()) {
        pData_->getPlotExpression().setExpressionName(expressionNameInput_.get());
        PlotFunction* oldData = pData_;
        pData_ = new PlotFunction(*oldData);
        pData_->setExpressionLength(static_cast<PlotFunction::ExpressionDescriptionLengthType>(expressionText_.getValue()),
            maxLength_.get(),selfDescription_.get());
        outPort_.setData(pData_);
        //oldData->getPlotExpression().deletePlotExpressionNodes();
        delete oldData;
    }
    else {
        PlotFunction* oldData = pData_;
        pData_ = new PlotFunction(*oldData);
        pData_->setExpressionLength(static_cast<PlotFunction::ExpressionDescriptionLengthType>(expressionText_.getValue()),
            maxLength_.get(),selfDescription_.get());
        outPort_.setData(pData_);
        //oldData->getPlotExpression().deletePlotExpressionNodes();
        delete oldData;
    }
}

PlotFunction* PlotFunctionSource::readInputString() {
    try {
        return new PlotFunction(PlotExpression(functionInput_.get(),expressionNameInput_.get()));
    }
    catch (VoreenException& ve) {
        LERROR("False Inputstring: " << ve.what());
        return new PlotFunction();
    }
}

void PlotFunctionSource::changeText() {
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

}



