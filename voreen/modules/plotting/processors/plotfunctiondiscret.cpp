/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany.                        *
 * Visualization and Computer Graphics Group <http://viscg.uni-muenster.de>        *
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

#include "plotfunctiondiscret.h"

#include "voreen/core/processors/processorwidget.h"
#include "voreen/core/processors/processorwidgetfactory.h"

#include "../datastructures/plotdata.h"
#include "../datastructures/plotfunction.h"

#include <limits>

namespace voreen {

const std::string PlotFunctionDiscret::loggerCat_("voreen.plotting.PlotFunctionDiscret");


PlotFunctionDiscret::PlotFunctionDiscret()
    :Processor()
    , outPort_(Port::OUTPORT,"PlotData.OutPort")
    , inPort_(Port::INPORT,"PlotPort.PlotFunction.InPort")
    , dimensionsList_("DimensionsList_","Dimension List",Processor::VALID)
    , leftInterval_("leftInterval_","Left Interval",-5,-std::numeric_limits<float>::max(),std::numeric_limits<float>::max(),Processor::VALID)
    , rightInterval_("rightInterval_","Right Interval",5,-std::numeric_limits<float>::max(),std::numeric_limits<float>::max(),Processor::VALID)
    , steps_("Steps_","Stepwidth", 0.1f,0.0001f,10,Processor::VALID)
    , startDiscretization_("startDiscretization_","Discretize again")
    , numberofCalculations_("NumberofCalculations","Number of Calculations",0,0,std::numeric_limits<float>::max(),Processor::VALID)
    , pData_(0)
    , pDataOut_(0)
{
    dimensionsList_.onChange(MemberFunctionCallback<PlotFunctionDiscret>(this, &PlotFunctionDiscret::choiceChange));
    leftInterval_.onChange(MemberFunctionCallback<PlotFunctionDiscret>(this, &PlotFunctionDiscret::leftIntervalChange));
    rightInterval_.onChange(MemberFunctionCallback<PlotFunctionDiscret>(this, &PlotFunctionDiscret::rightIntervalChange));
    steps_.onChange(MemberFunctionCallback<PlotFunctionDiscret>(this, &PlotFunctionDiscret::intervalStepChange));
    startDiscretization_.onClick(MemberFunctionCallback<PlotFunctionDiscret>(this,&PlotFunctionDiscret::calculate));
    addPort(inPort_);
    addPort(outPort_);

    addProperty(dimensionsList_);
    addProperty(leftInterval_);
    addProperty(rightInterval_);
    addProperty(steps_);
    addProperty(startDiscretization_);
    addProperty(numberofCalculations_);
    dimensionsList_.setGroupID("Discretization");
    leftInterval_.setGroupID("Discretization");
    rightInterval_.setGroupID("Discretization");
    steps_.setGroupID("Discretization");
    startDiscretization_.setGroupID("Discretization");
     setPropertyGroupGuiName("Discretization","Discretization");
    steps_.setNumDecimals(4);
    std::stringstream  optionText;
    std::stringstream  keyText;
    for (int i = 0; i < 26; ++i) {
        optionText.str("");
        optionText.clear();
        keyText.str("");
        keyText.clear();
        optionText << i+1 << " te Dimension";
        keyText << i+1 << " Column";
        dimensionsList_.addOption(keyText.str(),optionText.str(),i+1);
    }
    numberofCalculations_.setReadOnlyFlag(true);
    numberofCalculations_.setNumDecimals(0);
}


PlotFunctionDiscret::~PlotFunctionDiscret() {
}

Processor* PlotFunctionDiscret::create() const {
    return new PlotFunctionDiscret();
}

void PlotFunctionDiscret::initialize() {
    Processor::initialize();

    pDataOut_ = new PlotData(0,0);

    if (inPort_.hasData()) {
        readData();
        calculate();
    }
    setOutPortData();
}

void PlotFunctionDiscret::deinitialize() {
    outPort_.setData(0, false);
    delete pDataOut_;
    pData_ = 0;
    pDataOut_ = 0;

    Processor::deinitialize();
}

bool PlotFunctionDiscret::isReady() const {
    return true;
}

bool PlotFunctionDiscret::isEndProcessor() const {
    return (!outPort_.isConnected());
}

void PlotFunctionDiscret::serialize(Serializer &s) const {
    Processor::serialize(s);
    s.serialize("ResetList",resetList_);
}

void PlotFunctionDiscret::deserialize(Deserializer &s) {
    Processor::deserialize(s);
    try {
        s.deserialize("ResetList",resetList_);
    }
    catch (XmlSerializationNoSuchDataException&) {
        // superfluous or missing XML data item, just ignore
        s.removeLastError();
    }

}

void PlotFunctionDiscret::setOutPortData(){
    if (outPort_.isConnected()){
        if ((pDataOut_)){
            outPort_.setData(pDataOut_);
        }
        else
            outPort_.setData(0);
    }
}

void PlotFunctionDiscret::readData() {
    std::deque<Option<int> >  emptyDeque = std::deque<Option<int> >();
    dimensionsList_.setOptions(emptyDeque);
    if (dynamic_cast<const PlotFunction*>(inPort_.getData())) {
        pData_ = dynamic_cast<const PlotFunction*>(inPort_.getData());
    }
    else {
        dimensionsList_.addOption("1 Column","no key Column",1);
        resetList_.columns = std::pair<int,int>(0,0);
        LWARNING("PlotFunctionDiscret can only handle with PlotFunction-Objects.");
        return;
    }
    resetList_.columns = std::pair<int,int>(pData_->getKeyColumnCount(),pData_->getDataColumnCount());
    int length = -1;
    if (resetList_.columns.first > static_cast<int>(resetList_.interval.size())) {
        length = static_cast<int>(resetList_.interval.size());
        resetList_.interval.resize(resetList_.columns.first);
        resetList_.stepwidth.resize(resetList_.columns.first);
    }
    std::stringstream  optionText;
    std::stringstream  keyText;
    for (int i=0; i < resetList_.columns.first; ++i){
        optionText.str("");
        optionText.clear();
        keyText.str("");
        keyText.clear();
        optionText << i+1 << " te Dimension (" << pData_->getColumnLabel(i) << ")";
        if (length != -1 && i+1 > length) {
            resetList_.interval[i] = std::pair<int,int>(-5,5);
            resetList_.stepwidth[i] = 0.1;
        }
        keyText << i+1 << " Column";
        dimensionsList_.addOption(keyText.str(),optionText.str(),i+1);
    }
    dimensionsList_.updateWidgets();
    if (pData_->getKeyColumnCount() == 0) {
        dimensionsList_.addOption("1 Column","no key Column",1);
    }
    else {
        dimensionsList_.selectByValue(1);
    }
    calculateResultLength();
}



const PlotData* PlotFunctionDiscret::getPlotData() const {
    return pDataOut_;
}

void PlotFunctionDiscret::process() {
    tgtAssert(isInitialized(), "not initialized");
    if (inPort_.hasData() && inPort_.hasChanged()) {
        readData();
        calculate();
    }
    else if (!inPort_.hasData()){
        pData_ = 0;
        PlotData* oldPlotData;
        oldPlotData = pDataOut_;
        pDataOut_ = new PlotData(0,0);
        setOutPortData();
        if (getProcessorWidget()){
            getProcessorWidget()->updateFromProcessor();
        }
        delete oldPlotData;
    }
    else if (inPort_.hasData() && !inPort_.hasChanged()) {
        calculate();
    }
    else if (getProcessorWidget()){
        getProcessorWidget()->updateFromProcessor();
    }
}

void PlotFunctionDiscret::choiceChange() {
    int value = dimensionsList_.getValue();
    int length = static_cast<int>(resetList_.interval.size());
    if (value <= length) {
        leftInterval_.set(static_cast<float>(resetList_.interval.at(value-1).first));
        rightInterval_.set(static_cast<float>(resetList_.interval.at(value-1).second));
        steps_.set(static_cast<float>(resetList_.stepwidth.at(value-1)));
    }
    calculateResultLength();
}

void PlotFunctionDiscret::intervalStepChange() {
    int value = dimensionsList_.getValue();
    int length = static_cast<int>(resetList_.interval.size());
    if (value <= length) {
        resetList_.stepwidth[value-1] = steps_.get();
    }
    calculateResultLength();
}

void PlotFunctionDiscret::leftIntervalChange() {
    int value = dimensionsList_.getValue();
    int length = static_cast<int>(resetList_.interval.size());
    if (value <= length) {
        resetList_.interval[value-1].first = leftInterval_.get();
    }
    calculateResultLength();
}

void PlotFunctionDiscret::rightIntervalChange() {
    int value = dimensionsList_.getValue();
    int length = static_cast<int>(resetList_.interval.size());
    if (value <= length) {
        resetList_.interval[value-1].second = rightInterval_.get();
    }
    calculateResultLength();
}


void PlotFunctionDiscret::calculate() {
    if (!pData_)
        return;

    PlotData* newPlotData = new PlotData(0,0);

    if (!pData_->select(resetList_.interval,resetList_.stepwidth,*newPlotData)) {
        newPlotData->reset(0,0);
    }
    PlotData* oldData;
    oldData = pDataOut_;
    pDataOut_ = newPlotData;
    setOutPortData();
    if (getProcessorWidget()){
        getProcessorWidget()->updateFromProcessor();
    }
    delete oldData;
}

void PlotFunctionDiscret::calculateResultLength() {
    plot_t value = 1;
    plot_t plot_tValue;
    for (size_t i = 0; i < resetList_.interval.size(); ++i) {

        plot_tValue = (resetList_.interval[i].second - resetList_.interval[i].first)/resetList_.stepwidth[i];
        value *= std::floor(plot_tValue+0.500001);
    }
    numberofCalculations_.set(static_cast<float>(value));
    numberofCalculations_.updateWidgets();
}

void PlotFunctionDiscret::DiscretizeValues::serialize(Serializer& s) const {
    s.serialize("Interval",interval);
    s.serialize("Stepwidth",stepwidth);
    s.serialize("Columns",columns);
}

void PlotFunctionDiscret::DiscretizeValues::deserialize(Deserializer& s){
    s.deserialize("Interval",interval);
    s.deserialize("Stepwidth",stepwidth);
    s.deserialize("Columns",columns);
}

}
