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

#include "plotdatagroup.h"

#include "voreen/core/processors/processorwidget.h"
#include "voreen/core/processors/processorwidgetfactory.h"
#include "voreen/core/io/serialization/serialization.h"

namespace voreen {

const std::string PlotDataGroup::loggerCat_("voreen.PlotDataGroup");

PlotDataGroup::PlotDataGroup():
    Processor()
    , inPort_(Port::INPORT,"PlotData.InPort")
    , outPort_(Port::OUTPORT,"PlotData.OutPort")
    , pData_(0)
    , pDataOut_(0)
    , resetList_(0)
    , maxTableColumnWidth_("Max Column Width","Max Column Width",0,0,500,Processor::VALID)
    , keyColumnColor_("Key-Column Color","Key-Column Color",tgt::vec4(200.0/255.0,200.0/255.0,200.0/255.0,1.0),
        Processor::VALID)
    , dataColumnColor_("Data-Column Color","Data-Column Color",tgt::vec4(240.0/255.0,240.0/255.0,240.0/255.0,1.0),
        Processor::VALID)
{
    keyColumnColor_.onChange(MemberFunctionCallback<PlotDataGroup>(this, &PlotDataGroup::updateView));
    dataColumnColor_.onChange(MemberFunctionCallback<PlotDataGroup>(this, &PlotDataGroup::updateView));
    maxTableColumnWidth_.onChange(MemberFunctionCallback<PlotDataGroup>(this, &PlotDataGroup::updateView));

    addPort(inPort_);
    addPort(outPort_);

    addProperty(maxTableColumnWidth_);
    addProperty(keyColumnColor_);
    addProperty(dataColumnColor_);

    keyColumnColor_.setGroupID("Widget-Configuration");
    dataColumnColor_.setGroupID("Widget-Configuration");
    maxTableColumnWidth_.setGroupID("Widget-Configuration");

    setPropertyGroupGuiName("Widget-Configuration","Widget-Configuration");
}

PlotDataGroup::~PlotDataGroup() {
    deleteList();
}

Processor* PlotDataGroup::create() const {
    return new PlotDataGroup();
}

void PlotDataGroup::initialize() {
    Processor::initialize();

    pDataOut_ = new PlotData(0, 0);
    if (inPort_.hasData()) {
        readData();
        calculate();
    }
    setOutPortData();
}

void PlotDataGroup::deinitialize() {
    outPort_.setData(0, false);
    delete pDataOut_;
    pData_ = 0;
    pDataOut_ = 0;

    Processor::deinitialize();
}

void PlotDataGroup::process() {
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
    else {
        if (getProcessorWidget()){
            getProcessorWidget()->updateFromProcessor();
        }
    }
}

void PlotDataGroup::serialize(Serializer &s) const {
    Processor::serialize(s);
    s.serialize("ResetList",resetList_);
}

void PlotDataGroup::deserialize(Deserializer &s) {
    Processor::deserialize(s);
    try {
        s.deserialize("ResetList",resetList_);
    }
    catch (XmlSerializationNoSuchDataException&) {
        // superfluous or missing XML data item, just ignore
        s.removeLastError();
    }

}

void PlotDataGroup::fillResetList(int groupbyColumn, const AggregationVector& aggregationVector) {
    GroupValues groupValues;
    groupValues.aggregationVector = aggregationVector;
    groupValues.groupbyColumn = groupbyColumn;
    resetList_.push_back(groupValues);
}


const PlotData* PlotDataGroup::getPlotData() const {
    return pDataOut_;
}

void PlotDataGroup::updateView() {
    if (getProcessorWidget()){
        getProcessorWidget()->updateFromProcessor();
    }
}

void PlotDataGroup::deleteList() {
    for (size_t i = 0; i < resetList_.size(); ++i) {
        if (resetList_.at(i).aggregationVector.size() > 0) {
            for (size_t j = 0; j < resetList_.at(i).aggregationVector.size(); ++j) {
                if (resetList_.at(i).aggregationVector.at(j).second != NULL) {
                    delete resetList_.at(i).aggregationVector.at(j).second;
                }
            }
            resetList_.at(i).aggregationVector.clear();
        }
    }
    resetList_.clear();
}

void PlotDataGroup::calculate() {
    tgtAssert(isInitialized(), "not initialized");
    if (!pData_)
        return;

    PlotData* newPlotData = new PlotData(*pData_);
    PlotData* oldData;
    oldData = pDataOut_;
    pDataOut_ = PlotDataGroupby(newPlotData);
    setOutPortData();
    if (getProcessorWidget()){
        getProcessorWidget()->updateFromProcessor();
    }
    if (pData_ != oldData && oldData != 0) {
        delete oldData;
    }
}

void PlotDataGroup::groupbyWidget(int column, const std::vector< std::pair<int,AggregationFunction*> >& functions){
    PlotData* newData = new PlotData(1,1);
    PlotData* oldpData;
    pDataOut_->groupBy(column,functions,*newData);
    oldpData = pDataOut_;
    pDataOut_ = newData;
    setOutPortData();
    if (getProcessorWidget()){
        getProcessorWidget()->updateFromProcessor();
    }
    delete oldpData;
    fillResetList(column,functions);
}

bool PlotDataGroup::isReady() const {
    return true;
}

bool PlotDataGroup::isEndProcessor() const {
    return (!outPort_.isConnected());
}

PlotData* PlotDataGroup::PlotDataGroupby(PlotData* actPlotData){
    PlotData* newData = 0;
    PlotData* oldData = actPlotData;
    for (size_t i = 0; i < resetList_.size(); ++i) {
        newData = new PlotData(1,1);
        AggregationVector* aggvec = &resetList_.at(i).aggregationVector;

        if (aggvec != NULL) {
            oldData->groupBy(resetList_.at(i).groupbyColumn,*aggvec,*newData);
        }
        else {
            oldData->groupBy(resetList_.at(i).groupbyColumn,*aggvec,*newData);
        }
        delete oldData;
        oldData = newData;
    }
    if (resetList_.size() == 0)
        newData = actPlotData;
    return newData;
}

void PlotDataGroup::resetAllData() {
    deleteList();
    calculate();
}

void PlotDataGroup::resetLastData() {
    if (resetList_.size() > 0) {
        resetList_.pop_back();
        calculate();
    }
}

void PlotDataGroup::setOutPortData(){
    if (outPort_.isConnected()){
        if ((pDataOut_)){
            outPort_.setData(pDataOut_);
        }
        else
            outPort_.setData(0);
    }
}

void PlotDataGroup::readData() {
    if (dynamic_cast<const PlotData*>(inPort_.getData())) {
        pData_ = dynamic_cast<const PlotData*>(inPort_.getData());
    }
    else {
        pData_ = 0;
        LWARNING("PlotDataGroup can only handle with PlotData-Objects.");
    }
}


void PlotDataGroup::GroupValues::serialize(Serializer& s) const {
    s.serialize("GroupbyColumn",groupbyColumn);
    s.serialize("AggregationVector",aggregationVector);
}

void PlotDataGroup::GroupValues::deserialize(Deserializer& s){
    s.deserialize("GroupbyColumn",groupbyColumn);
    s.deserialize("AggregationVector",aggregationVector);
}



}// NAMESPACE VOREEN

