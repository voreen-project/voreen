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

#include "plotdataselect.h"

#include "voreen/core/processors/processorwidget.h"
#include "voreen/core/processors/processorwidgetfactory.h"
#include "voreen/core/io/serialization/serialization.h"
#include "../datastructures/plotdata.h"
#include "../datastructures/plotpredicate.h"

#include <limits>

namespace voreen {

const std::string PlotDataSelect::loggerCat_("voreen.PlotDataSelect");

PlotDataSelect::PlotDataSelect():
    Processor()
    , inPort_(Port::INPORT,"PlotData.InPort")
    , outPort_(Port::OUTPORT,"PlotData.OutPort")
    , pData_(0)
    , pDataOut_(0)
    , resetList_(0)
    , select_("select","Select",false,Processor::VALID)
    , newKeyColumnCount_("newKeyColumnCount","new Key Column Count",1,1,200,Processor::VALID)
    , plotPredicateVector_("PlotPredicateVector","PlotPredicateVector",Processor::VALID)
    , columnorder_("columnorder","Column Order",Processor::VALID)
    , filter_("filter","Filter",false,Processor::VALID)
    , maxTableColumnWidth_("Max Column Width","Max Column Width",0,0,500,Processor::VALID)
    , keyColumnColor_("Key-Column Color","Key-Column Color",tgt::vec4(200.0/255.0,200.0/255.0,200.0/255.0,1.0),
       Processor::VALID)
    , dataColumnColor_("Data-Column Color","Data-Column Color",tgt::vec4(240.0/255.0,240.0/255.0,240.0/255.0,1.0),
       Processor::VALID)
{
    select_.onChange(MemberFunctionCallback<PlotDataSelect>(this, &PlotDataSelect::checkBoxChange));
    columnorder_.onChange(MemberFunctionCallback<PlotDataSelect>(this, &PlotDataSelect::checkBoxChange));
    filter_.onChange(MemberFunctionCallback<PlotDataSelect>(this, &PlotDataSelect::checkBoxChange));
    newKeyColumnCount_.onChange(MemberFunctionCallback<PlotDataSelect>(this,&PlotDataSelect::keyColumnChange));
    keyColumnColor_.onChange(MemberFunctionCallback<PlotDataSelect>(this, &PlotDataSelect::updateView));
    dataColumnColor_.onChange(MemberFunctionCallback<PlotDataSelect>(this, &PlotDataSelect::updateView));
    maxTableColumnWidth_.onChange(MemberFunctionCallback<PlotDataSelect>(this, &PlotDataSelect::updateView));

    addProperty(select_);
    addProperty(newKeyColumnCount_);
    addProperty(plotPredicateVector_);
    addProperty(columnorder_);
    addProperty(filter_);

    addProperty(maxTableColumnWidth_);
    addProperty(keyColumnColor_);
    addProperty(dataColumnColor_);

    addPort(inPort_);
    addPort(outPort_);

    select_.setGroupID("Select");
    newKeyColumnCount_.setGroupID("Select");
    columnorder_.setGroupID("Column Order");
    filter_.setGroupID("Filter");
    keyColumnColor_.setGroupID("Widget-Configuration");
    dataColumnColor_.setGroupID("Widget-Configuration");
    maxTableColumnWidth_.setGroupID("Widget-Configuration");

    setPropertyGroupGuiName("Select","Select");
    setPropertyGroupGuiName("Column Order","Column Order");
    setPropertyGroupGuiName("Filter","Filter");
    setPropertyGroupGuiName("Widget-Configuration","Widget-Configuration");
    newKeyColumnCount_.setVisibleFlag(false);
    plotPredicateVector_.setVisibleFlag(false);
}

PlotDataSelect::~PlotDataSelect() {
    deleteList();
}

Processor* PlotDataSelect::create() const {
    return new PlotDataSelect();
}

void PlotDataSelect::initialize() {
    Processor::initialize();

    pDataOut_ = new PlotData(0,0);

    /*if (inPort_.hasData()) {
        readData();
        calculate();
    }*/
    setOutPortData();
}

void PlotDataSelect::deinitialize() {
    outPort_.setData(0, false);
    delete pDataOut_;
    pData_ = 0;
    pDataOut_ = 0;

    Processor::deinitialize();
}

void PlotDataSelect::process() {
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

void PlotDataSelect::serialize(Serializer &s) const {
    Processor::serialize(s);
    s.serialize("resetList",resetList_);
}

void PlotDataSelect::deserialize(Deserializer &s) {
    Processor::deserialize(s);
    try {
        s.deserialize("resetList",resetList_);
    }
    catch (SerializationNoSuchDataException&) {
        // superfluous or missing XML data item, just ignore
        s.removeLastError();
    }
}

void PlotDataSelect::deleteList() {
    for (size_t i = 0; i < resetList_.size(); ++i) {
        if (resetList_.at(i).predicateVector.size() != 0) {
            for (size_t j = 0; j < resetList_.at(i).predicateVector.size(); ++j) {
                if (resetList_.at(i).predicateVector.at(j).second != NULL) {
                    delete resetList_.at(i).predicateVector.at(j).second;
                }
            }
            resetList_.at(i).predicateVector.clear();
        }
    }
    resetList_.clear();
    plotPredicateVector_.clear();
}

const PlotData* PlotDataSelect::getPlotData() const {
    return pDataOut_;
}

void PlotDataSelect::updateView() {
    if (getProcessorWidget()){
        getProcessorWidget()->updateFromProcessor();
    }
}

void PlotDataSelect::calculate() {
    tgtAssert(isInitialized(), "not initialized");
    if (!pData_)
        return;

    if ((select_.get()&& !newKeyColumnCount_.isVisibleFlagSet()) ||
        (filter_.get() && !select_.isReadOnlyFlagSet()) ||
        (columnorder_.get() && !select_.isReadOnlyFlagSet()))
        checkBoxChange();

    PlotData* newPlotData;
    PlotData* oldPlotData;
    newPlotData = new PlotData(*pData_);
    if (!select_.get() && !columnorder_.get() && !filter_.get()) {
        deleteList();
    }
    else if (select_.get() && newPlotData->getColumnCount() > 0) {
        newPlotData = selectData(newPlotData);
        newKeyColumnCount_.setMaxValue(std::max<int>(newPlotData->getColumnCount(),1));
        if (newPlotData->getKeyColumnCount() != newKeyColumnCount_.get()) {
            newKeyColumnCount_.set(std::max<int>(newPlotData->getKeyColumnCount(),1));
        }
    }
    else if (columnorder_.get() && newPlotData->getColumnCount() > 0) {
        newPlotData = selectData(newPlotData);
    }
    else if (filter_.get() && newPlotData->getColumnCount() > 0) {
        filterData(newPlotData);
    }
    oldPlotData = pDataOut_;
    pDataOut_ = newPlotData;
    setOutPortData();
    if (getProcessorWidget()){
        getProcessorWidget()->updateFromProcessor();
    }

    if (pData_ != oldPlotData && oldPlotData != 0) {
        delete oldPlotData;
    }
}

bool PlotDataSelect::isReady() const {
    return true;
}

bool PlotDataSelect::isEndProcessor() const {
    return (!outPort_.isConnected());
}

void PlotDataSelect::checkBoxChange() {

    if (!isInitialized())
        return;

    if (select_.get()) {
        columnorder_.setReadOnlyFlag(true);
        filter_.setReadOnlyFlag(true);
        newKeyColumnCount_.setVisibleFlag(true);
    }
    else if (columnorder_.get()) {
        select_.setReadOnlyFlag(true);
        filter_.setReadOnlyFlag(true);
    }
    else if (filter_.get()) {
        columnorder_.setReadOnlyFlag(true);
        select_.setReadOnlyFlag(true);
    }
    else {
        select_.setReadOnlyFlag(false);
        columnorder_.setReadOnlyFlag(false);
        filter_.setReadOnlyFlag(false);
        newKeyColumnCount_.setVisibleFlag(false);
    }
    calculate();
    if (!pData_)
        updateView();
}

void PlotDataSelect::fillResetList(std::vector<int> columns, const std::vector<int>& rows,
                                   const PredicateVector& predVector, const LabelVector& columntexts,
                                  int keyColumnCount, int dataColumnCount) {
    SelectValues selectValues;
    selectValues.projectionVector = columns;
    selectValues.selectedRows = rows;
    selectValues.predicateVector = predVector;
    selectValues.labelVector = columntexts;
    selectValues.keyColumnCount = keyColumnCount;
    selectValues.dataColumnCount = dataColumnCount;
    resetList_.push_back(selectValues);
}

size_t PlotDataSelect::getResetListSize() const {
    return resetList_.size();
}

void PlotDataSelect::keyColumnChange() {

    if (!isInitialized() || !pData_ || !pDataOut_)
        return;
    if (pData_->getColumnCount() == pDataOut_->getColumnCount()) {
        int len = pDataOut_->getColumnCount();
        if (pDataOut_->getKeyColumnCount() == newKeyColumnCount_.get())
            return;
        std::vector<int> columns;
        columns.resize(len);
        for (int i = 0; i < len; ++i) {
            columns[i] = i;
        }
        PlotData* newData = new PlotData(1,1);
        PlotData* oldData;
        pDataOut_->select(columns,newKeyColumnCount_.get(),len-newKeyColumnCount_.get(),*newData);
        fillResetList(columns,std::vector<int>(),PredicateVector(),LabelVector(),newData->getKeyColumnCount(),newData->getDataColumnCount());
        oldData = pDataOut_;
        pDataOut_ = newData;
        setOutPortData();
        if (getProcessorWidget()){
            getProcessorWidget()->updateFromProcessor();
        }
        delete oldData;
    }
}

void PlotDataSelect::selectFunctionalityType(FunctionLibrary::ProcessorFunctionalityType type) {
    if (select_.get()) {
        select_.set(false);
        if (type != FunctionLibrary::NONE) {
            selectFunctionalityType(type);
        }
    }
    else if (columnorder_.get()) {
        columnorder_.set(false);
        if (type != FunctionLibrary::NONE) {
            selectFunctionalityType(type);
        }
    }
    else if (filter_.get()) {
        filter_.set(false);
        if (type != FunctionLibrary::NONE) {
            selectFunctionalityType(type);
        }
    }
    else {
        switch (type) {
            case FunctionLibrary::SELECT: {
                        select_.set(true);
                        break;
                    }
            case FunctionLibrary::COLUMNORDER: {
                        columnorder_.set(true);
                        break;
                    }
            case FunctionLibrary::FILTER: {
                        filter_.set(true);
                        break;
                    }
            default : break;
        }
    }

}

PlotData* PlotDataSelect::selectData(PlotData* plotData){
    PlotData* newPlotData;
    PlotData* oldPlotData;
    oldPlotData = plotData;
    newPlotData = plotData;
    bool falseValue = false;
    for (size_t i = 0; i < resetList_.size(); ++i) {
        PredicateVector predvec = resetList_.at(i).predicateVector;
        if (resetList_.at(i).labelVector.size() == 0) {
            newPlotData = new PlotData(1,1);
            for (size_t j = 0; j < resetList_.at(i).selectedRows.size() && !falseValue; ++j) {
                falseValue = resetList_.at(i).selectedRows.at(j) >= oldPlotData->getRowsCount();
            }
            for (size_t j = 0; j < resetList_.at(i).projectionVector.size() && !falseValue; ++j) {
                falseValue = resetList_.at(i).projectionVector.at(j) >= oldPlotData->getColumnCount();
            }
            if (predvec.size() != 0 && !falseValue) {
                if (resetList_.at(i).selectedRows.size() > 0)
                    oldPlotData->select(resetList_.at(i).projectionVector,resetList_.at(i).selectedRows,resetList_.at(i).keyColumnCount, resetList_.at(i).dataColumnCount,*newPlotData);
                else
                    oldPlotData->select(resetList_.at(i).projectionVector,resetList_.at(i).keyColumnCount, resetList_.at(i).dataColumnCount,predvec,*newPlotData);
            }
            else if (!falseValue) {
                if (resetList_.at(i).selectedRows.size() > 0)
                    oldPlotData->select(resetList_.at(i).projectionVector,resetList_.at(i).selectedRows,resetList_.at(i).keyColumnCount, resetList_.at(i).dataColumnCount,*newPlotData);
                else
                    oldPlotData->select(resetList_.at(i).projectionVector,resetList_.at(i).keyColumnCount, resetList_.at(i).dataColumnCount,*newPlotData);
            }
        }
        else {
            newPlotData = new PlotData(*oldPlotData);
            for (size_t j = 0; j < resetList_.at(i).labelVector.size(); ++j) {
                if (resetList_.at(i).labelVector.at(j).first < newPlotData->getColumnCount()) {
                    newPlotData->setColumnLabel(resetList_.at(i).labelVector.at(j).first,resetList_.at(i).labelVector.at(j).second);
                }
                else {
                    falseValue = true;
                    break;
                }
            }
        }
        delete oldPlotData;
        oldPlotData = newPlotData;
        if (falseValue) {
            delete newPlotData;
            newPlotData = new PlotData(0,0);
            break;
        }
    }
    if (resetList_.size() == 0 || falseValue) {
        plotPredicateVector_.clear();
    }
    else {
        PredicateVector predvec = resetList_.at(resetList_.size() - 1).predicateVector;
        if (predvec.size() != 0)
            plotPredicateVector_.set(predvec);
        else
            plotPredicateVector_.set(PredicateVector());

    }
    return newPlotData;
}

void PlotDataSelect::disSelectColumn(int column) {
    int len = pDataOut_->getColumnCount();
    if (column >= len)
        return;
    std::vector<int> columns;
    if (len-1 == 0) {
        columns.resize(0);
    }
    else {
        columns.resize(len - 1);
        for (int i = 0; i < len; ++i) {
            if (i == column) {
                continue;
            }
            else if (i < column) {
                columns[i] = i;
            }
            else {
                columns[i - 1] = i;
            }
        }
    }
    PlotData* newData = new PlotData(1,1);
    PlotData* oldData;
    int newkeycount = 0;
    int newdatacount;
    if (len-1-newKeyColumnCount_.get() < 0) {
        newkeycount = pDataOut_->getColumnCount() - 1;
        newdatacount = 0;
    }
    else {
        newdatacount = len - 1 - newKeyColumnCount_.get();
        newkeycount = newKeyColumnCount_.get();
    }
    std::vector< std::pair<int, PlotPredicate*> > liste = (plotPredicateVector_.getCloned());
    if (liste.size() > 0) {
        pDataOut_->select(columns,newkeycount,newdatacount, liste, *newData);
        fillResetList(columns,std::vector<int>(), liste, LabelVector(), newData->getKeyColumnCount(), newData->getDataColumnCount());
        plotPredicateVector_.select(columns);
    }
    else {
        pDataOut_->select(columns,newkeycount,newdatacount,*newData);
        fillResetList(columns,std::vector<int>(),PredicateVector(),LabelVector(),newData->getKeyColumnCount(),newData->getDataColumnCount());
    }
    oldData = pDataOut_;
    pDataOut_ = newData;
    setOutPortData();
    if (getProcessorWidget()){
        getProcessorWidget()->updateFromProcessor();
    }
    delete oldData;
    if (newkeycount != newKeyColumnCount_.get()) {
        newKeyColumnCount_.set(std::max<int>(newkeycount,1));
    }
}

void PlotDataSelect::disSelectColumns(const std::vector<int>& columns) {
    int len = pDataOut_->getColumnCount();
    std::vector<int> iColumns;
    if (len-columns.size() == 0) {
        iColumns.resize(0);
    }
    else {
        iColumns.resize(len-columns.size());
        bool match;
        int z = 0;
        for (int i = 0; i < len; ++i) {
            match = false;
            for (size_t j = 0; j < columns.size(); ++j) {
                if (columns.at(j) == i) {
                    match = true;
                    break;
                }
            }
            if (!match) {
                iColumns[z] = i;
                ++z;
            }
        }
    }
    PlotData* newData = new PlotData(1,1);
    PlotData* oldData;
    int newkeycount = 0;
    int newdatacount;
    if (len-static_cast<int>(columns.size())-newKeyColumnCount_.get() < 0) {
        newkeycount = len-static_cast<int>(columns.size());
        newdatacount = 0;
    }
    else {
        newdatacount = len-static_cast<int>(columns.size())-newKeyColumnCount_.get();
        newkeycount = newKeyColumnCount_.get();
    }
    std::vector< std::pair<int, PlotPredicate*> > liste = (plotPredicateVector_.getCloned());
    if (liste.size() > 0) {
        pDataOut_->select(iColumns,newkeycount,newdatacount, liste,*newData);
        fillResetList(iColumns,std::vector<int>(), liste,LabelVector(),newData->getKeyColumnCount(),newData->getDataColumnCount());
        plotPredicateVector_.select(iColumns);
    }
    else {
        pDataOut_->select(iColumns,newkeycount,newdatacount,*newData);
        fillResetList(iColumns,std::vector<int>(),PredicateVector(),LabelVector(),newData->getKeyColumnCount(),newData->getDataColumnCount());
    }
    oldData = pDataOut_;
    pDataOut_ = newData;
    setOutPortData();
    if (getProcessorWidget()){
        getProcessorWidget()->updateFromProcessor();
    }
    delete oldData;
    if (newkeycount != newKeyColumnCount_.get()) {
        newKeyColumnCount_.set(std::max<int>(newkeycount,1));
    }
}

void PlotDataSelect::selectColumns(const std::vector<int>& columns) {
    int len = static_cast<int>(columns.size());
    std::vector<int> iColumns;
    if (len == 0) {
        iColumns.resize(0);
    }
    else {
        iColumns.resize(len);
        for (int i = 0; i < len; ++i) {
            iColumns[i] = columns.at(i);
        }
    }
    PlotData* newData = new PlotData(1,1);
    PlotData* oldData;
    int newkeycount = 0;
    int newdatacount;
    if (len-newKeyColumnCount_.get() < 0) {
        newkeycount = len;
        newdatacount = 0;
    }
    else {
        newdatacount = len-newKeyColumnCount_.get();
        newkeycount = newKeyColumnCount_.get();
    }
    std::vector< std::pair<int, PlotPredicate*> > liste = (plotPredicateVector_.getCloned());
    if (liste.size() > 0) {
        pDataOut_->select(iColumns,newkeycount,newdatacount, liste,*newData);
        fillResetList(iColumns,std::vector<int>(), liste,LabelVector(),newData->getKeyColumnCount(),newData->getDataColumnCount());
        plotPredicateVector_.select(iColumns);
    }
    else {
        pDataOut_->select(iColumns,newkeycount,newdatacount,*newData);
        fillResetList(iColumns,std::vector<int>(),PredicateVector(),LabelVector(),newData->getKeyColumnCount(),newData->getDataColumnCount());
    }
    oldData = pDataOut_;
    pDataOut_ = newData;
    setOutPortData();
    if (getProcessorWidget()){
        getProcessorWidget()->updateFromProcessor();
    }
    delete oldData;
    if (newkeycount != newKeyColumnCount_.get()) {
        newKeyColumnCount_.set(std::max<int>(newkeycount,1));
    }
}

void PlotDataSelect::selectColumnsRows(const std::vector<int>& columns, const std::vector<int>& rows) {
    int len = static_cast<int>(columns.size());
    std::vector<int> iColumns;
    if (len == 0) {
        iColumns.resize(0);
    }
    else {
        iColumns.resize(len);
        for (int i = 0; i < len; ++i) {
            iColumns[i] = columns.at(i);
        }
    }
    PlotData* newData = new PlotData(1,1);
    PlotData* oldData;
    int newkeycount = 0;
    int newdatacount;
    if (len-newKeyColumnCount_.get() < 0) {
        newkeycount = len;
        newdatacount = 0;
    }
    else {
        newdatacount = len-  newKeyColumnCount_.get();
        newkeycount = newKeyColumnCount_.get();
    }
    std::vector<int> values = rows;
    pDataOut_->select(iColumns,values,newkeycount,newdatacount,*newData);
    fillResetList(iColumns,values,PredicateVector(),LabelVector(),newData->getKeyColumnCount(),newData->getDataColumnCount());
    plotPredicateVector_.select(iColumns);
    oldData = pDataOut_;
    pDataOut_ = newData;
    setOutPortData();
    if (getProcessorWidget()){
        getProcessorWidget()->updateFromProcessor();
    }
    delete oldData;
    if (newkeycount != newKeyColumnCount_.get()) {
        newKeyColumnCount_.set(std::max<int>(newkeycount,1));
    }
}

void PlotDataSelect::selectIncludeArea(const std::pair<int,int>& point1, const std::pair<int,int>& point2) {
    int len = point2.first - point1.first + 1;
    std::vector<int> columns;
    if (len == 0) {
        columns.resize(0);
    }
    else {
        columns.resize(len);
        for (int i = 0; i < len; ++i) {
            columns[i] = point1.first + i;
        }
    }
    PlotData* newData = new PlotData(1,1);
    PlotData* oldData;
    int newkeycount = 0;
    int newdatacount;
    if (len-newKeyColumnCount_.get() < 0) {
        newkeycount = len;
        newdatacount = 0;
    }
    else {
        newdatacount = len-newKeyColumnCount_.get();
        newkeycount = newKeyColumnCount_.get();
    }
    std::vector<int> rows;
    for (int j = point1.second; j <= point2.second; ++j) {
        rows.push_back(j);
    }
    pDataOut_->select(columns,rows,newkeycount,newdatacount,*newData);
    fillResetList(columns,rows,PredicateVector(),LabelVector(),newData->getKeyColumnCount(),newData->getDataColumnCount());
    plotPredicateVector_.select(columns);
    oldData = pDataOut_;
    pDataOut_ = newData;
    setOutPortData();
    if (getProcessorWidget()){
        getProcessorWidget()->updateFromProcessor();
    }
    delete oldData;
    if (newkeycount != newKeyColumnCount_.get()) {
        newKeyColumnCount_.set(std::max<int>(newkeycount,1));
    }
}

void PlotDataSelect::renameColumn(int column,std::string caption) {
    LabelVector labelVector;
    labelVector.push_back(std::pair< int, std::string>(column,caption));
    fillResetList(std::vector<int>(),std::vector<int>(),PredicateVector(),labelVector,pDataOut_->getKeyColumnCount(),pDataOut_->getDataColumnCount());
    pDataOut_->setColumnLabel(column,caption);
    setOutPortData();
    if (getProcessorWidget()){
        getProcessorWidget()->updateFromProcessor();
    }
}

void PlotDataSelect::resetAllData() {
    deleteList();
    calculate();
}

void PlotDataSelect::resetLastData() {
    if (resetList_.size() > 0) {
        if (resetList_.at(resetList_.size()-1).predicateVector.size() != 0) {
            for (size_t j = 0; j < resetList_.at(resetList_.size()-1).predicateVector.size(); ++j) {
                if (resetList_.at(resetList_.size()-1).predicateVector.at(j).second != NULL)
                    delete resetList_.at(resetList_.size()-1).predicateVector.at(j).second;
            }
        }
        resetList_.pop_back();
        calculate();
    }
}

void PlotDataSelect::selectPredicate() {
    int len = pDataOut_->getColumnCount();
    std::vector<int> columns;
    columns.resize(len);
    for (int i = 0; i < len; ++i) {
        columns[i] = i;
    }
    std::vector< std::pair<int, PlotPredicate*> > predicates = (plotPredicateVector_.getCloned());
    PlotData* newData = new PlotData(1,1);
    PlotData* oldData;
    pDataOut_->select(columns,pDataOut_->getKeyColumnCount(),pDataOut_->getDataColumnCount(), predicates,*newData);
    fillResetList(columns,std::vector<int>(), predicates,LabelVector(),newData->getKeyColumnCount(),newData->getDataColumnCount());
    oldData = pDataOut_;
    pDataOut_ = newData;
    setOutPortData();
    if (getProcessorWidget()){
        getProcessorWidget()->updateFromProcessor();
    }
    delete oldData;
}

void PlotDataSelect::copyColumn(int copycolumn, int targetcolumn) {
    int len = pDataOut_->getColumnCount() + 1;
    if (copycolumn >= len - 1 || targetcolumn > len - 1)
        return;
    if (targetcolumn < 0)
        targetcolumn = len - 1;
    std::vector<int> columns;
    columns.resize(len);
    for (int i = 0; i < len; ++i) {
        if (i < targetcolumn)
            columns[i] = i;
        else if (i == targetcolumn)
            columns[i] = copycolumn;
        else
            columns[i] = i-1;
    }
    PlotData* newData = new PlotData(1,1);
    PlotData* oldData;
    std::vector< std::pair<int, PlotPredicate*> > liste = (plotPredicateVector_.getCloned());
    if (liste.size() > 0) {
        pDataOut_->select(columns,pDataOut_->getKeyColumnCount(),len-pDataOut_->getKeyColumnCount(), liste,*newData);
        fillResetList(columns,std::vector<int>(), liste,LabelVector(),newData->getKeyColumnCount(),newData->getDataColumnCount());
        plotPredicateVector_.select(columns);
    }
    else {
        pDataOut_->select(columns,pDataOut_->getKeyColumnCount(),len-pDataOut_->getKeyColumnCount(),*newData);
        fillResetList(columns,std::vector<int>(),PredicateVector(),LabelVector(),newData->getKeyColumnCount(),newData->getDataColumnCount());
    }
    oldData = pDataOut_;
    pDataOut_ = newData;
    setOutPortData();
    if (getProcessorWidget()){
        getProcessorWidget()->updateFromProcessor();
    }
    delete oldData;
}

void PlotDataSelect::slideColumnBefore(int oldcolumn, int beforecolumn) {
    int len = pDataOut_->getColumnCount();
    if (oldcolumn >= len || beforecolumn >= len || oldcolumn == beforecolumn)
        return;
    std::vector<int> columns;
    columns.resize(len);
    for (int i = 0; i < len; ++i) {
        if (oldcolumn > beforecolumn) {
            if (i < beforecolumn || i > oldcolumn) {
                columns[i] = i;
            }
            else if (i == beforecolumn) {
                columns[i] = oldcolumn;
            }
            else {
                columns[i] = i-1;
            }
        }
        else {
            if (i >= beforecolumn || i < oldcolumn) {
                columns[i] = i;
            }
            else if (i == beforecolumn-1) {
                columns[i] = oldcolumn;
            }
            else {
                columns[i] = i+1;
            }

        }
    }
    PlotData* newData = new PlotData(1,1);
    PlotData* oldData;
    std::vector< std::pair<int, PlotPredicate*> > liste = (plotPredicateVector_.getCloned());
    if (liste.size() > 0) {
        pDataOut_->select(columns,pDataOut_->getKeyColumnCount(),pDataOut_->getDataColumnCount(), liste,*newData);
        fillResetList(columns,std::vector<int>(), liste,LabelVector(),newData->getKeyColumnCount(),newData->getDataColumnCount());
        plotPredicateVector_.select(columns);
    }
    else {
        pDataOut_->select(columns,pDataOut_->getKeyColumnCount(),pDataOut_->getDataColumnCount(),*newData);
        fillResetList(columns,std::vector<int>(),PredicateVector(),LabelVector(),newData->getKeyColumnCount(),newData->getDataColumnCount());
    }
    oldData = pDataOut_;
    pDataOut_ = newData;
    setOutPortData();
    if (getProcessorWidget()){
        getProcessorWidget()->updateFromProcessor();
    }
    delete oldData;
}

void PlotDataSelect::slideColumnAfter(int oldcolumn, int aftercolumn) {
    int len = pDataOut_->getColumnCount();
    if (oldcolumn >= len || aftercolumn >= len || oldcolumn == aftercolumn)
        return;
    std::vector<int> columns;
    columns.resize(len);
    for (int i = 0; i < len; ++i) {
        if (oldcolumn > aftercolumn) {
            if (i <= aftercolumn || i > oldcolumn) {
                columns[i] = i;
            }
            else if (i == aftercolumn+1) {
                columns[i] = oldcolumn;
            }
            else {
                columns[i] = i-1;
            }
        }
        else {
            if (i > aftercolumn || i < oldcolumn) {
                columns[i] = i;
            }
            else if (i == aftercolumn) {
                columns[i] = oldcolumn;
            }
            else {
                columns[i] = i+1;
            }

        }
    }
    PlotData* newData = new PlotData(1,1);
    PlotData* oldData;
    std::vector< std::pair<int, PlotPredicate*> > liste = (plotPredicateVector_.getCloned());
    if (liste.size() > 0) {
        pDataOut_->select(columns,pDataOut_->getKeyColumnCount(),pDataOut_->getDataColumnCount(), liste,*newData);
        fillResetList(columns,std::vector<int>(), liste,LabelVector(),newData->getKeyColumnCount(),newData->getDataColumnCount());
        plotPredicateVector_.select(columns);
    }
    else {
        pDataOut_->select(columns,pDataOut_->getKeyColumnCount(),pDataOut_->getDataColumnCount(),*newData);
        fillResetList(columns,std::vector<int>(),PredicateVector(),LabelVector(),newData->getKeyColumnCount(),newData->getDataColumnCount());
    }
    oldData = pDataOut_;
    pDataOut_ = newData;
    setOutPortData();
    if (getProcessorWidget()){
        getProcessorWidget()->updateFromProcessor();
    }
    delete oldData;
}

void PlotDataSelect::slideColumnSwitch(int firstcolumn, int secondcolumn) {
    int len = pDataOut_->getColumnCount();
    if (firstcolumn >= len || secondcolumn >= len || firstcolumn == secondcolumn)
        return;
    std::vector<int> columns;
    columns.resize(len);
    for (int i = 0; i < len; ++i) {
        if (i == firstcolumn) {
            columns[i] = secondcolumn;
        }
        else if (i == secondcolumn) {
            columns[i] = firstcolumn;
        }
        else {
            columns[i] = i;
        }
    }
    PlotData* newData = new PlotData(1,1);
    PlotData* oldData;
    std::vector< std::pair<int, PlotPredicate*> > liste = (plotPredicateVector_.getCloned());
    if (liste.size() >0 ) {
        pDataOut_->select(columns,pDataOut_->getKeyColumnCount(),pDataOut_->getDataColumnCount(), liste,*newData);
        fillResetList(columns,std::vector<int>(), liste,LabelVector(),newData->getKeyColumnCount(),newData->getDataColumnCount());
        plotPredicateVector_.select(columns);
    }
    else {
        pDataOut_->select(columns,pDataOut_->getKeyColumnCount(),pDataOut_->getDataColumnCount(),*newData);
        fillResetList(columns,std::vector<int>(),PredicateVector(),LabelVector(),newData->getKeyColumnCount(),newData->getDataColumnCount());
    }
    oldData = pDataOut_;
    pDataOut_ = newData;
    setOutPortData();
    if (getProcessorWidget()){
        getProcessorWidget()->updateFromProcessor();
    }
    delete oldData;
}

void PlotDataSelect::filterData(PlotData* plotData){
    for (size_t i = 0; i < resetList_.size(); ++i) {
        plotData->remove((resetList_.at(i).predicateVector));
    }
}

void PlotDataSelect::filterEmptyCells() {
    PlotData* newData = new PlotData(*pDataOut_);
    PlotData* oldData;
    std::vector< std::pair<int, PlotPredicate*> >* predicates = new std::vector< std::pair<int, PlotPredicate*> >();
    for (int i = 0; i < newData->getColumnCount(); ++i) {
        predicates->push_back(std::pair<int, PlotPredicateEmpty*>(i, new PlotPredicateEmpty()));
    }
    if (newData->remove(*predicates) > 0) {
        fillResetList(std::vector<int>(),std::vector<int>(),*predicates,LabelVector(),newData->getKeyColumnCount(),newData->getDataColumnCount());
    }
    oldData = pDataOut_;
    pDataOut_ = newData;
    setOutPortData();
    if (getProcessorWidget()){
        getProcessorWidget()->updateFromProcessor();
    }
    delete oldData;
}

void PlotDataSelect::filterAlphaNumeric() {
    PlotData* newData = new PlotData(*pDataOut_);
    PlotData* oldData;
    std::vector< std::pair<int, PlotPredicate*> >* predicates = new std::vector< std::pair<int, PlotPredicate*> >();
    for (int i = 0; i < newData->getColumnCount(); ++i) {
        predicates->push_back(std::pair<int, PlotPredicateNotAlphaNumeric*>(i, new PlotPredicateNotAlphaNumeric()));
    }
    if (newData->remove(*predicates) > 0) {
        fillResetList(std::vector<int>(),std::vector<int>(),*predicates,LabelVector(),newData->getKeyColumnCount(),newData->getDataColumnCount());
    }
    oldData = pDataOut_;
    pDataOut_ = newData;
    setOutPortData();
    if (getProcessorWidget()){
        getProcessorWidget()->updateFromProcessor();
    }
    delete oldData;
}

void PlotDataSelect::setOutPortData(){
    if (outPort_.isConnected()){
        if ((pDataOut_)){
            outPort_.setData(pDataOut_);
        }
        else
            outPort_.setData(0);
    }

}

void PlotDataSelect::readData() {

    if (inPort_.hasChanged()) {
        newKeyColumnCount_.setMaxValue(std::numeric_limits<int>::max());
        newKeyColumnCount_.set(std::max<int>(inPort_.getData()->getKeyColumnCount(),1));
        newKeyColumnCount_.setMaxValue(std::max<int>(inPort_.getData()->getColumnCount()-1,1));
    }
    if (dynamic_cast<const PlotData*>(inPort_.getData())) {
        pData_ = dynamic_cast<const PlotData*>(inPort_.getData());
        if (inPort_.hasChanged()) {
            newKeyColumnCount_.setMaxValue(std::numeric_limits<int>::max());
            newKeyColumnCount_.set(std::max<int>(pData_->getKeyColumnCount(),1));
            newKeyColumnCount_.setMaxValue(std::max<int>(pData_->getColumnCount()-1,1));
        }
    }
    else {
        pData_ = 0;
        LWARNING("PlotDataSelect can only handle with PlotData-Objects.");
    }
}

void PlotDataSelect::SelectValues::serialize(Serializer& s) const {
    s.serialize("KeyColumnCount",keyColumnCount);
    s.serialize("DataColumnCount",dataColumnCount);
    s.serialize("LabelVector",labelVector);
    s.serialize("PredicateVector",predicateVector);
    s.serialize("SelectedRows",selectedRows);
    s.serialize("ProjectionVector",projectionVector);
}

void PlotDataSelect::SelectValues::deserialize(Deserializer& s){
    s.deserialize("KeyColumnCount",keyColumnCount);
    s.deserialize("DataColumnCount",dataColumnCount);
    s.deserialize("LabelVector",labelVector);
    s.deserialize("PredicateVector",predicateVector);
    s.deserialize("SelectedRows",selectedRows);
    s.deserialize("ProjectionVector",projectionVector);
}

}// namespace voreen

