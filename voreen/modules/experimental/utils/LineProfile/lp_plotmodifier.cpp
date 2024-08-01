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

#include "lp_plotmodifier.h"

namespace voreen {

//******************************************************************************************
//                    find in plotdata
//******************************************************************************************
float LP_PlotModifier::findMax(const PlotData* pData, int column){
    if(!pData || (column >= pData->getColumnCount()))
        return 0.0f;
    return static_cast<float>(pData->getInterval(column).getRight());
}

float LP_PlotModifier::findMin(const PlotData* pData, int column){
    if(!pData || (column >= pData->getColumnCount()))
        return 0.0f;
    return static_cast<float>(pData->getInterval(column).getLeft());
}

float LP_PlotModifier::findXValue(const PlotData* pData, float percent){
    if(!pData || (percent > 100.0f) ) return 0.0f;
    if(percent > 1) percent /= 100.0f;
    Interval<plot_t> x = pData->getInterval(0);
    return static_cast<float>(x.getLeft()+((x.getRight() -x.getLeft())*percent));
}

float LP_PlotModifier::findYValue(const PlotData* pData, float percent, int column){
    if(!pData || (percent > 100.0f)) return 0.0f;
    if(percent > 1) percent /= 100.0f;
    Interval<plot_t> y = pData->getInterval(column);
    return static_cast<float>(y.getLeft()+((y.getRight() -y.getLeft())*percent));
}

//******************************************************************************************
//                    add/remove in plotdata
//******************************************************************************************
void LP_PlotModifier::addColumn(PlotData* pData, const std::vector<float>& vec, const std::string& label){
    if (static_cast<size_t>(pData->getRowsCount()) != vec.size())
        return;
    std::vector<PlotCellValue> cells;
    std::vector<PlotCellValue> row;
    PlotData* helpData = new PlotData(pData->getKeyColumnCount(),pData->getDataColumnCount()+1);
    int count = pData->getRowsCount();
    for (int i=0; i< count; ++i){ // number of rows
        cells.clear();
        row = pData->getRow(i).getCells();
        for( size_t j=0; j< row.size(); ++j){
             cells.push_back(row.at(j));
        }
        cells.push_back(vec[i]);
        helpData->insert(cells);
    }
    for (int i = 0; i < pData->getColumnCount(); ++i) {
        helpData->setColumnLabel(i,pData->getColumnLabel(i));
    }
    helpData->setColumnLabel(pData->getColumnCount(),label);
    pData->reset(helpData->getKeyColumnCount(),helpData->getDataColumnCount());
    for(int i = 0; i < helpData->getColumnCount(); i++)
        pData->setColumnLabel(i,helpData->getColumnLabel(i));
    for(int i = 0; i < count; i++)
        pData->insert(helpData->getRow(i).getCells());
    delete helpData;
}

void LP_PlotModifier::addFloatColumn(PlotData* pData, float value, const std::string& label){
    std::vector<float> vec(pData->getRowsCount(), value);
    addColumn(pData,vec,label);
}

void LP_PlotModifier::getColumn(const PlotData* pData, std::vector<float>& vec, int column) {
    if(column >= pData->getColumnCount()) return;
    vec.resize(pData->getRowsCount());
    for(int i = 0; i < pData->getRowsCount(); i++) {
        vec[i] = static_cast<float>(pData->getRow(i).getValueAt(column));
    }
}

void LP_PlotModifier::removeColumn(PlotData* pData, int pColumn){
    std::vector<PlotCellValue> cells;
    std::vector<PlotCellValue> row;
    PlotData* helpData = new PlotData(pData->getKeyColumnCount(),pData->getDataColumnCount()-1);
    int count = pData->getRowsCount();
    for (int i=0; i< count; ++i){ // number of rows
        cells.clear();
        row = pData->getRow(i).getCells();
        for( size_t j=0; j< row.size(); ++j){
            if(j != static_cast<size_t>(pColumn))
                cells.push_back(row.at(j));
        }
        helpData->insert(cells);
    }

    for (int i = 0; i < pData->getColumnCount(); ++i) {
        if(i < pColumn)
            helpData->setColumnLabel(i,pData->getColumnLabel(i));
        if(i > pColumn)
            helpData->setColumnLabel(i-1,pData->getColumnLabel(i));
    }

    pData->reset(helpData->getKeyColumnCount(),helpData->getDataColumnCount());
    for(int i = 0; i < helpData->getColumnCount(); i++)
        pData->setColumnLabel(i,helpData->getColumnLabel(i));
    for(int i = 0; i < count; i++)
        pData->insert(helpData->getRow(i).getCells());
    delete helpData;
}

void LP_PlotModifier::mergePlotData(PlotData* pInput1, const PlotData* pInput2, PlotData* pOutput){
    if(pOutput == 0)
        pOutput = pInput1;

    std::vector<PlotCellValue> cells;
    std::vector<PlotCellValue> row;
    PlotData* helpData = new PlotData(pInput1->getKeyColumnCount(),pInput1->getDataColumnCount()+pInput2->getDataColumnCount());
    int count = pInput1->getRowsCount();
    for (int i=0; i< count; ++i){ // number of rows
        cells.clear();
        row = pInput1->getRow(i).getCells();
        for( size_t j=0; j< row.size(); ++j){
             cells.push_back(row.at(j));
        }
        row = pInput2->getRow(i).getCells();
        for( size_t j=1; j< row.size(); ++j){
             cells.push_back(row.at(j));
        }
        helpData->insert(cells);
    }
    for (int i = 0; i < pInput1->getColumnCount(); ++i) {
        helpData->setColumnLabel(i,pInput1->getColumnLabel(i));
    }
    for (int i = 1; i < pInput2->getColumnCount(); ++i) {
        helpData->setColumnLabel(pInput1->getDataColumnCount()+i,pInput2->getColumnLabel(i));
    }
    pOutput->reset(helpData->getKeyColumnCount(),helpData->getDataColumnCount());
    for(int i = 0; i < helpData->getColumnCount(); i++)
        pOutput->setColumnLabel(i,helpData->getColumnLabel(i));
    for(int i = 0; i < count; i++)
        pOutput->insert(helpData->getRow(i).getCells());
    delete helpData;
}

void LP_PlotModifier::addRowVector(PlotData* pData, std::vector< std::vector<PlotCellValue> > & vec){
    if(static_cast<size_t>(pData->getRowsCount()) != vec.size())
        return;
    //save old label
    std::vector<std::string> sVec(pData->getColumnCount());
    for(int i = 0; i < pData->getColumnCount(); i++)
        sVec[i] = pData->getColumnLabel(i);
    //add old rows to the new ones
    for(size_t i = 0; i < vec.size(); i++)
        vec[i].insert(vec[i].begin(),pData->getRow(static_cast<int>(i)).getCells().begin(),pData->getRow(static_cast<int>(i)).getCells().end());
    //reset data and insert the rows
    pData->reset(1, static_cast<int>(vec[0].size()-1));
    for(size_t i = 0; i < vec.size(); i++)
        pData->insert(vec[i]);
    //retore old label
    for(size_t i = 0; i < sVec.size(); i++)
        pData->setColumnLabel(static_cast<int>(i),sVec[i]);
}

} // namespace voreen
