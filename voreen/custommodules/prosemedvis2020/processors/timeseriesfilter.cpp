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

#include "timeseriesfilter.h"

#include "../modules/plotting/datastructures/plotcell.h"
#include "../modules/plotting/datastructures/plotrow.h"
#include <vector>
#include <string>


namespace voreen {
    
TimeseriesFilter::TimeseriesFilter()
    : Processor()
    , inport_(Port::INPORT, "inport", "Inport")
    , outport_(Port::OUTPORT, "outport", "PlotData Outport")
    , regionProp_("region_prop", "Selected Region")
    , tracerProp_("tracer_prop", "Selected Tracer")
    , tracers_()
{
    // Ports
    addPort(inport_);
    addPort(outport_);

    // Properties
    addProperty(regionProp_);
    addProperty(tracerProp_);
}

Processor* TimeseriesFilter::create() const {
    return new TimeseriesFilter();
}

std::string TimeseriesFilter::getClassName() const {
    return "TimeseriesFilter";
}

std::string TimeseriesFilter::getCategory() const {
    return "Ensemble Processing";
}

void TimeseriesFilter::setDescriptions() {
    setDescription("Processor which passes all columns of a given PlotData-Object, that start with a given prefix.");
}

void TimeseriesFilter::process() {

    if (!inport_.isConnected() || inport_.getData() == nullptr) {
        outport_.setData(nullptr, true);
        return;
        
    }
    
    // If cannot be cast to PlotData, return null
    if (!dynamic_cast<const PlotData*>(inport_.getData())) {
        LWARNING("TimeseriesFilter: Cannot cast input to PlotData. Returning copy of input instead!");
        outport_.setData(nullptr, true);
        return;
    }
    const PlotData* oldPlot = dynamic_cast<const PlotData*>(inport_.getData());

    if(inport_.hasChanged()){
        // Update the tracer-Property
        tracers_ = getTracerNames();
        tracerProp_.reset();
        for (std::string tracer : tracers_) {
            tracerProp_.addRow(tracer);
        }
    }


    // If region = "" -> Return copy of input
    if (regionProp_.get() == "") {
        PlotData* newData = new PlotData(*oldPlot);
        outport_.setData(newData, true);
        return;
    }


    // Get indices of all columns, that start with the selected prefix:: or are a key-cloumn
    // Column-format: REGION::ID::TRACERNAME_MOUSEID
    std::vector<int> selectedCols;
    for (int c = 0; c < oldPlot->getColumnCount(); ++c) {
        const std::string& colName = oldPlot->getColumnLabel(c);
        if (c < oldPlot->getKeyColumnCount()) {
            // Column is a key-column
            selectedCols.push_back(c);
        } else if (colName.rfind(regionProp_.get() + "::", 0) == 0) {
            // Column matches selected region
            if (tracerProp_.getSelectedRowIndices().size() == 0) {
                // No tracer selected --> Show all tracers
                selectedCols.push_back(c);
            }else {
                // Check if column matches selected tracer
                bool found = false;
                for (int idx : tracerProp_.getSelectedRowIndices()) {
                    if (idx < 0 || idx >= tracers_.size()) continue;
                    const std::string tracer = tracers_.at(idx);
                    if (colName.find("::" + tracer + "_") != std::string::npos) {
                        found = true;
                        break;
                    }
                }
                if (found) selectedCols.push_back(c);
            }
        }
    }

    
    // Setup plot 
    const int numDataCols = selectedCols.size() - oldPlot->getKeyColumnCount();
    PlotData* newData = new PlotData(oldPlot->getKeyColumnCount(), numDataCols);

    // Set column-names
    for (int i = 0; i < selectedCols.size(); ++i) {
        newData->setColumnLabel(i, oldPlot->getColumnLabel(selectedCols.at(i)));
    }

    // Add Rows
    for (auto it = oldPlot->getRowsBegin(); it != oldPlot->getRowsEnd(); ++it) {
        std::vector<PlotCellValue> row;
        for (int idx : selectedCols) {
            const PlotRowValue& plotRow = *it;
            row.push_back(PlotCellValue(plotRow.getValueAt(idx)));
        }
        newData->insert(row);
    }
    
    outport_.setData(newData, true);
    
}

std::vector<std::string> TimeseriesFilter::getTracerNames() const {
    if (!inport_.hasData()) {
        return std::vector<std::string>();
    }

    const PlotBase* plotData = inport_.getData();
    std::set<std::string> names;

    for (int i = plotData->getKeyColumnCount(); i < plotData->getColumnCount(); ++i) {
        // Get tracername from columname (REGION::ID::TRACERNAME_MOUSEID)
        std::string colLabel = plotData->getColumnLabel(i);

        // Get begin pos (Pos after second occurence of ::)
        size_t beginPos = colLabel.find("::");
        if (beginPos == std::string::npos || beginPos + 2 >= colLabel.size()) {
            continue;
        }
        beginPos = colLabel.find("::", beginPos + 2);
        if (beginPos == std::string::npos || beginPos + 2 >= colLabel.size()) {
            continue;
        }
        beginPos += 2;

        // Get end-Pos (Position before last occurence of '_')
        size_t endPos = colLabel.find_last_of("_");
        if (endPos == std::string::npos || endPos <= beginPos) {
            continue;
        }

        std::string tracerName = colLabel.substr(beginPos, endPos - beginPos);
        names.insert(tracerName);
    }

    return std::vector<std::string>(names.begin(), names.end());
}

} // namespace
