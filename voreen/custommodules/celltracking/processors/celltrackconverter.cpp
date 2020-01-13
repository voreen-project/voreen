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

#include "celltrackconverter.h"

#include "modules/plotting/datastructures/plotdata.h"
#include "modules/plotting/datastructures/plotrow.h"
#include "modules/plotting/datastructures/plotcell.h"

namespace voreen {

const std::string CellTrackConverter::loggerCat_("voreen.celltracking.CellTrackConverter");

CellTrackConverter::CellTrackConverter()
    : Processor()
    , inport_(Port::INPORT, "volumecollection", "VolumeList Input", false)
    , outport_(Port::OUTPORT, "textport", "Elapsed Time")
    , minimumLength_("minLength", "Minimum Track Length (in steps)", 1, 1, 1, Processor::INVALID_RESULT, NumericProperty<int>::DYNAMIC)
    , columnRange_("timeInterval", "Time step interval")
{
    addPort(inport_);
    addPort(outport_);

    addProperty(minimumLength_);
    //addProperty(columnRange_);
}

Processor* CellTrackConverter::create() const {
    return new CellTrackConverter();
}

void CellTrackConverter::process() {
    outport_.clear();

    const PlotData* data = dynamic_cast<const PlotData*>(inport_.getData());
    if (!data) {
        LERROR("No valid plot data");
        return;
    }

    // no output for empty plot data
    if (data->rowsEmpty())
        return;

    // if this can be not be vec3 data: abort
    if ((data->getColumnCount() % 3) != 0) {
        LERROR("Plot data does not contain 3 coordinate dimensions for each time step!");
        return;
    }

    // create output geometry
    PointSegmentListGeometryVec3* segmentList = new PointSegmentListGeometryVec3();

    for (int i = 0; i < data->getRowsCount(); ++i) {
        // get the row (i.e., track)
        const PlotRowValue& row = data->getRow(i);
        const std::vector<PlotCellValue>& values = row.getCells();
        // we have to iterate over the coordinate dimensions
        size_t dim = 0;
        tgt::vec3 currentValue = tgt::vec3::zero;
        std::vector<tgt::vec3> track;
        for (auto it = values.begin(); it != values.end(); ++it) {
            currentValue.elem[dim] = static_cast<float>(it->getValue());

            if ((dim == 2) && (currentValue != tgt::vec3(-1.f))) {
                track.push_back(currentValue);
                currentValue = tgt::vec3::zero;
            }
            dim = (dim + 1 ) % 3;
        }

        if (!track.empty() && track.size() >= static_cast<size_t>(minimumLength_.get()))
            segmentList->addSegment(track);
    }

    outport_.setData(segmentList, true);
}

void CellTrackConverter::adjustPropertiesToInput() {
    if (!inport_.getData())
        return;
    const PlotData* data = dynamic_cast<const PlotData*>(inport_.getData());

    if (!data)
        return;

    // only adjust if tracks are available 
    if (data->getRowsCount() > 0) {
        columnRange_.setMaxValue(data->getColumnCount() / 3 - 1);
        minimumLength_.setMaxValue(data->getColumnCount() / 3 - 1);
    }
}

} // namespace
