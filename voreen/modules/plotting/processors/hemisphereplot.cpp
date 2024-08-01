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

#include "hemisphereplot.h"

#include "voreen/core/interaction/camerainteractionhandler.h"
#include "../interaction/plotcamerainteractionhandler.h"
#include "../datastructures/plotrow.h"

namespace voreen {

const std::string HemispherePlot::loggerCat_("voreen.HemispherePlot");

HemispherePlot::HemispherePlot()
    : SurfacePlot()
{
}

Processor* HemispherePlot::create() const {
    return new HemispherePlot();
}

void HemispherePlot::initialize() {
    SurfacePlot::initialize();
}

void HemispherePlot::deinitialize() {
    SurfacePlot::deinitialize();
}

void HemispherePlot::readFromInport() {
    if (dynamic_cast<const PlotData*>(inport_.getData())->getRowsCount() > 0) {
        ++omitDelaunayTriangulation_;
        const PlotData* tempData = dynamic_cast<const PlotData*>(inport_.getData());
        data_.reset(tempData->getKeyColumnCount(), tempData->getDataColumnCount() + 1);

        int numRowsPerColumn = 1;
        plot_t xValue = tempData->getRowsBegin()->getValueAt(0);
        for (std::vector<PlotRowValue>::const_iterator it = ++tempData->getRowsBegin(); it < tempData->getRowsEnd() && xValue == it->getValueAt(0); ++it) {
            ++numRowsPerColumn;
        }

        // Convert sphere coordinates to cartesian coordinates:
        // column   usage
        //   0      Phi (azimuth angle)
        //   1      Theta (inclination angle)
        //   2      length
        //   3...   passthrough
        if (tempData->getColumnCount() >= 3) {
            data_.setColumnLabel(0, "x");
            data_.setColumnLabel(1, "y");
            data_.setColumnLabel(2, "z");
            data_.setColumnLabel(3, "#Samples");
            for (int i = 3; i < tempData->getColumnCount(); ++i)
                data_.setColumnLabel(i+1, tempData->getColumnLabel(i));

            Interval<plot_t> range = tempData->getInterval(2);

            for (std::vector<PlotRowValue>::const_iterator it = tempData->getRowsBegin(); it != tempData->getRowsEnd(); ++it) {
                plot_t phi = it->getValueAt(0);
                plot_t theta = it->getValueAt(1);
                plot_t r = .25 + .25 * ((it->getValueAt(2) - range.getLeft()) / range.size());

                plot_t x = r * sin(theta) * cos(phi);
                plot_t y = r * sin(theta) * sin(phi);
                plot_t z = -r * cos(theta);


                std::vector<PlotCellValue> cells;
                cells.push_back(PlotCellValue(x));
                cells.push_back(PlotCellValue(y));
                cells.push_back(PlotCellValue(z));
                cells.push_back(it->getCellAt(2));
                for (int i = 3; i < tempData->getColumnCount(); ++i)
                    cells.push_back(it->getCellAt(i));

                data_.insert(cells);
            }

            // add pole to hemisphere
            int numSamples = 0;
            int currentRow = tempData->getRowsCount()-1;
            tgt::dvec4 sum(0.0);
            while (currentRow >= 0) {
                const PlotRowValue& row = data_.getRow(currentRow);
                sum.x += row.getValueAt(0);
                sum.y += row.getValueAt(1);
                sum.z += row.getValueAt(2);
                sum.w += row.getValueAt(3);
                ++numSamples;
                currentRow -= numRowsPerColumn;
            }
            sum /= static_cast<double>(numSamples);

            std::vector<PlotCellValue> cells;
            cells.push_back(PlotCellValue(sum.x));
            cells.push_back(PlotCellValue(sum.y));
            cells.push_back(PlotCellValue(sum.z));
            cells.push_back(PlotCellValue(sum.w));
            for (int i = 3; i < tempData->getColumnCount(); ++i)
                cells.push_back(0);

            data_.insert(cells);

        }

        dataProp_.set(&data_);
        plotEntitiesProp_.setPlotData(&data_);
        selectionProp_.setPlotData(&data_);
        inportHasPlotFunction_ = false;
        discreteStep_.setVisibleFlag(false);

        // enforce re-triangulation
        --omitDelaunayTriangulation_;
        currentIndexX_ = -1;
        currentIndexY_ = -1;
        generateDelaunay();
        plotPickingManager_.setColumnCount(data_.getColumnCount());
    }
/*    else if (dynamic_cast<PlotFunction*>(inport_.getData()) && inport_.getData()->getKeyColumnCount() == 2) {
        ++omitDelaunayTriangulation_;
        function_ = *dynamic_cast<PlotFunction*>(inport_.getData());
        inportHasPlotFunction_ = true;
        selectDataFromFunction();
        discreteStep_.setVisibleFlag(true);

        // enforce re-triangulation
        --omitDelaunayTriangulation_;
        currentIndexX_ = -1;
        currentIndexY_ = -1;
        generateDelaunay();
    }*/
    else {
        if (!dynamic_cast<const PlotData*>(inport_.getData()))
            LWARNINGC("HemispherePlot", "HemispherePlot can only handle PlotData objects with two key columns");
        ++omitDelaunayTriangulation_;
        data_ = PlotData(0,0);
        dataProp_.set(&data_);
        plotEntitiesProp_.setPlotData(&data_);
        selectionProp_.setPlotData(&data_);
        plotPickingManager_.setColumnCount(data_.getColumnCount());
        voronoiRegions_.clear();
        triangleEdgeIndices_.clear();
        --omitDelaunayTriangulation_;
    }
}

void HemispherePlot::generateDelaunay() {
    if (omitDelaunayTriangulation_  == 0 && plotEntitiesProp_.dataValid()) {
        if (plotEntitiesProp_.getXColumnIndex() != currentIndexX_ || plotEntitiesProp_.getYColumnIndex() != currentIndexY_) {
            currentIndexX_ = plotEntitiesProp_.getXColumnIndex();
            currentIndexY_ = plotEntitiesProp_.getYColumnIndex();

            // we assume we have an uniform rectliniar grid, calculate the number of rows per column:
            const PlotData* tempData = dynamic_cast<const PlotData*>(inport_.getData());
            if (tempData == 0)
                return;
            int numRowsPerColumn = 1;
            plot_t xValue = tempData->getRowsBegin()->getValueAt(currentIndexX_);
            for (std::vector<PlotRowValue>::const_iterator it = ++tempData->getRowsBegin(); it < tempData->getRowsEnd() && xValue == it->getValueAt(currentIndexX_); ++it) {
                ++numRowsPerColumn;
            }

            // generate triangles
            triangleEdgeIndices_.clear();
            int numCols = data_.getRowsCount() / numRowsPerColumn;
            for (int col = 0; col < numCols; ++col) {
                for (int row = 0; row < numRowsPerColumn - 1; ++row) {
                    triangleEdgeIndices_.push_back(col*numRowsPerColumn + row);
                    triangleEdgeIndices_.push_back(col*numRowsPerColumn + row + 1);
                    triangleEdgeIndices_.push_back(((col+1) % numCols)*numRowsPerColumn + row);

                    triangleEdgeIndices_.push_back(col*numRowsPerColumn + row + 1);
                    triangleEdgeIndices_.push_back(((col+1) % numCols)*numRowsPerColumn + row + 1);
                    triangleEdgeIndices_.push_back(((col+1) % numCols)*numRowsPerColumn + row);
                }
            }
            int row = numRowsPerColumn-1;
            for (int col = 0; col < numCols; ++col) {
                triangleEdgeIndices_.push_back(col*numRowsPerColumn + row);
                triangleEdgeIndices_.push_back(((col+1) % numCols)*numRowsPerColumn + row);
                triangleEdgeIndices_.push_back(numRowsPerColumn * numCols);
            }
            // no voronoi regions to generate
            voronoiRegions_.clear();
        }
    }
}

void HemispherePlot::calcDomains() {
    if (plotEntitiesProp_.dataValid() && !inportHasPlotFunction_) {
        Interval<plot_t> xDomain(-.5, .5);
        Interval<plot_t> yDomain(-.5, .5);
        Interval<plot_t> zDomain(-.25, .75);
        selectionProp_.setBaseZoomState(PlotZoomState(xDomain, yDomain, zDomain));
    }
}

void HemispherePlot::setPlotStatus() {
    voreen::SurfacePlot::setPlotStatus();
    //plotLib_->setCenterAxesFlag(true);
}

} // namespace voreen
