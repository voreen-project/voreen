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

#include "plotdatainserter.h"

#include "tgt/assert.h"

namespace voreen {

bool PlotDataInserter::insert(PlotData& pData, const PlotCellValue& arg1) {
    std::vector<PlotCellValue> cells;
    cells.push_back(arg1);
    return pData.insert(cells);
}

bool PlotDataInserter::insert(PlotData& pData, const PlotCellValue& arg1, const PlotCellValue& arg2) {
    std::vector<PlotCellValue> cells;
    cells.push_back(arg1);
    cells.push_back(arg2);
    return pData.insert(cells);
}

bool PlotDataInserter::insert(PlotData& pData,
                              const PlotCellValue& arg1,
                              const PlotCellValue& arg2,
                              const PlotCellValue& arg3) {
    std::vector<PlotCellValue> cells;
    cells.push_back(arg1);
    cells.push_back(arg2);
    cells.push_back(arg3);
    return pData.insert(cells);
}

bool PlotDataInserter::insert(PlotData& pData,
                              const PlotCellValue& arg1,
                              const PlotCellValue& arg2,
                              const PlotCellValue& arg3,
                              const PlotCellValue& arg4) {
    std::vector<PlotCellValue> cells;
    cells.push_back(arg1);
    cells.push_back(arg2);
    cells.push_back(arg3);
    cells.push_back(arg4);
    return pData.insert(cells);
}

bool PlotDataInserter::insert(PlotData& pData,
                              const PlotCellValue& arg1,
                              const PlotCellValue& arg2,
                              const PlotCellValue& arg3,
                              const PlotCellValue& arg4,
                              const PlotCellValue& arg5) {
    std::vector<PlotCellValue> cells;
    cells.push_back(arg1);
    cells.push_back(arg2);
    cells.push_back(arg3);
    cells.push_back(arg4);
    cells.push_back(arg5);
    return pData.insert(cells);
}

bool PlotDataInserter::insert(PlotData& pData,
                              const PlotCellValue& arg1,
                              const PlotCellValue& arg2,
                              const PlotCellValue& arg3,
                              const PlotCellValue& arg4,
                              const PlotCellValue& arg5,
                              const PlotCellValue& arg6) {
    std::vector<PlotCellValue> cells;
    cells.push_back(arg1);
    cells.push_back(arg2);
    cells.push_back(arg3);
    cells.push_back(arg4);
    cells.push_back(arg5);
    cells.push_back(arg6);
    return pData.insert(cells);
}

bool PlotDataInserter::insert(PlotData& pData,
                              const PlotCellValue& arg1,
                              const PlotCellValue& arg2,
                              const PlotCellValue& arg3,
                              const PlotCellValue& arg4,
                              const PlotCellValue& arg5,
                              const PlotCellValue& arg6,
                              const PlotCellValue& arg7) {
    std::vector<PlotCellValue> cells;
    cells.push_back(arg1);
    cells.push_back(arg2);
    cells.push_back(arg3);
    cells.push_back(arg4);
    cells.push_back(arg5);
    cells.push_back(arg6);
    cells.push_back(arg7);
    return pData.insert(cells);
}

bool PlotDataInserter::insert(PlotData& pData,
                              const PlotCellValue& arg1,
                              const PlotCellValue& arg2,
                              const PlotCellValue& arg3,
                              const PlotCellValue& arg4,
                              const PlotCellValue& arg5,
                              const PlotCellValue& arg6,
                              const PlotCellValue& arg7,
                              const PlotCellValue& arg8) {
    std::vector<PlotCellValue> cells;
    cells.push_back(arg1);
    cells.push_back(arg2);
    cells.push_back(arg3);
    cells.push_back(arg4);
    cells.push_back(arg5);
    cells.push_back(arg6);
    cells.push_back(arg7);
    cells.push_back(arg8);
    return pData.insert(cells);
}

bool PlotDataInserter::insert(PlotData& pData,
                              const PlotCellValue& arg1,
                              const PlotCellValue& arg2,
                              const PlotCellValue& arg3,
                              const PlotCellValue& arg4,
                              const PlotCellValue& arg5,
                              const PlotCellValue& arg6,
                              const PlotCellValue& arg7,
                              const PlotCellValue& arg8,
                              const PlotCellValue& arg9) {
    std::vector<PlotCellValue> cells;
    cells.push_back(arg1);
    cells.push_back(arg2);
    cells.push_back(arg3);
    cells.push_back(arg4);
    cells.push_back(arg5);
    cells.push_back(arg6);
    cells.push_back(arg7);
    cells.push_back(arg8);
    cells.push_back(arg9);
    return pData.insert(cells);
}

bool PlotDataInserter::insert(PlotData& pData,
                              const PlotCellValue& arg1,
                              const PlotCellValue& arg2,
                              const PlotCellValue& arg3,
                              const PlotCellValue& arg4,
                              const PlotCellValue& arg5,
                              const PlotCellValue& arg6,
                              const PlotCellValue& arg7,
                              const PlotCellValue& arg8,
                              const PlotCellValue& arg9,
                              const PlotCellValue& arg10) {
    std::vector<PlotCellValue> cells;
    cells.push_back(arg1);
    cells.push_back(arg2);
    cells.push_back(arg3);
    cells.push_back(arg4);
    cells.push_back(arg5);
    cells.push_back(arg6);
    cells.push_back(arg7);
    cells.push_back(arg8);
    cells.push_back(arg9);
    cells.push_back(arg10);
    return pData.insert(cells);
}

} // namespace
