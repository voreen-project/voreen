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

#include "cmplotdata.h"

#include <algorithm>

namespace voreen {

	CMPlotData::CMPlotData()
	{
		data_ = std::vector<CMPlotDataRow>();
	}

	CMPlotData::CMPlotData(CMPlotDataRow& dataRow)
	{
		data_ = std::vector<CMPlotDataRow>();
		data_.push_back(dataRow);
		
	}

	CMPlotData::CMPlotData(std::vector<CMPlotDataRow>& dataVector)
	{
		data_ = std::vector<CMPlotDataRow>(dataVector);
	}

	void CMPlotData::addDataRow(CMPlotDataRow& dataRow)
	{
		CMPlotDataRow newRow = { dataRow.name, dataRow.data, dataRow.min, dataRow.max };
		data_.push_back(newRow);

	}

	void CMPlotData::deleteDataRow(uint32_t index)
	{
		if (index < data_.size()) {
			data_.erase(data_.begin() + index);
		}
	}

	size_t CMPlotData::getDataRowCount()
	{
		return data_.size();
	}

	size_t CMPlotData::getDataRowLength()
	{
		size_t max = 0;
		for (CMPlotDataRow row : data_) {
			max = std::max(max, row.data.size());
		}
		return max;
	}

	CMPlotDataRow CMPlotData::getDataRow(uint32_t index)
	{
		if (index < data_.size()) {
			return data_[index];
		}
	}

}