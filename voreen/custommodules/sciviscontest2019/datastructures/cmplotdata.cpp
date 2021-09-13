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