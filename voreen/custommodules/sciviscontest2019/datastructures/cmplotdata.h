#ifndef VRN_CMPLOTDATA_H
#define VRN_CMPLOTDATA_H

#include <string>
#include <vector>

namespace voreen {

	struct CMPlotDataRow {
		std::string name;
		std::vector<float> data;
		float min;
		float max;
	};

	class CMPlotData {
	public:
		CMPlotData();
		CMPlotData(CMPlotDataRow& dataRow);
		CMPlotData(std::vector<CMPlotDataRow>& dataVector);

		void addDataRow(CMPlotDataRow& dataRow);
		void deleteDataRow(uint32_t index);

		size_t getDataRowCount();
		size_t getDataRowLength();
		CMPlotDataRow getDataRow(uint32_t index);

	private:

		std::vector<CMPlotDataRow> data_;
	};

}
#endif
