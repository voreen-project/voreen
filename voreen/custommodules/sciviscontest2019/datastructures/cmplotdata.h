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
