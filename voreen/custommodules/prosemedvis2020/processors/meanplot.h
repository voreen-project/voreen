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

#ifndef VRN_MINEPLOT_H
#define VRN_MEANPLOT_H

#include "modules/plotting/processors/plotprocessor.h"

#include "voreen/core/properties/floatproperty.h"

namespace voreen {

	class VRN_CORE_API MeanPlot : public PlotProcessor {

	public:
		MeanPlot();
		virtual Processor* create() const;

		virtual std::string getClassName() const { return "MeanPlot"; }
		virtual CodeState getCodeState() const { return CODE_STATE_STABLE; }

	protected:
		virtual void setDescriptions() {
			setDescription("This processor is able to plot lines showing mean and standard deviation");
		}

	private:
		// cache mean and std
		class MeanData
		{
		public:
			plot_t* mean_; // mean
			plot_t* std_; // standard deviation
			plot_t std_sum_; // range between + and - standard deviation
			std::string name_;
		};

		// sort MeanData by std_sum_ 
		struct meandata_sort
		{
			inline bool operator() (const MeanData& mean1, const MeanData& mean2) {
				return (mean1.std_sum_ > mean2.std_sum_);
			}
		};

		// sort MeanData
		bool sortMeanData(MeanData i, MeanData j) { return (i.std_sum_ > j.std_sum_); }

	private:
		const int COLUMNS_PER_TRACER_ = 5; // Five columns are used: -std, mean - std, mean, mean + std, +std 
		// inherited methods
		virtual void render();
		virtual void renderData();
		virtual void renderAxes();
		virtual void setPlotStatus();
		virtual void readFromInport();
		virtual void calcDomains();
		virtual void createPlotLabels();

		// create line labels
		void createLineLabels();
		// get color from colormap - every #COLUMNS_PER_TRACER_ columns share same color
		tgt::Color getColor(int index);
		// shorten tracername
		std::string getTracerName(const PlotData* data);

		// calculate mean and std
		MeanData readFromInport(const PlotData* p);
		bool isMeanColumn(int loop_index, PlotEntitySettings it, bool checkAlsoStdColumn) {

			int main_column_index = it.getMainColumnIndex();
			int tracer_column_start = getTracerInLoop(loop_index, true);
			return  main_column_index == tracer_column_start + 3 || (checkAlsoStdColumn ? (main_column_index == tracer_column_start + 2 || main_column_index == tracer_column_start + 4) : false);
		}
		int getTracerInLoop(int loop_index, bool multipleWithColumns_per_Tracer) {
			int tracer = loop_index / COLUMNS_PER_TRACER_;
			return multipleWithColumns_per_Tracer ? COLUMNS_PER_TRACER_ * tracer : tracer;
		}

		// properties
		IntProperty heatMap_;
		FloatProperty lineWidth_;
		FloatProperty pointSize_;
		BoolProperty renderLineLabel_;
		BoolProperty colorArea_; // Used to enable and disable areas
		BoolProperty showStdLabels_; // Used to enable and disable labels showing std
	};

}   //namespace

#endif // VRN_MEANPLOT_H
