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

#include "plottingmodule.h"

#include "processors/plotdatasource.h"
#include "processors/plotfunctionsource.h"
#include "processors/lineplot.h"
#include "processors/hemisphereplot.h"
#include "processors/barplot.h"
#include "processors/scatterplot.h"
#include "processors/surfaceplot.h"
#include "processors/plotdataselect.h"
#include "processors/plotdatagroup.h"
#include "processors/plotdataexport.h"
#include "processors/plotdataexporttext.h"
#include "processors/plotdatafitfunction.h"
#include "processors/plotdatamerge.h"
#include "processors/plotfunctiondiscret.h"
#include "processors/imageanalyzer.h"

#include "properties/link/linkevaluatorplotselection.h"
#include "properties/link/linkevaluatorcolormapid.h"
#include "properties/link/linkevaluatorplotentitiesid.h"

#include "datastructures/aggregationfunction.h"

#include "utils/plotpredicatefactory.h"

#include "voreen/core/voreenapplication.h"

namespace voreen {

PlottingModule::PlottingModule(const std::string& modulePath)
    : VoreenModule(modulePath)
{
    setID("Plotting");
    setGuiName("Plotting");

    addShaderPath(getModulePath("glsl"));

    registerSerializableType(new PlotDataSource());
    registerSerializableType(new PlotFunctionSource());
    registerSerializableType(new BarPlot());
    registerSerializableType(new HemispherePlot());
    registerSerializableType(new ImageAnalyzer());
    registerSerializableType(new LinePlot());
    registerSerializableType(new PlotDataExport());
    registerSerializableType(new PlotDataExportText());
    registerSerializableType(new PlotDataFitFunction());
    registerSerializableType(new PlotDataGroup());
    registerSerializableType(new PlotDataMerge());
    registerSerializableType(new PlotDataSelect());
    registerSerializableType(new PlotFunctionDiscret());
    registerSerializableType(new ScatterPlot());
#ifndef WIN32
    // Surface Plot uses the included triangle library which will crash on win32 systems.
    // As long as this isn't resolved (e.g. by using another library), we don't include this processor.
    registerSerializableType(new SurfacePlot());
#endif

    registerSerializableType(new LinkEvaluatorColorMapId());
    registerSerializableType(new LinkEvaluatorPlotEntitiesId());
    registerSerializableType(new LinkEvaluatorPlotSelection());

    registerSerializableType(new AggregationFunctionMin());
    registerSerializableType(new AggregationFunctionMinHistogram());
    registerSerializableType(new AggregationFunctionMax());
    registerSerializableType(new AggregationFunctionMaxHistogram());
    registerSerializableType(new AggregationFunctionSum());
    registerSerializableType(new AggregationFunctionSumHistogram());
    registerSerializableType(new AggregationFunctionCount());
    registerSerializableType(new AggregationFunctionCountHistogram());
    registerSerializableType(new AggregationFunctionMedian());
    registerSerializableType(new AggregationFunctionMedianHistogram());
    registerSerializableType(new AggregationFunctionAverage());
    registerSerializableType(new AggregationFunctionGeometricAverage());
    registerSerializableType(new AggregationFunctionHarmonicAverage());
    registerSerializableType(new AggregationFunctionStandardDeviation());
    registerSerializableType(new AggregationFunctionVariance());
    registerSerializableType(new AggregationFunctionMode());

    // TODO: convert PlotPredicateFactory into VoreenSerializableObjectFactory
    if (VoreenApplication::app()) {
        VoreenApplication::app()->registerSerializerFactory(PlotPredicateFactory::getInstance());
    }
}

} // namespace
