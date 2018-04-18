/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2015 University of Muenster, Germany.                        *
 * Visualization and Computer Graphics Group <http://viscg.uni-muenster.de>        *
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

#include "viscontest2018module.h"

#include "processors/ensembledatasource.h"
#include "processors/ensemblefilter.h"
#include "processors/ensemblevolumeextractor.h"
#include "processors/ensemblesimilarityplot.h"
#include "processors/fieldparallelplotcreator.h"
#include "processors/fieldparallelplotviewer.h"
#include "processors/fieldparallelplothistogram.h"
#include "processors/mdsplot.h"
#include "processors/volumeintensityfilter.h"
#include "processors/probabilityvolumecreator.h"
#include "processors/waveheightextractor.h"

#include "io/fieldplotsave.h"
#include "io/fieldplotsource.h"
#include "io/vtivolumereader.h"

#include "properties/link/viscontest2018linkevaluatorid.h"

namespace voreen {

VisContest2018Module::VisContest2018Module(const std::string& modulePath)
    : VoreenModule(modulePath)
{
    setID("VisContest2018");
    setGuiName("VisContest2018");

    addShaderPath(getModulePath("glsl"));

    // Processors
    registerProcessor(new EnsembleDataSource());
    registerProcessor(new EnsembleFilter);

    // Plotting
    registerProcessor(new FieldParallelPlotCreator());
    registerProcessor(new EnsembleSimilarityPlot());
    registerProcessor(new FieldParallelPlotViewer());
    registerProcessor(new FieldParallelPlotHistogram());
    registerProcessor(new MDSPlot());
    registerProcessor(new VolumeIntensityFilter());
    registerProcessor(new ProbabilityVolumeCreator());

    // IO
    registerProcessor(new FieldPlotSave());
    registerProcessor(new FieldPlotSource());
    registerVolumeReader(new VTIVolumeReader());
    
    // Properties
    registerSerializableType(new LinkEvaluatorIntListId());

    // Misc
    registerProcessor(new WaveHeightExtractor());
    registerProcessor(new EnsembleVolumeExtractor());
}

} // namespace
