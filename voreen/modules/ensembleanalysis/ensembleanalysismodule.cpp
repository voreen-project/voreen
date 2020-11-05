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

#include "ensembleanalysismodule.h"

#include "io/similaritymatrixsave.h"
#include "io/similaritymatrixsource.h"

#include "processors/connectedcomponentselector.h"
#include "processors/ensembledatasource.h"
#include "processors/ensemblefilter.h"
#include "processors/ensemblevolumeextractor.h"
#include "processors/localsimilarityanalysis.h"
#include "processors/metadataadder.h"
#include "processors/parallelcoordinatesaxescreator.h"
#include "processors/parallelcoordinatessource.h"
#include "processors/parallelcoordinatesviewer.h"
#include "processors/parallelcoordinatesvoxelselection.h"
#include "processors/referencevolumecreator.h"
#include "processors/similaritymatrixcombine.h"
#include "processors/similaritymatrixcreator.h"
#include "processors/similarityplot.h"
#include "processors/volumelistmerger.h"
#ifdef VRN_MODULE_HDF5
#include "processors/volumemerger.h"
#endif

#include "properties/parallelcoordinatessectionsproperty.h"

namespace voreen {

EnsembleAnalysisModule::EnsembleAnalysisModule(const std::string& modulePath)
    : VoreenModule(modulePath)

{
    setID("Ensemble Analysis");
    setGuiName("Ensemble Analysis");

    // Processors
    registerProcessor(new ConnectedComponentSelector());
    registerProcessor(new EnsembleDataSource());
    registerProcessor(new EnsembleFilter);
    registerProcessor(new ReferenceVolumeCreator());

    // Plotting
    registerProcessor(new LocalSimilarityAnalysis());
    registerProcessor(new SimilarityMatrixCombine());
    registerProcessor(new ParallelCoordinatesAxesCreator());
    registerProcessor(new ParallelCoordinatesSource());
    registerProcessor(new ParallelCoordinatesViewer());
    registerProcessor(new ParallelCoordinatesVoxelSelection());
    registerProcessor(new SimilarityMatrixCreator());
    registerProcessor(new SimilarityPlot());

    // IO
    registerProcessor(new SimilarityMatrixSave());
    registerProcessor(new SimilarityMatrixSource());

    // Misc
    registerProcessor(new EnsembleVolumeExtractor());
    registerProcessor(new MetaDataAdder());
    registerProcessor(new VolumeListMerger());
#ifdef VRN_MODULE_HDF5
    registerProcessor(new VolumeMerger());
#endif

    // Link evaluators
    registerSerializableType(new LinkEvaluatorParallelCoordinatesSectionsId());
}

} // namespace
