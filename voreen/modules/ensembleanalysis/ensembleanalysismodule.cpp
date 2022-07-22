/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2021 University of Muenster, Germany,                        *
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

#include "io/ensembledatasource.h"
#include "io/parallelcoordinatessave.h"
#include "io/parallelcoordinatessource.h"
#include "io/similaritymatrixsave.h"
#include "io/similaritymatrixsource.h"

#include "processors/ensemblechannelmerger.h"
#include "processors/ensemblecolor.h"
#include "processors/ensemblefilter.h"
#include "processors/ensembleinformation.h"
#include "processors/ensemblemeancreator.h"
#include "processors/ensembletimesteplinker.h"
#include "processors/ensemblevarianceanalysis.h"
#include "processors/ensemblevolumeextractor.h"
#include "processors/metadataadder.h"
#include "processors/parallelcoordinatesaxescreator.h"
#include "processors/parallelcoordinatesviewer.h"
#include "processors/parallelcoordinatesvoxelselection.h"
#include "processors/similaritymatrixcombine.h"
#include "processors/similaritymatrixcreator.h"
#include "processors/similarityplot.h"

#include "properties/parallelcoordinatessectionsproperty.h"

namespace voreen {

EnsembleAnalysisModule* EnsembleAnalysisModule::instance_ = nullptr;

EnsembleAnalysisModule::EnsembleAnalysisModule(const std::string& modulePath)
    : VoreenModule(modulePath)
    , forceDiskRepresentation_("forceDiskRepresentation", "Force Disk Representation (experimental)", false)
{
    setID("Ensemble Analysis");
    setGuiName("Ensemble Analysis");

    instance_ = this;
    addProperty(forceDiskRepresentation_);

    // Processors
    registerProcessor(new EnsembleChannelMerger());
    registerProcessor(new EnsembleColor());
    registerProcessor(new EnsembleDataSource());
    registerProcessor(new EnsembleFilter);
    registerProcessor(new EnsembleMeanCreator());

    // Plotting
    registerProcessor(new EnsembleVarianceAnalysis());
    registerProcessor(new SimilarityMatrixCombine());
    registerProcessor(new ParallelCoordinatesAxesCreator());
    registerProcessor(new ParallelCoordinatesViewer());
    registerProcessor(new ParallelCoordinatesVoxelSelection());
    registerProcessor(new SimilarityMatrixCreator());
    registerProcessor(new SimilarityPlot());

    // IO
    registerProcessor(new ParallelCoordinatesSave());
    registerProcessor(new ParallelCoordinatesSource());
    registerProcessor(new SimilarityMatrixSave());
    registerProcessor(new SimilarityMatrixSource());

    // Misc
    registerProcessor(new EnsembleInformation());
    registerProcessor(new EnsembleTimeStepLinker());
    registerProcessor(new EnsembleVolumeExtractor());
    registerProcessor(new MetaDataAdder());

    // Link evaluators
    registerSerializableType(new LinkEvaluatorParallelCoordinatesSectionsId());
}

void EnsembleAnalysisModule::setForceDiskRepresentation(bool enabled) {
    forceDiskRepresentation_.set(enabled);
}

bool EnsembleAnalysisModule::getForceDiskRepresentation() const {
    return forceDiskRepresentation_.get();
}

EnsembleAnalysisModule* EnsembleAnalysisModule::getInstance() {
    return instance_;
}

} // namespace
