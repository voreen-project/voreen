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

#include "processors/connectedcomponentselector.h"
#include "processors/ensembledatasource.h"
#include "processors/ensemblefilter.h"
#include "processors/ensemblevolumeextractor.h"
#include "processors/fieldparallelplotcreator.h"
#include "processors/fieldparallelplotviewer.h"
#include "processors/fieldparallelplothistogram.h"
#include "processors/localsimilarityanalysis.h"
#include "processors/physicalclippinglinker.h"
#include "processors/referencevolumecreator.h"
#include "processors/similaritydatavolume.h"
#include "processors/similaritymatrixcombine.h"
#include "processors/similaritymatrixcreator.h"
#include "processors/similarityplot.h"
#include "processors/volumelistmerger.h"
#include "processors/volumemerger.h"

#include "io/fieldplotsave.h"
#include "io/fieldplotsource.h"
#include "io/similaritymatrixsave.h"
#include "io/similaritymatrixsource.h"
#ifdef VRN_USE_VTK
#include "io/netcdfvolumereader.h"
#include "io/niftivolumewriter.h"
#include "io/vtivolumereader.h"
#include "io/vtivolumewriter.h"
#include "io/vtmvolumereader.h"
#endif

#include "custommodules/ensembleanalysis/properties/link/ensembleanalysislinkevaluatorid.h"

namespace voreen {

EnsembleAnalysisModule* EnsembleAnalysisModule::instance_ = nullptr;

EnsembleAnalysisModule::EnsembleAnalysisModule(const std::string& modulePath)
    : VoreenModule(modulePath)
    , forceDiskRepresentation_("forceDiskRepresentation", "Force Disk Representation", false)
{
    setID("EnsembleAnalysis");
    setGuiName("EnsembleAnalysis");
    instance_ = this;

    addShaderPath(getModulePath("glsl"));
    addProperty(forceDiskRepresentation_);
#ifndef VRN_MODULE_HDF5
    forceDiskRepresentation_.setVisibleFlag(false);
#endif

    // Processors
    registerProcessor(new ConnectedComponentSelector());
    registerProcessor(new EnsembleDataSource());
    registerProcessor(new EnsembleFilter);
    registerProcessor(new PhysicalClippingLinker());
    registerProcessor(new ReferenceVolumeCreator());

    // Plotting
    registerProcessor(new FieldParallelPlotCreator());
    registerProcessor(new FieldParallelPlotViewer());
    registerProcessor(new FieldParallelPlotHistogram());
    registerProcessor(new LocalSimilarityAnalysis());
    registerProcessor(new SimilartyDataVolume());
    registerProcessor(new SimilarityMatrixCombine());
    registerProcessor(new SimilarityMatrixCreator());
    registerProcessor(new SimilarityPlot());

    // IO
    registerProcessor(new FieldPlotSave());
    registerProcessor(new FieldPlotSource());
    registerProcessor(new SimilarityMatrixSave());
    registerProcessor(new SimilarityMatrixSource());
#ifdef VRN_USE_VTK
    registerVolumeReader(new NetCDFVolumeReader());
    registerVolumeReader(new VTIVolumeReader());
    registerVolumeReader(new VTMVolumeReader());
    registerVolumeWriter(new NiftiVolumeWriter());
    registerVolumeWriter(new VTIVolumeWriter());
#endif
    
    // Properties
    registerSerializableType(new LinkEvaluatorIntListId());

    // Misc
    registerProcessor(new EnsembleVolumeExtractor());
    registerProcessor(new VolumeListMerger());
    registerProcessor(new VolumeMerger());
}

void EnsembleAnalysisModule::setForceDiskRepresentation(bool enabled) {
    forceDiskRepresentation_.set(enabled);
}

bool EnsembleAnalysisModule::getForceDiskRepresentation() const {
#ifndef VRN_MODULE_HDF5
    return false;
#else
    return forceDiskRepresentation_.get();
#endif
}

EnsembleAnalysisModule* EnsembleAnalysisModule::getInstance() {
    return instance_;
}

} // namespace
