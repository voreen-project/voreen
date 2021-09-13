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

#include "sciviscontest2019module.h"

#include "processors/cosmologyvolumeconverter.h"
#include "processors/cosmologyarrowrenderer3d.h"
#include "processors/cosmologyparticlerenderer.h"
#include "processors/cosmologytracingpath.h"
#include "processors/cosmologyparticlefilter.h"
#include "processors/cosmologyhalorenderer.h"
#include "processors/cosmologydatasource.h"
#include "processors/cosmologyfakedatasource.h"
#include "processors/cosmologytimestepanimationprocessor.h"
#include "processors/cosmologyparticleboxfilter.h"
#include "processors/cosmologyparticlehalofilter.h"
#include "processors/cosmologyparticlecompactor.h"
#include "processors/cosmologytimestepoverlay.h"
#include "processors/cosmologyorientationoverlay.h"
#include "processors/cosmologyoctreevolumepreprocess.h"
#include "processors/cmhalodatasource.h"
#include "processors/cmmergertreerenderer.h"
#include "processors/cmhaloidtimesteplinker.h"
#include "processors/cmlazylinker.h"
#include "processors/cmhalodescriptor.h"
#include "processors/cmhalosubstructurerenderer.h"
#include "processors/cmhalodistanceoverlay.h"
#include "processors/cminfooverlay.h"
#include "processors/cmplotviewer.h"
#include "processors/cmplotcreator.h"
#include "processors/cosmologyyearprovider.h"
#include "voreen/core/voreenapplication.h"

namespace voreen {

SciVisContest2019Module::SciVisContest2019Module(const std::string& modulePath)
    : VoreenModule(modulePath)
{
    setID("SciVisContest2019");
    setGuiName("SciVis Contest 2019");

    addShaderPath(getModulePath("glsl"));

    registerProcessor(new CosmologyArrowRenderer3D());

    registerProcessor(new CosmologyDataSource());
    registerProcessor(new CosmologyFakeDataSource());
    registerProcessor(new CMHaloDataSource());

    registerProcessor(new CosmologyVolumeConverter());
    registerProcessor(new CosmologyParticleRenderer());
    registerProcessor(new CosmologyTracingPath());
    registerProcessor(new CosmologyParticleFilter());
    registerProcessor(new CosmologyTimeStepAnimationProcessor());
    registerProcessor(new CosmologyParticleBoxFilter());
    registerProcessor(new CosmologyParticleHaloFilter());
    registerProcessor(new CosmologyParticleCompactor());
    registerProcessor(new CosmologyHaloRenderer());
    registerProcessor(new CosmologyTimeStepOverlay());
    registerProcessor(new CosmologyOrientationOverlay());
    registerProcessor(new CosmologyOctreeVolumePreprocess());
    registerProcessor(new CMMergerTreeRenderer());
    registerProcessor(new CMHaloIDTimeStepLinker());
    registerProcessor(new CMLazyLinker());
    registerProcessor(new CMHaloDescriptor());
    registerProcessor(new CMHaloSubstructureRenderer());
    registerProcessor(new CMHaloDistanceOverlay());
    registerProcessor(new CMInfoOverlay());
    registerProcessor(new CMPlotViewer());
    registerProcessor(new CMPlotCreator());
    registerProcessor(new CosmologyYearProvider());
}

} // namespace
