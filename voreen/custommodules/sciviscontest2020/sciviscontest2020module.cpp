/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany,                        *
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

#include "sciviscontest2020module.h"

#include "processors/approximateparallelvectors.h"
#include "processors/curlprocessor.h"
#include "processors/corelinedensityvolumecreator.h"
#include "processors/flowmapprocessor.h"
#include "processors/ftvaprocessor.h"
#include "processors/vortexprocessor.h"
#include "processors/particlerenderer.h"
#include "processors/rotationaldirectionprocessor.h"
#include "processors/uncertainvectorfieldprocessor.h"
#include "processors/vortexcollectioncreator.h"
#include "processors/vortexcollectionsource.h"
#include "processors/vortexlistselector.h"
#include "processors/vortexmatchselector.h"
#include "processors/vortexprocessor.h"
#include "processors/vortexselector.h"
#include "processors/vortextracking.h"

namespace voreen {

SciVisContest2020Module::SciVisContest2020Module(const std::string& modulePath)
    : VoreenModule(modulePath)
{
    setID("SciVisContest2020");
    setGuiName("SciVis Contest 2020");

    addShaderPath(getModulePath("glsl"));

    // processors
    registerProcessor(new ApproximateParallelVectors());
    registerProcessor(new CorelineDensityVolumeCreator());
    registerProcessor(new CurlProcessor());
    registerProcessor(new FlowMapProcessor());
    registerProcessor(new FTVAProcessor());
    registerProcessor(new ParticleRenderer());
    registerProcessor(new RotationalDirectionProcessor());
    registerProcessor(new UncertainVectorFieldProcessor());
    registerProcessor(new VortexCollectionCreator());
    registerProcessor(new VortexCollectionSource());
    registerProcessor(new VortexListSelector());
    registerProcessor(new VortexMatchSelector());
    registerProcessor(new VortexProcessor());
    registerProcessor(new VortexSelector());
    registerProcessor(new VortexTracking());
}

} // namespace
