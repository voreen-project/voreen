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

#include "modules/flowanalysis/flowanalysismodule.h"

// processors
#include "processors/geometry/corelinecreator.h"
#include "processors/geometry/parallelvectors.h"
#include "processors/geometry/streamlinetoboundingbox.h"
#include "processors/geometry/streamlinetogeometry.h"
#include "processors/render/flowdirectionoverlay.h"
#include "processors/render/streamlinerenderer3d.h"
#include "processors/streamline/pathlinecreator.h"
#include "processors/streamline/streamlinebundledetector.h"
#include "processors/streamline/streamlinecombine.h"
#include "processors/streamline/streamlinecreator.h"
#include "processors/streamline/streamlinepredicates.h"
#include "processors/streamline/streamlinerotation.h"
#include "processors/streamline/streamlinesave.h"
#include "processors/streamline/streamlineselector.h"
#include "processors/streamline/streamlinesource.h"
#include "processors/volume/acceleration.h"
#include "processors/volume/helicitydensity.h"
#include "processors/volume/jacobian.h"
#include "processors/volume/vortexcriterion.h"

namespace voreen {

FlowAnalysisModule::FlowAnalysisModule(const std::string& modulePath)
    : VoreenModule(modulePath)
{
    setID("Flow Analysis");
    setGuiName("Flow Analysis");

    addShaderPath(getModulePath("glsl"));

    // processors
    registerProcessor(new Acceleration());
    registerProcessor(new CorelineCreator());
    registerProcessor(new FlowDirectionOverlay());
    registerProcessor(new HelicityDensity());
    registerProcessor(new Jacobian());
    registerProcessor(new ParallelVectors());
    registerProcessor(new PathlineCreator());
    registerProcessor(new StreamlineBundleDetector());
    registerProcessor(new StreamlineCombine());
    registerProcessor(new StreamlineCreator());
    registerProcessor(new StreamlinePredicates());
    registerProcessor(new StreamlineRenderer3D());
    registerProcessor(new StreamlineRotation());
    registerProcessor(new StreamlineSave());
    registerProcessor(new StreamlineSelector());
    registerProcessor(new StreamlineSource());
    registerProcessor(new StreamlineToBoundingBox());
    registerProcessor(new StreamlineToGeometry());
    registerProcessor(new VortexCriterion());
}

} // namespace
