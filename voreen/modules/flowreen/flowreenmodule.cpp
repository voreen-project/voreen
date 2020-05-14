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

#include "modules/flowreen/flowreenmodule.h"

// processors
#include "processors/geometry/streamlinetoboundingbox.h"
#include "processors/geometry/streamlinetogeometry.h"
#include "processors/render/flowdirectionoverlay.h"
#include "processors/render/streamlinerenderer3d.h"
#include "processors/streamline/pathlinecreator.h"
#include "processors/streamline/streamlinebundledetector.h"
#include "processors/streamline/streamlinecombine.h"
#include "processors/streamline/streamlinecreator.h"
#include "processors/streamline/streamlinefilter.h"
#include "processors/streamline/streamlinerotation.h"
#include "processors/streamline/streamlinesave.h"
#include "processors/streamline/streamlineselector.h"
#include "processors/streamline/streamlinesource.h"
#include "processors/volume/helicitydensity.h"

namespace voreen {

FlowreenModule::FlowreenModule(const std::string& modulePath)
    : VoreenModule(modulePath)
{
    setID("Flowreen");
    setGuiName("Flowreen");

    addShaderPath(getModulePath("glsl"));

    // processors
    registerSerializableType(new FlowDirectionOverlay());
    registerSerializableType(new HelicityDensity());
    registerSerializableType(new PathlineCreator());
    registerSerializableType(new StreamlineBundleDetector());
    registerSerializableType(new StreamlineCombine());
    registerSerializableType(new StreamlineCreator());
    registerSerializableType(new StreamlineFilter());
    registerSerializableType(new StreamlineRenderer3D());
    registerSerializableType(new StreamlineRotation());
    registerSerializableType(new StreamlineSave());
    registerSerializableType(new StreamlineSelector());
    registerSerializableType(new StreamlineSource());
    registerSerializableType(new StreamlineToBoundingBox());
    registerSerializableType(new StreamlineToGeometry());
}

} // namespace
