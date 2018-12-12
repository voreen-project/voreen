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
#include "processors/flowdirectionoverlay.h"
#include "processors/streamline/pathlinecreator.h"
#include "processors/streamline/streamlinecombine.h"
#include "processors/streamline/streamlinecreator.h"
#include "processors/streamlinerenderer3d.h"
#include "processors/streamline/streamlinerotation.h"
#include "processors/streamline/streamlinesave.h"
#include "processors/streamline/streamlineselector.h"
#include "processors/streamline/streamlinesource.h"
#include "processors/streamline/streamlinetoboundingbox.h"

#ifdef FLOWREEN_USE_OPENLB
#ifdef VRN_MODULE_OPENMP
#define PARALLEL_MODE_OMP
#endif
#include <olb3D.h>

#include "processors/flowsimulation.h"
#endif

#ifdef VRN_OPENGL_COMPATIBILITY_PROFILE
    #include "processors/flowarrowrenderer2d.h"
    #include "processors/flowarrowrenderer3d.h"
    #include "processors/flowmagnitudes3d.h"
    #include "processors/floworthogonalslicerenderer.h"
    #include "processors/flowreenadapter.h"
    #include "processors/flowslicerenderer2d.h"
    #include "processors/flowslicerenderer3d.h"
    #include "processors/flowstreamlinestexture3d.h"
    #include "processors/pathlinerenderer3d.h"

    // I/O
    #include "modules/flowreen/io/flowreader.h"
    // VolumeOperators
    #include "modules/flowreen/datastructures/deprecated/volumeoperatorintensitymask.h"
#endif

namespace voreen {

FlowreenModule::FlowreenModule(const std::string& modulePath)
    : VoreenModule(modulePath)
{
    setID("Flowreen");
    setGuiName("Flowreen");

    addShaderPath(getModulePath("glsl"));

    // processors
    registerSerializableType(new FlowDirectionOverlay());
    registerSerializableType(new PathlineCreator());
    registerSerializableType(new StreamlineCombine());
    registerSerializableType(new StreamlineRenderer3D());
    registerSerializableType(new StreamlineCreator());
    registerSerializableType(new StreamlineRotation());
    registerSerializableType(new StreamlineSave());
    registerSerializableType(new StreamlineSelector());
    registerSerializableType(new StreamlineSource());
    registerSerializableType(new StreamlineToBoundingBox());

#ifdef FLOWREEN_USE_OPENLB
    registerSerializableType(new FlowSimulation());
#endif

#ifdef VRN_OPENGL_COMPATIBILITY_PROFILE
    registerSerializableType(new FlowArrowRenderer2D);
    registerSerializableType(new FlowArrowRenderer3D);
    registerSerializableType(new FlowMagnitudes3D());
    registerSerializableType(new FlowOrthogonalSliceRenderer());
    registerSerializableType(new FlowSliceRenderer2D());
    registerSerializableType(new FlowSliceRenderer3D());
    registerSerializableType(new FlowStreamlinesTexture3D());
    registerSerializableType(new FlowreenAdapter());
    registerSerializableType(new PathlineRenderer3D());
    // I/O
    registerVolumeReader(new FlowReader());

    INST_SCALAR_TYPES(VolumeOperatorIntensityMask, VolumeOperatorIntensityMaskGeneric)
#endif

}

void FlowreenModule::initialize() {
    VoreenModule::initialize();

#ifdef FLOWREEN_USE_OPENLB
    olb::singleton::directories().setOutputDir(VoreenApplication::app()->getTemporaryPath("simulation")+"/");
#endif
}

void FlowreenModule::deinitialize() {
    VoreenModule::deinitialize();
}


} // namespace
