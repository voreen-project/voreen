/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2024 University of Muenster, Germany,                        *
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

#include "stagingmodule.h"

#include "processors/alignedsliceproxygeometry.h"
#include "processors/clipregiongeometrycreator.h"
#include "processors/geometryslicerenderer.h"
#include "processors/interactiveregistrationwidget.h"
#include "processors/multislicerenderer.h"
#include "processors/multisliceviewer.h"
#include "processors/planegeometrycreator.h"
#include "processors/preintegrationtablerenderer.h"
#include "processors/samplingpositiontransformation.h"
#include "processors/toucheventsimulator.h"
#include "processors/transfuncoverlay.h"
#include "processors/particles.h"
#include "processors/tabbedview.h"
#include "processors/transfuncalphachannelanimation.h"
#include "processors/arbitraryvolumeclipping.h"
#include "processors/pong.h"
#include "processors/screenspaceambientocclusion.h"
#include "processors/registrationinitializer.h"
#include "processors/volumerealworldmapping.h"
#include "processors/volumeuncertaintymeasure.h"
#include "processors/slicepoints/slicepointrenderer2d.h"
#include "processors/slicepoints/slicepointrenderer3d.h"
#include "processors/simdraycaster/simdraycaster.h"

namespace voreen {

StagingModule::StagingModule(const std::string& modulePath)
    : VoreenModule(modulePath)
{
    setID("Staging");
    setGuiName("Staging");

    addShaderPath(getModulePath("glsl"));

    registerSerializableType(new AlignedSliceProxyGeometry());
    registerSerializableType(new ClipRegionGeometryCreator());
    registerSerializableType(new GeometrySliceRenderer());
    registerSerializableType(new InteractiveRegistrationWidget());
    registerSerializableType(new MultiSliceRenderer());
    registerSerializableType(new MultiSliceViewer());
    registerSerializableType(new PlaneGeometryCreator());
    registerSerializableType(new PreIntegrationTableRenderer());
    registerSerializableType(new SamplingPositionTransformation());
    registerSerializableType(new TouchEventSimulator());
    registerSerializableType(new TransFuncOverlay());
    registerSerializableType(new ArbitraryVolumeClipping());
    registerSerializableType(new Pong());
    registerSerializableType(new ScreenSpaceAmbientOcclusion());
    registerSerializableType(new RegistrationInitializer());
    registerSerializableType(new VolumeRealWorldMapping());
    registerSerializableType(new VolumeUncertaintyMeasure());
    registerSerializableType(new SlicePointRenderer2D());
    registerSerializableType(new SlicePointRenderer3D());
    registerSerializableType(new SIMDRayCaster());
    #ifdef GL_COMPUTE_SHADER //disable compilation for old gl headers
        registerSerializableType(new Particles());
    #endif
    registerSerializableType(new TabbedView());


    registerProcessor(new TransFuncAlphaChannelAnimation());
}

} // namespace
