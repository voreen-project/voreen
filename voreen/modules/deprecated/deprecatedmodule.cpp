/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2019 University of Muenster, Germany,                        *
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

#include "deprecatedmodule.h"

#ifdef VRN_OPENGL_COMPATIBILITY_PROFILE
#include "processors/buttonoverlayprocessor.h"
#include "processors/canny.h"
#include "processors/cubeproxygeometry.h"
#include "processors/geometryclippingwidget.h"
#include "processors/multiview.h"
#include "processors/textseriessource.h"
#include "processors/rawtexturesource.h"
#include "processors/rawtexturesave.h"
#include "processors/targettotexture.h"
#include "processors/texturetotarget.h"
#include "processors/volumenormalization.h"
#include "processors/volumeseriessource.h"
#endif

#include "processors/volume/volumefiltering.h"
#include "processors/volume/volumemorphology.h"

#include "io/philipsusvolumereader.h"

#include "operators/volumeoperatornormalize.h"

namespace voreen {

DeprecatedModule::DeprecatedModule(const std::string& modulePath)
    : VoreenModule(modulePath)
{
    setID("Deprecated");
    setGuiName("Deprecated");

    addShaderPath(getModulePath("glsl"));

#ifdef VRN_OPENGL_COMPATIBILITY_PROFILE
    registerSerializableType(new ButtonOverlayProcessor());
    registerSerializableType(new Canny());
    registerSerializableType(new CubeProxyGeometry());
    registerSerializableType(new GeometryClippingWidget());
    registerSerializableType(new MultiView());
    registerSerializableType(new TextSeriesSource());
    registerSerializableType(new RawTextureSource());
    registerSerializableType(new RawTextureSave());
    registerSerializableType(new TargetToTexture);
    registerSerializableType(new TextureToTarget());
    registerSerializableType(new VolumeSeriesSource());
    registerSerializableType(new VolumeNormalization());
#endif

    registerSerializableType(new VolumeFiltering());
    registerSerializableType(new VolumeMorphology());
    registerVolumeReader(new PhilipsUSVolumeReader());

    INST_SCALAR_TYPES(VolumeOperatorNormalize, VolumeOperatorNormalizeGeneric)
    //INST_VECTOR_TYPES(VolumeOperatorNormalize, VolumeOperatorNormalizeGeneric)
}

} // namespace
