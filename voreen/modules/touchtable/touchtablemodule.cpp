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

//incude header file
#include "touchtablemodule.h"

// include classes to be registered
#include "processors/bodyparts3d/bodyparts3dsource.h"
#include "processors/bodyparts3d/bodyparts3drenderer.h"
#include "processors/bodyparts3d/bodyparts3dmanager.h"
#include "processors/touchinteraction/touchtableoverlay.h"
#include "processors/touchinteraction/widgets/touchtablewidget.h"
#include "processors/touchinteraction/widgets/touchtableboolwidget.h"
#include "processors/touchinteraction/widgets/touchtablemovewidget.h"
#include "processors/touchinteraction/widgets/touchtablebodypartswidget.h"
#include "processors/touchinteraction/widgets/touchtablebuttonwidget.h"
#include "processors/touchinteraction/widgets/touchtableclippingwidget.h"
#include "processors/touchinteraction/widgets/touchtabletransfuncwidget.h"
#include "processors/touchinteraction/widgets/touchtablevolumewidget.h"
#include "processors/touchinteraction/widgets/touchtablecamerawidget.h"
#include "processors/touchinteraction/widgets/touchtablesnapshotwidget.h"
#include "processors/touchinteraction/widgets/touchtablelightsourcewidget.h"
//#include "processors/touchinteraction/widgets/touchtableanimationwidget.h"

#include "processors/touchpainter.h"

//use voreen namespace
namespace voreen {

TouchtableModule::TouchtableModule(const std::string& modulePath)
    : VoreenModule(modulePath)
    , usePerspectiveTransformClipping_("touchtabletransformclipping", "Use perspective transformation for free clipping", true, Processor::INVALID_RESULT, Property::LOD_DEVELOPMENT)
{
    // module name to be used internally
    setID("Touchtable Module");

    // module name to be used in the GUI
    setGuiName("Touchtable Module");

    // each module processor needs to be registered
    registerSerializableType(new BodyParts3DSource());
    registerSerializableType(new BodyParts3DRenderer());
    registerSerializableType(new BodyParts3DManager());
    registerSerializableType(new TouchTableOverlay());
    registerSerializableType(new TouchTableBoolWidget());
    registerSerializableType(new TouchTableBodyPartsWidget());
    registerSerializableType(new TouchTableButtonWidget());
    registerSerializableType(new TouchTableClippingWidget());
    registerSerializableType(new TouchTableTransFuncWidget());
    registerSerializableType(new TouchTableCameraWidget());
    registerSerializableType(new TouchTableMoveWidget());
    registerSerializableType(new TouchTableVolumeWidget());
    registerSerializableType(new TouchTableSnapShotWidget());
    //registerSerializableType(new TouchTableAnimationWidget());
    registerSerializableType(new TouchTableLightSourceWidget());

    registerSerializableType(new TouchPainter());

    //add properties
    addProperty(usePerspectiveTransformClipping_);

    //add shader path
    addShaderPath(getModulePath("glsl"));
}

std::string TouchtableModule::getDescription() const {
    return "This module realizes functionality for using Voreen with a multi-touch table.";
}

bool TouchtableModule::usePerspectiveTransformClipping() const {
    return usePerspectiveTransformClipping_.get();
}

} // namespace
