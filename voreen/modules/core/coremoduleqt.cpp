/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany.                        *
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

#include "coremoduleqt.h"

#include "qt/processor/coreprocessorwidgetfactory.h"
#include "voreen/qt/widgets/property/corepropertywidgetfactory.h"

//meta data registration
#include "voreen/core/datastructures/meta/serializablevectormetadata.h"
#include "voreen/qt/networkeditor/meta/textboxmetadata.h"
#include "voreen/qt/networkeditor/meta/zoommetadata.h"

namespace voreen {

const std::string CoreModuleQt::loggerCat_("voreen.qt.CoreModuleQt");

CoreModuleQt::CoreModuleQt(const std::string& modulePath)
    : VoreenModuleQt(modulePath)
{
    setID("Core (Qt)");
    setGuiName("Core (Qt)");

    addShaderPath(getModulePath("glsl/qt"));

    registerProcessorWidgetFactory(new CoreProcessorWidgetFactory());

    registerPropertyWidgetFactory(new CorePropertyWidgetFactory());

    //meta data
    registerSerializableType(new TextBoxMetaData());
    registerSerializableType("SerializableVectorMetaData::TextBoxMetaData", new SerializableVectorMetaData<TextBoxMetaData*>());
    registerSerializableType(new ZoomMetaData());

}

} // namespace
