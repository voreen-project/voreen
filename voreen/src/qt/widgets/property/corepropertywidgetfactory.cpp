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

#include "voreen/qt/widgets/property/corepropertywidgetfactory.h"

#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/boundingboxproperty.h"
#include "voreen/core/properties/buttonproperty.h"
#include "voreen/core/properties/cameraproperty.h"
#include "voreen/core/properties/colorproperty.h"
#include "voreen/core/properties/filedialogproperty.h"
#include "voreen/core/properties/floatproperty.h"
#include "voreen/core/properties/fontproperty.h"
#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/lightsourceproperty.h"
#include "voreen/core/properties/matrixproperty.h"
#include "voreen/core/properties/optionproperty.h"
#include "voreen/core/properties/progressproperty.h"
#include "voreen/core/properties/propertyvector.h"
#include "voreen/core/properties/shaderproperty.h"
#include "voreen/core/properties/stringexpressionproperty.h"
#include "voreen/core/properties/string/stringtableproperty.h"
#include "voreen/core/properties/stringproperty.h"
#include "voreen/core/properties/transfunc/1d/1dkeys/transfunc1dkeysproperty.h"
#include "voreen/core/properties/transfunc/1d/1dgaussian/transfunc1dgaussianproperty.h"
#include "voreen/core/properties/transfunc/2d/2dprimitives/transfunc2dprimitivesproperty.h"
#include "voreen/core/properties/transfunc/transfunctypeproperty.h"
#include "voreen/core/properties/vectorproperty.h"
#include "voreen/core/properties/volumeurllistproperty.h"
#include "voreen/core/properties/volumeurlproperty.h"
#include "voreen/core/properties/voxeltypeproperty.h"
#include "voreen/core/properties/numeric/intervalproperty.h"
#include "voreen/core/properties/color/colorswitchproperty.h"
#include "voreen/core/ports/renderport.h"

#include "voreen/qt/widgets/property/boolpropertywidget.h"
#include "voreen/qt/widgets/property/floatboundingboxpropertywidget.h"
#include "voreen/qt/widgets/property/intboundingboxpropertywidget.h"
#include "voreen/qt/widgets/property/buttonpropertywidget.h"
#include "voreen/qt/widgets/property/camerapropertywidget.h"
#include "voreen/qt/widgets/property/colorpropertywidget.h"
#include "voreen/qt/widgets/property/filedialogpropertywidget.h"
#include "voreen/qt/widgets/property/floatmat2propertywidget.h"
#include "voreen/qt/widgets/property/floatmat3propertywidget.h"
#include "voreen/qt/widgets/property/floatmat4propertywidget.h"
#include "voreen/qt/widgets/property/floatpropertywidget.h"
#include "voreen/qt/widgets/property/floatvec2propertywidget.h"
#include "voreen/qt/widgets/property/floatvec3propertywidget.h"
#include "voreen/qt/widgets/property/floatvec4propertywidget.h"
#include "voreen/qt/widgets/property/fontpropertywidget.h"
#include "voreen/qt/widgets/property/intpropertywidget.h"
#include "voreen/qt/widgets/property/intvec2propertywidget.h"
#include "voreen/qt/widgets/property/intvec3propertywidget.h"
#include "voreen/qt/widgets/property/intvec4propertywidget.h"
#include "voreen/qt/widgets/property/lightpropertywidget.h"
#include "voreen/qt/widgets/property/optionpropertywidget.h"
#include "voreen/qt/widgets/property/progresspropertywidget.h"
#include "voreen/qt/widgets/property/propertyvectorwidget.h"
#include "voreen/qt/widgets/property/shaderpropertywidget.h"
#include "voreen/qt/widgets/property/stringexpressionpropertywidget.h"
#include "voreen/qt/widgets/property/string/stringtablepropertywidget.h"
#include "voreen/qt/widgets/property/stringpropertywidget.h"
#include "voreen/qt/widgets/property/grouppropertywidget.h"
#include "voreen/qt/widgets/property/transfunc/1d/1dkeys/transfunc1dkeyspropertywidget.h"
#include "voreen/qt/widgets/property/transfunc/1d/1dgaussian/transfunc1dgaussianpropertywidget.h"
#include "voreen/qt/widgets/property/transfunc/2d/2dprimitives/transfunc2dprimitivespropertywidget.h"
#include "voreen/qt/widgets/property/volumeinfopropertywidget.h"
#include "voreen/qt/widgets/property/volumeurllistpropertywidget.h"
#include "voreen/qt/widgets/property/volumeurlpropertywidget.h"
#include "voreen/qt/widgets/property/voxeltypepropertywidget.h"
#include "voreen/qt/widgets/property/numeric/intintervalpropertywidget.h"
#include "voreen/qt/widgets/property/numeric/floatintervalpropertywidget.h"
#include "voreen/qt/widgets/property/color/colorswitchpropertywidget.h"

namespace voreen {

PropertyWidget* CorePropertyWidgetFactory::createAssociatedWidget(Property* prop) const {

    if (!prop)
        return 0;

    if (typeid(*prop) == typeid(BoolProperty))
        return new BoolPropertyWidget(static_cast<BoolProperty*>(prop), 0);

    if (typeid(*prop) == typeid(ButtonProperty))
        return new ButtonPropertyWidget(static_cast<ButtonProperty*>(prop), 0);

    if (typeid(*prop) == typeid(CameraProperty))
        return new CameraPropertyWidget(static_cast<CameraProperty*>(prop), 0);

    if (typeid(*prop) == typeid(FileDialogProperty))
        return new FileDialogPropertyWidget(static_cast<FileDialogProperty*>(prop), 0);

    if (typeid(*prop) == typeid(FloatProperty))
        return new FloatPropertyWidget(static_cast<FloatProperty*>(prop), 0);

    if (typeid(*prop) == typeid(FloatVec2Property))
        return new FloatVec2PropertyWidget(static_cast<FloatVec2Property*>(prop), 0);

    if (typeid(*prop) == typeid(FloatVec3Property))
        return new FloatVec3PropertyWidget(static_cast<FloatVec3Property*>(prop), 0);

    if (typeid(*prop) == typeid(FloatVec4Property))
        return new FloatVec4PropertyWidget(static_cast<FloatVec4Property*>(prop), 0);

    if (typeid(*prop) == typeid(ColorProperty))
        return new ColorPropertyWidget(static_cast<ColorProperty*>(prop), 0);

    if (typeid(*prop) == typeid(FontProperty))
        return new FontPropertyWidget(static_cast<FontProperty*>(prop), 0);

    if (typeid(*prop) == typeid(IntProperty))
        return new IntPropertyWidget(static_cast<IntProperty*>(prop), 0);

    if (typeid(*prop) == typeid(IntVec2Property))
        return new IntVec2PropertyWidget(static_cast<IntVec2Property*>(prop), 0);
    if (typeid(*prop) == typeid(RenderSizeOriginProperty))
        return new IntVec2PropertyWidget(static_cast<IntVec2Property*>(prop), 0);
    if (typeid(*prop) == typeid(RenderSizeReceiveProperty))
        return new IntVec2PropertyWidget(static_cast<IntVec2Property*>(prop), 0);

    if (typeid(*prop) == typeid(IntVec3Property))
        return new IntVec3PropertyWidget(static_cast<IntVec3Property*>(prop), 0);

    if (typeid(*prop) == typeid(IntVec4Property))
        return new IntVec4PropertyWidget(static_cast<IntVec4Property*>(prop), 0);

    if (typeid(*prop) == typeid(FloatMat2Property))
        return new FloatMat2PropertyWidget(static_cast<FloatMat2Property*>(prop), 0);

    if (typeid(*prop) == typeid(FloatMat3Property))
        return new FloatMat3PropertyWidget(static_cast<FloatMat3Property*>(prop), 0);

    if (typeid(*prop) == typeid(FloatMat4Property))
        return new FloatMat4PropertyWidget(static_cast<FloatMat4Property*>(prop), 0);

    if (typeid(*prop) == typeid(IntIntervalProperty))
        return new IntIntervalPropertyWidget(static_cast<IntIntervalProperty*>(prop), 0);

    if (typeid(*prop) == typeid(FloatIntervalProperty))
        return new FloatIntervalPropertyWidget(static_cast<FloatIntervalProperty*>(prop), 0);

    if (typeid(*prop) == typeid(FloatBoundingBoxProperty))
        return new FloatBoundingBoxPropertyWidget(static_cast<FloatBoundingBoxProperty*>(prop), 0);

    if (typeid(*prop) == typeid(IntBoundingBoxProperty))
        return new IntBoundingBoxPropertyWidget(static_cast<IntBoundingBoxProperty*>(prop), 0);

    if (typeid(*prop) == typeid(LightSourceProperty)) {
        LightSourceProperty* lsp = static_cast<LightSourceProperty*>(prop);
        GroupPropertyWidget* tab = new GroupPropertyWidget(lsp, true, "");
        tab->addWidget(new LightPropertyWidget(lsp, 0), "Widget");
        tab->addWidget(new FloatVec4PropertyWidget(lsp, 0), "Vector");
        return tab;
    }

    // dynamic cast necessary, since we are dealing with an abstract base class
    if (dynamic_cast<OptionPropertyBase*>(prop))
        return new OptionPropertyWidget(static_cast<OptionPropertyBase*>(prop), 0);

    if (typeid(*prop) == typeid(ProgressProperty))
        return new ProgressPropertyWidget(static_cast<ProgressProperty*>(prop), 0);

    if (typeid(*prop) == typeid(PropertyVector))
        return new PropertyVectorWidget(static_cast<PropertyVector*>(prop), 0);

    if (typeid(*prop) == typeid(ShaderProperty))
        return new ShaderPropertyWidget(static_cast<ShaderProperty*>(prop), 0);

    if (typeid(*prop) == typeid(StringExpressionProperty))
        return new StringExpressionPropertyWidget(static_cast<StringExpressionProperty*>(prop), 0);

    if (typeid(*prop) == typeid(StringTableProperty))
        return new StringTablePropertyWidget(static_cast<StringTableProperty*>(prop), 0);

    if (typeid(*prop) == typeid(StringProperty))
        return new StringPropertyWidget(static_cast<StringProperty*>(prop), 0);

    if (typeid(*prop) == typeid(TransFunc1DKeysProperty))
        return new TransFunc1DKeysPropertyWidget(static_cast<TransFunc1DKeysProperty*>(prop), 0);

    if (typeid(*prop) == typeid(TransFunc1DGaussianProperty))
        return new TransFunc1DGaussianPropertyWidget(static_cast<TransFunc1DGaussianProperty*>(prop), 0);

    if (typeid(*prop) == typeid(TransFunc2DPrimitivesProperty))
        return new TransFunc2DPrimitivesPropertyWidget(static_cast<TransFunc2DPrimitivesProperty*>(prop), 0);

    if (typeid(*prop) == typeid(VolumeInfoProperty))
        return new VolumeInfoPropertyWidget(static_cast<VolumeInfoProperty*>(prop), 0);

    if (typeid(*prop) == typeid(VolumeURLProperty))
        return new VolumeURLPropertyWidget(static_cast<VolumeURLProperty*>(prop), 0);

    if (typeid(*prop) == typeid(VolumeURLListProperty))
        return new VolumeURLListPropertyWidget(static_cast<VolumeURLListProperty*>(prop), 0);

    if (typeid(*prop) == typeid(VoxelTypeProperty))
        return new VoxelTypePropertyWidget(static_cast<VoxelTypeProperty*>(prop), 0);
    if (typeid(*prop) == typeid(ColorSwitchProperty))
        return new ColorSwitchPropertyWidget(static_cast<ColorSwitchProperty*>(prop), 0);
    return 0;
}

} // namespace voreen
