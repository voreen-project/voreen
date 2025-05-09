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

#ifndef VRN_FLOATVEC4PROPERTYWIDGET_H
#define VRN_FLOATVEC4PROPERTYWIDGET_H

#include "voreen/core/properties/vectorproperty.h"
#include "voreen/qt/widgets/property/floatpropertywidget.h"
#include "voreen/qt/widgets/property/vecpropertywidget.h"

namespace voreen {

class FloatVec4PropertyWidget : public VecPropertyWidget<DoubleSliderSpinBoxWidget, FloatVec4Property, float> {
Q_OBJECT;
public:
    FloatVec4PropertyWidget(FloatVec4Property* prop, QWidget* parent = 0);

public slots:
    void setProperty(double value);

signals:
    void valueChanged(FloatVec4Property::ElemType);
};


} // namespace voreen

#endif // VRN_COMPACTFLOATVEC4PROPERTYWIDGET_H
