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

#include "voreen/qt/widgets/property/intvec3propertywidget.h"

namespace voreen {

IntVec3PropertyWidget::IntVec3PropertyWidget(IntVec3Property* prop, QWidget* parent)
    : VecPropertyWidget<SliderSpinBoxWidget, IntVec3Property, int>(prop, 3, parent)
{
    connect((const QObject*)widgets_[0], SIGNAL(valueChanged(int)), this, SLOT(setProperty(int)));
    connect((const QObject*)widgets_[1], SIGNAL(valueChanged(int)), this, SLOT(setProperty(int)));
    connect((const QObject*)widgets_[2], SIGNAL(valueChanged(int)), this, SLOT(setProperty(int)));
    updateFromProperty();
}

void IntVec3PropertyWidget::setProperty(int value) {
    if (disconnected_)
        return;

    IntVec3Property::ElemType newValue = setPropertyComponent(sender(), value);
    emit valueChanged(newValue);
    emit valueModifiedByUser();
}

} // namespace
