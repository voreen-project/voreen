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

#include "voreen/qt/widgets/property/intvec2propertywidget.h"

namespace voreen {

IntVec2PropertyWidget::IntVec2PropertyWidget(IntVec2Property* prop, QWidget* parent)
    : VecPropertyWidget<SliderSpinBoxWidget, IntVec2Property, int>(prop, 2, parent)
{
    connect((const QObject*)widgets_[0], SIGNAL(valueChanged(int)), this, SLOT(setProperty(int)));
    connect((const QObject*)widgets_[1], SIGNAL(valueChanged(int)), this, SLOT(setProperty(int)));
    updateFromProperty();
}

void IntVec2PropertyWidget::setProperty(int value) {
    if (disconnected_)
        return;

    IntVec2Property::ElemType newValue = setPropertyComponent(sender(), value);
    emit valueChanged(newValue);
    emit valueModifiedByUser();
}

} //namespace voreen
