/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2021 University of Muenster, Germany,                        *
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

#include "voreen/qt/widgets/property/colorpropertywidget.h"

namespace voreen {

ColorPropertyWidget::ColorPropertyWidget(ColorProperty* prop, QWidget* parent)
    : QPropertyWidget(prop, parent)
    , property_(prop)
    , colorSelectorWidget_(0)
{
    tgtAssert(prop, "no property !!!");

    //initialize anf set color selector
    colorSelectorWidget_ = new ColorSelectorWidget("",this,prop->useAlphaChannel_);
    addWidget(colorSelectorWidget_);
    colorSelectorWidget_->setColor(toQColor(property_->get()));
    connect(colorSelectorWidget_, SIGNAL(colorChangedSignal(QColor)), this, SLOT(colorChangedSlot(QColor)));

    //set magic number...
    setMinimumHeight(18);
}

void ColorPropertyWidget::updateFromPropertySlot() {
    colorSelectorWidget_->setColor(toQColor(property_->get()));
}

void ColorPropertyWidget::colorChangedSlot(QColor color) {
    property_->set(toTgtColor(color));
    emit valueModifiedByUser();
}

//-------------------------------------------------------------------------------------------------
//      Helpers
//-------------------------------------------------------------------------------------------------
tgt::Color ColorPropertyWidget::toTgtColor(const QColor& color) {
    return tgt::Color(color.redF(),color.greenF(),color.blueF(), color.alphaF());
}

QColor ColorPropertyWidget::toQColor(const tgt::Color& color) {
    return QColor(static_cast<int>(color.r * 255), static_cast<int>(color.g * 255),
                  static_cast<int>(color.b * 255), static_cast<int>(color.a * 255));
}



} // namespace
