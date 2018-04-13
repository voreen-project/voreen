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

#include "voreen/qt/widgets/property/color/colorswitchpropertywidget.h"

#include <QColorDialog>
#include <QLabel>
#include <QMouseEvent>
#include <QPainter>
#include <QVBoxLayout>
#include <QHBoxLayout>

namespace voreen {

ColorSwitchPropertyWidget::ColorSwitchPropertyWidget(ColorSwitchProperty* prop, QWidget* parent)
    : QPropertyWidget(prop, parent)
    , property_(prop)
    , activeColorLabel_(new ColorSelectorWidget(""))
    , inactiveColorLabel_(new ColorSelectorWidget(""))
    , useActiveColor_(new QCheckBox("Use active Color (used for linking!)"))
{
    updateFromPropertySlot();

    QVBoxLayout* box = new QVBoxLayout();
    addLayout(box);


    QHBoxLayout* activeLayout = new QHBoxLayout;
    QHBoxLayout* inactiveLayout = new QHBoxLayout;
    box->addLayout(activeLayout);
    box->addLayout(inactiveLayout);
    box->addWidget(useActiveColor_);

    useActiveColor_->setEnabled(false);
    activeLayout->addWidget(new QLabel("Active color:"));
    activeLayout->addWidget(activeColorLabel_);
    inactiveLayout->addWidget(new QLabel("Inactive color:"));
    inactiveLayout->addWidget(inactiveColorLabel_);

    connect(activeColorLabel_, SIGNAL(colorChangedSignal(QColor)), this, SLOT(changeActiveColor(QColor)));

    connect(inactiveColorLabel_, SIGNAL(colorChangedSignal(QColor)), this, SLOT(changeInactiveColor(QColor)));
    connect(useActiveColor_, SIGNAL(stateChanged(int)), this, SLOT(useActiveColorClicked(int)));
    setMinimumHeight(18);
}

void ColorSwitchPropertyWidget::updateFromPropertySlot() {
    activeColorLabel_->setColor(toQColor(property_->getActiveColor()));
    inactiveColorLabel_->setColor(toQColor(property_->getInactiveColor()));
    useActiveColor_->setChecked(property_->getUseActiveColor());
}

void ColorSwitchPropertyWidget::changeActiveColor(QColor color) {
    if (!disconnected_) {
        tgt::Color value = toTgtColor(color);
        property_->setActiveColor(value);
    }
}

void ColorSwitchPropertyWidget::changeInactiveColor(QColor color) {
    if (!disconnected_) {
        tgt::Color value = toTgtColor(color);
        property_->setInactiveColor(value);
    }
}


tgt::Color ColorSwitchPropertyWidget::toTgtColor(QColor color) {
    return tgt::Color(color.redF(),color.greenF(),color.blueF(), color.alphaF());
}

QColor ColorSwitchPropertyWidget::toQColor(tgt::Color color) {
    return QColor(static_cast<int>(color.r * 255), static_cast<int>(color.g * 255),
                  static_cast<int>(color.b * 255), static_cast<int>(color.a * 255));
}

void ColorSwitchPropertyWidget::useActiveColorClicked(int)
{
    if (!disconnected_) {
        property_->setUseActiveColor(useActiveColor_->isChecked());
    }
}

} // namespace
