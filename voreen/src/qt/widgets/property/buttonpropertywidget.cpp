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

#include "voreen/qt/widgets/property/buttonpropertywidget.h"

#include "voreen/core/properties/buttonproperty.h"
#include "voreen/qt/widgets/customlabel.h"

#include <QPushButton>
#include <QLabel>

namespace voreen {

ButtonPropertyWidget::ButtonPropertyWidget(ButtonProperty* prop, QWidget* parent)
    : QPropertyWidget(prop, parent, false)
    , property_(prop)
    , button_(new CustomButton(prop->getGuiName().c_str(), this, true))
{
    connect(button_, SIGNAL(clicked()), this, SLOT(clicked()));
    connect(button_, SIGNAL(clicked(void)), this, SIGNAL(widgetChanged(void)));
    addWidget(button_);
    QFontInfo fontInfo(font());
    button_->setFont(QFont(fontInfo.family(), QPropertyWidget::fontSize_));
}

void ButtonPropertyWidget::clicked() {
    if (!disconnected_)
        property_->clicked();
}

CustomLabel* ButtonPropertyWidget::getOrCreateNameLabel() const {
    if (!nameLabel_)
        nameLabel_ = new CustomLabel("", const_cast<ButtonPropertyWidget*>(this));
    return nameLabel_;
    //return QPropertyWidget::getOrCreateNameLabel();
}

void ButtonPropertyWidget::updateFromPropertySlot() {
    button_->setText(property_->getGuiName().c_str());
}

void ButtonPropertyWidget::setPropertyGuiName(std::string name) {
    if(button_ && button_->text().compare(QString(name.c_str())))
        button_->setText(QString(name.c_str()));
    if(prop_ && prop_->getGuiName().compare(name))
        prop_->setGuiName(name);
}


} // namespace
