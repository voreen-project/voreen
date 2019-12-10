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

#include "voreen/qt/widgets/property/boolpropertywidget.h"

#include "voreen/core/properties/boolproperty.h"

#include <QCheckBox>

namespace voreen {

BoolPropertyWidget::BoolPropertyWidget(BoolProperty* prop, QWidget* parent)
    : QPropertyWidget(prop, parent)
    , property_(prop)
    , chkBox_(new QCheckBox)
{
    chkBox_->setChecked(prop->get());
    addWidget(chkBox_);

    connect(chkBox_, SIGNAL(toggled(bool)), this, SLOT(setProperty(bool)));
    connect(chkBox_, SIGNAL(toggled(bool)), this, SIGNAL(widgetChanged()));
}

void BoolPropertyWidget::updateFromPropertySlot() {
    chkBox_->blockSignals(true);
    chkBox_->setChecked(property_->get());
    chkBox_->blockSignals(false);
}

void BoolPropertyWidget::setProperty(bool value) {
    if (!disconnected_) {
        property_->set(value);
        emit valueModifiedByUser();
    }
}

} // namespace
