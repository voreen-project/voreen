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

#include "voreen/qt/widgets/property/lightpropertywidget.h"

#include "voreen/core/processors/volumeraycaster.h"
#include "voreen/core/properties/cameraproperty.h"
#include "voreen/core/properties//property.h"
#include "voreen/core/properties/vectorproperty.h"
#include "voreen/qt/widgets/lightwidget.h"
#include "voreen/core/datastructures/meta/templatemetadata.h"

#include <QCheckBox>

namespace voreen {

LightPropertyWidget::LightPropertyWidget(LightSourceProperty* prop, QWidget* parent)
    : QPropertyWidget(prop, parent)
    , property_(prop)
{
    light_ = new LightWidget(this);
    followCam_ = new QCheckBox(tr("Follow Camera"), this);

    layout_->setDirection(QBoxLayout::TopToBottom);
    layout_->addWidget(light_);
    layout_->addWidget(followCam_);
    connect(light_, SIGNAL(lightWidgetChanged(tgt::vec4)), this, SLOT(changeWidgetLight(tgt::vec4)));
    connect(followCam_, SIGNAL(toggled(bool)), this, SLOT(changeFollowCam(bool)));

    updateFromPropertySlot();
}

void LightPropertyWidget::changeFollowCam(bool b) {
    property_->setFollowCam(b);
}

void LightPropertyWidget::changeWidgetLight(tgt::vec4 lightPos) {
    property_->setLightPos(lightPos);
}

void LightPropertyWidget::updateFromPropertySlot() {
    light_->setMinDist(0.01);
    light_->setMaxDist(property_->getMaxDist());
    light_->setLightPosition(property_->getLightPos());
    followCam_->setCheckState(property_->getFollowCam() ? Qt::Checked : Qt::Unchecked);
}

} // namespace voreen
