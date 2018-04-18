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

#include "volumestreampropertywidget.h"
#include "volumestreamwidget.h"

#include "../properties/volumestreamproperty.h"
#include "voreen/qt/widgets/voreentoolwindow.h"

#include <QPushButton>

namespace voreen {

VolumeStreamPropertyWidget::VolumeStreamPropertyWidget(VolumeStreamProperty* prop, QWidget* parent)
    : QPropertyWidgetWithToolWindow(prop, parent, true, true)
    , editBt_(new QPushButton(tr("launch")))
    , property_(prop)
    , volumeStreamWidget_(0)
{
    createToolWindow(Qt::LeftDockWidgetArea);

    addWidget(editBt_);
    connect(editBt_, SIGNAL(clicked()), this, SLOT(toggleWidgetVisibility()));
}


void VolumeStreamPropertyWidget::toggleWidgetVisibility() {
    toolWindow_->setVisible(!toolWindow_->isVisible());
}

void VolumeStreamPropertyWidget::customizeToolWindow() {
    toolWindow_->adjustSize();
}

QWidget* VolumeStreamPropertyWidget::createToolWindowWidget() {
    volumeStreamWidget_ = new VolumeStreamWidget(property_, parentWidget());
    return volumeStreamWidget_;
}

Property* VolumeStreamPropertyWidget::getProperty() {
    return property_;
}

} // namespace
