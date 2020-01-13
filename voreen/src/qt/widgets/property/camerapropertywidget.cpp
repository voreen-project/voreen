/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2020 University of Muenster, Germany,                        *
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

#include "voreen/qt/widgets/property/camerapropertywidget.h"

#include "voreen/core/properties/cameraproperty.h"

#include "voreen/qt/widgets/voreentoolwindow.h"
#include "voreen/qt/widgets/property/camerawidget.h"

#include "tgt/camera.h"

#include <QPushButton>

namespace voreen {

using tgt::Camera;

CameraPropertyWidget::CameraPropertyWidget(CameraProperty* prop, QWidget* parent)
    : QPropertyWidgetWithToolWindow(prop, parent, true, false)
    , editBt_(new QPushButton(tr("edit")))
    , property_(prop)
    , cameraWidget_(0)
{
    addWidget(editBt_);
    connect(editBt_, SIGNAL(clicked()), this, SLOT(buttonOnClick()));

    QFontInfo fontInfo(font());
    editBt_->setFont(QFont(fontInfo.family(), QPropertyWidget::fontSize_));
}

void CameraPropertyWidget::updateFromPropertySlot() {
    if (cameraWidget_ && cameraWidget_->isVisible())
        cameraWidget_->updateFromCamera();
}

void CameraPropertyWidget::setProperty(Camera* value) {
    if (!disconnected_)
        property_->set(*value);
}

void CameraPropertyWidget::buttonOnClick() {
    if(!toolWindow_) {
        createToolWindow(Qt::LeftDockWidgetArea); //< @see QPropertyWindgetWithToolWindow
        toolWindow_->setVisible(true);
    }else
        toggleToolWindow(); //< @see QPropertyWindgetWithToolWindow
}

void CameraPropertyWidget::customizeToolWindow() {
    toolWindow_->adjustSize();
    //editorWindow_->setFixedSize(editorWindow_->size());
}

QWidget* CameraPropertyWidget::createToolWindowWidget() {
    cameraWidget_ = new CameraWidget(property_, parentWidget());
    cameraWidget_->createWidgets();
    cameraWidget_->createConnections();

    return cameraWidget_;
}

} // namespace
