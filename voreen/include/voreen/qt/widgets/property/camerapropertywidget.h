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

#ifndef VRN_CAMERAPROPERTYWIDGET_H
#define VRN_CAMERAPROPERTYWIDGET_H

#include "voreen/qt/widgets/property/qpropertywidgetwithtoolwindow.h"

class QPushButton;

namespace tgt {
    class Camera;
}
namespace voreen {

class CameraProperty;
class VoreenToolWindow;
class CameraWidget;


/**
 * The property widget class associated with the CameraProperty.
 * It consists of a single button, which opens a tool window conatining
 * all adjustable properties of the camera.
 * @see CameraWidget
 */
class VRN_QT_API CameraPropertyWidget : public QPropertyWidgetWithToolWindow {
    Q_OBJECT
public:
    CameraPropertyWidget(CameraProperty* prop, QWidget* parent = 0);

public slots:
    void setProperty(tgt::Camera* value);

private slots:
    void buttonOnClick();

protected:
    virtual QWidget* createToolWindowWidget();
    virtual void customizeToolWindow();

    QPushButton* editBt_;

    CameraProperty* property_;
    CameraWidget* cameraWidget_;

protected slots:
    virtual void updateFromPropertySlot();
};

} // namespace

#endif // VRN_CAMERAPROPERTYWIDGET_H
