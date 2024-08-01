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

#ifndef VRN_EVENTPROPERTYWIDGET_H
#define VRN_EVENTPROPERTYWIDGET_H

class QBoxLayout;

#include "voreen/core/properties/propertywidget.h"
#include <QWidget>

class QCheckBox;
class QComboBox;
class QCheckBox;

namespace voreen {

class EventPropertyBase;
class ModifierDetectorWidget;
class KeyDetectorWidget;

class EventPropertyWidget : public QWidget, public PropertyWidget {
Q_OBJECT
public:
    EventPropertyWidget(EventPropertyBase* property, QWidget* parent = 0);
    ~EventPropertyWidget();

    virtual void setEnabled(bool enabled);
    virtual void setVisible(bool state);

    virtual void setPropertyGuiName(std::string);
    virtual std::string getPropertyGuiName();

    virtual void disconnect();
    virtual void updateFromProperty();
    virtual void updateViewFlags(Property::ViewFlags flags) {}; //TODO

public slots:
    void modifierChanged(Qt::KeyboardModifiers modifier);
    void keyChanged(int key);
    void buttonChanged(int button);
    void enabledChanged(bool enabled);
    void sharingChanged(bool shared);

protected:
    void createEnabledBox();
    void createSharingBox();
    void createMouseWidgets();
    void createKeyWidgets();
    void adjustWidgetState();

    QBoxLayout* layout_;
    EventPropertyBase* property_;

    QCheckBox* checkEnabled_;
    ModifierDetectorWidget* modifierWidget_;
    KeyDetectorWidget* keyWidget_;
    QComboBox* buttonBox_;
    QCheckBox* checkSharing_;

    bool disconnected_;
};

} // namespace

#endif // VRN_EVENTPROPERTYWIDGET_H
