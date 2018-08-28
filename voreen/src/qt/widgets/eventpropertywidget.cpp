/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany,                        *
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

#include "voreen/qt/widgets/eventpropertywidget.h"
#include "voreen/core/properties/eventproperty.h"
#include "voreen/core/processors/processor.h"
#include "voreen/qt/widgets/keydetectorwidget.h"

#include <QComboBox>
#include <QHBoxLayout>
#include <QLabel>
#include <QCheckBox>

namespace voreen {

EventPropertyWidget::EventPropertyWidget(EventPropertyBase* property, QWidget* parent)
    : QWidget(parent)
    , PropertyWidget(property)
    , property_(property)
    , checkEnabled_(0)
    , modifierWidget_(0)
    , keyWidget_(0)
    , buttonBox_(0)
    , checkSharing_(0)
    , disconnected_(false)
{
    tgtAssert(property, "No property");

    layout_ = new QHBoxLayout(this);
    layout_->setContentsMargins(2,2,2,2);

    createEnabledBox();

    if (property->receivesMouseEvents())
        createMouseWidgets();
    if (property->receivesKeyEvents())
        createKeyWidgets();

    layout_->addSpacing(5);
    createSharingBox();

    enabledChanged(property->isEnabled());
    if (!property->isVisibleFlagSet())
        setVisible(false);
}

EventPropertyWidget::~EventPropertyWidget() {
    if (!disconnected_ && property_) {
        property_->removeWidget(this);
        property_ = 0;
    }
}

void EventPropertyWidget::setPropertyGuiName(std::string name) {
    checkEnabled_->setText(QString(name.c_str()));
    property_->setGuiName(name);
}

std::string EventPropertyWidget::getPropertyGuiName() {
    tgtAssert(checkEnabled_, "no widget");
    return checkEnabled_->text().toStdString();
}

void EventPropertyWidget::createMouseWidgets() {
    modifierWidget_ = new ModifierDetectorWidget;
    modifierWidget_->setText(ModifierDetectorWidget::getStringForModifier(ModifierDetectorWidget::getQtModifierFromTGT(property_->getModifier())));
    modifierWidget_->setFixedWidth(60);
    connect(modifierWidget_, SIGNAL(modifierChanged(Qt::KeyboardModifiers)), this, SLOT(modifierChanged(Qt::KeyboardModifiers)));
    layout_->addWidget(modifierWidget_);

    std::string suffix = "";
    if(property_->getMouseAction() == tgt::MouseEvent::DOUBLECLICK) {
        suffix = " (Doubleclick)";
    }
    buttonBox_ = new QComboBox;
    buttonBox_->addItem(("Left mouse button" + suffix).c_str());
    buttonBox_->addItem(("Middle mouse button" + suffix).c_str());
    buttonBox_->addItem(("Right mouse button" + suffix).c_str());
    buttonBox_->addItem(("Mouse wheel" + suffix).c_str());
    buttonBox_->addItem(("Any mouse button" + suffix).c_str());
    buttonBox_->addItem(("No mouse button" + suffix).c_str());
    switch (property_->getMouseButtons()) {
        case tgt::MouseEvent::MOUSE_BUTTON_LEFT:
            buttonBox_->setCurrentIndex(0);
            break;
        case tgt::MouseEvent::MOUSE_BUTTON_MIDDLE:
            buttonBox_->setCurrentIndex(1);
            break;
        case tgt::MouseEvent::MOUSE_BUTTON_RIGHT:
            buttonBox_->setCurrentIndex(2);
            break;
        case tgt::MouseEvent::MOUSE_WHEEL_UP:
            buttonBox_->setCurrentIndex(3);
            break;
        case tgt::MouseEvent::MOUSE_WHEEL_DOWN:
            buttonBox_->setCurrentIndex(3);
            break;
        case tgt::MouseEvent::MOUSE_WHEEL:
            buttonBox_->setCurrentIndex(3);
            break;
        case tgt::MouseEvent::MOUSE_BUTTON_ALL:
            buttonBox_->setCurrentIndex(4);
            break;
        case tgt::MouseEvent::MOUSE_BUTTON_NONE:
        default:
            buttonBox_->setCurrentIndex(5);
            break;
    }
    connect(buttonBox_, SIGNAL(currentIndexChanged(int)), this, SLOT(buttonChanged(int)));
    layout_->addWidget(buttonBox_);
    layout_->addStretch();

    buttonBox_->setFixedWidth(125);
}

void EventPropertyWidget::createKeyWidgets() {
    modifierWidget_ = new ModifierDetectorWidget;
    modifierWidget_->setText(ModifierDetectorWidget::getStringForModifier(ModifierDetectorWidget::getQtModifierFromTGT(property_->getModifier())));
    modifierWidget_->setFixedWidth(60);
    connect(modifierWidget_, SIGNAL(modifierChanged(Qt::KeyboardModifiers)), this, SLOT(modifierChanged(Qt::KeyboardModifiers)));
    layout_->addWidget(modifierWidget_);

    keyWidget_ = new KeyDetectorWidget;
    keyWidget_->setText(KeyDetectorWidget::getStringForKey(KeyDetectorWidget::getQtKeyFromTGT(property_->getKeyCode())));
    keyWidget_->setFixedWidth(60);
    connect(keyWidget_, SIGNAL(keyChanged(int)), this, SLOT(keyChanged(int)));
    layout_->addWidget(keyWidget_);
    layout_->addStretch();
}

void EventPropertyWidget::createEnabledBox() {
    tgtAssert(property_, "No property");
    checkEnabled_ = new QCheckBox(QString::fromStdString(property_->getGuiName()));
    checkEnabled_->setChecked(property_->isEnabled());
    checkEnabled_->setFixedWidth(150);
    connect(checkEnabled_, SIGNAL(toggled(bool)), this, SLOT(enabledChanged(bool)));
    layout_->addWidget(checkEnabled_);
}

void EventPropertyWidget::createSharingBox() {
    tgtAssert(property_, "No property");
    checkSharing_ = new QCheckBox("Sharing");
    checkSharing_->setChecked(property_->isSharing());
    connect(checkSharing_, SIGNAL(toggled(bool)), this, SLOT(sharingChanged(bool)));
    layout_->addWidget(checkSharing_);
}

void EventPropertyWidget::modifierChanged(Qt::KeyboardModifiers modifier) {
    if (disconnected_)
        return;

    tgt::Event::Modifier m = ModifierDetectorWidget::getTGTModifierFromQt(modifier);
    property_->setModifier(m);
}

void EventPropertyWidget::keyChanged(int key) {
    if (disconnected_)
        return;

    property_->setKeyCode(KeyDetectorWidget::getTGTKeyFromQt(key));
}

void EventPropertyWidget::buttonChanged(int button) {
    if (disconnected_)
        return;

    switch (button) {
    case 0:
        property_->setMouseButtons(tgt::MouseEvent::MOUSE_BUTTON_LEFT);
        break;
    case 1:
        property_->setMouseButtons(tgt::MouseEvent::MOUSE_BUTTON_MIDDLE);
        break;
    case 2:
        property_->setMouseButtons(tgt::MouseEvent::MOUSE_BUTTON_RIGHT);
        break;
    case 3:
        property_->setMouseButtons(tgt::MouseEvent::MOUSE_WHEEL);
        break;
    case 4:
        property_->setMouseButtons(tgt::MouseEvent::MOUSE_BUTTON_ALL);
        break;
    case 5:
        property_->setMouseButtons(tgt::MouseEvent::MOUSE_BUTTON_NONE);
        break;
    }
}

void EventPropertyWidget::enabledChanged(bool enabled) {
    if (disconnected_)
        return;

    property_->setEnabled(enabled);
    if (modifierWidget_)
        modifierWidget_->setEnabled(enabled);
    if (keyWidget_)
        keyWidget_->setEnabled(enabled);
    if (buttonBox_)
        buttonBox_->setEnabled(enabled);
    if (checkSharing_)
        checkSharing_->setEnabled(enabled);
}

void EventPropertyWidget::sharingChanged(bool shared) {
    if (disconnected_)
        return;

    property_->setSharing(shared);
}

void EventPropertyWidget::adjustWidgetState() {

}

void EventPropertyWidget::setEnabled(bool enabled) {
    QWidget::setEnabled(enabled);
}

void EventPropertyWidget::setVisible(bool state) {
    QWidget::setVisible(state);
}

void EventPropertyWidget::disconnect() {
    disconnected_ = true;
    property_ = 0;
}

void EventPropertyWidget::updateFromProperty() {
    if (disconnected_)
        return;
    bool enabled = property_->isEnabled();
    //set widget's enabled state
    if (modifierWidget_)
        modifierWidget_->setEnabled(enabled);
    if (keyWidget_)
        keyWidget_->setEnabled(enabled);
    if (buttonBox_)
        buttonBox_->setEnabled(enabled);
    if (checkSharing_)
        checkSharing_->setEnabled(enabled);
    //set checkbox state
    if (checkEnabled_)
        checkEnabled_->setChecked(enabled);
    //set visibility
    setVisible(property_->isVisibleFlagSet());
}

} // namespace
