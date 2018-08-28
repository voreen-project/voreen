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

#include "voreen/qt/widgets/property/voxeltypepropertywidget.h"

#include "voreen/core/properties/voxeltypeproperty.h"
#include "voreen/qt/widgets/sliderspinboxwidget.h"

#include "voreen/core/datastructures/volume/volumeram.h"
#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"

#include <QHBoxLayout>
#include <QMouseEvent>
#include <QMenu>
#include <QWidget>

namespace voreen {

VoxelTypePropertyWidget::VoxelTypePropertyWidget(VoxelTypeProperty* prop, QWidget* parent, bool addVisibilityControl)
    : QPropertyWidget(prop, parent)
    , property_(prop)
    , volume_(0)
    , widget_(new QWidget(this))
{
    tgtAssert(property_, "no property");

    widget_->setLayout(new QVBoxLayout());
    addWidget(widget_);
    createWidgets();

    updateFromPropertySlot();

    /*connect(widget_, SIGNAL(valueChanged(int)), this, SLOT(setProperty(int)));
    connect(widget_, SIGNAL(sliderPressedChanged(bool)), this, SLOT(toggleInteractionMode(bool)));
    connect(widget_, SIGNAL(valueChanged(int)), this, SIGNAL(widgetChanged())); */

    instantValueChangeMenu_ = new QMenu(this);
    instantValueChangeAction_ = instantValueChangeMenu_->addAction("Tracking Mode");
    instantValueChangeAction_->setCheckable(true);
/*    if (property_->hasTracking())
        instantValueChangeAction_->toggle(); */
}

VoxelTypePropertyWidget::~VoxelTypePropertyWidget() {
    widget_->disconnect();
    delete widget_;
}

void VoxelTypePropertyWidget::updateFromPropertySlot() {
    tgtAssert(property_, "no property");
    setVolume(property_->getVolume());
}

void VoxelTypePropertyWidget::mousePressEvent(QMouseEvent* event) {
    if(event->button() == Qt::RightButton) {
        /*QAction* prec = instantValueChangeMenu_->exec(mapToGlobal(event->pos()));
        if(prec == instantValueChangeAction_) {
            property_->setTracking(!property_->hasTracking());
        } */
        updateFromPropertySlot();
    }
    QWidget::mousePressEvent(event);
}

void VoxelTypePropertyWidget::setVolume(const VolumeBase* volumeHandle) {
    volume_ = volumeHandle;
    updateWidgetFromProperty();
}

void VoxelTypePropertyWidget::createWidgets() {
    for (int i=0; i<3; i++) {
        intWidgets_.push_back(new SliderSpinBoxWidget(this));
        intWidgets_.back()->setVisible(false);
        connect(intWidgets_.back(), SIGNAL(valueChanged(int)), this, SLOT(updatePropertyFromWidget()));
    }
    for (int i=0; i<3; i++) {
        floatWidgets_.push_back(new DoubleSliderSpinBoxWidget(this));
        floatWidgets_.back()->setSingleStep(0.01f);
        floatWidgets_.back()->setVisible(false);
        connect(floatWidgets_.back(), SIGNAL(valueChanged(double)), this, SLOT(updatePropertyFromWidget()));
    }
}

void VoxelTypePropertyWidget::updateWidgetFromProperty() {
    if (widget_->layout()) {
        while (widget_->layout()->count() > 0)
            widget_->layout()->takeAt(0);
    }
    else {
        widget_->setLayout(new QHBoxLayout());
    }

    for (size_t i=0; i<intWidgets_.size(); i++)
        intWidgets_.at(i)->setVisible(false);
    for (size_t i=0; i<floatWidgets_.size(); i++)
        floatWidgets_.at(i)->setVisible(false);

    if (!volume_)
        return;

    tgtAssert(widget_->layout(), "widget has no layout");
    const VolumeRAM* volume = volume_->getRepresentation<VolumeRAM>();
    if (dynamic_cast<const VolumeRAM_UInt8*>(volume) || dynamic_cast<const VolumeRAM_Int8*>(volume) ||
        dynamic_cast<const VolumeRAM_UInt16*>(volume) || dynamic_cast<const VolumeRAM_Int16*>(volume) ||
        dynamic_cast<const VolumeRAM_UInt32*>(volume) || dynamic_cast<const VolumeRAM_Int32*>(volume) ) {
            // voxel type: integer (scalar)
            tgtAssert(intWidgets_.size() > 0, "no int widget");
            SliderSpinBoxWidget* intWidget = intWidgets_.at(0);
            widget_->layout()->addWidget(intWidget);
            intWidget->setVisible(true);

            intWidget->blockSignals(true);
            intWidget->setMinValue(property_->getMinValue<int>());
            intWidget->setMaxValue(property_->getMaxValue<int>());
            intWidget->setValue(property_->getValue<int>());
            intWidget->blockSignals(false);
    }
    else if (dynamic_cast<const VolumeRAM_3xUInt8*>(volume) || dynamic_cast<const VolumeRAM_3xInt8*>(volume) ||
            dynamic_cast<const VolumeRAM_3xUInt16*>(volume) || dynamic_cast<const VolumeRAM_3xInt16*>(volume) ) {
        // voxel type: RGB integer
        tgtAssert(intWidgets_.size() >= 3, "less than three int widgets");
        tgt::ivec3 value = property_->getValue<tgt::ivec3>();
        tgt::ivec3 minValue = property_->getMinValue<tgt::ivec3>();
        tgt::ivec3 maxValue = property_->getMaxValue<tgt::ivec3>();

        for (int i=0; i<3; i++) {
            widget_->layout()->addWidget(intWidgets_.at(i));
            intWidgets_.at(i)->setVisible(true);
            intWidgets_.at(i)->blockSignals(true);
            intWidgets_.at(i)->setMinValue(minValue[i]);
            intWidgets_.at(i)->setMaxValue(maxValue[i]);
            intWidgets_.at(i)->setValue(value[i]);
            intWidgets_.at(i)->blockSignals(false);
        }
    }
    else if (dynamic_cast<const VolumeRAM_Float*>(volume) || dynamic_cast<const VolumeRAM_Double*>(volume)) {
        // voxel type: float (scalar)
        tgtAssert(floatWidgets_.size() > 0, "no float widget");
        DoubleSliderSpinBoxWidget* floatWidget = floatWidgets_.at(0);
        widget_->layout()->addWidget(floatWidget);
        floatWidget->setVisible(true);

        floatWidget->blockSignals(true);
        floatWidget->setMinValue(property_->getMinValue<float>());
        floatWidget->setMaxValue(property_->getMaxValue<float>());
        floatWidget->setSingleStep((property_->getMaxValue<float>() - property_->getMinValue<float>()) / 100.0);
        floatWidget->setValue(property_->getValue<float>());
        floatWidget->blockSignals(false);
    }
    else if (dynamic_cast<const VolumeRAM_3xFloat*>(volume) || dynamic_cast<const VolumeRAM_3xDouble*>(volume)) {
        // voxel type: RGB float
        tgtAssert(floatWidgets_.size() >= 3, "less than three float widgets");
        tgt::vec3 value = property_->getValue<tgt::vec3>();
        tgt::vec3 minValue = property_->getMinValue<tgt::vec3>();
        tgt::vec3 maxValue = property_->getMaxValue<tgt::vec3>();

        for (int i=0; i<3; i++) {
            widget_->layout()->addWidget(floatWidgets_.at(i));
            floatWidgets_.at(i)->setVisible(true);
            floatWidgets_.at(i)->blockSignals(true);
            floatWidgets_.at(i)->setMinValue(minValue[i]);
            floatWidgets_.at(i)->setMaxValue(maxValue[i]);
            floatWidgets_.at(i)->setSingleStep((property_->getMaxValue<tgt::vec3>()[i] - property_->getMinValue<tgt::vec3>()[i]) / 100.0);
            floatWidgets_.at(i)->setValue(value[i]);
            floatWidgets_.at(i)->blockSignals(false);
        }
    }
    else {
        LWARNINGC("voreen.voreenqt.VoxelTypePropertyWidget",
            "updateWidgetFromProperty() unknown volume type: " << typeid(volume).name());
    }
}

void VoxelTypePropertyWidget::updatePropertyFromWidget() const {
    const VolumeRAM* volume = volume_->getRepresentation<VolumeRAM>();
    if (dynamic_cast<const VolumeRAM_UInt8*>(volume) || dynamic_cast<const VolumeRAM_Int8*>(volume) ||
            dynamic_cast<const VolumeRAM_UInt16*>(volume) || dynamic_cast<const VolumeRAM_Int16*>(volume) ||
            dynamic_cast<const VolumeRAM_UInt32*>(volume) || dynamic_cast<const VolumeRAM_Int32*>(volume) ) {
        // voxel type: integer (scalar)
        tgtAssert(property_, "no property");
        tgtAssert(intWidgets_.size() > 0, "no int widget");
        SliderSpinBoxWidget* intWidget = intWidgets_.at(0);
        property_->setValue<int>(intWidget->getValue());
    }
    else if (dynamic_cast<const VolumeRAM_3xUInt8*>(volume) || dynamic_cast<const VolumeRAM_3xInt8*>(volume) ||
        dynamic_cast<const VolumeRAM_3xUInt16*>(volume) || dynamic_cast<const VolumeRAM_3xInt16*>(volume) ) {
            // voxel type: RGB integer
            tgtAssert(intWidgets_.size() >= 3, "less than three int widgets");
            tgt::ivec3 value;
            for (int i=0; i<3; i++) {
                value[i] = intWidgets_.at(i)->getValue();
            }
            property_->setValue<tgt::ivec3>(value);
    }
    else if (dynamic_cast<const VolumeRAM_Float*>(volume) || dynamic_cast<const VolumeRAM_Double*>(volume)) {
        // voxel type: float (scalar)
        tgtAssert(property_, "no property");
        tgtAssert(floatWidgets_.size() > 0, "no float widget");
        DoubleSliderSpinBoxWidget* floatWidget = floatWidgets_.at(0);
        property_->setValue<float>(floatWidget->getValue());
    }
    else if (dynamic_cast<const VolumeRAM_3xFloat*>(volume) || dynamic_cast<const VolumeRAM_3xDouble*>(volume)) {
        // voxel type: RGB float
        tgtAssert(floatWidgets_.size() >= 3, "less than three int widgets");
        tgt::vec3 value;
        for (int i=0; i<3; i++) {
            value[i] = static_cast<float>(floatWidgets_.at(i)->getValue());
        }
        property_->setValue<tgt::vec3>(value);
    }
    else {
        LWARNINGC("voreen.voreenqt.VoxelTypePropertyWidget",
            "updatePropertyFromWidget(): unknown volume type: " << typeid(volume).name());
    }
}

} // namespace
