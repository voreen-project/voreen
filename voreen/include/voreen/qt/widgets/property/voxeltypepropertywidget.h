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

#ifndef VRN_VOXELTYPEPROPERTYWIDGET_H
#define VRN_VOXELTYPEPROPERTYWIDGET_H

#include "voreen/qt/widgets/property/qpropertywidget.h"
#include <QWidget>

namespace voreen {

class VolumeBase;
class VoxelTypeProperty;

class SliderSpinBoxWidget;
class DoubleSliderSpinBoxWidget;

class VoxelTypePropertyWidget : public QPropertyWidget {
Q_OBJECT;
public:
    VoxelTypePropertyWidget(VoxelTypeProperty* prop, QWidget* parent = 0, bool addVisibilityControl = true);
    virtual ~VoxelTypePropertyWidget();

    void setVolume(const VolumeBase* volume);

signals:
    void valueChanged(int);

protected:
    void createWidgets();
    void updateWidgetFromProperty();

    virtual void mousePressEvent(QMouseEvent*);

    VoxelTypeProperty* property_;
    const VolumeBase* volume_;

    QWidget* widget_;
    std::vector<SliderSpinBoxWidget*> intWidgets_;
    std::vector<DoubleSliderSpinBoxWidget*> floatWidgets_;

    QMenu* instantValueChangeMenu_;
    QAction* instantValueChangeAction_;

protected slots:
    void updatePropertyFromWidget() const;
    void updateFromPropertySlot();

};

} // namespace

#endif // VRN_VOXELTYPEPROPERTYWIDGET_H
