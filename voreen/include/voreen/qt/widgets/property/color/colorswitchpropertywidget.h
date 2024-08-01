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

#ifndef VRN_COLORSWITCHPROPERTYWIDGET_H
#define VRN_COLORSWITCHPROPERTYWIDGET_H

#include "voreen/core/properties/color/colorswitchproperty.h"

#include "voreen/qt/widgets/property/qpropertywidget.h"
#include "voreen/qt/widgets/colorselectorwidget.h"

#include "tgt/vector.h"

#include <QLabel>
#include <QCheckBox>

namespace voreen {

class ColorSwitchPropertyWidget : public QPropertyWidget {
    Q_OBJECT
public:
    ColorSwitchPropertyWidget(ColorSwitchProperty* prop, QWidget* parent = 0);

public slots:
    void changeActiveColor(QColor color);
    void changeInactiveColor(QColor color);

protected:
    tgt::Color toTgtColor(QColor color);
    QColor toQColor(tgt::Color color);

protected slots:
    virtual void updateFromPropertySlot();
    void useActiveColorClicked(int);

private:
    ColorSwitchProperty* property_;
    ColorSelectorWidget* activeColorLabel_;
    ColorSelectorWidget* inactiveColorLabel_;
    QCheckBox* useActiveColor_;
};

} // namespace

#endif // VRN_COLORPROPERTYWIDGET_H
