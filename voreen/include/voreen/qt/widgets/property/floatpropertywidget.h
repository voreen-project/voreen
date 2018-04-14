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

#ifndef VRN_FLOATPROPERTYWIDGET_H
#define VRN_FLOATPROPERTYWIDGET_H

#include "voreen/qt/widgets/property/qpropertywidget.h"

class QWidget;

namespace voreen {

class FloatProperty;
class DoubleSliderSpinBoxWidget;

/**
 * @see SliderSpinBoxWidget
 */
class FloatPropertyWidget : public QPropertyWidget {
Q_OBJECT;
public:
    FloatPropertyWidget(FloatProperty* prop, QWidget* parent = 0, bool addVisibilityControl = true);
    virtual ~FloatPropertyWidget();

    void setWidget(const float value, const float minValue, const float maxValue, const float stepping);

public slots:
    void setProperty(double value);

signals:
    void valueChanged(double);

protected:
    FloatProperty* property_;
    DoubleSliderSpinBoxWidget* widget_;
    void mousePressEvent(QMouseEvent*);
    QMenu* precisionMenu_;
    QAction* highAction_;
    QAction* instantValueChangeAction_;

protected slots:
    virtual void updateFromPropertySlot();

};

} // namespace voreen

#endif // VRN_FLOATPROPERTYWIDGET_H
