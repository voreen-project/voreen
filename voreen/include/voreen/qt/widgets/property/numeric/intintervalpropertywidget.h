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

#ifndef VRN_INTINTERVALPROPERTYWIDGET_H
#define VRN_INTINTERVALPROPERTYWIDGET_H

#include "voreen/qt/widgets/property/qpropertywidget.h"
#include "voreen/core/properties/numeric/intervalproperty.h"
namespace voreen {


class SliderSpinBoxWidget;

class IntIntervalPropertyWidget : public QPropertyWidget {
Q_OBJECT;
public:
    IntIntervalPropertyWidget(IntIntervalProperty* prop, QWidget* parent = 0, bool addVisibilityControl = true);
    virtual ~IntIntervalPropertyWidget();
    //void setWidget(const int value, const int minValue, const int maxValue, const int stepping);

public slots:
    void setPropertyMin(int value);
    void setPropertyMax(int value);

signals:
    void valueChanged(tgt::ivec2);

protected:
    IntIntervalProperty* property_;
    SliderSpinBoxWidget* minWidget_;
    SliderSpinBoxWidget* maxWidget_;
    QMenu* instantValueChangeMenu_;
    QAction* instantValueChangeAction_;
    QVBoxLayout* vboxLayout_;

protected slots:
    virtual void updateFromPropertySlot();

private:


};

} // namespace

#endif // VRN_INTPROPERTYWIDGET_H
