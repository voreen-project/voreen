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

#ifndef VRN_INTPROPERTYWIDGET_H
#define VRN_INTPROPERTYWIDGET_H

#include "voreen/qt/widgets/property/qpropertywidget.h"

namespace voreen {

class IntProperty;
class SliderSpinBoxWidget;

class IntPropertyWidget : public QPropertyWidget {
Q_OBJECT;
public:
    IntPropertyWidget(IntProperty* prop, QWidget* parent = 0, bool addVisibilityControl = true);
    virtual ~IntPropertyWidget();
    void setWidget(const int value, const int minValue, const int maxValue, const int stepping);

public slots:
    void setProperty(int value);

signals:
    void valueChanged(int);

protected:
    IntProperty* property_;
    SliderSpinBoxWidget* widget_;
    virtual void mousePressEvent(QMouseEvent*);
    QMenu* instantValueChangeMenu_;
    QAction* instantValueChangeAction_;

protected slots:
    virtual void updateFromPropertySlot();

};

} // namespace

#endif // VRN_INTPROPERTYWIDGET_H
