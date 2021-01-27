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

#ifndef VRN_INTBOUNDINGBOXPROPERTYWIDGET_H
#define VRN_INTBOUNDINGBOXPROPERTYWIDGET_H

#include "voreen/qt/widgets/property/qpropertywidget.h"
#include "voreen/qt/widgets/property/numeric/intintervalpropertywidget.h"
#include "voreen/core/properties/numeric/intervalproperty.h"
#include "voreen/core/properties/boundingboxproperty.h"


#include <QBoxLayout>
#include <QGroupBox>
#include <QLabel>

class QMenu;

namespace voreen {

class IntBoundingBoxPropertyWidget : public QPropertyWidget {
    Q_OBJECT
public:
    IntBoundingBoxPropertyWidget(IntBoundingBoxProperty* const prop, QWidget* parent = 0);
    virtual ~IntBoundingBoxPropertyWidget();
protected:
    IntIntervalProperty xProperty_;
    IntIntervalProperty yProperty_;
    IntIntervalProperty zProperty_;

    IntIntervalPropertyWidget* xWidget_;
    IntIntervalPropertyWidget* yWidget_;
    IntIntervalPropertyWidget* zWidget_;

    QLabel* xLabel_;
    QLabel* yLabel_;
    QLabel* zLabel_;

    QVBoxLayout* vbox_;
    QVBoxLayout* vboxX_;
    QVBoxLayout* vboxY_;
    QVBoxLayout* vboxZ_;
    QGroupBox*   groupBoxX_;
    QGroupBox*   groupBoxY_;
    QGroupBox*   groupBoxZ_;
protected:

    IntBoundingBoxProperty* prop_;

protected slots:
    virtual void updateFromPropertySlot();
    void valueChanged(tgt::ivec2);

};

}   // namespace

#endif
