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

#ifndef VRN_PLOTZOOMPROPERTYWIDGET_H
#define VRN_PLOTZOOMPROPERTYWIDGET_H

#include "voreen/qt/widgets/property/qpropertywidget.h"

class QDoubleSpinBox;
class QLabel;
class QPushButton;

namespace voreen {

class PlotSelectionProperty;

class PlotZoomPropertyWidget : public QPropertyWidget {
Q_OBJECT
public:
    PlotZoomPropertyWidget(PlotSelectionProperty* prop, QWidget* parent = 0);

    void updateFromProperty();



protected:
    PlotSelectionProperty* property_;
    QLabel* lblX_;
    QLabel* lblY_;
    QLabel* lblZ_;

    QPushButton* btnPush_;
    QPushButton* btnPop_;

    QDoubleSpinBox* spXlow_;
    QDoubleSpinBox* spXhigh_;
    QDoubleSpinBox* spYlow_;
    QDoubleSpinBox* spYhigh_;
    QDoubleSpinBox* spZlow_;
    QDoubleSpinBox* spZhigh_;

private slots:
    void btnPushClicked();
    void btnPopClicked();
    void setSpinBoxRanges();

private:
    bool checkSanityOfInput();
};

} // namespace

#endif // VRN_PLOTZOOMPROPERTYWIDGET_H

