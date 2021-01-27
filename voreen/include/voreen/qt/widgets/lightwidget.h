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

#ifndef VRN_LIGHTWIDGET_H
#define VRN_LIGHTWIDGET_H

#include <QWidget>
#include <tgt/vector.h>

#include "tgt/vector.h"

class DoubleSliderSpinBoxWidget;
class QGridLayout;
class QPushButton;

namespace voreen {

class LightWidget : public QWidget {
Q_OBJECT
public:
    LightWidget(QWidget* parent = 0);
    QPointF direction();

    bool getHemisphereStatus() const;

public slots:
    void setLightPosition(const tgt::vec4&);
    tgt::vec4 getLightPosition() const;
    void setMinDist(double);
    void setMaxDist(double);
    float getMaxDist() const;

protected:
    void mousePressEvent(QMouseEvent*);
    void mouseMoveEvent(QMouseEvent*);
    void paintEvent(QPaintEvent*);

private:
    int radius_;
    bool hemisphere_;
    QPointF lightPosition_;
    float distance_;

    tgt::vec2 offset_;

    QPushButton* hemisphereButton_;
    DoubleSliderSpinBoxWidget* distanceSlider_;

private slots:
    void setHemisphere(bool); //true is an open sphere, false is a closed sphere
    void updateDistance(double);

signals:
    void lightWidgetChanged(tgt::vec4);
};

} // namespace voreen

#endif // VRN_LIGHTWIDGET_H
