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

#include "voreen/qt/widgets/property/transfunc/utils/colorluminancepicker.h"

#include <QColor>
#include <QMouseEvent>
#include <QPainter>
#include <QPaintEvent>
#include <qdrawutil.h>

namespace voreen {

ColorLuminancePicker::ColorLuminancePicker(QWidget* parent)
    : QWidget(parent)
    , hue_(0)
    , sat_(0)
    , val_(0)
    , pix_(0)
{
    setSizePolicy(QSizePolicy(QSizePolicy::Fixed, QSizePolicy::Preferred));
    setMinimumSize(10, 50);
}

ColorLuminancePicker::~ColorLuminancePicker() {
    delete pix_;
}

//-------------------------------------------------------------------------------
//  Slots Functions
//-------------------------------------------------------------------------------
void ColorLuminancePicker::updateHSSlot(int h, int s) {
    updateAllMembers(h, s, val_);
}

void ColorLuminancePicker::setHSVByColorSlot(const QColor& c) {
    updateAllMembers(c.hue(), c.saturation(), c.value());
}

//-------------------------------------------------------------------------------
//  Converter Functions
//-------------------------------------------------------------------------------
int ColorLuminancePicker::yPositionToValue(int y) {
    int d = height() - 2*C_OFF_ - 1;
    if (d != 0)
        return 255 - (y - C_OFF_)*255/d;
    else
        return 0;
}

int ColorLuminancePicker::valueToYPosition(int v) {
    int d = height() - 2*C_OFF_ - 1;
    return C_OFF_ + (255-v)*d/255;
}

void ColorLuminancePicker::updateValue(int v) {
    //check, if a update is needed
    int nVal = qMax(0, qMin(v,255));
    if (val_ == nVal)
        return;

    val_ = nVal;
    repaint();
    emit newHSVSignal(hue_, sat_, val_);
}

void ColorLuminancePicker::updateAllMembers(int h, int s , int v) {
    bool changed = false;
    int nhue = qMin(qMax(0,h), 359);
    int nsat = qMin(qMax(0,s), 255);
    int nval = qMin(qMax(0,v), 255);

    if(hue_ != nhue || sat_ != nsat) {
        hue_ = nhue;
        sat_ = nsat;
        delete pix_;
        pix_ = 0;
        changed = true;
    }
    if(val_ != nval) {
        val_ = nval;
        changed = true;
    }
    if(changed) {
        repaint();
        emit newHSVSignal(hue_, sat_, val_);
    }
}

//-------------------------------------------------------------------------------
//  Event Functions
//-------------------------------------------------------------------------------
void ColorLuminancePicker::mousePressEvent(QMouseEvent* event) {
    emit toggleInteractionModeSignal(true);
    updateValue(yPositionToValue(event->y()));
}

void ColorLuminancePicker::mouseMoveEvent(QMouseEvent* event) {
    updateValue(yPositionToValue(event->y()));
}

void ColorLuminancePicker::mouseReleaseEvent(QMouseEvent* event) {
    event->accept();
    emit toggleInteractionModeSignal(false);
}

void ColorLuminancePicker::paintEvent(QPaintEvent* /*event*/) {
    //calc size of current picker
    int w = width() - 5;
    QRect r(0, F_OFF_, w, height() - 2*F_OFF_);
    int wi = r.width() - 2;
    int hi = r.height() - 2;
    if (wi <= 1 || hi <= 1)
        return;
    //update texture if needed
    if ((pix_ == 0) || (pix_->height() != hi) || (pix_->width() != wi)) {
        QImage img(wi, hi, QImage::Format_RGB32);
        for (int y = 0 ; y < hi ; ++y) {
            QColor c;
            c.setHsv(hue_, sat_, yPositionToValue(y+C_OFF_));
            QRgb r = c.rgb();
            for (int x = 0 ; x < wi ; ++x)
                img.setPixel(x, y, r);
        }
        delete pix_;
        pix_ = new QPixmap(QPixmap::fromImage(img));
    }

    // color bar
    QPainter p(this);
    p.drawPixmap(1, C_OFF_, *pix_);
    const QPalette &g = palette();
    qDrawShadePanel(&p, r, g, true);
    p.setPen(g.foreground().color());
    p.setBrush(g.foreground());

    // arrow
    QPolygon a;
    int y = valueToYPosition(val_);
    a.setPoints(3, w, y, w+5, y+5, w+5, y-5);
    p.eraseRect(w, 0, 5, height());
    p.drawPolygon(a);
}

} // namespace voreen
