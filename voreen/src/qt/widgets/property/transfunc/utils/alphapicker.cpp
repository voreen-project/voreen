/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2020 University of Muenster, Germany,                        *
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

#include "voreen/qt/widgets/property/transfunc/utils/alphapicker.h"

#include <QColor>
#include <QMouseEvent>
#include <QPainter>
#include <QPaintEvent>
#include <qdrawutil.h>

namespace voreen {

AlphaPicker::AlphaPicker(QWidget* parent)
    : QWidget(parent)
    , hue_(0), sat_(0), val_(0), alpha_(0), pix_(0)
{
    setSizePolicy(QSizePolicy(QSizePolicy::Fixed, QSizePolicy::Preferred));
    setMinimumSize(10, 50);
}

AlphaPicker::~AlphaPicker() {
    delete pix_;
}
//-------------------------------------------------------------------------------
//  Slots Functions
//-------------------------------------------------------------------------------
void AlphaPicker::updateHSVSlot(int h, int s, int v) {
    updateAllMembers(h, s, v, alpha_);
}

void AlphaPicker::setHSVAByColorSlot(const QColor& c) {
    updateAllMembers(c.hue(), c.saturation(), c.value(),c.alpha());
}

//-------------------------------------------------------------------------------
//  Converter Functions
//-------------------------------------------------------------------------------
int AlphaPicker::yPositionToAlpha(int a) {
    int d = height() - 2*C_OFF_ - 1;
    if (d != 0)
        return 255 - (a - C_OFF_)*255/d;
    else
        return 0;
}

int AlphaPicker::alphaToYPosition(int a) {
    int d = height() - 2*C_OFF_ - 1;
    return C_OFF_ + (255-a)*d/255;
}

void AlphaPicker::updateAlpha(int a) {
    //check, if a update is needed
    int nAlpha = qMax(0, qMin(a,255));
    if (alpha_ == nAlpha)
        return;

    alpha_ = nAlpha;
    repaint();
    emit newHSVASignal(hue_, sat_, val_, alpha_);
}

void AlphaPicker::updateAllMembers(int h, int s , int v, int a) {
    bool changed = false;
    int nhue   = qMin(qMax(0,h), 359);
    int nsat   = qMin(qMax(0,s), 255);
    int nval   = qMin(qMax(0,v), 255);
    int nalpha = qMin(qMax(0,a), 255);

    if(hue_ != nhue || sat_ != nsat || val_ != nval) {
        hue_ = nhue; sat_ = nsat; val_ = nval;
        delete pix_; pix_ = 0;
        changed = true;
    }
    if(alpha_ != nalpha) {
        alpha_ = nalpha;
        changed = true;
    }
    if(changed) {
        repaint();
        emit newHSVASignal(hue_, sat_, val_, alpha_);
    }
}

//-------------------------------------------------------------------------------
//  Event Functions
//-------------------------------------------------------------------------------
void AlphaPicker::mousePressEvent(QMouseEvent* event) {
    emit toggleInteractionModeSignal(true);
    updateAlpha(yPositionToAlpha(event->y()));
}

void AlphaPicker::mouseMoveEvent(QMouseEvent* event) {
    updateAlpha(yPositionToAlpha(event->y()));
}

void AlphaPicker::mouseReleaseEvent(QMouseEvent* event) {
    event->accept();
    emit toggleInteractionModeSignal(false);
}

void AlphaPicker::paintEvent(QPaintEvent* /*event*/) {
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
        img.fill(Qt::white);
        //draw checker board
        QPainter painter(&img);
        painter.setPen(QColor(220, 220, 220));
        painter.setBrush(QColor(220, 220, 220));
        painter.setRenderHint(QPainter::Antialiasing, false);
        int checkerWi = wi/2, checkerHi = hi/7;
        for(int i = 0; i < 7; i++) {
            painter.drawRect((i%2)*checkerWi, i*checkerHi,checkerWi,checkerHi);
        }
        //draw alpha
        painter.setRenderHint(QPainter::Antialiasing, true);
        QColor c; c.setHsv(hue_, sat_, val_);
        for (int y = 0 ; y < hi ; ++y) {
            c.setAlpha(yPositionToAlpha(y+C_OFF_));
            painter.setPen(c);
            painter.drawLine(0,y,wi-1,y);
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
    int y = alphaToYPosition(alpha_);
    a.setPoints(3, w, y, w+5, y+5, w+5, y-5);
    p.eraseRect(w, 0, 5, height());
    p.drawPolygon(a);
}

} // namespace voreen
