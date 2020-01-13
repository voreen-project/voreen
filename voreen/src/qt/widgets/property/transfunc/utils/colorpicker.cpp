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

#include "voreen/qt/widgets/property/transfunc/utils/colorpicker.h"

#include <QColor>
#include <QMouseEvent>
#include <QPainter>
#include <QPaintEvent>

namespace voreen {

ColorPicker::ColorPicker(QWidget* parent)
    : QFrame(parent)
    , hue_(0)
    , sat_(0)
    , pix_(0)
{
    setAttribute(Qt::WA_NoSystemBackground);
    setSizePolicy(QSizePolicy(QSizePolicy::Expanding, QSizePolicy::Preferred));
    setMinimumSize(128, 50);
}

ColorPicker::~ColorPicker() {
    delete pix_;
}

//-------------------------------------------------------------------------------
//  Slots Functions
//-------------------------------------------------------------------------------
void ColorPicker::setHSSlot(int h, int s) {
    //check, if values must be updated
    int nhue = qMin(qMax(0,h), 359);
    int nsat = qMin(qMax(0,s), 255);
    if ((nhue == hue_) && (nsat == sat_))
        return;
    //set values
    hue_ = nhue; sat_ = nsat;
    //update area
    repaint();
}

void ColorPicker::setHSByColorSlot(const QColor& c) {
    setHSSlot(c.hue(), c.saturation());
}

//-------------------------------------------------------------------------------
//  Converter Functions
//-------------------------------------------------------------------------------
QPoint ColorPicker::getCurrentPosition() {
    return QPoint((360-hue_)*(contentsRect().width()-1)/360, (255-sat_)*(contentsRect().height()-1)/255);
}

int ColorPicker::getHueFromPosition(const QPoint &pt) {
    if (contentsRect().width()-1 != 0)
        return 360 - pt.x()*360/(contentsRect().width()-1);
    else
        return 0;
}

int ColorPicker::getSaturationFromPosition(const QPoint &pt) {
    if (contentsRect().height()-1 != 0)
        return 255 - pt.y()*255/(contentsRect().height()-1) ;
    else
        return 0;
}

void ColorPicker::setColorByPosition(const QPoint &pt) {
    setHSSlot(getHueFromPosition(pt), getSaturationFromPosition(pt));
}

//-------------------------------------------------------------------------------
//  Event Functions
//-------------------------------------------------------------------------------
void ColorPicker::mousePressEvent(QMouseEvent* event) {
    QPoint p = event->pos() - contentsRect().topLeft();
    setColorByPosition(p);
    emit toggleInteractionModeSignal(true);
    emit newHSSignal(hue_, sat_);
}

void ColorPicker::mouseMoveEvent(QMouseEvent* event) {
    QPoint p = event->pos() - contentsRect().topLeft();
    setColorByPosition(p);
    emit newHSSignal(hue_, sat_);
}

void ColorPicker::mouseReleaseEvent(QMouseEvent* event) {
    event->accept();
    emit toggleInteractionModeSignal(false);
}

void ColorPicker::paintEvent(QPaintEvent* /*event*/) {
    QPainter p(this);
    drawFrame(&p);
    QRect r = contentsRect();
    if ((r.width() <= 1) || (r.height() <= 1))
        return;

    int w = r.width();
    int h = r.height();
    if (!pix_ || (pix_->width() != w) || (pix_->height() != h)) {
        QImage img(w, h, QImage::Format_RGB32);
        img.fill(Qt::blue);
        for (int y = 0 ; y < h ; ++y) {
            for (int x = 0 ; x < w ; ++x) {
                QPoint p(x, y);
                QColor c;
                c.setHsv(getHueFromPosition(p), getSaturationFromPosition(p), VAL_);
                img.setPixel(x, y, c.rgb());
            }
        }
        delete pix_;
        pix_ = new QPixmap(QPixmap::fromImage(img));
    }
    p.drawPixmap(r.topLeft(), *pix_);

    // black cross
    QPoint pt = getCurrentPosition() + r.topLeft();
    p.setPen(Qt::black);
    p.fillRect(pt.x()-9, pt.y(), 20, 2, Qt::black);
    p.fillRect(pt.x(), pt.y()-9, 2, 20, Qt::black);
}

} // namespace voreen
