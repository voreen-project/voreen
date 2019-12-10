/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2019 University of Muenster, Germany,                        *
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

#include "voreen/qt/widgets/colorselectorwidget.h"

#include <QColorDialog>
#include <QLabel>
#include <QMouseEvent>
#include <QPainter>

namespace voreen {

ColorSelectorWidget::ColorSelectorWidget(const QString& text, QWidget* parent, bool useAlphaChannel)
    : QLabel(text, parent)
    , useAlphaChannel_(useAlphaChannel)
{}

void ColorSelectorWidget::setColor(const QColor& color) {
    currentColor_ = color;
    update();
}

QColor ColorSelectorWidget::getColor() const {
    return currentColor_;
}

void ColorSelectorWidget::openColorDialog() {
    QColor col = QColorDialog::getColor(currentColor_, this, "Color",
                                       (useAlphaChannel_ ? QColorDialog::ShowAlphaChannel : static_cast<QColorDialog::ColorDialogOptions>(0)));
    if(col.isValid()){
        currentColor_ = col;
        emit colorChangedSignal(col);
        update();
    }
}


//-------------------------------------------------------------------------------------------------
//      Events
//-------------------------------------------------------------------------------------------------
void ColorSelectorWidget::mousePressEvent(QMouseEvent* e) {
    if (e->button() == Qt::LeftButton)
        openColorDialog();
}

void ColorSelectorWidget::paintEvent(QPaintEvent* /*event*/) {
    QPainter painter(this);
    painter.setBrush(currentColor_);
    painter.drawRect(1, 1, rect().width() - 2, rect().height() - 2);
}


} // namespace
