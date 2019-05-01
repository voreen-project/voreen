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

#include "colormappropertywidget.h"

#include "../../properties/colormapproperty.h"

#include <QColorDialog>
#include <QComboBox>
#include <QLabel>
#include <QMouseEvent>
#include <QPainter>
#include <QGradient>
#include <QPushButton>
namespace voreen {

ClickableColorMapLabel::ClickableColorMapLabel(const QString& text, QWidget* parent, Qt::WindowFlags f)
    : QLabel(text, parent, f)
    , cm_(ColorMap::createColdHot())
{}

void ClickableColorMapLabel::setColorMap(const ColorMap& cm) {
    cm_ = cm;
    update();
}

void ClickableColorMapLabel::mousePressEvent(QMouseEvent* e) {
    if (e->button() == Qt::LeftButton)
        emit(clicked());
}

void ClickableColorMapLabel::paintEvent(QPaintEvent* /*event*/) {
    QPainter painter(this);

    QLinearGradient backGrad(QPointF(1, 1), QPointF(rect().width() - 2, 1));
    int count = cm_.getColorCount()-1;
    for (int i=0; i <= count; ++i) {
        backGrad.setColorAt(qreal(i)/qreal(count), toQColor(cm_.getColorAtIndex(i)));
    }

    painter.setBrush(backGrad);
    painter.drawRect(1, 1, rect().width() - 2, rect().height() - 2);
}

tgt::Color ClickableColorMapLabel::toTgtColor(QColor color) {
    return tgt::Color(color.redF(),color.greenF(),color.blueF(), color.alphaF());
}

QColor ClickableColorMapLabel::toQColor(tgt::Color color) {
    return QColor(static_cast<int>(color.r * 255), static_cast<int>(color.g * 255),
                  static_cast<int>(color.b * 255), static_cast<int>(color.a * 255));
}




// ------------------------------------------------------------------------------------------------

ColorMapPropertyWidget::ColorMapPropertyWidget(ColorMapProperty* prop, QWidget* parent)
    : QPropertyWidget(prop, parent)
    , property_(prop)
    , cbMaps_(new QComboBox(this))
    , colorLbl_(new ClickableColorMapLabel(""))
{
    updateColorLabel();

    addWidget(colorLbl_);

    if (cbMaps_) {
        customColorMap_ = property_->get();
        cbMaps_->addItem("custom");
        for(const std::string& label : ColorMap::getColorMapLabels()) {
            cbMaps_->addItem(label.c_str());
        }
    }
    addWidget(cbMaps_);
    connect(cbMaps_, SIGNAL(currentIndexChanged(int)), this, SLOT(selectColormap(int)));

}

void ColorMapPropertyWidget::updateFromPropertySlot() {
    updateColorLabel();
}

void ColorMapPropertyWidget::updateColorLabel() {
    colorLbl_->setColorMap(property_->get());
}

void ColorMapPropertyWidget::selectColormap(int index) {
    if(index == 0) {
        property_->set(customColorMap_);
    }
    else {
        property_->set(ColorMap::createColorMap(index - 1));
    }
}


tgt::Color ColorMapPropertyWidget::toTgtColor(QColor color) {
    return tgt::Color(color.redF(),color.greenF(),color.blueF(), color.alphaF());
}

QColor ColorMapPropertyWidget::toQColor(tgt::Color color) {
    return QColor(static_cast<int>(color.r * 255), static_cast<int>(color.g * 255),
                  static_cast<int>(color.b * 255), static_cast<int>(color.a * 255));
}

} // namespace
