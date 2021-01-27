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

#include "voreen/core/datastructures/callback/memberfunctioncallback.h"
#include "voreen/core/properties/fontproperty.h"
#include "voreen/core/processors/volumeraycaster.h"
#include "voreen/core/voreenapplication.h"
#include "voreen/qt/widgets/property/fontpropertywidget.h"
#include "voreen/qt/widgets/sliderspinboxwidget.h"

#include <QMouseEvent>
#include <QPainter>
#include <QPaintEvent>
#include <QPushButton>
#include <QSlider>
#include <QStyle>
#include <QTabWidget>
#include <QWidget>
#include <QLCDNumber>
//#include <QColorGroup>
#include <QLabel>

namespace voreen {

FontPropertyWidget::FontPropertyWidget(FontProperty* prop, QWidget* parent)
    : QPropertyWidget(prop, parent)
    , property_(prop)
{
    colorPropertyWidget_ = new ColorPropertyWidget(&colorProperty_);
    colorProperty_.addWidget(colorPropertyWidget_);

    groupBox_ = new QGroupBox(this);
    addWidget(groupBox_);

    QVBoxLayout* tgtFontLayout = new QVBoxLayout(groupBox_);
    tgtFontLayout->setContentsMargins(2,2,2,2);

    tgtFontLayout->addWidget(colorPropertyWidget_);

    tgtFontName_ = new QComboBox();
    tgtFontName_->addItem("VeraMono", "VeraMono.ttf");
    tgtFontName_->addItem("Vera", "Vera.ttf");
    tgtFontLayout->addWidget(tgtFontName_);

    QLabel* fontSizeLabel = new QLabel("Fontsize:");
    tgtFontLayout->addWidget(fontSizeLabel);

    QWidget* tgtFontSize = new QWidget();
    QHBoxLayout* tgtFontSizeLayout = new QHBoxLayout(tgtFontSize);
    tgtFontSizeLayout->setContentsMargins(0,0,0,0);
    fontSizeSlider_ = new SliderSpinBoxWidget();
    tgtFontSizeLayout->addWidget(fontSizeSlider_);
    tgtFontLayout->addWidget(tgtFontSize);



    fontSizeSlider_->setMinValue(10);
    fontSizeSlider_->setMaxValue(100);

    updateFromPropertySlot();

    connect(tgtFontName_,    SIGNAL(activated(int)),    this, SLOT(updateProperty()));
    connect(fontSizeSlider_, SIGNAL(valueChanged(int)), this, SLOT(updateProperty()));
    colorProperty_.onChange(MemberFunctionCallback<FontPropertyWidget>(this, &FontPropertyWidget::updateProperty));


}

void FontPropertyWidget::updateFromPropertySlot() {

    std::string fontName = property_->get()->getFontName();
    if(fontName.length() > 0)
        fontName = fontName.substr(fontName.find_last_of("/\\"));
    for(int i=0; i<tgtFontName_->count(); i++)
        if(fontName.compare(tgtFontName_->itemData(i).toString().toStdString()) == 0)
            tgtFontName_->setCurrentIndex(i);

    fontSizeSlider_->blockSignals(true);
    fontSizeSlider_->setValue(property_->get()->getFontSize());
    fontSizeSlider_->blockSignals(false);

    if (colorProperty_.get() != property_->get()->getFontColor()){
        colorProperty_.set(property_->get()->getFontColor());
        colorProperty_.invalidate();
    }
}

void FontPropertyWidget::updateProperty() {
    property_->get()->setFontSize(fontSizeSlider_->getValue());
    property_->get()->setFontName(VoreenApplication::app()->getFontPath(tgtFontName_->itemData(tgtFontName_->currentIndex()).toString().toStdString()));
    property_->get()->setFontColor(colorProperty_.get());
    property_->invalidate();
}

} // namespace voreen
