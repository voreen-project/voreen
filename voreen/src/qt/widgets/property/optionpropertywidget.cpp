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

#include "voreen/qt/widgets/property/optionpropertywidget.h"

#include "voreen/core/properties/optionproperty.h"

#include <QComboBox>

namespace voreen {

OptionPropertyWidget::OptionPropertyWidget(OptionPropertyBase* prop, QWidget* parent)
    : QPropertyWidget(prop, parent)
    , property_(prop)
    , cBox_(new QComboBox)
{
    updateFromPropertySlot();
    QFontInfo fontInfo(font());
    cBox_->setFont(QFont(fontInfo.family(), QPropertyWidget::fontSize_));
    addWidget(cBox_);
    connect(cBox_, SIGNAL(currentIndexChanged(int)), this, SLOT(setProperty(int)));
    connect(cBox_, SIGNAL(currentIndexChanged(int)), this, SIGNAL(widgetChanged()));
}

void OptionPropertyWidget::updateFromPropertySlot() {
    cBox_->blockSignals(true);
    cBox_->clear();

    // build combo box from descriptions
    std::vector<std::string> keys = property_->getKeys();
    std::vector<std::string> descriptions = property_->getDescriptions();
    std::vector<tgt::col4>  colors = property_->getGUIColors();
    tgtAssert(keys.size() == descriptions.size(), "Dimension mismatch"); // can not happen!!
    for (size_t i=0; i<keys.size(); i++) {
        cBox_->addItem(QString::fromStdString(descriptions.at(i)), QString::fromStdString(keys.at(i)));
        if(colors[i].a != 0)
            cBox_->setItemData(i,QColor(colors[i].r,colors[i].g,colors[i].b), Qt::ForegroundRole);
    }

    // set selected options
    if (!(property_->get() == "")) {
        int itemIndex = cBox_->findData(QString::fromStdString(property_->get()));
        if (itemIndex == -1) {
            LWARNINGC("OptionPropertyWidget", std::string("Data item ") + property_->get() + " not found");
        }
        else {
            cBox_->setCurrentIndex(itemIndex);
        }
    }

    cBox_->blockSignals(false);
}

void OptionPropertyWidget::setProperty(int index) {
    if (!disconnected_) {
        property_->set(cBox_->itemData(index).toString().toStdString());
        emit valueModifiedByUser();
    }
}

} // namespace
