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

#include "voreen/qt/widgets/property/appmodepropertygroupwidget.h"

#include "voreen/qt/widgets/customlabel.h"

#include "voreen/core/voreenapplication.h"
#include "voreen/core/voreenmodule.h"
#include "voreen/core/interaction/interactionhandler.h"
#include "voreen/core/processors/processor.h"
#include "voreen/core/ports/port.h"

#include "voreen/qt/widgets/expandableheaderbutton.h"
#include "voreen/qt/widgets/property/lightpropertywidget.h"
#include "voreen/qt/widgets/property/propertyvectorwidget.h"
#include "voreen/qt/widgets/property/qpropertywidget.h"
#include "voreen/qt/widgets/property/grouppropertywidget.h"
#include "voreen/core/ports/renderport.h" //for size properties

#include <QVBoxLayout>
#include <QGridLayout>

namespace voreen {

AppModePropertyGroupWidget::AppModePropertyGroupWidget(const std::string& groupName, QWidget* parent, bool expanded)
    : QWidget(parent)
    , propertyFrame_(0)
    , expanded_(expanded)
{
    setObjectName("AppModePropertyGroupWidget");

    mainLayout_ = new QVBoxLayout();
    mainLayout_->setContentsMargins(1, 1, 1, 1);
    mainLayout_->setSpacing(0);
    //mainLayout_->setMargin(0);
    setLayout(mainLayout_);

    header_ = new ExpandableHeaderButton(groupName.c_str(), this, expanded_);
    connect(header_, SIGNAL(toggled(bool)), this, SLOT(updateState()));
    header_->setFrameShape(QFrame::StyledPanel);
    header_->setBackgroundRole(QPalette::ButtonText);
    QFont headerFont = header_->font();
    headerFont.setBold(true);
    //headerFont.setItalic(true);
    //headerFont.setPointSize(headerFont.pointSize()+1);
    header_->setFontStyle(headerFont);
    header_->setStyleSheet("QFrame { background-color: #d7d7d7; }");
    mainLayout_->addWidget(header_);
    mainLayout_->addSpacing(-1);

    propertyFrame_ = new QFrame();
    propertyLayout_ = new QGridLayout(propertyFrame_);
    propertyLayout_->setContentsMargins(2, 5, 5, 5);
    propertyLayout_->setSpacing(3);
    propertyLayout_->setColumnStretch(0, 1);
    propertyLayout_->setColumnStretch(1, 2);
    propertyFrame_->setFrameShape(QFrame::StyledPanel);
    propertyFrame_->setStyleSheet("QFrame { background-color: white; }");
    mainLayout_->addWidget(propertyFrame_);

}

void AppModePropertyGroupWidget::updateState() {
    tgtAssert(propertyFrame_, "no property frame");
    propertyFrame_->setVisible(header_->isExpanded());
    updateGeometry();  // prevent flicker when hiding property widgets
}

bool AppModePropertyGroupWidget::isExpanded() const {
    return header_->isExpanded();
}

void AppModePropertyGroupWidget::setExpanded(bool expanded) {
    header_->setExpanded(expanded);
    updateState();
}

void AppModePropertyGroupWidget::toggleExpansionState() {
    header_->setExpanded(!header_->isExpanded());
    updateState();
}

void AppModePropertyGroupWidget::instantiateWidgets() {
    tgtAssert(propertyLayout_, "no property layout");

    if (!VoreenApplication::app()) {
        LERRORC("voreen.qt.ProcessorPropertiesWidget", "VoreenApplication not instantiated");
        return;
    }

    // create widget for every property and put them into a vertical layout
    int rows = 0;

    // create non-lazy widgets (such as widgets with editor windows)
    for (std::vector<Property*>::iterator iter = properties_.begin(); iter != properties_.end(); ++iter) {
        Property* prop = *iter;
        PropertyWidget* propWidget = VoreenApplication::app()->createPropertyWidget(prop);
        if (propWidget)
            prop->addWidget(propWidget);

        QPropertyWidget* w = dynamic_cast<QPropertyWidget*>(propWidget);
        if (w) {
            connect(w, SIGNAL(checkVisibility()), this, SLOT(updatePropertyVisibility()));
            widgets_.push_back(w);
            CustomLabel* nameLabel = w->getOrCreateNameLabel();
            if(nameLabel) {
                nameLabel->setContentsMargins(3,0,0,0);
                propertyLayout_->addWidget(w, rows, 1, 1, 1);
                propertyLayout_->addWidget(nameLabel, rows, 0, 1, 1);
            }
            else {
                propertyLayout_->addWidget(w, rows, 0, 1, 2);
            }

            w->setVisible(prop->isVisibleFlagSet());
        }
        ++rows;
    }

    updateState();
}

void AppModePropertyGroupWidget::showEvent(QShowEvent* event) {
    //instantiateWidgets();
    QWidget::showEvent(event);
}

void AppModePropertyGroupWidget::setProperties(const std::vector<Property*> properties) {
    properties_ = properties;
}

void AppModePropertyGroupWidget::addProperty(Property* property) {
    tgtAssert(property, "null pointer passed");
    properties_.push_back(property);
}

void AppModePropertyGroupWidget::updatePropertyVisibility() {
    QObject* obj = QObject::sender();
    if (QPropertyWidget* widget = dynamic_cast<QPropertyWidget*>(obj)) {
        if (Property* prop = widget->getProperty()) {
            widget->setVisible(prop->isVisibleFlagSet());
        }
    }
}

} // namespace
