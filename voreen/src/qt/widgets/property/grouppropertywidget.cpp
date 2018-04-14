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

#include "voreen/qt/widgets/property/grouppropertywidget.h"

#include <QTabWidget>
#include <QPushButton>

namespace voreen {

    GroupPropertyWidget::GroupPropertyWidget(Property* prop, bool tabbed, std::string guiName, QWidget* parent)
    : QPropertyWidget(prop, parent)
    , prop_(prop)
    , tabbed_(tabbed)
    , groupName_(guiName)
{
    if(tabbed) {
        tabWidget_ = new QTabWidget(this);
        layout_->addWidget(tabWidget_);
        groupBox_ = 0;
    }
    else {
        tabWidget_ = 0;
        groupBox_ = new QGroupBox(QString::fromStdString(guiName), this);
        gridLayout_ = new QGridLayout();
        gridLayout_->setSpacing(2);
        gridLayout_->setColumnStretch(0, 1);
        gridLayout_->setColumnStretch(1, 2);
        gridLayout_->setContentsMargins(5, 2, 2, 5);
        groupBox_->setLayout(gridLayout_);
        layout_->addWidget(groupBox_);
        /*QIcon icon = QIcon(":/qt/icons/expand-minus.png");
        hideWidgetButton_ = new QPushButton(icon, "", this);
        hideWidgetButton_->setGeometry(QRect(10, 10, 10, 10));
        hideWidgetButton_->setFixedSize(10,10);
        hideWidgetButton_->setFlat(true);
        hideWidgetButton_->setCheckable(true);
        hideWidgetButton_->setStyleSheet("QToolButton { border: none; padding: 1px; }");
        connect(hideWidgetButton_, SIGNAL(toggled(bool)), this, SLOT(hideGroup(bool)));
        gridLayout_->addWidget(hideWidgetButton_, 0, 1, Qt::AlignRight);
        // disable hide button for now
        hideWidgetButton_->setVisible(false);*/
    }
}

void GroupPropertyWidget::addWidget(QPropertyWidget* widget, const QString labelName) {
    if(tabbed_) {
        tabWidget_->addTab(widget, labelName);
    }
    else {
        gridLayout_->addWidget(widget);
    }
    propertyWidgets_.push_back(widget);
}

void GroupPropertyWidget::addWidget(QPropertyWidget* widget, QWidget* label, const QString labelName) {
    if(tabbed_) {
        QWidget* container = new QWidget(this);
        QGridLayout* layout = new QGridLayout(container);
        if(label) { //handle widgets without label (e.g. volumeurlproperty)
            layout->addWidget(label, 0, 0);
            layout->addWidget(widget, 0, 1);
        } else {
            layout->addWidget(widget, 0, 0, 1, 2);
        }
        tabWidget_->addTab(container, labelName);
    }
    else {
        if(label) { //handle widgets without label (e.g. volumeurlproperty)
            gridLayout_->addWidget(label, gridLayout_->rowCount(), 0);
            gridLayout_->addWidget(widget, gridLayout_->rowCount()-1, 1);
        } else {
           gridLayout_->addWidget(widget, gridLayout_->rowCount(), 0, 1, 2);
        }
    }
    propertyWidgets_.push_back(widget);
}
void GroupPropertyWidget::hideGroup(bool toggled){
    std::vector<QPropertyWidget*>::iterator it = propertyWidgets_.begin();
    if(toggled) {
        hideWidgetButton_ ->setIcon(QIcon(":/qt/icons/expand-plus.png"));

        while(it != propertyWidgets_.end()) {
            (*it)->hide();
            ++it;
        }
    }
    else {
        hideWidgetButton_ ->setIcon(QIcon(":/qt/icons/expand-minus.png"));

        while(it != propertyWidgets_.end()) {
            (*it)->show();
            ++it;
        }
    }
}

void GroupPropertyWidget::updateFromPropertySlot() {
    std::vector<QPropertyWidget*>::iterator it = propertyWidgets_.begin();
    while(it != propertyWidgets_.end()) {
        (*it)->updateFromProperty();
        ++it;
    }
}

bool GroupPropertyWidget::isAnyPropertyVisible(Property::LevelOfDetail lod) {
    for(size_t i=0; i<propertyWidgets_.size(); i++) {
        QPropertyWidget* propWidget = propertyWidgets_[i];
        if(propWidget) {
            Property* prop = propWidget->getProperty();
            if(prop && prop->isVisibleFlagSet() && (prop->getLevelOfDetail() >= lod))
                return true;
        }
    }
    return false;
}

Property* GroupPropertyWidget::getProperty() {
    return prop_;
}

std::string GroupPropertyWidget::getGroupName() const {
    return groupName_;
}

CustomLabel* GroupPropertyWidget::getOrCreateNameLabel() const {
    return 0;
}

} // namespace
