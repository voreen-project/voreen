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

#include "voreen/qt/widgets/property/propertyownerwidget.h"

#include "voreen/qt/widgets/customlabel.h"

#include "voreen/core/voreenapplication.h"
#include "voreen/core/properties/propertyowner.h"
#include "voreen/qt/widgets/expandableheaderbutton.h"
#include "voreen/qt/widgets/property/qpropertywidget.h"
#include "voreen/qt/widgets/property/grouppropertywidget.h"
#include "voreen/qt/widgets/property/lightpropertywidget.h"
#include "voreen/qt/widgets/property/propertyvectorwidget.h"

#include <QPushButton>
#include <QGridLayout>

namespace voreen {

PropertyOwnerWidget::PropertyOwnerWidget(PropertyOwner* propertyOwner, QWidget* parent, std::string title,
                                                     bool expanded, bool userExpandable, bool addResetButton)
    : QWidget(parent)
    , propertyOwner_(propertyOwner)
    , mainLayout_(0)
    , header_(0)
    , propertyWidget_(0)
    , hasBeenInitialized_(false)
    , addResetButton_(addResetButton)
    , currentLOD_(Property::LOD_DEFAULT)
{
    tgtAssert(propertyOwner,"null pointer passed");

    setObjectName("PropertyOwnerTitleWidget"); //?

    mainLayout_ = new QVBoxLayout(this);
    mainLayout_->setContentsMargins(0, 0, 0, 0);
    mainLayout_->setSpacing(0);
    mainLayout_->setMargin(0);

    //setUpdatesEnabled(false); //?

    header_ = new ExpandableHeaderButton(title != "" ? title.c_str() : propertyOwner_->getGuiName().c_str(), this,
            expanded, userExpandable);
    connect(header_, SIGNAL(toggled(bool)), this, SLOT(updateState()));
    connect(header_, SIGNAL(toggled(bool)), this, SLOT(instantiateWidgets()));
    mainLayout_->addWidget(header_);

    dynamic_cast<const Observable<PropertyOwnerObserver>* >(propertyOwner_)->addObserver(dynamic_cast<PropertyOwnerObserver*>(this));
}

PropertyOwnerWidget::~PropertyOwnerWidget() {
    std::vector<QPropertyWidget*>::iterator it = widgets_.begin();
    for ( ; it != widgets_.end(); ++it) {
        layout()->removeWidget(*it);
        delete (*it);
    }
}

//----------------------------------------------------------------------------------------------------
//                   getter and setter
//----------------------------------------------------------------------------------------------------
PropertyOwner* PropertyOwnerWidget::getPropertyOwner() const {
    return propertyOwner_;
}

bool PropertyOwnerWidget::isExpanded() const {
    return header_->isExpanded();
}

void PropertyOwnerWidget::setExpanded(bool expanded) {
    header_->setExpanded(expanded);
    updateState();
}

void PropertyOwnerWidget::toggleExpansionState() {
    header_->setExpanded(!header_->isExpanded());
    updateState();
}

bool PropertyOwnerWidget::isUserExpandable() const {
    return header_->userExpandable();
}

void PropertyOwnerWidget::setUserExpandable(bool expandable) {
    header_->setUserExpandable(expandable);
    updateState();
}

void PropertyOwnerWidget::updateHeaderTitle() {
    header_->updateNameLabel(propertyOwner_->getGuiName());
}

void PropertyOwnerWidget::resetAllProperties() {
    propertyOwner_->resetAllProperties();
}

std::vector<Property*>* PropertyOwnerWidget::createPropertyList() {
    return new std::vector<Property*>(propertyOwner_->getProperties());
}

//----------------------------------------------------------------------------------------------------
//                   widget creation
//----------------------------------------------------------------------------------------------------
void PropertyOwnerWidget::instantiateWidgets() {
    setUpdatesEnabled(false);
    if (!VoreenApplication::app()) {
        LERRORC("voreen.qt.PropertyOwnerWidget", "VoreenApplication not instantiated");
        return;
    }
    if (!hasBeenInitialized_) {
        propertyWidget_ = new QWidget;
        QGridLayout* gridLayout = new QGridLayout(propertyWidget_);
        gridLayout->setContentsMargins(0, 4, 3, 2); //3,4,0,2
        gridLayout->setSpacing(2);//2
        gridLayout->setColumnStretch(0, 1);
        gridLayout->setColumnStretch(1, 2);
        gridLayout->setEnabled(false);
        std::vector<Property*>* propertyList = createPropertyList();

        // create widget for every property and put them into a vertical layout
        int rows = 0;

        for (std::vector<Property*>::iterator iter = propertyList->begin(); iter != propertyList->end(); ++iter) {
            Property* prop = *iter;
            PropertyWidget* propWidget = VoreenApplication::app()->createPropertyWidget(prop);
            if (propWidget)
                prop->addWidget(propWidget);

            QPropertyWidget* w = dynamic_cast<QPropertyWidget*>(propWidget);
            if (w != 0) {
                widgets_.push_back(w);
                connect(w,SIGNAL(checkVisibility()),this,SLOT(updatePropertyVisibility()));
                connect(w,SIGNAL(checkGroupVisibility()),this,SLOT(updateGroupVisibility()));
                // we are dealing with a propertygroup
                if (prop->getGroupID() != "") {
                    std::map<std::string, GroupPropertyWidget*>::iterator it = propertyGroupsMap_.find(prop->getGroupID());
                    // the group does not exist
                    if (it == propertyGroupsMap_.end()) {
                        std::string guiName = prop->getOwner()->getPropertyGroupGuiName(prop->getGroupID());
                        propertyGroupsMap_[prop->getGroupID()] = new GroupPropertyWidget(prop, false, guiName, this);
                    }
                    propertyGroupsMap_[prop->getGroupID()]->addWidget(w, w->getOrCreateNameLabel(), QString::fromStdString(w->getPropertyGuiName()));
                    gridLayout->addWidget(propertyGroupsMap_[prop->getGroupID()], rows, 0, 1, 2);
                }
                else {
                    CustomLabel* nameLabel = w->getOrCreateNameLabel();

                    if(dynamic_cast<LightPropertyWidget*>(w)) {     // HACK: this prevents a cut off gui element e.g. seen in the clipping plane widget
                        gridLayout->setRowStretch(rows, 1);
                    }

                    if (nameLabel) {
                        nameLabel->setContentsMargins(4,0,0,0); //HACK: otherwise, the group widgets and TF are not aligned
                        gridLayout->addWidget(nameLabel, rows, 0, 1, 1);
                        if(dynamic_cast<PropertyVectorWidget*>(w))
                            ++rows;
                        gridLayout->addWidget(w, rows, 1, 1, 1);
                    }
                    else {
                        gridLayout->addWidget(w, rows, 0, 1, 2);
                    }
                }
            }
            ++rows;
        }
        //adding a button to reset all properties if needed.
        if(addResetButton_){
            QPushButton* pb = new QPushButton(QString("Reset All Properties"),this);
            connect(pb, SIGNAL(clicked()), this, SLOT(resetAllProperties()));
            gridLayout->addWidget(pb , rows, 0, 1, 2);
        }

        gridLayout->setEnabled(true);

        mainLayout_->addWidget(propertyWidget_);

        //check visibility of property group
        std::map<std::string, GroupPropertyWidget*>::iterator it = propertyGroupsMap_.begin();
        for(; it != propertyGroupsMap_.end(); it++)
            it->second->setVisible(it->second->isAnyPropertyVisible(currentLOD_));

        setUpdatesEnabled(true);
        //updateState();
        hasBeenInitialized_ = true;
        delete propertyList;
    }
    setUpdatesEnabled(true);
    //forced update
    Property::LevelOfDetail tmp = currentLOD_;
    currentLOD_ = static_cast<Property::LevelOfDetail>(currentLOD_ + 1);
    setLevelOfDetail(tmp);
}

void PropertyOwnerWidget::updateState() {
    if(propertyWidget_!=0) {
        propertyWidget_->setVisible(header_->isExpanded());
    }
    updateGeometry();  // prevent flicker when hiding property widgets
}

void PropertyOwnerWidget::showEvent(QShowEvent* event) {
    if(!hasBeenInitialized_) instantiateWidgets();
    QWidget::showEvent(event);
}

//----------------------------------------------------------------------------------------------------
//                   propertyowner observer function
//----------------------------------------------------------------------------------------------------
void PropertyOwnerWidget::propertiesChanged(const PropertyOwner*) {

    if (!VoreenApplication::app()) {
        LERRORC("voreen.qt.PropertyOwnerWidget", "VoreenApplication not instantiated");
        return;
    }

    if(!propertyWidget_) {
        LERRORC("voreen.qt.PropertyOwnerWidget", "Property widget not instantiated");
        return;
    }

    bool newWidgetAdded = false;
    QGridLayout* gridLayout = dynamic_cast<QGridLayout*>(propertyWidget_->layout());
    std::vector<Property*>* properties = createPropertyList();
    for (unsigned int i=0; i<properties->size(); i++) {
        Property* prop = properties->at(i);
        const std::set<PropertyWidget*> propWidgets = prop->getPropertyWidgets();
        if (propWidgets.empty()) {
            PropertyWidget* propWidget = VoreenApplication::app()->createPropertyWidget(prop);
            if (propWidget)
                prop->addWidget(propWidget);
            QPropertyWidget* w = dynamic_cast<QPropertyWidget*>(propWidget);
            if (w) {
                widgets_.push_back(w);
                connect(w,SIGNAL(checkVisibility()),this,SLOT(updatePropertyVisibility()));
                connect(w,SIGNAL(checkGroupVisibility()),this,SLOT(updateGroupVisibility()));
                CustomLabel* nameLabel = w->getOrCreateNameLabel();
                int curRow = gridLayout->rowCount();
                if(nameLabel) {
                    gridLayout->addWidget(nameLabel, curRow, 0, 1, 1);
                    gridLayout->addWidget(w, curRow, 1, 1, 1);
                }
                else {
                    gridLayout->addWidget(w, curRow, 0, 1, 2);
                }
                newWidgetAdded = true;
            }
        }
    }
    delete properties;
    //forced update
    if(newWidgetAdded) {
        Property::LevelOfDetail tmp = currentLOD_;
        currentLOD_ = static_cast<Property::LevelOfDetail>(currentLOD_ + 1);
        setLevelOfDetail(tmp);
    }
}

void PropertyOwnerWidget::updatePropertyVisibility() {
    QObject* obj = QObject::sender();
    if(QPropertyWidget* widget = dynamic_cast<QPropertyWidget*>(obj)) {
        if(Property* prop = widget->getProperty()) {
            widget->setVisible(prop->isVisibleFlagSet() && (prop->getLevelOfDetail() >= currentLOD_));
        }
    }
}

void PropertyOwnerWidget::updateGroupVisibility() {
    QObject* obj = QObject::sender();
    if(QPropertyWidget* widget = dynamic_cast<QPropertyWidget*>(obj)) {
        if(Property* prop = widget->getProperty()) {
            std::map<std::string, GroupPropertyWidget*>::iterator it = propertyGroupsMap_.find(prop->getGroupID());
            if(it != propertyGroupsMap_.end()) //groups do not exist at the beginning
                it->second->setVisible(it->second->isAnyPropertyVisible(currentLOD_));
        }
    }
}

void PropertyOwnerWidget::setLevelOfDetail(Property::LevelOfDetail lod) {
    //return, if lod havent chnaged
    if(currentLOD_ == lod) return;
    //update lod
    currentLOD_ = lod;
    //return if not initialized. initWidgets must call setLevelOfDetail again
    if(!hasBeenInitialized_) return;
    //get all properties
    std::vector<Property*> propertyList(propertyOwner_->getProperties());

    bool atLeastOnVisibleWidget_ = false;

    // update property widgets visibility and LOD controls
    for (size_t i = 0; i < propertyList.size(); ++i) {
        std::set<PropertyWidget*> set = propertyList[i]->getPropertyWidgets(); // This set usually contains only one property
        for (std::set<PropertyWidget*>::iterator iter = set.begin(); iter != set.end(); ++iter) {
            QPropertyWidget* wdt = dynamic_cast<QPropertyWidget*> (*iter);
            for (std::vector<QPropertyWidget*>::iterator innerIter = widgets_.begin(); innerIter != widgets_.end(); ++innerIter) {
                if (wdt == (*innerIter)) {
                    bool visible = propertyList[i]->isVisibleFlagSet() && (propertyList[i]->getLevelOfDetail() >= currentLOD_);
                    wdt->setVisible(visible);
                    //wdt->showNameLabel(visible); // TODO: needed?
                    if(visible)
                        atLeastOnVisibleWidget_ = true;
                }
            }
        }
    }

    std::map<std::string, GroupPropertyWidget*>::iterator it;
    for(it = propertyGroupsMap_.begin(); it != propertyGroupsMap_.end(); it++) {
        GroupPropertyWidget* gpw = it->second;
        if(gpw) {
            if(gpw->isAnyPropertyVisible(lod))
                gpw->setVisible(true);
            else
                gpw->setVisible(false);
        }
    }

    setVisible(atLeastOnVisibleWidget_);
}

} // namespace
