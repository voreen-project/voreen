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

#include <string>

#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QWidget>

#include "voreen/qt/animation/animationtimelinewidget.h"
#include "tgt/camera.h"
#include "tgt/vector.h"
#include "voreen/core/properties/cameraproperty.h"

#include "voreen/qt/animation/animationpropertylistwidget.h"
#include "voreen/core/voreenapplication.h"
#include "voreen/core/voreenmodule.h"
#include "voreen/qt/widgets/property/qpropertywidget.h"
#include "voreen/core/properties/property.h"
#include "voreen/core/properties/templateproperty.h"

using tgt::vec2;
using tgt::vec3;
using tgt::vec4;

using tgt::ivec2;
using tgt::ivec3;
using tgt::ivec4;

using tgt::Camera;

using std::string;

namespace voreen {

//---------------------------------------------------------------------
// AnimationPropertyListLineWidget
//---------------------------------------------------------------------

AnimationPropertyListLineWidget::AnimationPropertyListLineWidget(AnimationPropertyListWidget* parent, PropertyTimeline* tl) :
    AnimationLineWidgetBase(parent, tl)
{  }

AnimationPropertyListLineWidget::~AnimationPropertyListLineWidget() {
    lineLayout_->removeWidget(deleteButton_);
    delete deleteButton_;
}

void AnimationPropertyListLineWidget::setFixedWidthSlot(int width) {
    setFixedWidth(width-30);

    int componentWidth = std::max(10, tgt::iround(static_cast<float>(width - 166) / 5.f));

    lineLayout_->setColumnMinimumWidth(2, 5*componentWidth);

    nameLabel_->setFixedWidth(4*componentWidth);

    //lineLayout_->setColumnMinimumWidth(2, componentWidth);
    //deleteButton_->setFixedwidth(...);
}

//-----------------------------------------------------------------------
// AnimationPropertyListWidget
//-----------------------------------------------------------------------

AnimationPropertyListWidget::AnimationPropertyListWidget(AnimationTimelineWidget* parent) : AnimationInfoWidgetBase(parent) {}

void AnimationPropertyListWidget::createHeader() {

    headerWidget_ = new QWidget(this);
    headerLayout_ = new QGridLayout(headerWidget_);
    headerLayout_->setMargin(0);
    headerLayout_->setSpacing(30);

    //QSpacerItem* spacer = new QSpacerItem(10, 0);

    QLabel* nameLabel = new QLabel(this);
    nameLabel->setText(QString("Animated Properties:"));
    nameLabel->setDisabled(true);

    //headerLayout_->addItem(spacer, 0, 0);
    //headerLayout_->addItem(spacer, 0, 1);
    headerLayout_->addWidget(nameLabel, 0, 2);
}

void AnimationPropertyListWidget::resizeHeader(int width) {

    if (!headerLayout_)
        return;

    int componentWidth = std::max(10, tgt::iround(static_cast<float>(width - 166) / 5.f));

    headerLayout_->setColumnMinimumWidth(0, 16);
    headerLayout_->setColumnMinimumWidth(1, 135);
    headerLayout_->setColumnMinimumWidth(2, 5 * componentWidth);
}

void AnimationPropertyListWidget::clear() {

    //remove all group boxes from layout, and delete the corresponding line widgets
    for (std::map<std::string, QGroupBox*>::iterator i = applicationModeGroups_.begin(); i != applicationModeGroups_.end(); ++i) {
        mainLayout_->removeWidget(i->second);
        std::map<std::string, std::vector<std::pair<int, PropertyTimeline*> > >::iterator groupTimelineIter = groupTimelines_.find(i->first);
        std::vector<std::pair<int, PropertyTimeline*> >::iterator timelineIter = groupTimelineIter->second.begin();
        for ( ; timelineIter != groupTimelineIter->second.end(); ++timelineIter) {
            std::map<PropertyTimeline*, AnimationPropertyListLineWidget*>::iterator lineWidgetIter = timelineWidgets_.find(timelineIter->second);
            i->second->layout()->removeWidget(lineWidgetIter->second);
            delete lineWidgetIter->second;
        }
        delete i->second;
    }

    // clear maps
    applicationModeGroups_.clear();
    groupTimelines_.clear();
    timelineWidgets_.clear();

    // remove header widget so that it is not destroyed
    mainLayout_->removeWidget(headerWidget_);

    // remove and delete everything still present on the layout
    removeWidgetsFromLayout(mainLayout_);

    // add the header again
    if (headerWidget_) {
        mainLayout_->addWidget(headerWidget_);
        mainLayout_->addSpacerItem(new QSpacerItem(0, 15));
    }
}

void AnimationPropertyListWidget::addPropertyTimeline(PropertyTimeline* prop) {

    if (!appConfig_) {
        LERRORC("voreen.qt.AnimationPropertyListWidget", "Could not add property timeline. No application mode configuration found.");
        return;
    }

    // create line widget for the property timeline
    AnimationPropertyListLineWidget* newLine = new AnimationPropertyListLineWidget(this, prop);
    if (!newLine->isValid()) {
        LERRORC("voreen.qt.AnimationPropertyListWidget", "Could not create PropertyWidget for " + prop->getProperty()->getGuiName());
        delete newLine;
        return;
    }

    connect(newLine, SIGNAL(deleteTimelineSignal(PropertyTimeline*)), this, SIGNAL(deleteTimelineSignal(PropertyTimeline*)));

    //if this is the first property timeline: create header, but don't add it to the layout, else: remove from layout
    if (!headerWidget_ && !headerLayout_)
        createHeader();
    else
        mainLayout_->removeWidget(headerWidget_);

    //remove all group boxes from layout
    for (std::map<std::string, QGroupBox*>::iterator i = applicationModeGroups_.begin(); i != applicationModeGroups_.end(); ++i)
        mainLayout_->removeWidget(i->second);

    // get group and priority for the property timeline
    std::pair<std::string, int> groupMembership = appConfig_->getPropertyGroupMembership(prop->getProperty());

    //if group does not exist already: create a new QGroupBox for it
    if (applicationModeGroups_.find(groupMembership.first) == applicationModeGroups_.end()) {
        QGroupBox* groupBox = new QGroupBox(QString::fromStdString(groupMembership.first), this);

        groupBox->setStyleSheet(" QGroupBox { border: 1px solid gray; border-radius: 9px; margin-top: 0.5em; } \
            QGroupBox::title { subcontrol-origin: margin; left: 10px; padding: 0 3px 0 3px; }");

        QVBoxLayout* boxLayout = new QVBoxLayout(/*this*/);
        boxLayout->setMargin(4);
        boxLayout->setSpacing(4);
        boxLayout->setAlignment(Qt::AlignTop);
        groupBox->setLayout(boxLayout);
        applicationModeGroups_.insert(std::make_pair(groupMembership.first, groupBox));
        groupTimelines_.insert(std::make_pair(groupMembership.first, std::vector<std::pair<int, PropertyTimeline*> >()));
    }

    //find the correct group and remove all its line widgets
    QGroupBox* box = applicationModeGroups_.find(groupMembership.first)->second;
    std::vector<std::pair<int, PropertyTimeline*> >& timelines = groupTimelines_.find(groupMembership.first)->second;
    for (std::vector<std::pair<int, PropertyTimeline*> >::iterator i = timelines.begin(); i != timelines.end(); ++i)
        if (timelineWidgets_.find(i->second) != timelineWidgets_.end())
            box->layout()->removeWidget(timelineWidgets_.find(i->second)->second);

    //add the new property to the group and sort the group properties
    timelines.push_back(std::make_pair(groupMembership.second, prop));
    std::sort(timelines.begin(), timelines.end(), AnimationInfoWidgetBase::sortTimelinesByPriority);
    timelineWidgets_.insert(std::make_pair(prop, newLine));
    connect(this, SIGNAL(setFixedWidthSignal(int)), newLine, SLOT(setFixedWidthSlot(int)));

    while (box->layout()->count() != 0) {
        QLayoutItem* item = box->layout()->takeAt(0);
        box->layout()->removeItem(item);
        delete item;
    }

    while(mainLayout_->count() != 0) {
        QLayoutItem* item = mainLayout_->takeAt(0);
        mainLayout_->removeItem(item);
        delete item;
    }

    //add the line widgets in the correct order to the box
    for (std::vector<std::pair<int, PropertyTimeline*> >::iterator i = timelines.begin(); i != timelines.end(); ++i) {
        box->layout()->addWidget(timelineWidgets_.find(i->second)->second);
        dynamic_cast<QBoxLayout*>(box->layout())->addSpacerItem(new QSpacerItem(0, 10));
    }

    // add the header again to the main layout
    mainLayout_->addWidget(headerWidget_);
    mainLayout_->addSpacerItem(new QSpacerItem(0, 15));

    //add the boxes in the correct order to the main layout
    std::vector<std::string> groupsInOrder = appConfig_->getPropertyGroups();
    for (std::vector<std::string>::iterator i = groupsInOrder.begin(); i != groupsInOrder.end(); ++i) {
        std::map<std::string, QGroupBox*>::iterator findIter = applicationModeGroups_.find(*i);
        if (findIter != applicationModeGroups_.end()) {
            mainLayout_->addWidget(findIter->second);
            mainLayout_->addSpacerItem(new QSpacerItem(0, 10));
        }
    }

    setFixedWidthSlot(width() - 20);

    mainLayout_->update();
}

void AnimationPropertyListWidget::updateActivationStatus() {

    std::map<PropertyTimeline*, AnimationPropertyListLineWidget*>::iterator it = timelineWidgets_.begin();
    for ( ; it != timelineWidgets_.end(); ++it) {
        it->second->updateActivationStatus();
    }
}

} // namespace voreen
