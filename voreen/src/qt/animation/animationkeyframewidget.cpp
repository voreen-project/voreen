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

#include <string>

#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QWidget>

#include "voreen/qt/animation/animationtimelinewidget.h"
#include "tgt/camera.h"
#include "tgt/vector.h"
#include "voreen/core/properties/cameraproperty.h"

#include "voreen/qt/animation/animationkeyframewidget.h"
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
// AnimationKeyframeLineWidget
//---------------------------------------------------------------------

AnimationKeyframeLineWidget::AnimationKeyframeLineWidget(AnimationKeyframeWidget* parent, PropertyTimeline* tl)
        : AnimationLineWidgetBase(parent, tl)
        , property_(0)
        , propertyWidget_(0)
        , currentKey_(0)
{
    //create 'dummy' widget and corresponding property for the current key
    property_ = dynamic_cast<Property*>(tl->getProperty()->create());
    propertyWidget_ = VoreenApplication::app()->createPropertyWidget(property_);

    if (!tl->getKeyValuesAbstract().empty())
        updateProperty(tl->getKeyValuesAbstract().begin()->second);

    QPropertyWidget* qWidget = dynamic_cast<QPropertyWidget*>(propertyWidget_);

    if (qWidget) {

       lineLayout_->addWidget(qWidget, 0, 3);

        TemplateProperty<Camera>* testCamera = dynamic_cast<TemplateProperty<Camera>*>(property_);
        TemplateProperty<TransFunc1DKeys*>* test1dTf = dynamic_cast<TemplateProperty<TransFunc1DKeys*>*>(property_);
        TemplateProperty<TransFunc2DPrimitives*>* test2dTf = dynamic_cast<TemplateProperty<TransFunc2DPrimitives*>*>(property_);
        TemplateProperty<tgt::IntBounds>* testIntBounds = dynamic_cast<TemplateProperty<tgt::IntBounds>*>(property_);
        if (testCamera) {
            testCamera->onChange(MemberFunctionCallback<PropertyWidget>(propertyWidget_, &PropertyWidget::updateFromProperty));
            testCamera->onChange(MemberFunctionCallback<AnimationKeyframeLineWidget>(this, &AnimationKeyframeLineWidget::updateKeyframe));
        }
        else if (test1dTf) {
            test1dTf->onChange(MemberFunctionCallback<PropertyWidget>(propertyWidget_, &PropertyWidget::updateFromProperty));
            test1dTf->onChange(MemberFunctionCallback<AnimationKeyframeLineWidget>(this, &AnimationKeyframeLineWidget::updateKeyframe));

            //if (TransFuncPropertyWidget* tfWidget = dynamic_cast<TransFuncPropertyWidget*>(propertyWidget_))
            //    connect(tfWidget, SIGNAL(zoomLevelChanged()), this, SLOT(tfWidgetZoomLevelChanged()));

            /*for (std::set<PropertyWidget*>::const_iterator i = tl->getProperty()->getPropertyWidgets().begin();
                    i != tl->getProperty()->getPropertyWidgets().end(); ++i) {

                if (QPropertyWidget* qWidget = dynamic_cast<QPropertyWidget*>(*i));
                    connect(qWidget, SIGNAL(zoomLevelChanged()), this, SLOT(adjustTfWidgetZoomLevel()));
            }*/
        }
        else if (test2dTf) {
            test2dTf->onChange(MemberFunctionCallback<PropertyWidget>(propertyWidget_, &PropertyWidget::updateFromProperty));
            test2dTf->onChange(MemberFunctionCallback<AnimationKeyframeLineWidget>(this, &AnimationKeyframeLineWidget::updateKeyframe));

            //if (TransFuncPropertyWidget* tfWidget = dynamic_cast<TransFuncPropertyWidget*>(propertyWidget_))
            //    connect(tfWidget, SIGNAL(zoomLevelChanged()), this, SLOT(tfWidgetZoomLevelChanged()));

            /*for (std::set<PropertyWidget*>::const_iterator i = tl->getProperty()->getPropertyWidgets().begin();
                    i != tl->getProperty()->getPropertyWidgets().end(); ++i) {

                if (QPropertyWidget* qWidget = dynamic_cast<QPropertyWidget*>(*i));
                    connect(qWidget, SIGNAL(zoomLevelChanged()), this, SLOT(adjustTfWidgetZoomLevel()));
            }*/
        }
        else if (testIntBounds){
            testIntBounds->onChange(MemberFunctionCallback<PropertyWidget>(propertyWidget_, &PropertyWidget::updateFromProperty));
            testIntBounds->onChange(MemberFunctionCallback<AnimationKeyframeLineWidget>(this, &AnimationKeyframeLineWidget::updateKeyframe));
        }

        else {
            connect(qWidget, SIGNAL(widgetChanged()), this, SLOT(updateKeyframe()));
            connect(qWidget, SIGNAL(valueModifiedByUser()), this, SLOT(updateKeyframe()));
        }

    }
    else {
        delete property_; property_ = 0;
        delete propertyWidget_; propertyWidget_ = 0;
    }

    updateActivationStatus();
    lineLayout_->update();
}

AnimationKeyframeLineWidget::~AnimationKeyframeLineWidget() {

    QPropertyWidget* qWidget = dynamic_cast<QPropertyWidget*>(propertyWidget_);
    lineLayout_->removeWidget(qWidget);
    delete propertyWidget_;

    delete property_;
}

bool AnimationKeyframeLineWidget::isValid() const {
    return (property_ && propertyWidget_);
}

void AnimationKeyframeLineWidget::updateActivationStatus() {
    bool active = propertyTimeline_->getActiveOnRendering();
    propertyWidget_->setEnabled(active);
    AnimationLineWidgetBase::updateActivationStatus();
}

void AnimationKeyframeLineWidget::updateProperty(PropertyKeyValueBase* key) {

    //remember the left and right keys
    currentKey_ = key;

    //update properties and interpolation
    if (informationWidget_)
        informationWidget_->updatePropertyFromKey(property_, propertyWidget_, key, propertyTimeline_);
}

void AnimationKeyframeLineWidget::updateKeyframe() {
    if (informationWidget_)
        informationWidget_->updateKeyFromWidget(propertyTimeline_, currentKey_, property_);
}

/*void AnimationKeyframeLineWidget::tfWidgetZoomLevelChanged() {
    if (property_->getMetaDataContainer().hasMetaData("TransfuncPropertyWidgetPainterZoom")) {
        propertyTimeline_->getProperty()->getMetaDataContainer().addMetaData("TransfuncPropertyWidgetPainterZoom",
                                      property_->getMetaDataContainer().getMetaData("TransfuncPropertyWidgetPainterZoom")->clone());
        propertyTimeline_->getProperty()->updateWidgets();
        std::cout << "Zoom" << std::endl;
    }
}

void AnimationKeyframeLineWidget::adjustTfWidgetZoomLevel() {
    std::cout << "Network Property Zoom Level changed" << std::endl;
}*/

//-----------------------------------------------------------------------
// AnimationKeyframeWidget
//-----------------------------------------------------------------------

AnimationKeyframeWidget::AnimationKeyframeWidget(AnimationTimelineWidget* parent) : AnimationInfoWidgetBase(parent) {}

void AnimationKeyframeWidget::createHeader() {

    headerWidget_ = new QWidget(this);
    headerLayout_ = new QGridLayout(headerWidget_);
    headerLayout_->setMargin(0);
    headerLayout_->setSpacing(30);

    //QSpacerItem* spacer = new QSpacerItem(46, 0);

    QLabel* nameLabel = new QLabel(this);
    nameLabel->setText(QString("Property Name"));
    nameLabel->setDisabled(true);

    QLabel* valueLabel = new QLabel(this);
    valueLabel->setText(QString("Key Value"));
    valueLabel->setDisabled(true);

    //headerLayout_->addItem(spacer, 0, 0);
    headerLayout_->addWidget(nameLabel, 0, 2);
    headerLayout_->addWidget(valueLabel, 0, 3);
}

void AnimationKeyframeWidget::resizeHeader(int width) {

    if (!headerLayout_)
        return;

    int componentWidth = std::max(10, tgt::iround(static_cast<float>(width - 166) / 5.f));

    headerLayout_->setColumnMinimumWidth(0, 16);
    headerLayout_->setColumnMinimumWidth(1, 135);
    headerLayout_->setColumnMinimumWidth(2, tgt::iround(1.42 * static_cast<float>(componentWidth)));
    //headerLayout_->setColumnMinimumWidth(3, 3 * componentWidth);
}

void AnimationKeyframeWidget::clear() {

    //remove all group boxes from layout, and delete the corresponding line widgets
    for (std::map<std::string, QGroupBox*>::iterator i = applicationModeGroups_.begin(); i != applicationModeGroups_.end(); ++i) {
        mainLayout_->removeWidget(i->second);
        std::map<std::string, std::vector<std::pair<int, PropertyTimeline*> > >::iterator groupTimelineIter = groupTimelines_.find(i->first);
        std::vector<std::pair<int, PropertyTimeline*> >::iterator timelineIter = groupTimelineIter->second.begin();
        for ( ; timelineIter != groupTimelineIter->second.end(); ++timelineIter) {
            std::map<PropertyTimeline*, AnimationKeyframeLineWidget*>::iterator lineWidgetIter = timelineWidgets_.find(timelineIter->second);
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

void AnimationKeyframeWidget::addPropertyTimeline(PropertyTimeline* prop) {

    if (!appConfig_) {
        LERRORC("voreen.qt.AnimationKeyframeWidget", "Could not add property timeline. No application mode configuration found.");
        return;
    }

    // create line widget for the property timeline
    AnimationKeyframeLineWidget* newLine = new AnimationKeyframeLineWidget(this, prop);
    if (!newLine->isValid()) {
        LERRORC("voreen.qt.AnimationKeyframeWidget", "Could not create PropertyWidget for " + prop->getProperty()->getGuiName());
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

void AnimationKeyframeWidget::updateProperties(const std::vector<std::pair<PropertyKeyValueBase*, PropertyTimeline*> >& propertyValues) {
    //iterate over property timelines
    for (size_t i = 0; i < propertyValues.size(); ++i) {
        //for this property timeline: geht the corresponding entry in the map and update property and widget
        std::map<PropertyTimeline*, AnimationKeyframeLineWidget*>::iterator propIt = timelineWidgets_.find(propertyValues.at(i).second);

        if (propIt != timelineWidgets_.end()) {
            propIt->second->updateProperty(propertyValues.at(i).first);
        }
        else
            LERRORC("voreen.qt.AnimationKeyframeWidget", "PropertyWidget not found for " + propertyValues.at(i).second->getProperty()->getGuiName());
    }
}

void AnimationKeyframeWidget::updateActivationStatus() {

    std::map<PropertyTimeline*, AnimationKeyframeLineWidget*>::iterator it = timelineWidgets_.begin();
    for ( ; it != timelineWidgets_.end(); ++it) {
        it->second->updateActivationStatus();
    }
}

} // namespace voreen
