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

#include <string>

#include "tgt/camera.h"
#include "tgt/vector.h"
#include "voreen/core/properties/transfunc/1d/1dkeys/transfunc1dkeysproperty.h"
#include "voreen/core/properties/transfunc/2d/2dprimitives/transfunc2dprimitivesproperty.h"
#include "voreen/core/properties/cameraproperty.h"

#include "voreen/qt/animation/animationintervalwidget.h"
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
// AnimationIntervalLineWidget
//---------------------------------------------------------------------

AnimationIntervalLineWidget::AnimationIntervalLineWidget(AnimationIntervalWidget* parent, PropertyTimeline* tl)
        : AnimationLineWidgetBase(parent,tl)
        , lProperty_(0)
        , rProperty_(0)
        , lPropertyWidget_(0)
        , rPropertyWidget_(0)
        , interpolationMenu_(0)
        , leftKey_(0)
        , rightKey_(0)
{
    interpolationFunction_ = new QPushButton("< unknown >", this);

    //create two widgets and corresponding properties (for left and right key)
    lProperty_ = dynamic_cast<Property*>(tl->getProperty()->create());
    lPropertyWidget_ = VoreenApplication::app()->createPropertyWidget(lProperty_);

    rProperty_ = dynamic_cast<Property*>(tl->getProperty()->create());
    rPropertyWidget_ = VoreenApplication::app()->createPropertyWidget(rProperty_);

    if (!tl->getKeyValuesAbstract().empty())
        updateProperties(tl->getKeyValuesAbstract().begin()->second, tl->getKeyValuesAbstract().begin()->second);

    QPropertyWidget* lW = dynamic_cast<QPropertyWidget*>(lPropertyWidget_);
    QPropertyWidget* rW = dynamic_cast<QPropertyWidget*>(rPropertyWidget_);

    if (lW && rW) {

        lineLayout_->addWidget(lW, 0, 3);
        lineLayout_->addWidget(interpolationFunction_, 0, 4);
        lineLayout_->addWidget(rW, 0, 5);

        TemplateProperty<Camera>* lTestCamera = dynamic_cast<TemplateProperty<Camera>*>(lProperty_);
        TemplateProperty<TransFunc1DKeys*>* lTest1dTf = dynamic_cast<TemplateProperty<TransFunc1DKeys*>*>(lProperty_);
        TemplateProperty<TransFunc2DPrimitives*>* lTest2dTf = dynamic_cast<TemplateProperty<TransFunc2DPrimitives*>*>(lProperty_);
        TemplateProperty<Camera>* rTestCamera = dynamic_cast<TemplateProperty<Camera>*>(rProperty_);
        TemplateProperty<TransFunc1DKeys*>* rTest1dTf = dynamic_cast<TemplateProperty<TransFunc1DKeys*>*>(rProperty_);
        TemplateProperty<TransFunc2DPrimitives*>* rTest2dTf = dynamic_cast<TemplateProperty<TransFunc2DPrimitives*>*>(rProperty_);
        if (lTestCamera && rTestCamera) {
            lTestCamera->onChange(MemberFunctionCallback<PropertyWidget>(lPropertyWidget_, &PropertyWidget::updateFromProperty));
            lTestCamera->onChange(MemberFunctionCallback<AnimationIntervalLineWidget>(this, &AnimationIntervalLineWidget::updateLeftKeyframe));
            rTestCamera->onChange(MemberFunctionCallback<PropertyWidget>(rPropertyWidget_, &PropertyWidget::updateFromProperty));
            rTestCamera->onChange(MemberFunctionCallback<AnimationIntervalLineWidget>(this, &AnimationIntervalLineWidget::updateRightKeyframe));
        }
        else if (lTest1dTf && rTest1dTf) {
            lTest1dTf->onChange(MemberFunctionCallback<PropertyWidget>(lPropertyWidget_, &PropertyWidget::updateFromProperty));
            lTest1dTf->onChange(MemberFunctionCallback<AnimationIntervalLineWidget>(this, &AnimationIntervalLineWidget::updateLeftKeyframe));
            rTest1dTf->onChange(MemberFunctionCallback<PropertyWidget>(rPropertyWidget_, &PropertyWidget::updateFromProperty));
            rTest1dTf->onChange(MemberFunctionCallback<AnimationIntervalLineWidget>(this, &AnimationIntervalLineWidget::updateRightKeyframe));
        }
        else if (lTest2dTf && rTest2dTf) {
            lTest2dTf->onChange(MemberFunctionCallback<PropertyWidget>(lPropertyWidget_, &PropertyWidget::updateFromProperty));
            lTest2dTf->onChange(MemberFunctionCallback<AnimationIntervalLineWidget>(this, &AnimationIntervalLineWidget::updateLeftKeyframe));
            rTest2dTf->onChange(MemberFunctionCallback<PropertyWidget>(rPropertyWidget_, &PropertyWidget::updateFromProperty));
            rTest2dTf->onChange(MemberFunctionCallback<AnimationIntervalLineWidget>(this, &AnimationIntervalLineWidget::updateRightKeyframe));
        }
        else {
            connect(lW, SIGNAL(widgetChanged()), this, SLOT(updateLeftKeyframe()));
            connect(lW, SIGNAL(valueModifiedByUser()), this, SLOT(updateLeftKeyframe()));
            connect(rW, SIGNAL(widgetChanged()), this, SLOT(updateRightKeyframe()));
            connect(rW, SIGNAL(valueModifiedByUser()), this, SLOT(updateRightKeyframe()));
        }

    }
    else {
        delete lProperty_; lProperty_ = 0;
        delete rProperty_; rProperty_ = 0;
        delete lPropertyWidget_; lPropertyWidget_ = 0;
        delete rPropertyWidget_; rPropertyWidget_ = 0;
        delete interpolationFunction_; interpolationFunction_ = 0;
    }

    updateActivationStatus();

    interpolationMenu_ = new QMenu(this);
    populateInterpolationMenu();
    interpolationFunction_->setMenu(interpolationMenu_);

    connect(interpolationMenu_, SIGNAL(triggered(QAction*)), this, SLOT(changeInterpolation(QAction*)));

    lineLayout_->update();
}

AnimationIntervalLineWidget::~AnimationIntervalLineWidget() {

    QPropertyWidget* lW = dynamic_cast<QPropertyWidget*>(lPropertyWidget_);
    QPropertyWidget* rW = dynamic_cast<QPropertyWidget*>(rPropertyWidget_);
    lineLayout_->removeWidget(lW);
    lineLayout_->removeWidget(rW);
    delete lPropertyWidget_;
    delete rPropertyWidget_;

    delete lProperty_;
    delete rProperty_;

    lineLayout_->removeWidget(interpolationFunction_);
    delete interpolationFunction_;

    interpolationMenu_->clear();
    delete interpolationMenu_;
}

bool AnimationIntervalLineWidget::isValid() const {
    return (lProperty_ && rProperty_ && lPropertyWidget_ && rPropertyWidget_);
}

void AnimationIntervalLineWidget::updateActivationStatus() {

    bool active = propertyTimeline_->getActiveOnRendering();

    lPropertyWidget_->setEnabled(active);
    rPropertyWidget_->setEnabled(active);
    interpolationFunction_->setEnabled(active);

    AnimationLineWidgetBase::updateActivationStatus();
}

void AnimationIntervalLineWidget::setFixedWidthSlot(int width) {

    AnimationLineWidgetBase::setFixedWidthSlot(width);

    int componentWidth = std::max(10, tgt::iround(static_cast<float>(width - 166) / 5.f));

    lineLayout_->setColumnMinimumWidth(3, componentWidth);
    lineLayout_->setColumnMinimumWidth(4, componentWidth);
    lineLayout_->setColumnMinimumWidth(5, componentWidth);

    interpolationFunction_->setFixedWidth(componentWidth);

    QPropertyWidget* lW = dynamic_cast<QPropertyWidget*>(lPropertyWidget_);
    QPropertyWidget* rW = dynamic_cast<QPropertyWidget*>(rPropertyWidget_);

    lW->setFixedWidth(componentWidth);
    rW->setFixedWidth(componentWidth);
}

void AnimationIntervalLineWidget::updateProperties(PropertyKeyValueBase* lKey, PropertyKeyValueBase* rKey) {

    //remember the left and right keys
    leftKey_ = lKey;
    rightKey_ = rKey;

    //update properties and interpolation
    if (informationWidget_) {
        informationWidget_->updatePropertyFromKey(lProperty_, lPropertyWidget_, lKey, propertyTimeline_);
        informationWidget_->updatePropertyFromKey(rProperty_, rPropertyWidget_, rKey, propertyTimeline_);
    }

    const InterpolationFunctionBase* func = 0;

    //update interpolation function
    if (TemplatePropertyTimeline<float>* tl = dynamic_cast<TemplatePropertyTimeline<float>* >(propertyTimeline_)) {
        func = dynamic_cast<PropertyKeyValue<float>*>(lKey)->getFollowingInterpolationFunction();
    }
    else if (TemplatePropertyTimeline<int>* tl = dynamic_cast<TemplatePropertyTimeline<int>* >(propertyTimeline_)) {
        func = dynamic_cast<PropertyKeyValue<int>*>(lKey)->getFollowingInterpolationFunction();
    }
    else if (TemplatePropertyTimeline<bool>* tl = dynamic_cast<TemplatePropertyTimeline<bool>* >(propertyTimeline_)) {
        func = dynamic_cast<PropertyKeyValue<bool>*>(lKey)->getFollowingInterpolationFunction();
    }
    else if (TemplatePropertyTimeline<ivec2>* tl = dynamic_cast<TemplatePropertyTimeline<ivec2>* >(propertyTimeline_)) {
        func = dynamic_cast<PropertyKeyValue<ivec2>*>(lKey)->getFollowingInterpolationFunction();
    }
    else if (TemplatePropertyTimeline<ivec3>* tl = dynamic_cast<TemplatePropertyTimeline<ivec3>* >(propertyTimeline_)) {
        func = dynamic_cast<PropertyKeyValue<ivec3>*>(lKey)->getFollowingInterpolationFunction();
    }
    else if (TemplatePropertyTimeline<ivec4>* tl = dynamic_cast<TemplatePropertyTimeline<ivec4>* >(propertyTimeline_)) {
        func = dynamic_cast<PropertyKeyValue<ivec4>*>(lKey)->getFollowingInterpolationFunction();
    }
    else if (TemplatePropertyTimeline<vec2>* tl = dynamic_cast<TemplatePropertyTimeline<vec2>* >(propertyTimeline_)) {
        func = dynamic_cast<PropertyKeyValue<vec2>*>(lKey)->getFollowingInterpolationFunction();
    }
    else if (TemplatePropertyTimeline<vec3>* tl = dynamic_cast<TemplatePropertyTimeline<vec3>* >(propertyTimeline_)) {
        func = dynamic_cast<PropertyKeyValue<vec3>*>(lKey)->getFollowingInterpolationFunction();
    }
    else if (TemplatePropertyTimeline<vec4>* tl = dynamic_cast<TemplatePropertyTimeline<vec4>* >(propertyTimeline_)) {
        func = dynamic_cast<PropertyKeyValue<vec4>*>(lKey)->getFollowingInterpolationFunction();
    }
    else if (TemplatePropertyTimeline<string>* tl = dynamic_cast<TemplatePropertyTimeline<string>* >(propertyTimeline_)) {
        func = dynamic_cast<PropertyKeyValue<string>*>(lKey)->getFollowingInterpolationFunction();
    }
    else if (TemplatePropertyTimeline<ShaderSource>* tl = dynamic_cast<TemplatePropertyTimeline<ShaderSource>* >(propertyTimeline_)) {
        func = dynamic_cast<PropertyKeyValue<ShaderSource>*>(lKey)->getFollowingInterpolationFunction();
    }
    else if (TemplatePropertyTimeline<TransFunc1DKeys*>* tl = dynamic_cast<TemplatePropertyTimeline<TransFunc1DKeys*>* >(propertyTimeline_)) {
        func = dynamic_cast<PropertyKeyValue<TransFunc1DKeys*>*>(lKey)->getFollowingInterpolationFunction();
    }
    else if (TemplatePropertyTimeline<TransFunc2DPrimitives*>* tl = dynamic_cast<TemplatePropertyTimeline<TransFunc2DPrimitives*>* >(propertyTimeline_)) {
        func = dynamic_cast<PropertyKeyValue<TransFunc2DPrimitives*>*>(lKey)->getFollowingInterpolationFunction();
    }
    else if (TemplatePropertyTimeline<Camera>* tl = dynamic_cast<TemplatePropertyTimeline<Camera>* >(propertyTimeline_)) {
        func = dynamic_cast<PropertyKeyValue<Camera>*>(lKey)->getFollowingInterpolationFunction();
    }
    else if (TemplatePropertyTimeline<tgt::IntBounds>* tl = dynamic_cast<TemplatePropertyTimeline<tgt::IntBounds>* >(propertyTimeline_)){
        func = dynamic_cast<PropertyKeyValue<tgt::IntBounds>*>(lKey)->getFollowingInterpolationFunction();
    }
    else if (TemplatePropertyTimeline<tgt::Bounds>* tl = dynamic_cast<TemplatePropertyTimeline<tgt::Bounds>* >(propertyTimeline_)){
        func = dynamic_cast<PropertyKeyValue<tgt::Bounds>*>(lKey)->getFollowingInterpolationFunction();
    }

    if (func)
        interpolationFunction_->setText(QString::fromStdString(func->getGuiName()));
}

void AnimationIntervalLineWidget::updateLeftKeyframe() {
    if (informationWidget_)
        informationWidget_->updateKeyFromWidget(propertyTimeline_, leftKey_, lProperty_);
}

void AnimationIntervalLineWidget::updateRightKeyframe() {
    if (informationWidget_)
        informationWidget_->updateKeyFromWidget(propertyTimeline_, rightKey_, rProperty_);
}

void AnimationIntervalLineWidget::populateInterpolationMenu() {

    if (TemplatePropertyTimeline<float>* tl = dynamic_cast<TemplatePropertyTimeline<float>* >(propertyTimeline_)) {
        populateTemplateInterpolationMenu<float>();
    }
    else if (TemplatePropertyTimeline<int>* tl = dynamic_cast<TemplatePropertyTimeline<int>* >(propertyTimeline_)) {
        populateTemplateInterpolationMenu<int>();
    }
    else if (TemplatePropertyTimeline<bool>* tl = dynamic_cast<TemplatePropertyTimeline<bool>* >(propertyTimeline_)) {
        populateTemplateInterpolationMenu<bool>();
    }
    else if (TemplatePropertyTimeline<ivec2>* tl = dynamic_cast<TemplatePropertyTimeline<ivec2>* >(propertyTimeline_)) {
        populateTemplateInterpolationMenu<ivec2>();
    }
    else if (TemplatePropertyTimeline<ivec3>* tl = dynamic_cast<TemplatePropertyTimeline<ivec3>* >(propertyTimeline_)) {
        populateTemplateInterpolationMenu<ivec3>();
    }
    else if (TemplatePropertyTimeline<ivec4>* tl = dynamic_cast<TemplatePropertyTimeline<ivec4>* >(propertyTimeline_)) {
        populateTemplateInterpolationMenu<ivec4>();
    }
    else if (TemplatePropertyTimeline<vec2>* tl = dynamic_cast<TemplatePropertyTimeline<vec2>* >(propertyTimeline_)) {
        populateTemplateInterpolationMenu<vec2>();
    }
    else if (TemplatePropertyTimeline<vec3>* tl = dynamic_cast<TemplatePropertyTimeline<vec3>* >(propertyTimeline_)) {
        populateTemplateInterpolationMenu<vec3>();
    }
    else if (TemplatePropertyTimeline<vec4>* tl = dynamic_cast<TemplatePropertyTimeline<vec4>* >(propertyTimeline_)) {
        populateTemplateInterpolationMenu<vec4>();
    }
    else if (TemplatePropertyTimeline<string>* tl = dynamic_cast<TemplatePropertyTimeline<string>* >(propertyTimeline_)) {
        populateTemplateInterpolationMenu<string>();
    }
    else if (TemplatePropertyTimeline<ShaderSource>* tl = dynamic_cast<TemplatePropertyTimeline<ShaderSource>* >(propertyTimeline_)) {
        populateTemplateInterpolationMenu<ShaderSource>();
    }
    else if (TemplatePropertyTimeline<TransFunc1DKeys*>* tl = dynamic_cast<TemplatePropertyTimeline<TransFunc1DKeys*>* >(propertyTimeline_)) {
        populateTemplateInterpolationMenu<TransFunc1DKeys*>();
    }
    else if (TemplatePropertyTimeline<TransFunc2DPrimitives*>* tl = dynamic_cast<TemplatePropertyTimeline<TransFunc2DPrimitives*>* >(propertyTimeline_)) {
        populateTemplateInterpolationMenu<TransFunc2DPrimitives*>();
    }
    else if (TemplatePropertyTimeline<Camera>* tl = dynamic_cast<TemplatePropertyTimeline<Camera>* >(propertyTimeline_)) {
        populateTemplateInterpolationMenu<Camera>();
    }
    else if (TemplatePropertyTimeline<tgt::IntBounds>* tl = dynamic_cast<TemplatePropertyTimeline<tgt::IntBounds>* >(propertyTimeline_)){
        populateTemplateInterpolationMenu<tgt::IntBounds>();
    }
    else if (TemplatePropertyTimeline<tgt::Bounds>* tl = dynamic_cast<TemplatePropertyTimeline<tgt::Bounds>* >(propertyTimeline_)){
        populateTemplateInterpolationMenu<tgt::Bounds>();
    }
    else
        LERRORC("voreen.qt.AnimationIntervalLineWidget", "Could not set interpolation menu, unknown type of property");
}

void AnimationIntervalLineWidget::changeInterpolation(QAction* action) {

    if (TemplatePropertyTimeline<float>* tl = dynamic_cast<TemplatePropertyTimeline<float>* >(propertyTimeline_)) {
        changeInterpolationTemplate<float>(action);
    }
    else if (TemplatePropertyTimeline<int>* tl = dynamic_cast<TemplatePropertyTimeline<int>* >(propertyTimeline_)) {
        changeInterpolationTemplate<int>(action);
    }
    else if (TemplatePropertyTimeline<bool>* tl = dynamic_cast<TemplatePropertyTimeline<bool>* >(propertyTimeline_)) {
        changeInterpolationTemplate<bool>(action);
    }
    else if (TemplatePropertyTimeline<ivec2>* tl = dynamic_cast<TemplatePropertyTimeline<ivec2>* >(propertyTimeline_)) {
        changeInterpolationTemplate<ivec2>(action);
    }
    else if (TemplatePropertyTimeline<ivec3>* tl = dynamic_cast<TemplatePropertyTimeline<ivec3>* >(propertyTimeline_)) {
        changeInterpolationTemplate<ivec3>(action);
    }
    else if (TemplatePropertyTimeline<ivec4>* tl = dynamic_cast<TemplatePropertyTimeline<ivec4>* >(propertyTimeline_)) {
        changeInterpolationTemplate<ivec4>(action);
    }
    else if (TemplatePropertyTimeline<vec2>* tl = dynamic_cast<TemplatePropertyTimeline<vec2>* >(propertyTimeline_)) {
        changeInterpolationTemplate<vec2>(action);
    }
    else if (TemplatePropertyTimeline<vec3>* tl = dynamic_cast<TemplatePropertyTimeline<vec3>* >(propertyTimeline_)) {
        changeInterpolationTemplate<vec3>(action);
    }
    else if (TemplatePropertyTimeline<vec4>* tl = dynamic_cast<TemplatePropertyTimeline<vec4>* >(propertyTimeline_)) {
        changeInterpolationTemplate<vec4>(action);
    }
    else if (TemplatePropertyTimeline<string>* tl = dynamic_cast<TemplatePropertyTimeline<string>* >(propertyTimeline_)) {
        changeInterpolationTemplate<string>(action);
    }
    else if (TemplatePropertyTimeline<ShaderSource>* tl = dynamic_cast<TemplatePropertyTimeline<ShaderSource>* >(propertyTimeline_)) {
        changeInterpolationTemplate<ShaderSource>(action);
    }
    else if (TemplatePropertyTimeline<TransFunc1DKeys*>* tl = dynamic_cast<TemplatePropertyTimeline<TransFunc1DKeys*>* >(propertyTimeline_)) {
        changeInterpolationTemplate<TransFunc1DKeys*>(action);
    }
    else if (TemplatePropertyTimeline<TransFunc2DPrimitives*>* tl = dynamic_cast<TemplatePropertyTimeline<TransFunc2DPrimitives*>* >(propertyTimeline_)) {
        changeInterpolationTemplate<TransFunc2DPrimitives*>(action);
    }
    else if (TemplatePropertyTimeline<Camera>* tl = dynamic_cast<TemplatePropertyTimeline<Camera>* >(propertyTimeline_)) {
        changeInterpolationTemplate<Camera>(action);
    }
    else if (TemplatePropertyTimeline<tgt::IntBounds>* tl = dynamic_cast<TemplatePropertyTimeline<tgt::IntBounds>* >(propertyTimeline_)){
        changeInterpolationTemplate<tgt::IntBounds>(action);
    }
    else if (TemplatePropertyTimeline<tgt::Bounds>* tl = dynamic_cast<TemplatePropertyTimeline<tgt::Bounds>* >(propertyTimeline_)){
        changeInterpolationTemplate<tgt::Bounds>(action);
    }
}


//-----------------------------------------------------------------------
// AnimationIntervalWidget
//-----------------------------------------------------------------------


AnimationIntervalWidget::AnimationIntervalWidget(QWidget* parent) : AnimationInfoWidgetBase(parent) {}

void AnimationIntervalWidget::createHeader() {

    headerWidget_ = new QWidget(this);
    headerLayout_ = new QGridLayout(headerWidget_);
    headerLayout_->setMargin(0);
    headerLayout_->setSpacing(30);

    //QSpacerItem* spacer = new QSpacerItem(46, 0);

    QLabel* nameLabel = new QLabel(this);
    nameLabel->setText(QString("Property Name"));
    nameLabel->setDisabled(true);

    QLabel* lValueLabel = new QLabel(this);
    lValueLabel->setText(QString("Left Key Value"));
    lValueLabel->setDisabled(true);

    QLabel* intervalLabel = new QLabel(this);
    intervalLabel->setText(QString("Interpolation Function"));
    intervalLabel->setDisabled(true);

    QLabel* rValueLabel = new QLabel(this);
    rValueLabel->setText(QString("Right Key Value"));
    rValueLabel->setDisabled(true);

    //headerLayout_->addItem(spacer, 0, 0);
    headerLayout_->addWidget(nameLabel, 0, 2);
    headerLayout_->addWidget(lValueLabel, 0, 3);
    headerLayout_->addWidget(intervalLabel, 0, 4);
    headerLayout_->addWidget(rValueLabel, 0, 5);
}

void AnimationIntervalWidget::resizeHeader(int width) {

    if (!headerLayout_)
        return;

    int componentWidth = std::max(10, tgt::iround(static_cast<float>(width - 166) / 5.f));

    headerLayout_->setColumnMinimumWidth(0, 16);
    headerLayout_->setColumnMinimumWidth(1, 135);
    headerLayout_->setColumnMinimumWidth(2, tgt::iround(1.42 * static_cast<float>(componentWidth)));
    headerLayout_->setColumnMinimumWidth(3, componentWidth);
    headerLayout_->setColumnMinimumWidth(4, componentWidth);
    headerLayout_->setColumnMinimumWidth(5, componentWidth);
}

void AnimationIntervalWidget::clear() {

    //remove all group boxes from layout, and delete the corresponding line widgets
    for (std::map<std::string, QGroupBox*>::iterator i = applicationModeGroups_.begin(); i != applicationModeGroups_.end(); ++i) {
        mainLayout_->removeWidget(i->second);
        std::map<std::string, std::vector<std::pair<int, PropertyTimeline*> > >::iterator groupTimelineIter = groupTimelines_.find(i->first);
        std::vector<std::pair<int, PropertyTimeline*> >::iterator timelineIter = groupTimelineIter->second.begin();
        for ( ; timelineIter != groupTimelineIter->second.end(); ++timelineIter) {
            std::map<PropertyTimeline*, AnimationIntervalLineWidget*>::iterator lineWidgetIter = timelineWidgets_.find(timelineIter->second);
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

void AnimationIntervalWidget::addPropertyTimeline(PropertyTimeline* prop) {

    if (!appConfig_) {
        LERRORC("voreen.qt.SimpleKeyframeWidget", "Could not add property timeline. No application mode configuration found.");
        return;
    }

    AnimationIntervalLineWidget* newLine = new AnimationIntervalLineWidget(this, prop);

    if (!newLine->isValid()) {
        LERRORC("voreen.qt.AnimationIntervalWidget", "Could not create PropertyWidget for " + prop->getProperty()->getGuiName());
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

void AnimationIntervalWidget::updateProperties(const std::vector<std::pair<PropertyKeyValueBase*, PropertyTimeline*> >& leftPropertyValues,
                                            const std::vector<std::pair<PropertyKeyValueBase*, PropertyTimeline*> >& rightPropertyValues) {


    //iterate over property timelines
    for (size_t i = 0; i < leftPropertyValues.size(); ++i) {
        //for this property timeline: geht the corresponding entry in the map and update property and widget
        std::map<PropertyTimeline*, AnimationIntervalLineWidget*>::iterator propIt = timelineWidgets_.find(leftPropertyValues.at(i).second);

        if (propIt != timelineWidgets_.end()) {
            propIt->second->updateProperties(leftPropertyValues.at(i).first, rightPropertyValues.at(i).first);
        }
        else
            LERRORC("voreen.qt.AnimationIntervalWidget", "PropertyWidget not found for " + leftPropertyValues.at(i).second->getProperty()->getGuiName());
    }
}

void AnimationIntervalWidget::updateActivationStatus() {

    std::map<PropertyTimeline*, AnimationIntervalLineWidget*>::iterator it = timelineWidgets_.begin();
    for ( ; it != timelineWidgets_.end(); ++it) {
        it->second->updateActivationStatus();
    }
}

} // namespace voreen
