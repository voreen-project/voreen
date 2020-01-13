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
#include "voreen/core/properties/transfunc/1d/1dkeys/transfunc1dkeysproperty.h"
#include "voreen/core/properties/transfunc/2d/2dprimitives/transfunc2dprimitivesproperty.h"
#include "voreen/core/properties/cameraproperty.h"

#include "voreen/qt/animation/animationinfowidgetbase.h"
#include "voreen/core/voreenapplication.h"
#include "voreen/core/voreenmodule.h"
#include "voreen/qt/widgets/property/qpropertywidget.h"
#include "voreen/core/properties/property.h"
#include "voreen/core/properties/templateproperty.h"
#include "voreen/core/properties/lightsourceproperty.h"
#include "voreen/core/properties/boundingboxproperty.h"
#include "voreen/core/properties/numeric/intervalproperty.h"

#include "voreen/core/datastructures/volume/slice/slicehelper.h"

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
// AnimationLineWidgetBase
//---------------------------------------------------------------------

AnimationLineWidgetBase::AnimationLineWidgetBase(AnimationInfoWidgetBase* parent, PropertyTimeline* tl)
        : QWidget(parent)
        , propertyTimeline_(tl)
        , informationWidget_(parent)
{
    setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);

    lineLayout_ = new QGridLayout();
    lineLayout_->setMargin(0);
    lineLayout_->setSpacing(30);
    setLayout(lineLayout_);

    activateTimelineButton_ = new QPushButton(QIcon(":/qt/icons/apply.png"), "", this);
    activateTimelineButton_->setToolTip(tr("Activate or deactivate timeline"));
    activateTimelineButton_->setStyleSheet("QToolButton { border: none; padding: 1px; }");
    activateTimelineButton_->setFlat(true);
    activateTimelineButton_->setFocusPolicy(Qt::NoFocus);
    activateTimelineButton_->setCheckable(true);
    activateTimelineButton_->setMaximumWidth(16);
    activateTimelineButton_->setMaximumHeight(16);
    activateTimelineButton_->setChecked(tl->getActiveOnRendering());

    connect(activateTimelineButton_, SIGNAL(clicked(bool)), this, SLOT(timelineActivated()));

    nameLabel_ = new QLabel(this);
    nameLabel_->setFixedWidth(300);
    nameLabel_->setText(QString(tl->getProperty()->getGuiName().c_str()));

    deleteButton_ = new QPushButton("Remove", this);
    connect(deleteButton_, SIGNAL(clicked(bool)), this, SLOT(deleteButtonClicked()));

    lineLayout_->addWidget(activateTimelineButton_, 0, 0);
    lineLayout_->addWidget(deleteButton_, 0, 1);
    lineLayout_->addWidget(nameLabel_, 0, 2);

    lineLayout_->update();
}

AnimationLineWidgetBase::~AnimationLineWidgetBase() {

    lineLayout_->removeWidget(activateTimelineButton_);
    delete activateTimelineButton_;

    lineLayout_->removeWidget(nameLabel_);
    delete nameLabel_;
}

void AnimationLineWidgetBase::updateActivationStatus() {

    bool active = propertyTimeline_->getActiveOnRendering();
    activateTimelineButton_->setChecked(active);

    if (!active)
        activateTimelineButton_->setIcon(QIcon(":/qt/icons/button_cancel.png"));
    else
        activateTimelineButton_->setIcon(QIcon(":/qt/icons/apply.png"));

    lineLayout_->update();
}

void AnimationLineWidgetBase::setFixedWidthSlot(int width) {
    setFixedWidth(width-30);

    int componentWidth = std::max(10, tgt::iround(static_cast<float>(width - 166) / 5.f));

    lineLayout_->setColumnMinimumWidth(2, componentWidth);

    nameLabel_->setFixedWidth(2*componentWidth);
}

void AnimationLineWidgetBase::timelineActivated() {
    //update property timeline in core and then update gui
    propertyTimeline_->setActiveOnRendering(activateTimelineButton_->isChecked());
    updateActivationStatus();
}

void AnimationLineWidgetBase::deleteButtonClicked() {
    emit deleteTimelineSignal(propertyTimeline_);
}

//-----------------------------------------------------------------------
// AnimationInfoWidgetBase
//-----------------------------------------------------------------------

AnimationInfoWidgetBase::AnimationInfoWidgetBase(QWidget* parent)
        : QWidget(parent)
        , appConfig_(0)
        , mainLayout_(0)
        , headerWidget_(0)
        , headerLayout_(0)
{

    setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
    mainLayout_ = new QVBoxLayout(this);
    mainLayout_->setMargin(0);
    mainLayout_->setSpacing(0);
    mainLayout_->setAlignment(Qt::AlignTop);
}

QSize AnimationInfoWidgetBase::sizeHint() const {
    return QSize(800, 100);
}

void AnimationInfoWidgetBase::setFixedWidthSlot(int width) {
    setFixedWidth(std::max(width-30, 0));
    resizeHeader(width);
    emit setFixedWidthSignal(width - 20);
}

void AnimationInfoWidgetBase::removeWidgetsFromLayout(QLayout* layout) {
    QLayoutItem *child;
    while ((child = layout->takeAt(0)) != 0) {
        if (child->layout() != 0)
            removeWidgetsFromLayout(child->layout());
        else if (child->widget() != 0)
            delete child->widget();

        delete child;
    }
}

void AnimationInfoWidgetBase::updateKeyFromWidget(PropertyTimeline* timeline, const PropertyKeyValueBase* key, Property* prop) {
    //check type and update key
    if (TemplatePropertyTimeline<float>* tl = dynamic_cast<TemplatePropertyTimeline<float>* >(timeline)) {
        const PropertyKeyValue<float>* k = dynamic_cast<const PropertyKeyValue<float>* >(key);
        TemplateProperty<float>* p = dynamic_cast<TemplateProperty<float>* >(prop);

        if (!k || !p) {
            LERRORC("voreen.qt.SimpleKeyframeWidget", "types of propertytimeline, property and key value do not match - cannot update keyframe");
            return;
        }
        bool set = tl->changeValueOfKeyValue(p->get(), k);
        //bool set = tl->changeValueOfKeyValue(p->get().clone(), k);
        if (!set)
            LWARNINGC("voreen.qt.SimpleKeyframeWidget", "Could not set key value for property " + prop->getGuiName());
    }
    else if (TemplatePropertyTimeline<int>* tl = dynamic_cast<TemplatePropertyTimeline<int>* >(timeline)) {
        const PropertyKeyValue<int>* k = dynamic_cast<const PropertyKeyValue<int>* >(key);
        TemplateProperty<int>* p = dynamic_cast<TemplateProperty<int>* >(prop);

        if (!k || !p) {
            LERRORC("voreen.qt.SimpleKeyframeWidget", "types of propertytimeline, property and key value do not match - cannot update keyframe");
            return;
        }
        bool set = tl->changeValueOfKeyValue(p->get(), k);
        //bool set = tl->changeValueOfKeyValue(p->get().clone(), k);
        if (!set)
            LWARNINGC("voreen.qt.SimpleKeyframeWidget", "Could not set key value for property " + prop->getGuiName());
    }
    else if (TemplatePropertyTimeline<bool>* tl = dynamic_cast<TemplatePropertyTimeline<bool>* >(timeline)) {
        const PropertyKeyValue<bool>* k = dynamic_cast<const PropertyKeyValue<bool>* >(key);
        TemplateProperty<bool>* p = dynamic_cast<TemplateProperty<bool>* >(prop);

        if (!k || !p) {
            LERRORC("voreen.qt.SimpleKeyframeWidget", "types of propertytimeline, property and key value do not match - cannot update keyframe");
            return;
        }
        bool set = tl->changeValueOfKeyValue(p->get(), k);
        //bool set = tl->changeValueOfKeyValue(p->get().clone(), k);
        if (!set)
            LWARNINGC("voreen.qt.SimpleKeyframeWidget", "Could not set key value for property " + prop->getGuiName());
    }
    else if (TemplatePropertyTimeline<ivec2>* tl = dynamic_cast<TemplatePropertyTimeline<ivec2>* >(timeline)) {
        const PropertyKeyValue<ivec2>* k = dynamic_cast<const PropertyKeyValue<ivec2>* >(key);
        TemplateProperty<ivec2>* p = dynamic_cast<TemplateProperty<ivec2>* >(prop);

        if (!k || !p) {
            LERRORC("voreen.qt.SimpleKeyframeWidget", "types of propertytimeline, property and key value do not match - cannot update keyframe");
            return;
        }
        bool set = tl->changeValueOfKeyValue(p->get(), k);
        //bool set = tl->changeValueOfKeyValue(p->get().clone(), k);
        if (!set)
            LWARNINGC("voreen.qt.SimpleKeyframeWidget", "Could not set key value for property " + prop->getGuiName());
    }
    else if (TemplatePropertyTimeline<ivec3>* tl = dynamic_cast<TemplatePropertyTimeline<ivec3>* >(timeline)) {
        const PropertyKeyValue<ivec3>* k = dynamic_cast<const PropertyKeyValue<ivec3>* >(key);
        TemplateProperty<ivec3>* p = dynamic_cast<TemplateProperty<ivec3>* >(prop);

        if (!k || !p) {
            LERRORC("voreen.qt.SimpleKeyframeWidget", "types of propertytimeline, property and key value do not match - cannot update keyframe");
            return;
        }
        bool set = tl->changeValueOfKeyValue(p->get(), k);
        //bool set = tl->changeValueOfKeyValue(p->get().clone(), k);
        if (!set)
            LWARNINGC("voreen.qt.SimpleKeyframeWidget", "Could not set key value for property " + prop->getGuiName());
    }
    else if (TemplatePropertyTimeline<ivec4>* tl = dynamic_cast<TemplatePropertyTimeline<ivec4>* >(timeline)) {
        const PropertyKeyValue<ivec4>* k = dynamic_cast<const PropertyKeyValue<ivec4>* >(key);
        TemplateProperty<ivec4>* p = dynamic_cast<TemplateProperty<ivec4>* >(prop);

        if (!k || !p) {
            LERRORC("voreen.qt.SimpleKeyframeWidget", "types of propertytimeline, property and key value do not match - cannot update keyframe");
            return;
        }
        bool set = tl->changeValueOfKeyValue(p->get(), k);
        //bool set = tl->changeValueOfKeyValue(p->get().clone(), k);
        if (!set)
            LWARNINGC("voreen.qt.SimpleKeyframeWidget", "Could not set key value for property " + prop->getGuiName());
    }
    else if (TemplatePropertyTimeline<vec2>* tl = dynamic_cast<TemplatePropertyTimeline<vec2>* >(timeline)) {
        const PropertyKeyValue<vec2>* k = dynamic_cast<const PropertyKeyValue<vec2>* >(key);
        TemplateProperty<vec2>* p = dynamic_cast<TemplateProperty<vec2>* >(prop);

        if (!k || !p) {
            LERRORC("voreen.qt.SimpleKeyframeWidget", "types of propertytimeline, property and key value do not match - cannot update keyframe");
            return;
        }
        bool set = tl->changeValueOfKeyValue(p->get(), k);
        //bool set = tl->changeValueOfKeyValue(p->get().clone(), k);
        if (!set)
            LWARNINGC("voreen.qt.SimpleKeyframeWidget", "Could not set key value for property " + prop->getGuiName());
    }
    else if (TemplatePropertyTimeline<vec3>* tl = dynamic_cast<TemplatePropertyTimeline<vec3>* >(timeline)) {
        const PropertyKeyValue<vec3>* k = dynamic_cast<const PropertyKeyValue<vec3>* >(key);
        TemplateProperty<vec3>* p = dynamic_cast<TemplateProperty<vec3>* >(prop);

        if (!k || !p) {
            LERRORC("voreen.qt.SimpleKeyframeWidget", "types of propertytimeline, property and key value do not match - cannot update keyframe");
            return;
        }
        bool set = tl->changeValueOfKeyValue(p->get(), k);
        //bool set = tl->changeValueOfKeyValue(p->get().clone(), k);
        if (!set)
            LWARNINGC("voreen.qt.SimpleKeyframeWidget", "Could not set key value for property " + prop->getGuiName());
    }
    else if (TemplatePropertyTimeline<vec4>* tl = dynamic_cast<TemplatePropertyTimeline<vec4>* >(timeline)) {
        const PropertyKeyValue<vec4>* k = dynamic_cast<const PropertyKeyValue<vec4>* >(key);
        TemplateProperty<vec4>* p = dynamic_cast<TemplateProperty<vec4>* >(prop);

        if (!k || !p) {
            LERRORC("voreen.qt.SimpleKeyframeWidget", "types of propertytimeline, property and key value do not match - cannot update keyframe");
            return;
        }
        bool set = tl->changeValueOfKeyValue(p->get(), k);
        //bool set = tl->changeValueOfKeyValue(p->get().clone(), k);
        if (!set)
            LWARNINGC("voreen.qt.SimpleKeyframeWidget", "Could not set key value for property " + prop->getGuiName());
    }
    else if (TemplatePropertyTimeline<string>* tl = dynamic_cast<TemplatePropertyTimeline<string>* >(timeline)) {
        const PropertyKeyValue<string>* k = dynamic_cast<const PropertyKeyValue<string>* >(key);
        TemplateProperty<string>* p = dynamic_cast<TemplateProperty<string>* >(prop);

        if (!k || !p) {
            LERRORC("voreen.qt.SimpleKeyframeWidget", "types of propertytimeline, property and key value do not match - cannot update keyframe");
            return;
        }
        bool set = tl->changeValueOfKeyValue(p->get(), k);
        //bool set = tl->changeValueOfKeyValue(p->get().clone(), k);
        if (!set)
            LWARNINGC("voreen.qt.SimpleKeyframeWidget", "Could not set key value for property " + prop->getGuiName());
    }
    else if (TemplatePropertyTimeline<ShaderSource>* tl = dynamic_cast<TemplatePropertyTimeline<ShaderSource>* >(timeline)) {
        const PropertyKeyValue<ShaderSource>* k = dynamic_cast<const PropertyKeyValue<ShaderSource>* >(key);
        TemplateProperty<ShaderSource>* p = dynamic_cast<TemplateProperty<ShaderSource>* >(prop);

        if (!k || !p) {
            LERRORC("voreen.qt.SimpleKeyframeWidget", "types of propertytimeline, property and key value do not match - cannot update keyframe");
            return;
        }
        bool set = tl->changeValueOfKeyValue(p->get(), k);
        //bool set = tl->changeValueOfKeyValue(p->get().clone(), k);
        if (!set)
            LWARNINGC("voreen.qt.SimpleKeyframeWidget", "Could not set key value for property " + prop->getGuiName());
    }
    else if (TemplatePropertyTimeline<TransFunc1DKeys*>* tl = dynamic_cast<TemplatePropertyTimeline<TransFunc1DKeys*>* >(timeline)) {
        const PropertyKeyValue<TransFunc1DKeys*>* k = dynamic_cast<const PropertyKeyValue<TransFunc1DKeys*>* >(key);
        TemplateProperty<TransFunc1DKeys*>* p = dynamic_cast<TemplateProperty<TransFunc1DKeys*>* >(prop);

        if (!k || !p) {
            LERRORC("voreen.qt.SimpleKeyframeWidget", "types of propertytimeline, property and key value do not match - cannot update keyframe");
            return;
        }
        bool set = tl->changeValueOfKeyValue(p->get()->clone(), k);
        //bool set = tl->changeValueOfKeyValue(p->get(), k);
        if (!set)
            LWARNINGC("voreen.qt.SimpleKeyframeWidget", "Could not set key value for property " + prop->getGuiName());
    }
    else if (TemplatePropertyTimeline<TransFunc2DPrimitives*>* tl = dynamic_cast<TemplatePropertyTimeline<TransFunc2DPrimitives*>* >(timeline)) {
        const PropertyKeyValue<TransFunc2DPrimitives*>* k = dynamic_cast<const PropertyKeyValue<TransFunc2DPrimitives*>* >(key);
        TemplateProperty<TransFunc2DPrimitives*>* p = dynamic_cast<TemplateProperty<TransFunc2DPrimitives*>* >(prop);

        if (!k || !p) {
            LERRORC("voreen.qt.SimpleKeyframeWidget", "types of propertytimeline, property and key value do not match - cannot update keyframe");
            return;
        }
        bool set = tl->changeValueOfKeyValue(p->get()->clone(), k);
        //bool set = tl->changeValueOfKeyValue(p->get(), k);
        if (!set)
            LWARNINGC("voreen.qt.SimpleKeyframeWidget", "Could not set key value for property " + prop->getGuiName());
    }
    else if (TemplatePropertyTimeline<Camera>* tl = dynamic_cast<TemplatePropertyTimeline<Camera>* >(timeline)) {
        const PropertyKeyValue<Camera>* k = dynamic_cast<const PropertyKeyValue<Camera>* >(key);
        TemplateProperty<Camera>* p = dynamic_cast<TemplateProperty<Camera>* >(prop);

        if (!k || !p) {
            LERRORC("voreen.qt.SimpleKeyframeWidget", "types of propertytimeline, property and key value do not match - cannot update keyframe");
            return;
        }
        bool set = tl->changeValueOfKeyValue(p->get(), k);
        //bool set = tl->changeValueOfKeyValue(p->get().clone(), k);
        if (!set)
            LWARNINGC("voreen.qt.SimpleKeyframeWidget", "Could not set key value for property " + prop->getGuiName());
    }
    else if (TemplatePropertyTimeline<tgt::IntBounds>* tl = dynamic_cast<TemplatePropertyTimeline<tgt::IntBounds>* >(timeline)) {
        const PropertyKeyValue<tgt::IntBounds>* k = dynamic_cast<const PropertyKeyValue<tgt::IntBounds>* >(key);
        TemplateProperty<tgt::IntBounds>* p = dynamic_cast<TemplateProperty<tgt::IntBounds>* >(prop);

        if (!k || !p) {
            LERRORC("voreen.qt.SimpleKeyframeWidget", "types of propertytimeline, property and key value do not match - cannot update keyframe");
            return;
        }
        bool set = tl->changeValueOfKeyValue(p->get(), k);
        //bool set = tl->changeValueOfKeyValue(p->get().clone(), k);
        if (!set)
            LWARNINGC("voreen.qt.SimpleKeyframeWidget", "Could not set key value for property " + prop->getGuiName());
    }
    else if (TemplatePropertyTimeline<tgt::Bounds>* tl = dynamic_cast<TemplatePropertyTimeline<tgt::Bounds>* >(timeline)) {
        const PropertyKeyValue<tgt::Bounds>* k = dynamic_cast<const PropertyKeyValue<tgt::Bounds>* >(key);
        TemplateProperty<tgt::Bounds>* p = dynamic_cast<TemplateProperty<tgt::Bounds>* >(prop);

        if (!k || !p) {
            LERRORC("voreen.qt.SimpleKeyframeWidget", "types of propertytimeline, property and key value do not match - cannot update keyframe");
            return;
        }
        bool set = tl->changeValueOfKeyValue(p->get(), k);
        //bool set = tl->changeValueOfKeyValue(p->get().clone(), k);
        if (!set)
            LWARNINGC("voreen.qt.SimpleKeyframeWidget", "Could not set key value for property " + prop->getGuiName());
    }
    else
        LERRORC("voreen.qt.SimpleKeyframeWidget", "Could not update key, unknown type of property");
}

void AnimationInfoWidgetBase::updatePropertyFromKey(Property* prop, PropertyWidget* widget, const PropertyKeyValueBase* key, PropertyTimeline* tl) {

    //check type and update
    if (NumericProperty<float>* p = dynamic_cast<NumericProperty<float>* >(prop)) {
        const PropertyKeyValue<float>* k = dynamic_cast<const PropertyKeyValue<float>* >(key);
        NumericProperty<float>* origProp = dynamic_cast<NumericProperty<float>* >(tl->getProperty());

        p->setMinValue(origProp->getMinValue());
        p->setMaxValue(origProp->getMaxValue());

        p->set(k->getValue());
        widget->updateFromProperty();
        p->updateWidgets();
    }
    else if (NumericProperty<int>* p = dynamic_cast<NumericProperty<int>* >(prop)) {
        const PropertyKeyValue<int>* k = dynamic_cast<const PropertyKeyValue<int>* >(key);
        NumericProperty<int>* origProp = dynamic_cast<NumericProperty<int>* >(tl->getProperty());

        p->setMinValue(origProp->getMinValue());
        p->setMaxValue(origProp->getMaxValue());

        p->set(k->getValue());
        widget->updateFromProperty();
        p->updateWidgets();
    }
    else if (TemplateProperty<bool>* p = dynamic_cast<TemplateProperty<bool>* >(prop)) {
        const PropertyKeyValue<bool>* k = dynamic_cast<const PropertyKeyValue<bool>* >(key);
        //p->set(k->getValue().clone());
        p->set(k->getValue());
        widget->updateFromProperty();
        p->updateWidgets();
    }
    else if (NumericProperty<ivec2>* p = dynamic_cast<NumericProperty<ivec2>* >(prop)) {
        const PropertyKeyValue<ivec2>* k = dynamic_cast<const PropertyKeyValue<ivec2>* >(key);
        NumericProperty<ivec2>* origProp = dynamic_cast<NumericProperty<ivec2>* >(tl->getProperty());

        p->setMinValue(origProp->getMinValue());
        p->setMaxValue(origProp->getMaxValue());

        p->set(k->getValue());
        widget->updateFromProperty();
        p->updateWidgets();
    }
    else if (NumericProperty<ivec3>* p = dynamic_cast<NumericProperty<ivec3>* >(prop)) {
        const PropertyKeyValue<ivec3>* k = dynamic_cast<const PropertyKeyValue<ivec3>* >(key);
        NumericProperty<ivec3>* origProp = dynamic_cast<NumericProperty<ivec3>* >(tl->getProperty());

        p->setMinValue(origProp->getMinValue());
        p->setMaxValue(origProp->getMaxValue());

        p->set(k->getValue());
        widget->updateFromProperty();
        p->updateWidgets();
    }
    else if (NumericProperty<ivec4>* p = dynamic_cast<NumericProperty<ivec4>* >(prop)) {
        const PropertyKeyValue<ivec4>* k = dynamic_cast<const PropertyKeyValue<ivec4>* >(key);
        NumericProperty<ivec4>* origProp = dynamic_cast<NumericProperty<ivec4>* >(tl->getProperty());

        p->setMinValue(origProp->getMinValue());
        p->setMaxValue(origProp->getMaxValue());

        p->set(k->getValue());
        widget->updateFromProperty();
        p->updateWidgets();
    }
    else if (NumericProperty<vec2>* p = dynamic_cast<NumericProperty<vec2>* >(prop)) {
        const PropertyKeyValue<vec2>* k = dynamic_cast<const PropertyKeyValue<vec2>* >(key);
        NumericProperty<vec2>* origProp = dynamic_cast<NumericProperty<vec2>* >(tl->getProperty());

        p->setMinValue(origProp->getMinValue());
        p->setMaxValue(origProp->getMaxValue());

        p->set(k->getValue());
        widget->updateFromProperty();
        p->updateWidgets();
    }
    else if (NumericProperty<vec3>* p = dynamic_cast<NumericProperty<vec3>* >(prop)) {
        const PropertyKeyValue<vec3>* k = dynamic_cast<const PropertyKeyValue<vec3>* >(key);
        NumericProperty<vec3>* origProp = dynamic_cast<NumericProperty<vec3>* >(tl->getProperty());

        p->setMinValue(origProp->getMinValue());
        p->setMaxValue(origProp->getMaxValue());

        p->set(k->getValue());
        widget->updateFromProperty();
        p->updateWidgets();
    }
    else if (NumericProperty<vec4>* p = dynamic_cast<NumericProperty<vec4>* >(prop)) {
        const PropertyKeyValue<vec4>* k = dynamic_cast<const PropertyKeyValue<vec4>* >(key);
        NumericProperty<vec4>* origProp = dynamic_cast<NumericProperty<vec4>* >(tl->getProperty());

        p->setMinValue(origProp->getMinValue());
        p->setMaxValue(origProp->getMaxValue());

        if (LightSourceProperty* lsProp = dynamic_cast<LightSourceProperty*>(prop)) {
            LightSourceProperty* origProp = dynamic_cast<LightSourceProperty*>(tl->getProperty());
            lsProp->setFollowCam(origProp->getFollowCam());
            lsProp->setMaxDist(origProp->getMaxDist());
        }

        p->set(k->getValue());
        widget->updateFromProperty();
        p->updateWidgets();
    }
    else if (TemplateProperty<string>* p = dynamic_cast<TemplateProperty<string>* >(prop)) {
        const PropertyKeyValue<string>* k = dynamic_cast<const PropertyKeyValue<string>* >(key);

        //check for option property types and set the options before setting the value
        if (OptionProperty<string>* optProp = dynamic_cast<OptionProperty<string>*>(prop)) {
            OptionProperty<string>* origProp = dynamic_cast<OptionProperty<string>*>(tl->getProperty());
            optProp->setOptions(origProp->getOptions());
        }
        else if (OptionProperty<int>* optProp = dynamic_cast<OptionProperty<int>*>(prop)) {
            OptionProperty<int>* origProp = dynamic_cast<OptionProperty<int>*>(tl->getProperty());
            optProp->setOptions(origProp->getOptions());
        }
        else if (OptionProperty<float>* optProp = dynamic_cast<OptionProperty<float>*>(prop)) {
            OptionProperty<float>* origProp = dynamic_cast<OptionProperty<float>*>(tl->getProperty());
            optProp->setOptions(origProp->getOptions());
        }
        else if (OptionProperty<bool>* optProp = dynamic_cast<OptionProperty<bool>*>(prop)) {
            OptionProperty<bool>* origProp = dynamic_cast<OptionProperty<bool>*>(tl->getProperty());
            optProp->setOptions(origProp->getOptions());
        }
        else if (OptionProperty<GLenum>* optProp = dynamic_cast<OptionProperty<GLenum>*>(prop)) {
            OptionProperty<GLenum>* origProp = dynamic_cast<OptionProperty<GLenum>*>(tl->getProperty());
            optProp->setOptions(origProp->getOptions());
        }
        else if (OptionProperty<SliceAlignment>* optProp = dynamic_cast<OptionProperty<SliceAlignment>*>(prop)) {
            OptionProperty<SliceAlignment>* origProp = dynamic_cast<OptionProperty<SliceAlignment>*>(tl->getProperty());
            optProp->setOptions(origProp->getOptions());
        }
        else if (OptionPropertyBase* optProp = dynamic_cast<OptionPropertyBase*>(prop)) {
            //type not supported -> cannot set options
            LERRORC("voreen.qt.InfoWidget", "OptionProperty type not supported: " + tl->getProperty()->getGuiName());
            return;
        }

        p->set(k->getValue());
        widget->updateFromProperty();
        p->updateWidgets();
    }
    else if (TemplateProperty<ShaderSource>* p = dynamic_cast<TemplateProperty<ShaderSource>* >(prop)) {
        const PropertyKeyValue<ShaderSource>* k = dynamic_cast<const PropertyKeyValue<ShaderSource>* >(key);
        //p->set(k->getValue().clone());
        p->set(k->getValue());
        widget->updateFromProperty();
        p->updateWidgets();
    }
    //else if (/*TemplateProperty<TransFuncBase*>* p = dynamic_cast<TemplateProperty<TransFuncBase*>* >(prop)*/) {
    else if (TransFunc1DKeysProperty* p = dynamic_cast<TransFunc1DKeysProperty*>(prop)) {
        const PropertyKeyValue<TransFunc1DKeys*>* k = dynamic_cast<const PropertyKeyValue<TransFunc1DKeys*>* >(key);
        p->set1DKeys(k->getValue()->clone());
        //p->set(k->getValue());
        TransFunc1DKeysProperty* refProp = dynamic_cast<TransFunc1DKeysProperty*>(tl->getProperty());
        if (refProp) {
            p->setVolume(refProp->getVolume());
            if (refProp->getMetaDataContainer().hasMetaData("TransfuncPropertyWidgetPainterZoom")) {
                p->getMetaDataContainer().addMetaData("TransfuncPropertyWidgetPainterZoom",
                        refProp->getMetaDataContainer().getMetaData("TransfuncPropertyWidgetPainterZoom")->clone());
            }
        }
        //p->set(k->getValue());
        widget->updateFromProperty();
        p->updateWidgets();
    }
    else if (TransFunc2DPrimitivesProperty* p = dynamic_cast<TransFunc2DPrimitivesProperty*>(prop)) {
        const PropertyKeyValue<TransFunc2DPrimitives*>* k = dynamic_cast<const PropertyKeyValue<TransFunc2DPrimitives*>* >(key);
        p->set2DPrimitives(k->getValue()->clone());
        //p->set(k->getValue());
        TransFunc2DPrimitivesProperty* refProp = dynamic_cast<TransFunc2DPrimitivesProperty*>(tl->getProperty());
        if (refProp) {
            p->setVolume(refProp->getVolume());
            if (refProp->getMetaDataContainer().hasMetaData("TransfuncPropertyWidgetPainterZoom")) {
                p->getMetaDataContainer().addMetaData("TransfuncPropertyWidgetPainterZoom",
                        refProp->getMetaDataContainer().getMetaData("TransfuncPropertyWidgetPainterZoom")->clone());
            }
        }
        //p->set(k->getValue());
        widget->updateFromProperty();
        p->updateWidgets();
    }
    else if (CameraProperty* p = dynamic_cast<CameraProperty*>(prop)) {
        const PropertyKeyValue<Camera>* k = dynamic_cast<const PropertyKeyValue<Camera>* >(key);
        //CameraProperty* refProp = dynamic_cast<CameraProperty*>(tl->getProperty());
        //p->setFrustum(k->getValue().getFrustum());
        p->set(k->getValue());
        widget->updateFromProperty();
        p->updateWidgets();
    }else if(IntBoundingBoxProperty* p = dynamic_cast<IntBoundingBoxProperty*>(prop)){
        const PropertyKeyValue<tgt::IntBounds>* k =
            dynamic_cast<const PropertyKeyValue<tgt::IntBounds >*>(key);
        IntBoundingBoxProperty* origProp = dynamic_cast<IntBoundingBoxProperty* >(tl->getProperty());
        p->setMinValue(origProp->getMinValue());
        p->setMaxValue(origProp->getMaxValue());
        p->setMinRange(origProp->getMinRange());
        p->setMaxRange(origProp->getMaxRange());
        p->setStepping(origProp->getStepping());
        p->set(k->getValue());
        widget->updateFromProperty();
        p->updateWidgets();

    }
    else if(FloatBoundingBoxProperty* p = dynamic_cast<FloatBoundingBoxProperty*>(prop)){
        const PropertyKeyValue<tgt::Bounds>* k =
            dynamic_cast<const PropertyKeyValue<tgt::Bounds >*>(key);
        FloatBoundingBoxProperty* origProp = dynamic_cast<FloatBoundingBoxProperty* >(tl->getProperty());
        p->setMinValue(origProp->getMinValue());
        p->setMaxValue(origProp->getMaxValue());
        p->setMinRange(origProp->getMinRange());
        p->setMaxRange(origProp->getMaxRange());
        p->setStepping(origProp->getStepping());
        p->set(k->getValue());
        widget->updateFromProperty();
        p->updateWidgets();

    }
    else if(FloatIntervalProperty* p = dynamic_cast<FloatIntervalProperty*>(prop)){
        const PropertyKeyValue<tgt::vec2>* k =
            dynamic_cast<const PropertyKeyValue<tgt::vec2 >*>(key);
        FloatIntervalProperty* origProp = dynamic_cast<FloatIntervalProperty* >(tl->getProperty());
        p->setMinValue(origProp->getMinValue());
        p->setMaxValue(origProp->getMaxValue());
        p->setMinRange(origProp->getMinRange());
        p->setMaxRange(origProp->getMaxRange());
        p->setStepping(origProp->getStepping());
        p->set(k->getValue());
        widget->updateFromProperty();
        p->updateWidgets();
    }
    else if(IntIntervalProperty* p = dynamic_cast<IntIntervalProperty*>(prop)){
        const PropertyKeyValue<tgt::ivec2>* k =
            dynamic_cast<const PropertyKeyValue<tgt::ivec2 >*>(key);
        IntIntervalProperty* origProp = dynamic_cast<IntIntervalProperty* >(tl->getProperty());
        p->setMinValue(origProp->getMinValue());
        p->setMaxValue(origProp->getMaxValue());
        p->setMinRange(origProp->getMinRange());
        p->setMaxRange(origProp->getMaxRange());
        p->setStepping(origProp->getStepping());
        p->set(k->getValue());
        widget->updateFromProperty();
        p->updateWidgets();
    }
    else
        LERRORC("voreen.qt.SimpleKeyframeWidget", "Could not update widget, unknown type of property");
}

void AnimationInfoWidgetBase::setApplicationModeConfig(ApplicationModeConfiguration* appConfig) {
    appConfig_ = appConfig;
}

bool AnimationInfoWidgetBase::sortTimelinesByPriority(const std::pair<int, PropertyTimeline*>& a, const std::pair<int, PropertyTimeline*>& b) {
    return (a.first < b.first);
}

} // namespace voreen
