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

#include <QGraphicsScene>
#include <QGraphicsView>
#include <QGroupBox>
#include <QMessageBox>
#include <QPushButton>
#include <QSlider>
#include <QSpinBox>
#include <QKeyEvent>
#include <QScrollBar>

#include "voreen/core/voreenapplication.h"
#include "voreen/core/voreenmodule.h"
#include "voreen/qt/widgets/property/qpropertywidget.h"
#include "voreen/core/animation/animatedprocessor.h"
#include "voreen/core/animation/templatepropertytimeline.h"
#include "voreen/core/animation/animatedprocessor.h"
#include "voreen/core/animation/interpolationfunctionfactory.h"
#include "voreen/core/animation/interpolation/intinterpolationfunctions.h"

#include "voreen/qt/animation/animationoverviewwidget.h"
#include "voreen/qt/animation/animationtimelinewidget.h"


namespace voreen {

// helper function for sorting the keyframe graphics items
bool compareKeyframeGraphicsItems(KeyframeGraphicsItem* a, KeyframeGraphicsItem* b) {
    return (a->pos().x() < b->pos().x());
}

AnimationTimelineWidget::AnimationTimelineWidget(Animation* animation, AnimationEditor* parent, NetworkEvaluator* networkEval)
        : QWidget(parent)
        , animation_(animation)
        , currentTime_(0)
        , changed_(false)
        , keyframeSelected_(false)
        , currentKfgi_(0)
        , intervalSelected_(false)
        , leftKfgi_(0)
        , rightKfgi_(0)
        , keyframeWidget_(0)
        , intervalWidget_(0)
        , listWidget_(0)
        , currentItemLabel_(0)
        , currentItemTime_(0)
{

    setMinimumWidth(550);

    mainLayout_ = new QVBoxLayout(this);
    QHBoxLayout* controlLayout = new QHBoxLayout();

    timeBox_ = new QGroupBox(tr("Time"), this);
    timeBox_->setStyleSheet("QGroupBox { font-weight: bold; } ");
    QVBoxLayout* timeLayout = new QVBoxLayout();
    timeBox_->setLayout(timeLayout);
    timeBox_->setFixedWidth(timeBox_->width() * 1.5f);

    QLabel* timeLabel = new QLabel(this);
    timeLabel->setText("Current Frame:");
    timeLabel->setFixedWidth(timeBox_->width());
    timeLayout->addWidget(timeLabel);

    timeCounter_ = new QLabel(this);
    timeCounter_->setText("00:00");
    timeCounter_->setFixedWidth(timeBox_->width());
    timeLayout->addWidget(timeCounter_);

    currentItemLabel_ = new QLabel(this);
    currentItemLabel_->setText(tr(" ")); // this is changed according to the type of current item (key / interval)
    currentItemLabel_->setFixedWidth(timeBox_->width());
    timeLayout->addWidget(currentItemLabel_);

    currentItemTime_ = new QLabel(this);
    currentItemTime_->setText(" "); // this is changed according to the time of the current item (key / interval)
    currentItemTime_->setFixedWidth(timeBox_->width());
    timeLayout->addWidget(currentItemTime_);

    mainLayout_->addLayout(controlLayout);
    controlLayout->addWidget(timeBox_);
    //controlLayout->addSpacerItem(new QSpacerItem(20,0));
    overviewTimeline_ = new AnimationOverviewWidget(this, networkEval);
    controlLayout->addWidget(overviewTimeline_);

    connect(this, SIGNAL(viewResizeSignal(int)), overviewTimeline_, SLOT(updateViewportRect(int)));
    connect(this, SIGNAL(autoPreview(bool)), overviewTimeline_, SLOT(autoPreview(bool)));
    connect(this, SIGNAL(updatePreviews()), overviewTimeline_, SLOT(updatePreviews()));
    connect(overviewTimeline_, SIGNAL(currentFrameChanged(int)), this, SIGNAL(currentFrameChanged(int)));
    connect(overviewTimeline_, SIGNAL(intervalSelectedAt(QPointF)), this, SLOT(intervalSelectedAt(QPointF)));
    connect(overviewTimeline_, SIGNAL(clearSelection()), this, SLOT(clearGuiSelection()));

    connect(parent, SIGNAL(currentFrameChanged(int)), overviewTimeline_, SLOT(setCurrentFrame(int)));
    connect(overviewTimeline_, SIGNAL(currentFrameChanged(int)), parent, SLOT(setCurrentFrame(int)));
    connect(this, SIGNAL(scrollBarOrder(int)), overviewTimeline_, SLOT(scrollBarOrder(int)));
    connect(this, SIGNAL(sceneOrder(QMatrix)), overviewTimeline_, SLOT(sceneOrder(QMatrix)));

    connect(overviewTimeline_, SIGNAL(addKeyframe(QPointF)), this, SLOT(addKeyframe(QPointF)));
    connect(overviewTimeline_, SIGNAL(removeKeyframe(KeyframeGraphicsItem*)), this, SLOT(removeKeyframe(KeyframeGraphicsItem*)));
    connect(overviewTimeline_, SIGNAL(itemSelected(KeyframeGraphicsItem*)), this, SLOT(itemClicked(KeyframeGraphicsItem*)));
    connect(overviewTimeline_, SIGNAL(recordToKeyframe(KeyframeGraphicsItem*)), this, SLOT(recordToKeyframe(KeyframeGraphicsItem*)));
    connect(overviewTimeline_, SIGNAL(recordAt(int)), this, SLOT(recordAt(int)));
    connect(this, SIGNAL(playingSignal(bool)), overviewTimeline_, SLOT(playingStateChanged(bool)));

    mainLayout_->setMargin(1);
    mainLayout_->setSpacing(1);
    mainLayout_->setAlignment(Qt::AlignTop);

    scrollArea_ = new QScrollArea();
    containerWidget_ = new QWidget(scrollArea_);
    scrollArea_->setWidget(containerWidget_);
    scrollAreaLayout_ = new QVBoxLayout(containerWidget_);;
    scrollAreaLayout_->setSizeConstraint(QLayout::SetMinimumSize);
    mainLayout_->addWidget(scrollArea_);
    scrollArea_->setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
    scrollArea_->verticalScrollBar()->setContextMenuPolicy(Qt::NoContextMenu);
    scrollArea_->setContextMenuPolicy(Qt::NoContextMenu);
    scrollArea_->setAlignment(Qt::AlignTop);
    scrollAreaLayout_->setAlignment(Qt::AlignTop);
    scrollArea_->setWidgetResizable(true);

    keyframeWidget_ = new AnimationKeyframeWidget(this);
    keyframeWidget_->setVisible(false);
    connect(this, SIGNAL(resizeSignal(int)), keyframeWidget_, SLOT(setFixedWidthSlot(int)));
    connect(keyframeWidget_, SIGNAL(deleteTimelineSignal(PropertyTimeline*)), this, SIGNAL(deleteTimelineSignal(PropertyTimeline*)));

    intervalWidget_ = new AnimationIntervalWidget(this);
    intervalWidget_->setVisible(false);
    connect(this, SIGNAL(resizeSignal(int)), intervalWidget_, SLOT(setFixedWidthSlot(int)));
    connect(intervalWidget_, SIGNAL(deleteTimelineSignal(PropertyTimeline*)), this, SIGNAL(deleteTimelineSignal(PropertyTimeline*)));

    listWidget_ = new AnimationPropertyListWidget(this);
    listWidget_->setVisible(false);
    connect(this, SIGNAL(resizeSignal(int)), listWidget_, SLOT(setFixedWidthSlot(int)));
    connect(listWidget_, SIGNAL(deleteTimelineSignal(PropertyTimeline*)), this, SIGNAL(deleteTimelineSignal(PropertyTimeline*)));

    rebuildAnimation(animation_);
}

void AnimationTimelineWidget::rebuildCurrentAnimation() {
    rebuildAnimation(animation_);
}

void AnimationTimelineWidget::rebuildAnimation(Animation* animation) {

    animation_ = animation;

    //clear property timelines and keyframes
    propertyTimelines_.clear();
    for (std::vector<KeyframeGraphicsItem*>::iterator i = keyframeGraphicsItems_.begin(); i != keyframeGraphicsItems_.end(); ++i) {
        //remove item
        overviewTimeline_->removeKeyframeFromScene(*i);
        delete (*i);
    }
    keyframeGraphicsItems_.clear();

    //clear association map
    keyframeAssociations_.clear();

    if (!animation_)
        return;

    animation_->addObserver(this);

    if ((animation_->getDuration() * 30) != overviewTimeline_->getDuration())
        overviewTimeline_->setDuration(animation->getDuration() * 30);

    //add animated processors and corresponding property timelines
    populateProcessors();

    //build union of keys of all property timelines and build corresponding graphical key value items
    insertAllKeyValuesAsGraphicalItems();

    //insert with empty association vectors into map
    for (std::vector<KeyframeGraphicsItem*>::iterator i = keyframeGraphicsItems_.begin(); i != keyframeGraphicsItems_.end(); ++i) {
        std::vector<std::pair<PropertyKeyValueBase*, PropertyTimeline*> > associations;
        keyframeAssociations_.insert(std::make_pair(*i, associations));
    }

    //put missing key frames into property timelines and put everything into map
    putKeyframesIntoTimelines();

    overviewTimeline_->setStandardRenderport();
    overviewTimeline_->reset();
    //emit setAnimationEditorDuration(static_cast<int>(animation_->getDuration() * 30.0f));

    if (!animation_->isEmpty()) {
        emit showActiveTimelines();
        overviewTimeline_->renderPreviews();
    }

    resetSelectionOfItems();

    keyframeWidget_->clear();
    intervalWidget_->clear();
    listWidget_->clear();
    for (size_t i = 0; i < propertyTimelines_.size(); ++i) {
        keyframeWidget_->addPropertyTimeline(propertyTimelines_.at(i));
        intervalWidget_->addPropertyTimeline(propertyTimelines_.at(i));
        listWidget_->addPropertyTimeline(propertyTimelines_.at(i));
    }

    emit(resizeSignal(scrollArea_->width()));

    updateScrollArea();
}

void AnimationTimelineWidget::clearGuiSelection() {
    resetSelectionOfItems();
    updateScrollArea();
}

void AnimationTimelineWidget::resetSelectionOfItems() {

    if (intervalSelected_) {
        intervalSelected_ = false;
        leftKfgi_ = 0;
        rightKfgi_ = 0;

        overviewTimeline_->removeIntervalFromScene(intervalGraphicsItem_);
        delete intervalGraphicsItem_;
        intervalGraphicsItem_ = 0;
    }
    else if (keyframeSelected_) {
        currentKfgi_ = 0;
        keyframeSelected_ = false;
    }

}

void AnimationTimelineWidget::updateScrollArea() {
    bool keyframePresent = (scrollAreaLayout_->indexOf(keyframeWidget_) != -1) ? true : false;
    bool intervalPresent = (scrollAreaLayout_->indexOf(intervalWidget_) != -1) ? true : false;
    bool listPresent = (scrollAreaLayout_->indexOf(listWidget_) != -1) ? true : false;

    //remove widgets
    if (!keyframeSelected_ && keyframePresent) {
        scrollAreaLayout_->removeWidget(keyframeWidget_);
        keyframeWidget_->setVisible(false);
    }
    else if (!intervalSelected_ && intervalPresent) {
        scrollAreaLayout_->removeWidget(intervalWidget_);
        intervalWidget_->setVisible(false);
    }
    else if ((intervalSelected_ || keyframeSelected_) && listPresent) {
        scrollAreaLayout_->removeWidget(listWidget_);
        listWidget_->setVisible(false);
    }

    if (!keyframeSelected_ && !intervalSelected_) {
        //nothing selected -> set list of animated properties to scroll area
        listWidget_->updateActivationStatus();
        if (!listPresent) {
            scrollAreaLayout_->addWidget(listWidget_);
            listWidget_->show();
        }

        // do not show the current item labels
        currentItemLabel_->setText(" ");
        currentItemTime_->setText(" ");
    }
    else if (keyframeSelected_) {
        keyframeWidget_->updateActivationStatus();
        if (!keyframePresent) {
            scrollAreaLayout_->addWidget(keyframeWidget_);
            keyframeWidget_->show();
        }

        // set text and show the current item time
        currentItemLabel_->setText("Current Keyframe:");
        currentItemTime_->setText(getTimeString(currentKfgi_->pos().x()));
    }
    else if (intervalSelected_) {
        intervalWidget_->updateActivationStatus();
        if (!intervalPresent) {
            scrollAreaLayout_->addWidget(intervalWidget_);
            intervalWidget_->show();
        }

        // set text and show the current item time
        currentItemLabel_->setText("Current Interval:");
        currentItemTime_->setText(getTimeString(leftKfgi_->pos().x()) + " - " + getTimeString(rightKfgi_->pos().x()));
    }

    scrollAreaLayout_->update();
}

bool AnimationTimelineWidget::hasAnimatedProperties() const {
    return (!propertyTimelines_.empty());
}

void AnimationTimelineWidget::removeWidgetsFromLayout(QLayout* layout) {
    QLayoutItem *child;
    while ((child = layout->takeAt(0)) != 0) {
        if (child->layout() != 0)
            removeWidgetsFromLayout(child->layout());
        else if (child->widget() != 0)
            delete child->widget();

        delete child;
    }
}

void AnimationTimelineWidget::insertAllKeyValuesAsGraphicalItems() {
    //for each property timeline
    for (std::vector<PropertyTimeline*>::iterator i = propertyTimelines_.begin(); i != propertyTimelines_.end(); ++i) {
        //get the time values of the keys
        std::vector<float> times = (*i)->getTimesteps();
        //iterate over timesteps and check if they are already present as graphics items, else build a new graphics item
        std::vector<KeyframeGraphicsItem*>::iterator currentItem = keyframeGraphicsItems_.begin();
        for (std::vector<float>::iterator currentTime = times.begin(); currentTime != times.end(); ++currentTime) {
            while ((currentItem != keyframeGraphicsItems_.end()) && ((*currentItem)->pos().x() < ((*currentTime) * 30)))
                ++currentItem;
            if ((currentItem == keyframeGraphicsItems_.end()) || ((*currentItem)->pos().x() > ((*currentTime) * 30))) {
                KeyframeGraphicsItem* kfgi = new KeyframeGraphicsItem(QRectF(-9, -24, 18, 48), static_cast<float>(animation_->getDuration()) * 30.f);

                kfgi->setToolTip(tr("Time: ") + QString::number(*currentTime));

                kfgi->moveBy((*currentTime) * 30, 45);
                kfgi->setZValue(2.1);

                overviewTimeline_->addKeyframeToScene(kfgi);
                connect(kfgi, SIGNAL(itemClicked(KeyframeGraphicsItem*)), this, SLOT(itemClicked(KeyframeGraphicsItem*)));
                connect(kfgi, SIGNAL(itemReleased(KeyframeGraphicsItem*, bool)), this, SLOT(itemReleased(KeyframeGraphicsItem*, bool)));
                connect(kfgi, SIGNAL(itemMoving(KeyframeGraphicsItem*)), this, SLOT(itemMoving(KeyframeGraphicsItem*)));

                currentItem = keyframeGraphicsItems_.insert(currentItem, kfgi);
            }
        }
    }
}

void AnimationTimelineWidget::putKeyframesIntoTimelines() {
    //for each property timeline
    for (std::vector<PropertyTimeline*>::iterator currentTimeline = propertyTimelines_.begin(); currentTimeline != propertyTimelines_.end(); ++currentTimeline) {
        //get key values
        std::map<float, PropertyKeyValueBase*> keyValues = (*currentTimeline)->getKeyValuesAbstract();
        //for each graphical key item
        for (std::vector<KeyframeGraphicsItem*>::iterator currentItem = keyframeGraphicsItems_.begin(); currentItem != keyframeGraphicsItems_.end(); ++currentItem) {
            //look for a key value with the current time in the property timeline
            float time = static_cast<float>((*currentItem)->pos().x()) / 30.f;
            std::map<float, PropertyKeyValueBase*>::iterator findIterator = keyValues.find(time);
            if (findIterator == keyValues.end()) {
                //not found -> put new key value into property timeline
                (*currentTimeline)->setCurrentSettingAsKeyvalue(time, true);
                //get new map and find newly inserted key value
                keyValues = (*currentTimeline)->getKeyValuesAbstract();
                findIterator = keyValues.find(time);
            }
            //put the property key value into the association map
            std::map<KeyframeGraphicsItem*, std::vector<std::pair<PropertyKeyValueBase*, PropertyTimeline*> > >::iterator associationIterator = keyframeAssociations_.find(*currentItem);
            associationIterator->second.push_back(std::make_pair(findIterator->second, *currentTimeline));
        }
    }
}

void AnimationTimelineWidget::populateProcessors() {
    const std::vector<AnimatedProcessor*> procs = animation_->getAnimatedProcessors();

    for(size_t i=0; i<procs.size();i++) {
        if (procs.at(i)->getPropertyTimelines().size() != 0) {
            animatedProcessorAdded(procs.at(i));
        }
    }
}

void AnimationTimelineWidget::reloadAnimation() {
    rebuildAnimation(animation_);
}

void AnimationTimelineWidget::addPropertyTimeline(PropertyTimeline* tl) {
    propertyTimelines_.push_back(tl);

    if (!keyframeGraphicsItems_.empty() && (keyframeGraphicsItems_.at(0)->pos().x() != 0)) {
        //move the only key to the first keyframe to avoid new keys at time = 0
        std::map<float, PropertyKeyValueBase*> keyValues = tl->getKeyValuesAbstract();
        std::map<float, PropertyKeyValueBase*>::iterator zeroItem = keyValues.find(0);
        if (zeroItem != keyValues.end()) {
            tl->setCurrentSettingAsKeyvalue(static_cast<float>(keyframeGraphicsItems_.at(0)->pos().x()) / 30.f, true);
            tl->removeKeyValueAbstract(zeroItem->second);
        }
    }

    keyframeWidget_->addPropertyTimeline(tl);
    intervalWidget_->addPropertyTimeline(tl);
    listWidget_->addPropertyTimeline(tl);

    //rebuild the make sure everything is updated correctly
    //TODO: necessary?!
    rebuildAnimation(animation_);
}

void AnimationTimelineWidget::animatedProcessorAdded(const AnimatedProcessor* processor) {
    AnimatedProcessor* proc = const_cast<AnimatedProcessor*>(processor);    // Maybe we should do this another way, but today we feel quick and dirty

    //push back all property timelines for that animated processor available in simple animation mode so that we can insert and delete keyframes
    Processor* correspondingProcessor = proc->getCorrespondingProcessor();
    for (size_t i = 0; i < proc->getPropertyTimelines().size(); ++i) {
        PropertyTimeline* currentTimeline = proc->getPropertyTimelines().at(i);
        Property* correspondingProperty = currentTimeline->getProperty();
        bool available = dynamic_cast<AnimationEditor*>(parent())->isPropertyAvailable(correspondingProcessor, correspondingProperty);
        if (available)
            if (std::find(propertyTimelines_.begin(), propertyTimelines_.end(), currentTimeline) == propertyTimelines_.end())
                propertyTimelines_.push_back(currentTimeline);
    }
}

void AnimationTimelineWidget::animatedProcessorRemoved(const AnimatedProcessor* processor) {

    //get the property timelines of this processor and remove them one by one
    AnimatedProcessor* proc = const_cast<AnimatedProcessor*>(processor);    // Maybe we should do this another way, but today we feel quick and dirty
    Processor* correspondingProcessor = proc->getCorrespondingProcessor();
    for (size_t i = 0; i < proc->getPropertyTimelines().size(); ++i) {
        PropertyTimeline* currentTimeline = proc->getPropertyTimelines().at(i);

        //remove from:
        // 1. association map
        std::map<KeyframeGraphicsItem*, std::vector<std::pair<PropertyKeyValueBase*, PropertyTimeline*> > >::iterator keyframeIter = keyframeAssociations_.begin();
        while(keyframeIter != keyframeAssociations_.end()) {
            for (std::vector<std::pair<PropertyKeyValueBase*, PropertyTimeline*> >::iterator timelineIter = keyframeIter->second.begin(); timelineIter != keyframeIter->second.end(); ++timelineIter) {
                //if the right timeline has been found: remove the timeline from the vector
                if (currentTimeline == timelineIter->second) {
                    keyframeIter->second.erase(timelineIter);
                    break;
                }
            }
            //check if no associations are left and remove keyframe in that case
            if (keyframeIter->second.empty()) {
                KeyframeGraphicsItem* kfgi = keyframeIter->first;
                overviewTimeline_->removeKeyframeFromScene(kfgi);
                delete kfgi;
                keyframeAssociations_.erase(kfgi);
                std::vector<KeyframeGraphicsItem*>::iterator iter = std::find(keyframeGraphicsItems_.begin(), keyframeGraphicsItems_.end(), kfgi);
                if (iter != keyframeGraphicsItems_.end())
                    keyframeGraphicsItems_.erase(iter);

                keyframeIter = keyframeAssociations_.begin();
            }
            else
                ++keyframeIter;
        }

        // 2. list of property timelines
        std::vector<PropertyTimeline*>::iterator vecIter = std::find(propertyTimelines_.begin(), propertyTimelines_.end(), currentTimeline);
        if (vecIter != propertyTimelines_.end())
            propertyTimelines_.erase(vecIter);

        // 3. scroll area (e.g. list)
        resetSelectionOfItems();
        updateScrollArea();

        // 4. keyframe widget and interval widget (quick and dirty: just clear and rebuild them both)
        keyframeWidget_->clear();
        intervalWidget_->clear();
        listWidget_->clear();
        for (size_t i = 0; i < propertyTimelines_.size(); ++i) {
            keyframeWidget_->addPropertyTimeline(propertyTimelines_.at(i));
            intervalWidget_->addPropertyTimeline(propertyTimelines_.at(i));
            listWidget_->addPropertyTimeline(propertyTimelines_.at(i));
        }
    }
}

bool AnimationTimelineWidget::keysLaterThan(int time) const {
    if (keyframeGraphicsItems_.empty())
        return false;
    else return (keyframeGraphicsItems_.at(keyframeGraphicsItems_.size() - 1)->pos().x() > time);
}

void AnimationTimelineWidget::durationChanged(int duration) {
    overviewTimeline_->setDuration(duration);
    rebuildAnimation(animation_);
}

void AnimationTimelineWidget::removeKeysAboveDuration(int duration) {

    //TODO: this assumes that the items are not ordered. If this has been changed, revise this method, as it is extremely inefficient
    bool finished = false;
    while(!finished) {
        for (std::vector<KeyframeGraphicsItem*>::iterator i = keyframeGraphicsItems_.begin(); ; ++i) {
            if (i == keyframeGraphicsItems_.end()) {
                finished = true;
                break;
            }
            else if ((*i)->pos().x() < duration)
                continue;
            else {
                //not at end and the item should be removed
                overviewTimeline_->removeKeyframeFromScene(*i);
                delete (*i);
                //delete from association map ->
                //note that we do not have to remove the key values from the property timelines as the animation does this automatically
                keyframeAssociations_.erase(*i);

                keyframeGraphicsItems_.erase(i);
                break;
            }
        }
    }
}


void AnimationTimelineWidget::keyPressEvent(QKeyEvent* e) {
    if(e->key() == Qt::Key_Delete) {
        //if a keyframe had been selected before: remove it in core and gui
        if (keyframeSelected_) {
            removeKeyframe(currentKfgi_);
        }
        else if (intervalSelected_) {
            //interval: reset the selection
            resetSelectionOfItems();
            updateScrollArea();
        }
    }
}

void AnimationTimelineWidget::removeKeyframe(KeyframeGraphicsItem* kfgi) {
    //delete the keyframes in the core
    if (keyframeAssociations_.find(kfgi) == keyframeAssociations_.end()) {
        LERRORC("voreen.qt.AnimationTimelineWidget", "KeyframeGraphicsItem to be deleted not found!");
        return;
    }
    std::vector<std::pair<PropertyKeyValueBase*, PropertyTimeline*> >& assocs = keyframeAssociations_[kfgi];
    for (size_t i = 0; i < assocs.size(); ++i) {
        DeleteKeyValueReturn dkvr = assocs.at(i).second->removeKeyValueAbstract(assocs.at(i).first);
        if (dkvr == KV_IS_THE_ONLY_ONE) {
            if (i == 0) {
                LINFOC("voreen.qt.AnimationTimelineWidget", "Removing the only keyframe from animation is not possible.");
                return;
            }
            else {
                LERRORC("voreen.qt.AnimationTimelineWidget",
                        "Could not delete keyframe correctly, number of keyframes differs in property timelines. Rebuilding Animation.");
                rebuildAnimation(animation_);
                return;
            }
        }
    }

    //delete keyframe graphics item
    if (keyframeSelected_ && kfgi == currentKfgi_) {
        resetSelectionOfItems();
        updateScrollArea();
    }

    overviewTimeline_->removeKeyframeFromScene(kfgi);
    delete kfgi;
    keyframeAssociations_.erase(kfgi);
    std::vector<KeyframeGraphicsItem*>::iterator iter = std::find(keyframeGraphicsItems_.begin(), keyframeGraphicsItems_.end(), kfgi);
    if (iter != keyframeGraphicsItems_.end())
        keyframeGraphicsItems_.erase(iter);
}

KeyframeGraphicsItem* AnimationTimelineWidget::addKeyframe(QPointF pos) {    //this adds also to the scene

    // check if there are actual timelines in the animation
    if (propertyTimelines_.empty()) {
        QString text;
        text.append(tr("No properties are animated - adding a keyframe is not possible."));
        QMessageBox::information(this, tr("No properties animated"), text);

        return 0;
    }

    // check if the key is inside the time frame (i.e. <= duration)
    if (static_cast<float>(pos.x()) / 30.f > animation_->getDuration()) {
        LWARNINGC("voreen.qt.AnimationTimelineWidget", "Cannot insert keyframe outside of the animation timeframe");
        return 0;
    }

    KeyframeGraphicsItem* kfgi = new KeyframeGraphicsItem(QRectF(-9, -24, 18, 48), static_cast<float>(animation_->getDuration()) * 30.f);
    //insert in correct place
    std::vector<KeyframeGraphicsItem*>::iterator graphicsIterator = keyframeGraphicsItems_.begin();
    while((graphicsIterator != keyframeGraphicsItems_.end()) && ((*graphicsIterator)->pos().x() < pos.x()))
       ++graphicsIterator;
    keyframeGraphicsItems_.insert(graphicsIterator, kfgi);

    kfgi->moveBy(pos.x(), 45);
    kfgi->setZValue(2.1);

    kfgi->setToolTip(tr("Time: ") + QString::number(pos.x() / 30.f));
    //propertyTimelineScene_->addItem(kfgi);
    overviewTimeline_->addKeyframeToScene(kfgi);
    connect(kfgi, SIGNAL(itemClicked(KeyframeGraphicsItem*)), this, SLOT(itemClicked(KeyframeGraphicsItem*)));
    connect(kfgi, SIGNAL(itemReleased(KeyframeGraphicsItem*, bool)), this, SLOT(itemReleased(KeyframeGraphicsItem*, bool)));
    connect(kfgi, SIGNAL(itemMoving(KeyframeGraphicsItem*)), this, SLOT(itemMoving(KeyframeGraphicsItem*)));

    //put entry into association map
    std::vector<std::pair<PropertyKeyValueBase*, PropertyTimeline*> > associations;
    keyframeAssociations_.insert(std::make_pair(kfgi, associations));

    //put key values into property timelines, add to association map
    for (size_t i = 0; i < propertyTimelines_.size(); ++i) {
        float time = pos.x() / 30.f;
        PropertyKeyValueBase* newKeyValue = const_cast<PropertyKeyValueBase*>(propertyTimelines_.at(i)->newKeyValueAbstract(time));
        if (!newKeyValue)
            LERRORC("voreen.qt.AnimationTimelineWidget", "Could not insert key value at time " << time << " into timeline: " << propertyTimelines_.at(i)->getProperty()->getGuiName());
        else {
            //put the property key value into the association map
            std::map<KeyframeGraphicsItem*, std::vector<std::pair<PropertyKeyValueBase*, PropertyTimeline*> > >::iterator associationIterator = keyframeAssociations_.find(kfgi);
            associationIterator->second.push_back(std::make_pair(newKeyValue, propertyTimelines_.at(i)));
        }
    }

    return kfgi;
}

void AnimationTimelineWidget::recordAt(int frame) {

    if (hasAnimatedProperties()) {

        // check if the key is inside the time frame (i.e. <= duration)
        if (static_cast<float>(frame) / 30.f  > animation_->getDuration()) {
            LWARNINGC("voreen.qt.AnimationTimelineWidget", "Cannot insert keyframe outside of the animation timeframe");
            return;
        }

        animation_->setActualNetworkAsKeyvalues(static_cast<float>(frame)/30.f, true);

        //force rebuilding the gui components
        rebuildAnimation(animation_);
    }
    else {
        // notify user
        QString text;
        text.append(tr("No properties are animated - adding a keyframe is not possible."));
        QMessageBox::information(this, tr("No properties animated"), text);
    }
}

void AnimationTimelineWidget::recordToKeyframe(KeyframeGraphicsItem* kfgi) {

    std::map<KeyframeGraphicsItem*, std::vector<std::pair<PropertyKeyValueBase*, PropertyTimeline*> > >::iterator associationIterator = keyframeAssociations_.find(kfgi);
    if (associationIterator != keyframeAssociations_.end()) {
        for (size_t i = 0; i < associationIterator->second.size(); ++i) {
            associationIterator->second.at(i).second->setCurrentSettingToKeyValue(associationIterator->second.at(i).first);
        }
        keyframeWidget_->updateProperties(associationIterator->second);
        updateScrollArea();
    }
    else {
        LERRORC("voreen.qt.AnimationTimelineWidget", "KeyframeGraphicsItem not found!");
    }
}

void AnimationTimelineWidget::intervalSelectedAt(QPointF pos) {
    if (keyframeSelected_) {
        keyframeSelected_ = false;
        currentKfgi_ = 0;
    }
    else if (intervalSelected_) {
        //remove graphical item if an interval has been selected before
        overviewTimeline_->removeIntervalFromScene(intervalGraphicsItem_);
        delete intervalGraphicsItem_;
        intervalGraphicsItem_ = 0;
        intervalSelected_ = false;
        leftKfgi_ = 0;
        rightKfgi_ = 0;
    }

    //find keyframes left and right (if existing)
    KeyframeGraphicsItem* currentItem = 0;
    KeyframeGraphicsItem* lastItem = 0;
    for (size_t i = 0; i < keyframeGraphicsItems_.size(); ++i) {
        lastItem = currentItem;
        currentItem = keyframeGraphicsItems_.at(i);

        if (currentItem->pos().x() > pos.x()) {
            //found first key right to the position
            leftKfgi_ = lastItem;
            rightKfgi_ = currentItem;
            break;
        }
        else if (i == (keyframeGraphicsItems_.size()-1)) {
            //position is right to the last key
            leftKfgi_ = currentItem;
            rightKfgi_ = 0;
        }
    }

    if (leftKfgi_ == 0 || rightKfgi_ == 0)
        intervalSelected_ = false;
    else {
        intervalSelected_ = true;

        //create graphics item for this interval
        float width = rightKfgi_->pos().x() - leftKfgi_->pos().x();
        float xPos = (rightKfgi_->pos().x() + leftKfgi_->pos().x()) / 2;
        intervalGraphicsItem_ = new QGraphicsRectItem(QRectF(-width/2 + 1, -23, width - 1, 47));
        intervalGraphicsItem_->moveBy(xPos, 45);
        intervalGraphicsItem_->setZValue(2.2);
        intervalGraphicsItem_->setOpacity(0.5);

        QLinearGradient gradient(0,0,0, 10);
        gradient.setSpread(QGradient::PadSpread);

        gradient.setColorAt(0.0, QColor(155, 255, 155, 255));
        gradient.setColorAt(1.0, QColor(25, 255, 25, 255));

        QBrush brush(gradient);
        intervalGraphicsItem_->setBrush(gradient);
        intervalGraphicsItem_->setVisible(true);

        overviewTimeline_->addIntervalToScene(intervalGraphicsItem_);

        intervalWidget_->updateProperties(keyframeAssociations_[leftKfgi_], keyframeAssociations_[rightKfgi_]);
    }

    // display information at scrollable area
    updateScrollArea();
}

void AnimationTimelineWidget::itemClicked(KeyframeGraphicsItem* kfgi) {

    KeyframeGraphicsItem* sKfgi = dynamic_cast<KeyframeGraphicsItem*>(kfgi);
    if (!sKfgi) {
        LERRORC("voreen.qt.AnimationTimelineWidget", "not a valid keyframe");
        return;
    }

    if (intervalSelected_) {
        overviewTimeline_->removeIntervalFromScene(intervalGraphicsItem_);
        delete intervalGraphicsItem_;
        intervalGraphicsItem_ = 0;
        intervalSelected_ = false;
        leftKfgi_ = 0;
        rightKfgi_ = 0;
    }
    else if (keyframeSelected_ && sKfgi == currentKfgi_)
        return;

    keyframeSelected_ = true;

    leftKfgi_ = 0;
    rightKfgi_ = 0;
    currentKfgi_ = sKfgi;

    keyframeWidget_->updateProperties(keyframeAssociations_[sKfgi]);

    // display information at scrollable area
    updateScrollArea();
}

void AnimationTimelineWidget::applyTimeStretch(int nDuration) {
    if (nDuration > (animation_->getDuration() * 30)) {

        //compute the new stretch ratio
        float stretch = static_cast<float>(nDuration) / (static_cast<float>(animation_->getDuration()) * 30.f);

        //first set new duration to have enough space for moving the keys
        animation_->setDuration(static_cast<float>(nDuration) / 30.f);

        //property timelines in core copy the keyframes -> pointers are invalid -> rebuild animation to fix this
        rebuildAnimation(animation_);

        //for all keyframegraphicsitems: shift keys in all property timelines
        std::map<KeyframeGraphicsItem*, std::vector<std::pair<PropertyKeyValueBase*, PropertyTimeline*> > >::iterator kfgiIter = keyframeAssociations_.begin();
        for (; kfgiIter != keyframeAssociations_.end(); ++kfgiIter) {
            float time = kfgiIter->second.at(0).first->getTime() * stretch;
            for (size_t i = 0; i < kfgiIter->second.size(); ++i)
                kfgiIter->second.at(i).second->changeTimeOfKeyValueAbstract(time, kfgiIter->second.at(i).first);
            //kfgiIter->first->setX(time * 30.f);
        }

        rebuildAnimation(animation_);
    }
    else if (nDuration < (animation_->getDuration() * 30)) {

        //compute the new stretch ratio
        float stretch = static_cast<float>(nDuration) / (static_cast<float>(animation_->getDuration()) * 30.f);

        //first change the time of all key values
        std::map<KeyframeGraphicsItem*, std::vector<std::pair<PropertyKeyValueBase*, PropertyTimeline*> > >::iterator kfgiIter = keyframeAssociations_.begin();
        for (; kfgiIter != keyframeAssociations_.end(); ++kfgiIter) {
            float time = kfgiIter->second.at(0).first->getTime() * stretch;
            for (size_t i = 0; i < kfgiIter->second.size(); ++i)
                kfgiIter->second.at(i).second->changeTimeOfKeyValueAbstract(time, kfgiIter->second.at(i).first);
            //kfgiIter->first->setX(time * 30.f);
        }

        //set new duration afterwards so that no keys get lost
        animation_->setDuration(static_cast<float>(nDuration) / 30.f);

        //rebuild animation
        rebuildAnimation(animation_);
    }
}

void AnimationTimelineWidget::itemReleased(KeyframeGraphicsItem* kfgi, bool shift) {

    KeyframeGraphicsItem* sKfgi = dynamic_cast<KeyframeGraphicsItem*>(kfgi);
    if (!sKfgi) {
        LERRORC("voreen.qt.AnimationTimelineWidget", "not a valid keyframe");
        return;
    }

    // if the item has not been moved: do nothing
    if (!sKfgi->moving())
        return;

    //else: compute new position of the keyframes in the animation core
    float newTime = static_cast<float>(kfgi->pos().x()) / 30.f;

    //check if there is another key at this position in the core
    //there is only need to check one timeline since the time values are the same in all property timelines
    if (!propertyTimelines_.empty()) {
        std::vector<float> timesteps = propertyTimelines_.at(0)->getTimesteps();
        float ownTime = keyframeAssociations_[sKfgi].at(0).first->getTime();
        if ((ownTime != newTime) && (std::find(timesteps.begin(), timesteps.end(), newTime) != timesteps.end())) {
            LWARNINGC("voreen.qt.AnimationTimelineWidget", "Positioning two keyvalues at the same time is not allowed");
            sKfgi->restoreOldPosition();
            return;
        }
    }
    else {
        //no property timelines to check -> there should be no keys to move
        LERRORC("voreen.qt.AnimationTimelineWidget", "No property timelines found!");
        return;
    }

    //get the map entry and shift keys in all property timelines
    std::vector<std::pair<PropertyKeyValueBase*, PropertyTimeline*> > keys = keyframeAssociations_[sKfgi];
    for (size_t i = 0; i < keys.size(); ++i) {
        if (shift) {
            keys.at(i).second->shiftKeyValueAbstract(newTime, keys.at(i).first);
            //also change position of following keyframes
            std::vector<KeyframeGraphicsItem*>::iterator iter = std::find(keyframeGraphicsItems_.begin(), keyframeGraphicsItems_.end(), sKfgi);
            for ( ; iter != keyframeGraphicsItems_.end(); ++iter) {
                float coreTime = keyframeAssociations_[*iter].at(0).first->getTime();
                (*iter)->setPos(coreTime * 30.f, (*iter)->pos().y());
                (*iter)->setToolTip(tr("Time: ") + QString::number(coreTime));
            }
        }
        else
            keys.at(i).second->changeTimeOfKeyValueAbstract(newTime, keys.at(i).first);
    }

    //change the tool tip of the graphics item that was actually moved
    kfgi->setToolTip(tr("Time: ") + QString::number(newTime));

    //sort vector of graphics items, because the order might have changed
    std::sort(keyframeGraphicsItems_.begin(), keyframeGraphicsItems_.end(), compareKeyframeGraphicsItems);
}

void AnimationTimelineWidget::itemMoving(KeyframeGraphicsItem* kfgi) {
    //std::cout << "Item moving" << std::endl;
    /*moveFrameHUD(kfgi->pos().x()-40);
    currentFrameCounter_->setHtml(getTimeString(kfgi->pos().x()));*/

    // change time string of current key
    if (currentKfgi_ == kfgi) {
        currentItemTime_->setText(getTimeString(kfgi->pos().x()));
    }

}

float AnimationTimelineWidget::getCurrentTime() {
    return currentTime_;
}

void AnimationTimelineWidget::setCurrentTime(float time) {
    currentTime_ = time;
}

void AnimationTimelineWidget::renderAt(float time) {
        animation_->renderAt(time);
}

void AnimationTimelineWidget::currentFrame(int frame) {
    timeCounter_->setText(getTimeString(frame));
}

QString AnimationTimelineWidget::getTimeString(int frame) {
    int min = 0;
    //int sec = frame / 30;
    float sec = static_cast<float>(frame) / 30.f;

    while(sec >= 60.f) {
        sec -= 60.f;
        min++;
    }
    QString minutes;
    QString seconds;
    if(min >= 10)
        minutes = QString::number(min);
    else
        minutes = "0"+QString::number(min);
    if(sec >= 10)
        seconds = QString::number(sec, 'f', 4);
    else
        seconds = "0"+QString::number(sec, 'f', 4);

    QString ret = minutes+":"+seconds;
    return ret;
}

void AnimationTimelineWidget::resizeEvent(QResizeEvent* event) {
    emit(resizeSignal(scrollArea_->width()));
    QWidget::resizeEvent(event);
}

void AnimationTimelineWidget::setApplicationModeConfig(ApplicationModeConfiguration* appConfig) {
    keyframeWidget_->setApplicationModeConfig(appConfig);
    intervalWidget_->setApplicationModeConfig(appConfig);
    listWidget_->setApplicationModeConfig(appConfig);
}

} // namespace voreen

