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

#include "voreen/qt/animation/animationeditor.h"

#include "voreen/qt/animation/animationexportwidget.h"
#include "voreen/qt/animation/animationtimelinewidget.h"
#include "voreen/qt/widgets/sliderspinboxwidget.h"

#include "voreen/core/animation/animatedprocessor.h"
#include "voreen/core/animation/serializationfactories.h"

#include <QCheckBox>
#include <QVBoxLayout>
#include <QLabel>
#include <QLCDNumber>
#include <QMenuBar>
#include <QMenu>
#include <QPushButton>
#include <QProgressBar>
#include <QTimeLine>
#include <QSlider>
#include <QTimeLine>
#include <QSlider>
#include <QKeyEvent>
#include <iostream>
#include <QGraphicsView>
#include <QGraphicsScene>
#include <QScrollArea>
#include <QTimer>
#include <QToolBar>
#include <QToolButton>
#include <QMessageBox>
#include <QString>

#include <QImageReader>

namespace voreen {

AnimationEditor::AnimationEditor(NetworkEvaluator* eval, Workspace* workspace, QWidget* parent, Qt::WindowFlags flags)
    : QWidget(parent, flags | Qt::Tool | Qt::Window)
    , WorkspaceObserver()
    , UsesWorkspace()
    , evaluator_(eval)
    , mainMenu_(new QMenuBar(this))
    , currentFrame_(0)
    , frameSkip_(1)
    , timeStretch_(1.0f)
    , interactionMode_(false)
    , workspace_(workspace)
    , applyTimeStretch_(false)
    , animation_(0)
{

    setFocusPolicy(Qt::StrongFocus);

    if (workspace->getAnimation()) {
        animation_ = workspace->getAnimation();
        duration_ = (int)animation_->getDuration()*30.0f;
    } else {
        animation_ = new Animation(workspace->getProcessorNetwork());
        animation_->setDuration(20);
        duration_ = 600;
        workspace->setAnimation(animation_);
    }

    QVBoxLayout* mainLayout = new QVBoxLayout(this);
    mainLayout->setMargin(0);
    mainLayout->setSpacing(0);

    QToolBar* toolbar = new QToolBar(this);
    QToolBar* controlToolbar = new QToolBar(toolbar);
    toolbar->setMovable(true);
    controlToolbar->setMovable(true);

    toolbar->setStyleSheet("* {background: qlineargradient(x1:0 y1:0, x2:0 y2:1, stop:0 rgba( 155, 155, 155, 255 ), stop:1 rgba( 225, 225, 225, 255 ) );}");

    QToolButton* newButton = new QToolButton(this);
    QToolButton* settingsButton = new QToolButton(this);
    QToolButton* videoExportButton = new QToolButton(this);
    newButton->setIcon(QIcon(":/qt/icons/clear.png"));
    settingsButton->setIcon(QIcon(":/qt/icons/properties.png"));
    videoExportButton->setIcon(QIcon(":/qt/icons/video_export.png"));

    //QToolButton* undoButton = new QToolButton(this);
    //QToolButton* redoButton = new QToolButton(this);
    //undoButton->setIcon(QIcon(":/qt/icons/edit_undo.png"));
    //redoButton->setIcon(QIcon(":/qt/icons/edit_redo.png"));

    //TODO: re-enable undo / redo functionality if the simple editor is functional
    //undoButton->setEnabled(false);
    //redoButton->setEnabled(false);

    QToolButton* startButton = new QToolButton(this);

    playerControlGroup_ = new QActionGroup(this);
    QAction* rewind = new QAction(QIcon(":/qt/icons/player_rew.png"), tr("rewind"), this);
    QAction* play = new QAction(QIcon(":/qt/icons/player_play.png"), tr("play"), this);
    QAction* pause = new QAction(QIcon(":/qt/icons/player_pause.png"), tr("pause"), this);
    QAction* stop = new QAction(QIcon(":/qt/icons/player_stop.png"), tr("stop"), this);
    QAction* forward = new QAction(QIcon(":/qt/icons/player_fwd.png"), tr("forward"), this);
    playerControlGroup_->addAction(rewind);
    playerControlGroup_->addAction(play);
    playerControlGroup_->addAction(pause);
    playerControlGroup_->addAction(stop);
    playerControlGroup_->addAction(forward);
    rewind->setCheckable(true);
    play->setCheckable(true);
    pause->setCheckable(true);
    stop->setCheckable(true);
    forward->setCheckable(true);
    connect(playerControlGroup_, SIGNAL(triggered(QAction*)), this, SLOT(playerControl(QAction*)));
    QToolButton* endButton = new QToolButton(this);
    //QToolButton* recordButton = new QToolButton(this);
    addTimelineButton_ = new QToolButton(this);

    QAction* interaction = new QAction("Use Interaction Mode", this);
    interaction->setCheckable(true);
    connect(interaction, SIGNAL(toggled(bool)), this, SLOT(setInteractionMode(bool)));
    QAction* autoPreview = new QAction("Render Previews", this);
    autoPreview->setCheckable(true);
    connect(autoPreview, SIGNAL(toggled(bool)), this, SIGNAL(autoPreview(bool)));
    QAction* updatePreviews = new QAction("Update Previews", this);

    addTimelineMenu_ = new QMenu();
    connect(addTimelineMenu_, SIGNAL(aboutToShow()), this, SLOT(populateAddTimelineMenu()));
    connect(addTimelineMenu_ , SIGNAL(aboutToShow()), this, SLOT(addButtonClicked()));
    connect(addTimelineMenu_, SIGNAL(triggered(QAction*)), this, SLOT(addTimeline(QAction*)));

    addTimelineButton_->setMenu(addTimelineMenu_);
    addTimelineButton_->setPopupMode(QToolButton::InstantPopup);
    //addTimelineButton_->setAutoRaise(true);

    startButton->setIcon(QIcon(":/qt/icons/player_start.png"));
    endButton->setIcon(QIcon(":/qt/icons/player_end.png"));
    //recordButton->setIcon(QIcon(":/qt/icons/player_record.png"));
    addTimelineButton_->setIcon(QIcon(":/qt/icons/edit_add.png"));

    toolbar->addWidget(newButton);
    toolbar->addWidget(settingsButton);
    toolbar->addWidget(videoExportButton);
    newButton->setToolTip(tr("New Animation"));
    settingsButton->setToolTip(tr("Animation Settings"));
    videoExportButton->setToolTip(tr("Export Animation"));

    toolbar->addSeparator();

    //toolbar->addWidget(undoButton);
    //toolbar->addWidget(redoButton);
    //TODO: set proper tool tips as soon as undo and redo are available again
    //undoButton->setToolTip(tr("Undo (Currently not available)"));
    //redoButton->setToolTip(tr("Redo (currently not available)"));

    //TODO: check if this button is necessary
    /*QToolButton* rebuildButton = new QToolButton(this);
    rebuildButton->setIcon(QIcon(":/qt/icons/edit_undo.png"));
    toolbar->addWidget(rebuildButton);
    rebuildButton->setToolTip(tr("Rebuild Animation (use if Application Mode Configuration was changed"));*/

    toolbar->addSeparator();

    controlToolbar->addWidget(startButton);
    controlToolbar->addAction(rewind);
    controlToolbar->addAction(pause);
    controlToolbar->addAction(play);
    controlToolbar->addAction(stop);
    controlToolbar->addAction(forward);
    controlToolbar->addWidget(endButton);

    startButton->setToolTip(tr("Start"));
    rewind->setToolTip(tr("Rewind"));
    pause->setToolTip(tr("Pause"));
    play->setToolTip(tr("Play"));
    stop->setToolTip(tr("Stop"));
    forward->setToolTip(tr("Forward"));
    endButton->setToolTip(tr("End"));
    //recordButton->setToolTip(tr("Take Snapshot"));
    addTimelineButton_->setToolTip(tr("Add Timeline"));

    toolbar->addWidget(controlToolbar);

    toolbar->addSeparator();
    //toolbar->addWidget(recordButton);
    toolbar->addWidget(addTimelineButton_);

    toolbar->addSeparator();
    toolbar->addAction(interaction);
    toolbar->addAction(autoPreview);
    toolbar->addAction(updatePreviews);

    mainLayout->addWidget(toolbar);

    timelineWidget_ = new AnimationTimelineWidget(animation_, this, evaluator_);
    mainLayout->addWidget(timelineWidget_);

    connect(this, SIGNAL(currentFrameChanged(int)), timelineWidget_, SLOT(currentFrame(int)));
    connect(this, SIGNAL(currentFrameChanged(int)), timelineWidget_, SIGNAL(currentFrameChanged(int)));
    connect(this, SIGNAL(durationChanged(int)), timelineWidget_, SLOT(durationChanged(int)));
    connect(this, SIGNAL(autoPreview(bool)), timelineWidget_, SIGNAL(autoPreview(bool)));
    connect(updatePreviews, SIGNAL(triggered(bool)), timelineWidget_, SIGNAL(updatePreviews()));

    connect(this, SIGNAL(newAnimation(Animation*)), timelineWidget_, SLOT(rebuildAnimation(Animation*)));
    connect(this, SIGNAL(playingSignal(bool)), timelineWidget_, SIGNAL(playingSignal(bool)));
    //connect(timelineWidget_, SIGNAL(setAnimationEditorDuration(int)), this, SLOT(setDuration(int)));
    connect(timelineWidget_, SIGNAL(deleteTimelineSignal(PropertyTimeline*)), this, SLOT(removeTimelineFromAnimation(PropertyTimeline*)));


    connect(newButton, SIGNAL(clicked()), this, SLOT(newAnimation()));
    connect(settingsButton, SIGNAL(clicked()), this, SLOT(settings()));
    connect(videoExportButton, SIGNAL(clicked()), this, SLOT(videoExport()));

    connect(startButton, SIGNAL(clicked()), this, SLOT(start()));
    connect(endButton, SIGNAL(clicked()), this, SLOT(end()));
    //connect(recordButton, SIGNAL(clicked()), this, SLOT(record()));
    //connect(this, SIGNAL(recordSignal()), timelineWidget_, SIGNAL(recordSignal()));
    //connect(undoButton, SIGNAL(clicked()), this, SLOT(undo()));
    //connect(redoButton, SIGNAL(clicked()), this, SLOT(redo()));

    //connect(rebuildButton, SIGNAL(clicked()), timelineWidget_, SLOT(rebuildCurrentAnimation()));

    //connect(this, SIGNAL(deactivatedWindowSignal()), timelineWidget_, SLOT(clearGuiSelection()));

    init();

    workspace_->addObserver(this);
}

void AnimationEditor::addButtonClicked() {
    if (groups_.empty() || addTimelineButton_->menu()->isEmpty()) {
        QString text;
        if (groups_.empty()) {
            //No properties available in application mode -> message box
            text.append(tr("There are no properties visible in application mode.\n"));
            text.append(tr("Please configure property visibility first.\n"));
        }
        else {
            //all properties already animated
            text.append(tr("All available properties have already been selected.\n"));
            text.append(tr("To add more properties, use the property visibility configuration.\n"));
        }
        /*int answer = */QMessageBox::information(this, tr("No properties available"), text);
        //addTimelineButton_->menu()->close();
    }
}

void AnimationEditor::analyzeNetwork() {

    groups_.clear();

    //TODO: revise this method
    if (!evaluator_->getProcessorNetwork())
        return;

    const std::vector<Processor*> processors = evaluator_->getProcessorNetwork()->getProcessors();

    //1. find properties to animate
    std::vector<std::pair<Processor*, Property*> > pList;
    for (std::vector<Processor*>::const_iterator i = processors.begin(); i != processors.end(); ++i) {
        const std::vector<Property*>& properties = (*i)->getProperties();
        for (std::vector<Property*>::const_iterator p = properties.begin(); p != properties.end(); ++p) {

            //only properties that are visible in application mode may be animated
            if (workspace_->getApplicationModeConfig().isPropertyVisible(*p))
                pList.push_back(std::make_pair(*i, *p));
        }
    }

    //2. find groups of linked properties

    while (!pList.empty()) {
        Processor* currentProcessor = pList.begin()->first;
        Property* currentProperty = pList.begin()->second;
        pList.erase(pList.begin());

        std::deque<std::pair<Processor*, Property*> > currentSet;
        currentSet.push_back(std::make_pair(currentProcessor, currentProperty));

        //now: find linked properties
        size_t currentIndex = 0;
        while(currentIndex < currentSet.size()) {
            size_t checkIndex = 0;
            while(!pList.empty() && (checkIndex < pList.size())) {
                if (currentSet.at(currentIndex).second->isLinkedWith(pList.at(checkIndex).second, true)) {
                    currentSet.push_back(pList.at(checkIndex));
                    pList.erase(pList.begin() + checkIndex);
                }
                else if (pList.at(checkIndex).second->isLinkedWith(currentSet.at(currentIndex).second, true)) {
                    //use element at checkIndex as master property
                    currentSet.push_front(pList.at(checkIndex));
                    pList.erase(pList.begin() + checkIndex);
                }
                else
                    checkIndex++;
            }
            currentIndex++;
        }

        groups_.push_back(currentSet);
    }
}

bool AnimationEditor::isPropertyAvailable(const Processor* processor, const Property* property) const {

    //iterate over our property groups
    bool found = false;
    for (std::vector<std::deque<std::pair<Processor*, Property*> > >::const_iterator group = groups_.begin(); group != groups_.end(); ++group) {
        //iterate over pairs in this group and look for the processor / property combination
        for (std::deque<std::pair<Processor*, Property*> >::const_iterator currentPair = group->begin(); currentPair != group->end(); ++currentPair) {
            if (currentPair->first == processor && currentPair->second == property) {
                found = true;
                break;
            }
        }
        //already found: no need to check other groups
        if (found)
        break;
    }

    return found;

}

bool AnimationEditor::isAnimationValid(Animation* animation, std::vector<std::pair<const Processor*, PropertyTimeline*> >& invalidPropertyTimelines) const {
    bool valid = true;

    // this vector should contain timelines that have to be tested later for group membership
    std::deque<std::pair<const Processor*, PropertyTimeline*> > timelinesToTest;

    //iterate over processors
    const std::vector<AnimatedProcessor*>& animatedProcessors = animation->getAnimatedProcessors();
    for (std::vector<AnimatedProcessor*>::const_iterator i = animatedProcessors.begin(); i != animatedProcessors.end(); ++i) {
        //get processor and property timelines
        const Processor* currentProcessor = (*i)->getCorrespondingProcessor();
        const std::vector<PropertyTimeline*>& propertyTimelines = (*i)->getPropertyTimelines();

        //iterate over property timelines
        for (std::vector<PropertyTimeline*>::const_iterator j = propertyTimelines.begin(); j != propertyTimelines.end(); ++j) {
            const Property* currentProperty = (*j)->getProperty();
            //if property is animated we have to check if this is allowed in the simple animation
            if ((*i)->isPropertyAnimated(currentProperty) && (*j)->getActiveOnRendering()) {
                bool found = isPropertyAvailable(currentProcessor, currentProperty);
                if (!found) {
                    valid = false;
                    invalidPropertyTimelines.push_back(std::make_pair(currentProcessor, *j));
                }
                else    //these property timelines have to be tested since some may be in the same group of linked properties
                    timelinesToTest.push_back(std::make_pair(currentProcessor, *j));
            }
        }
    }

    //now test properties for group membership -> only one property of each group may be animated
    while (!timelinesToTest.empty()) {
        //take the first property timeline and find its group (both exist since this has been tested before)
        std::pair<const Processor*, PropertyTimeline*> currentProperty = timelinesToTest.at(0);
        timelinesToTest.pop_front();
        std::deque<std::pair<Processor*,Property*> > currentGroup = groups_.at(0);
        size_t currentPosition = 0;
        bool foundGroup = false;
        for (size_t i = 0; i < groups_.size(); ++i) {
            currentGroup = groups_.at(i);
            for (size_t j = 0; j < currentGroup.size(); ++j) {
                if ((currentGroup.at(j).first == currentProperty.first) && (currentGroup.at(j).second == currentProperty.second->getProperty())) {
                    foundGroup = true;
                    currentPosition = j;
                    break;
                }
            }
            if (foundGroup)
                break;
        }
        //now check for all properties if they are in the same group and only keep the one which is the first of them in the group
        std::deque<std::pair<const Processor*, PropertyTimeline*> >::iterator it = timelinesToTest.begin();
        while (it != timelinesToTest.end()) {
            size_t testPosition;
            bool inGroup = false;
            for (size_t j = 0; j < currentGroup.size(); ++j) {
                if ((currentGroup.at(j).first == it->first) && (currentGroup.at(j).second == it->second->getProperty())) {
                    inGroup = true;
                    testPosition = j;
                    break;
                }
            }
            if (inGroup) {
                //this property obviously is in the same group. Depending on the position, mark one of the two as non-valid
                if (currentPosition < testPosition) {
                    // the current reference property is first, so the other one has to be removed
                    invalidPropertyTimelines.push_back(std::make_pair(it->first, it->second));
                    timelinesToTest.erase(it);
                }
                else {
                    // the current reference property has to be removed, the test property becomes the new reference
                    invalidPropertyTimelines.push_back(std::make_pair(currentProperty.first, currentProperty.second));
                    currentPosition = testPosition;
                    currentProperty = *it;
                    timelinesToTest.erase(timelinesToTest.begin() + testPosition);
                }
                it = timelinesToTest.begin();
                valid = false;
            }
            else
                ++it;
        }
    }

    return valid;
}

void AnimationEditor::applicationModeConfigurationChanged() {

    stop();

    if (!evaluator_->getProcessorNetwork()) {
        timelineWidget_->setApplicationModeConfig(0);
        return;
    }

    analyzeNetwork();
    //check if the animation is valid under the changed config and if not, display a message box to the user
    std::vector<std::pair<const Processor*, PropertyTimeline*> > invalidTimelines;

    if (!isAnimationValid(animation_, invalidTimelines)) {
        QString text(tr("The following animated properties are not available anymore and will be deactivated:\n"));
        text.append(tr("\n"));
        //display all invalid property timelines
        for (size_t i = 0; i < invalidTimelines.size(); ++i) {
           text.append(trim(invalidTimelines.at(i).first->getGuiName()).c_str());
           text.append(tr("."));
           text.append(trim(invalidTimelines.at(i).second->getProperty()->getGuiName()).c_str());
           text.append(tr("\n"));
        }

        QMessageBox::information(this, tr("Application Mode Configuration Changed"), text);

        for (size_t i = 0; i < invalidTimelines.size(); ++i)
           invalidTimelines.at(i).second->setActiveOnRendering(false);
    }

    if (workspace_)
        timelineWidget_->setApplicationModeConfig(&(workspace_->getApplicationModeConfig()));
    else
        timelineWidget_->setApplicationModeConfig(0);

    timelineWidget_->rebuildAnimation(animation_);
}

float AnimationEditor::getDuration() {
    return duration_;
}

void AnimationEditor::setWorkspace(Workspace* workspace) {

    stop();

    animation_ = 0;

    workspace_ = workspace;
    if (workspace_) {
        workspace_->addObserver(this);
        animation_ = workspace->getAnimation();
        timelineWidget_->setApplicationModeConfig(&(workspace_->getApplicationModeConfig()));
    }

    populateAddTimelineMenu();

    if (!workspace_) {
        timelineWidget_->setApplicationModeConfig(0);
        emit newAnimation(0);
        return;
    }
    else if(animation_) {
        duration_ = (int)animation_->getDuration()*30.0f;

        //check if this animation is valid and if not, display a message box to the user
        std::vector<std::pair<const Processor*, PropertyTimeline*> > invalidTimelines;

        if (!isAnimationValid(animation_, invalidTimelines)) {
            QString text(tr("The following animated properties are not available in the simple animation editor and will be deactivated:\n"));
            text.append(tr("\n"));
            //display all invalid property timelines
            for (size_t i = 0; i < invalidTimelines.size(); ++i) {
                text.append(trim(invalidTimelines.at(i).first->getGuiName()).c_str());
                text.append(tr("."));
                text.append(trim(invalidTimelines.at(i).second->getProperty()->getGuiName()).c_str());
                text.append(tr("\n"));
            }

            QMessageBox::information(this, tr("Deactivated Properties"), text);

            for (size_t i = 0; i < invalidTimelines.size(); ++i)
                invalidTimelines.at(i).second->setActiveOnRendering(false);
        }

        emit newAnimation(animation_);
    }
    else if(evaluator_ && evaluator_->getProcessorNetwork())
        newAnimation();
}

void AnimationEditor::init() {

    timer_ = new QTimer();
    connect(timer_, SIGNAL(timeout()), this, SLOT(update()));
    emit durationChanged(static_cast<int>(duration_));

}

void AnimationEditor::playerControl(QAction*action) {
    if(action->text() == "rewind") {
        rewind();
    }
    else if(action->text() == "play") {
        play();
    }
    else if(action->text() == "pause") {
        pause();
    }
    else if(action->text() == "stop") {
        stop();
    }
    else if(action->text() == "forward") {
        forward();
    }
}

void AnimationEditor::newAnimation() {
    stop();
    //temporarily store pointer to old animation to free memory later
    Animation* oldAnim = animation_;
    animation_ = new Animation(const_cast<ProcessorNetwork*>(evaluator_->getProcessorNetwork()));
    workspace_->setAnimation(animation_);
    duration_ = 600;
    animation_->setDuration(20);
    emit newAnimation(animation_);
    //delete old animation
    delete oldAnim;
}

void AnimationEditor::videoExport() {
    AnimationExportWidget* animationExportWidget
        = new AnimationExportWidget(this, animation_, evaluator_, duration_, 1, int(duration_), timeStretch_);
    animationExportWidget->resize(200,150);
    animationExportWidget->networkChanged();
    animationExportWidget->exec();
}

void AnimationEditor::undo() {
    animation_->undoLastChange();
}

void AnimationEditor::redo() {
    animation_->redoLastUndo();
}

void AnimationEditor::start() {
    currentFrame_  = 0;
    emit currentFrameChanged(currentFrame_);
}

void AnimationEditor::rewind() {
    emit playingSignal(true);
    if(interactionMode_)
        animation_->setInteractionMode(true);
    frameSkip_ = -4;
    timer_->start(30);
}

void AnimationEditor::pause(bool uncheck) {
    emit playingSignal(false);
    if(interactionMode_)
        animation_->setInteractionMode(false);
    animation_->renderAt((float)currentFrame_/30.0f);
    timer_->stop();
    if (uncheck)
        uncheckControls();
}

void AnimationEditor::play() {
    emit playingSignal(true);
    if(interactionMode_)
        animation_->setInteractionMode(true);
    frameSkip_ = 1;
    timer_->start(30);      // renders every 3/100 second
}

void AnimationEditor::stop() {
    emit playingSignal(false);
    if (!animation_)
        return;

    if(interactionMode_)
        animation_->setInteractionMode(false);
    timer_->stop();

    // Reset animation to first frame IF we have animated something already.
    if(currentFrame_ != 0) {
        currentFrame_ = 0;
        animation_->renderAt((float)currentFrame_/30.0f);
        emit currentFrameChanged(currentFrame_);
    }

    uncheckControls();
}

void AnimationEditor::forward() {
    emit playingSignal(true);
    if(interactionMode_)
        animation_->setInteractionMode(true);
    frameSkip_ = 4;
    timer_->start(30);
}
void AnimationEditor::end() {
    if(interactionMode_)
        animation_->setInteractionMode(false);
    currentFrame_ = static_cast<int>(duration_);
    timer_->stop();
    animation_->renderAt((float)currentFrame_/30.0f);
    emit currentFrameChanged(currentFrame_);

    uncheckControls();
}

void AnimationEditor::uncheckControls() {
    QList<QAction*> actions = playerControlGroup_->actions();
    for (int i = 0; i < actions.size(); ++i)
        actions.at(i)->setChecked(false);
}

void AnimationEditor::update() {
    if(interactionMode_)
        animation_->setInteractionMode(true);
    else
        animation_->setInteractionMode(false);
    if(currentFrame_ + frameSkip_ < 0) {
        currentFrame_ = 0;
        timer_->stop();
        uncheckControls();
        emit playingSignal(false);
    }
    else if(currentFrame_ + frameSkip_ > duration_) {
        currentFrame_ = static_cast<int>(duration_);
        timer_->stop();
        uncheckControls();
        emit playingSignal(false);
    }
    else
        currentFrame_ +=frameSkip_* timeStretch_;
    emit currentFrameChanged(currentFrame_);
    animation_->renderAt((float)currentFrame_/30.0f);
}

void AnimationEditor::setCurrentFrame(int currentFrame) {
    currentFrame_ = currentFrame;
    animation_->renderAt((float)currentFrame_/30.0f);
    emit currentFrameChanged(currentFrame_);
    timer_->stop();
    uncheckControls();
    emit playingSignal(false);
}

int AnimationEditor::getCurrentFrame() {
    return currentFrame_;
}

void AnimationEditor::setInteractionMode(bool on) {
    interactionMode_ = on;
}

void AnimationEditor::setDuration(int duration) {
    duration_ = duration;
    emit durationChanged(duration);
    animation_->setDuration((float)duration /30.0f);
}

void AnimationEditor::settings() {
    QDialog* dialog = new QDialog(this);
    dialog->setWindowTitle(tr("Animation Settings"));
    QVBoxLayout* lay = new QVBoxLayout(dialog);
    QHBoxLayout* buttonLayout = new QHBoxLayout();

    QGroupBox* durationBox = new QGroupBox(dialog);
    durationBox->setTitle("Duration (seconds)");
    QHBoxLayout* durationLayout = new QHBoxLayout(durationBox);
    QSpinBox dsp(dialog);
    dsp.setMinimum(1);
    dsp.setMaximum(200000);
    dsp.setValue(static_cast<int>(duration_/ 30.0f));
    durationLayout->addWidget(&dsp);
    QPushButton ok("Ok", dialog);
    QPushButton cancel("Cancel", dialog);
    connect (&ok, SIGNAL(clicked()), dialog, SLOT(accept()));
    connect (&cancel, SIGNAL(clicked()), dialog, SLOT(reject()));

    QLabel* tsLabel = new QLabel("Apply Time Stretch");
    tsLabel->setToolTip("Moves the keyframes to fit the new animation duration.");
    QHBoxLayout* timeStretchLayout = new QHBoxLayout();
    QCheckBox* applyTimeStretch = new QCheckBox(this);
    applyTimeStretch->setChecked(applyTimeStretch_);
    applyTimeStretch->setToolTip("Moves the keyframes to fit the new animation duration.");
    timeStretchLayout->addWidget(tsLabel);
    timeStretchLayout->addWidget(applyTimeStretch);
    //timeStretchLayout->addWidget(timeStretchSlider);
    timeStretchLayout->addStretch();

    buttonLayout->addWidget(&ok);
    buttonLayout->addWidget(&cancel);

    lay->addWidget(durationBox);
    lay->addLayout(timeStretchLayout);
    lay->addLayout(buttonLayout);
    //dialog->exec();
    if (dialog->exec() == QDialog::Accepted) {
        applyTimeStretch_ = applyTimeStretch->isChecked(); //remember the setting
        applyDurationSetting(static_cast<int>(dsp.value()* 30.0f), applyTimeStretch->isChecked());
    }
}

void AnimationEditor::applyDurationSetting(int duration, bool applyTimeStretch) {
    if (duration == duration_)
        return;

    //check if changing the duration will delete keys and notify user
    if (!applyTimeStretch && (duration < duration_)) {
        bool keysOverDuration = timelineWidget_->keysLaterThan(duration);
        if (keysOverDuration) {
            QString text(tr("This will delete all keys above the time selected as the new animation duration.\n"));
            text.append(tr("\n"));
            text.append(tr("Do you want to proceed?"));
            int answer = QMessageBox::question(this, tr("Proceed?"), text, QMessageBox::Yes, QMessageBox::No);

            //no action if user does not confirm
            if (!(answer == QMessageBox::Yes))
                return;
        }
    }

    duration_ = duration;

    if (!applyTimeStretch) {
        animation_->setDuration((float)duration /30.0f);
        emit durationChanged(duration);
    }
    else {
        timelineWidget_->applyTimeStretch(duration_);
    }
}
ApplicationModeConfiguration& AnimationEditor::getApplicationModeConfig() const {
    return workspace_->getApplicationModeConfig();
}

void AnimationEditor::populateAddTimelineMenu() {
    //clear the old menu and delete all submenus
    addTimelineMenu_->clear();

    for (std::map<std::string, QMenu*>::iterator i = applicationModeGroups_.begin(); i != applicationModeGroups_.end(); ++i) {
        i->second->clear();
        delete i->second;
    }
    applicationModeGroups_.clear();

    if (!workspace_)
       return;
    else
        analyzeNetwork();     // find property groups

    if (!animation_)
        return;

    //get application mode config groups for the correct order
    std::vector<std::string> groupOrderVector = workspace_->getApplicationModeConfig().getPropertyGroups();

    // store the actions for each group to sort them later
    std::map<std::string, std::map<Property*, QAction*> > groupActions;

    // iterate over animated processors and for each processor over the available properties
    const std::vector<AnimatedProcessor*> procs = animation_->getAnimatedProcessors();

    for(size_t i=0; i<procs.size();i++) {
        AnimatedProcessor* animProc = procs[i];
        const Processor* proc = animProc->getCorrespondingProcessor();

        const std::vector<Property*>& props = proc->getProperties();
        for(size_t i=0; i<props.size(); i++) {
            Property* prop = props[i];

            //check if this processor / property combination is master for a property group
            bool found = false;
            for (size_t group = 0; group < groups_.size(); ++group) {
                if ((groups_.at(group).at(0).first == proc) && (groups_.at(group).at(0).second == prop)) {
                    found = true;
                    break;
                }
            }
            if (!found)
                continue;

            // check if property is not animated yet, but can be animated
            if(!animProc->isPropertyAnimated(prop) && PropertyTimelineFactory::getInstance()->canPropertyBeAnimated(prop)) {

                //get the application mode group and add it to the corresponding submenu
                std::string currentApplicationModeGroup = workspace_->getApplicationModeConfig().getPropertyGroupMembership(prop).first;
                if ((currentApplicationModeGroup != "") && (applicationModeGroups_.find(currentApplicationModeGroup) == applicationModeGroups_.end())) {
                    // build a new submenu
                    applicationModeGroups_.insert(std::make_pair(currentApplicationModeGroup,
                                new QMenu(QString::fromStdString(currentApplicationModeGroup), this)));

                    groupActions.insert(std::make_pair(currentApplicationModeGroup, std::map<Property*, QAction*>()));
                }

                if (currentApplicationModeGroup == "") {
                    LERRORC("voreen.qt.AnimationEditor", "Found no application mode configuration group for property: " + prop->getGuiName());
                    continue;
                }
                else {
                    QMenu* currentSubmenu = applicationModeGroups_[currentApplicationModeGroup];
                    std::map<std::string, std::map<Property*, QAction*> >::iterator actionIterator = groupActions.find(currentApplicationModeGroup);
                    if (actionIterator != groupActions.end()) {
                        QAction* propAction = new QAction(QString::fromStdString(prop->getGuiName()), currentSubmenu);
                        QStringList sl;
                        sl << QString::fromStdString(proc->getID()) << QString::fromStdString(prop->getID());
                        propAction->setData(QVariant(sl));

                        actionIterator->second.insert(std::make_pair(prop, propAction));
                    }
                }
            }
        }
    }

    //add application mode groups in correct order
    for (std::vector<std::string>::iterator i = groupOrderVector.begin(); i != groupOrderVector.end(); ++i) {
        // if this group has a submenu: add it to the higher level menu
        std::map<std::string, QMenu*>::iterator mapIt = applicationModeGroups_.find(*i);
        if (mapIt != applicationModeGroups_.end()) {
            //get the QAction map for this group
            std::map<Property*, QAction*> actions = groupActions[*i];

            //get the properties for this group in correct order
            std::vector<Property*> currentGroupProperties = workspace_->getApplicationModeConfig().getGroupProperties(*i);

            //add all QActions that are animated to the current menu
            for (std::vector<Property*>::iterator propIt = currentGroupProperties.begin(); propIt != currentGroupProperties.end(); ++propIt) {
                if (actions.find(*propIt) != actions.end())
                    mapIt->second->addAction(actions[*propIt]);
            }

            addTimelineMenu_->addMenu(mapIt->second);
        }
    }

    if (addTimelineMenu_->isEmpty()) {
        // no items present -> add dummy item
        QAction* propAction = new QAction(QString::fromStdString("No properties available"), addTimelineMenu_);
        propAction->setEnabled(false);
        addTimelineMenu_->addAction(propAction);
    }
}

void AnimationEditor::removeTimelineFromAnimation(PropertyTimeline* tl) {
    //LINFOC("voreen.qt.AnimationEditor", "Delete " << tl);

    stop();

    //display message to the user
    QString text(tr("Removing this property from the animation will clear all of its animation state.\n"));
    text.append(tr("\n"));
    text.append(tr("Do you want to proceed?"));
    int answer = QMessageBox::question(this, tr("Proceed?"), text, QMessageBox::Yes, QMessageBox::No);

    //no action if user does not confirm
    if (!(answer == QMessageBox::Yes))
        return;

    AnimatedProcessor* animProc = 0;
    const std::vector<AnimatedProcessor*> procs = animation_->getAnimatedProcessors();
    for (std::vector<AnimatedProcessor*>::const_iterator i = procs.begin(); i != procs.end(); ++i) {
        if ((*i)->isPropertyAnimated(tl->getProperty())) {
            animProc = *i;
            break;
        }
    }

    if (animProc)
        animProc->removeTimeline(tl->getProperty());
    else
        LERRORC("voreen.qt.AnimationEditor", "Could not remove property timeline: found no animated processor");

    timelineWidget_->rebuildAnimation(animation_);
}

void AnimationEditor::addTimeline(QAction* action) {
    QStringList data = action->data().toStringList();
    if(data.size() == 2) {
        std::string procID = data.at(0).toStdString();
        std::string propID = data.at(1).toStdString();

        const std::vector<AnimatedProcessor*> procs = animation_->getAnimatedProcessors();
        for(size_t i=0; i<procs.size();i++) {
            AnimatedProcessor* animProc = procs[i];
            if(animProc->getProcessorName() == procID) {
                Processor* proc = animProc->getCorrespondingProcessor();
                if(proc) {
                    Property* prop = proc->getProperty(propID);
                    if(prop) {
                        PropertyTimeline* tl = animProc->addTimeline(prop);
                        if(tl) {
                            tl->registerUndoObserver(animation_);
                            //timelineWidget_->rebuildAnimation(animation_); //oldTODO call ProcessorTimelineWidget::showAnimatedProperties() instead
                            timelineWidget_->addPropertyTimeline(tl);
                        }
                    }
                }
            }
        }
    }
}

} // namespace voreen
