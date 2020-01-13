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

#ifndef ANIMATIONEDITOR_H
#define ANIMATIONEDITOR_H

#include <QWidget>
#include <QString>
#include <QToolButton>

#include "voreen/qt/voreenqtapi.h"
#include "voreen/core/animation/animation.h"
#include "voreen/core/network/networkevaluator.h"
#include "voreen/core/network/workspace.h"

class QLineEdit;
class QMenuBar;
class QMenu;
class QPushButton;
class QTimer;
class QLCDNumber;
class QActionGroup;

namespace voreen {

class AnimationTimelineWidget;

/**
* Global Container for all Animation specific Widgets.
* This is the place of animation control and it is the
* global interface node to the animation core
*/
class VRN_QT_API AnimationEditor : public QWidget, public WorkspaceObserver, public UsesWorkspace {
Q_OBJECT

public:
    AnimationEditor(NetworkEvaluator* eval, Workspace*, QWidget* = 0, Qt::WindowFlags = 0);

    virtual void setWorkspace(Workspace*);

    /// rebuilds the animation
    virtual void applicationModeConfigurationChanged();

    /// returns the current application mode configuration of the workspace
    ApplicationModeConfiguration& getApplicationModeConfig() const;

public slots:
    /// sets the current Frame
    void setCurrentFrame(int);
    float getDuration();

    /**
     * Checks if a processor / property combination is available in simple animation mode.
     */
    bool isPropertyAvailable(const Processor* processor, const Property* property) const;

protected slots:
    void newAnimation();
    void videoExport();
    void undo();
    void redo();
    void start();
    void rewind();
    void pause(bool uncheck = false);
    void play();
    void stop();
    void forward();
    void end();
    //void record();
    void setInteractionMode(bool);
    void settings();
    void playerControl(QAction*);
    void update();
    void setDuration(int);

    //void timeStretchChanged(double);
    void populateAddTimelineMenu();
    void addTimeline(QAction* action);

    void addButtonClicked();

    void removeTimelineFromAnimation(PropertyTimeline* tl);

protected:

    void applyDurationSetting(int duration, bool applyTimeStretch);

    //for calling stop() when selecting a different animation editor
    friend class AnimationTabView;

    NetworkEvaluator* evaluator_;
    QMenuBar* mainMenu_;
    QMenu* addTimelineMenu_;
    QActionGroup* playerControlGroup_;
    AnimationTimelineWidget* timelineWidget_;

    /// button for adding property timelines to the animation
    QToolButton* addTimelineButton_;

    bool applyTimeStretch_;  ///< remember the last user choice in the animation settings

    float currentFrame_;
    /// frameSkip determining the fastforward and backward speed
    int frameSkip_;
    float timeStretch_;
    bool interactionMode_;
    Workspace* workspace_;

    Animation* animation_;
    QTimer* timer_;
    float currentTime_;
    float duration_;        //in seconds
    /// String Representation for usage in LCD Displays
    QString getTimeString();
    //void clearGui();

    void init();

    /**
     * Analyzes the network to find groups of linked properties visible in application mode.
     * TODO: write proper comment and revise implementation
     */
    void analyzeNetwork();

    int getCurrentFrame();

    /** Checks if a given animation is valid for this animation editor (ie. if all animated properties are supported by this simple animation editor)
     * @param animation the animation to check
     * @param invalidPropertyTimelines if the animation is not valid, the timelines that are not allowed will be added to this vector (along with their corresponding processors)
     */
    bool isAnimationValid(Animation* animation, std::vector<std::pair<const Processor*, PropertyTimeline*> >& invalidPropertyTimelines) const;

    /**
     * Unchecks all player control buttons.
     */
    void uncheckControls();

    std::vector<std::deque<std::pair<Processor*, Property*> > > groups_; ///< groups of linked processors properties

    //remember the application mode config groups for submenus to the add timeline menu
    std::map<std::string, QMenu*> applicationModeGroups_;

signals:
    void currentFrameChanged(int);
    //void recordSignal();
    void durationChanged(int);
    void newAnimation(Animation*);
    void autoPreview(bool);
    //void deactivatedWindowSignal();
    void playingSignal(bool);
};

} // namespace voreen

#endif

