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

#ifndef ANIMATIONTIMELINEWIDGET_H
#define ANIMATIONTIMELINEWIDGET_H

#include "voreen/core/animation/animation.h"
#include "voreen/core/animation/animationobserver.h"
#include "voreen/core/network/networkevaluator.h"
#include "voreen/qt/animation/animationeditor.h"
#include "voreen/qt/animation/animationkeyframewidget.h"
#include "voreen/qt/animation/animationintervalwidget.h"
#include "voreen/qt/animation/animationpropertylistwidget.h"
#include "keyframegraphicsitem.h"

#include <QWidget>
#include <QGraphicsRectItem>

class AnimatedProcessor;
class AnimationOverviewWidget;
class QGroupBox;
class QVBoxLayout;
class QKeyEvent;


namespace voreen {

class AnimationTimelineWidget : public QWidget, public AnimationObserver {
Q_OBJECT
public:
    AnimationTimelineWidget(Animation*, AnimationEditor*, NetworkEvaluator*);
    float getCurrentTime();
    void setCurrentTime(float /*time*/);

    /**
    * this function is called if there was an undo or a redo
    * all current ProcessorsTimelineWidgets and Propertytimelinewidgets are not valid anymore and have to be rebuild
    * this function is not called after creation
    */
    void reloadAnimation();

    /**
    * this function is called if a new animatedProcessors is added
    */
    void animatedProcessorAdded(const AnimatedProcessor* processor);

    /**
    * this function is called by the animation if an animatedProcessor gets deleted (ie. a processor is removed from the network).
    * rebuilds the animation.
    */
    void animatedProcessorRemoved(const AnimatedProcessor* /*processor*/);

    Animation* getAnimation() { return animation_; }

    /// called by simple animation editor if a property is added, puts all keyframes in it and rebuilds the animation
    void addPropertyTimeline(PropertyTimeline* tl);

    /// applies a time stretch to the animation, ie. changes duration and moves keys according to the new duration
    void applyTimeStretch(int duration);

    /// checks if keys exist that are positioned at a later time position than time
    bool keysLaterThan(int time) const;

    /// returns if any property timelines are currently animated (even if currently deactivated!)
    bool hasAnimatedProperties() const;

    void setApplicationModeConfig(ApplicationModeConfiguration* appConfig);

public slots:

    /// setWorkspace
    void rebuildAnimation(Animation*);

    /// rebuild the current animation
    void rebuildCurrentAnimation();

    //void timelineActivated();

    void durationChanged(int);

    void clearGuiSelection();

protected:

    /// takes all present timelines from the core and initializes templatepropertytimelinewidgets for them
    void populateProcessors();

    /// the Animation core, usually obtained from the workspace if present
    Animation* animation_;
    /// the current time
    float currentTime_;
    /// Mainlayout
    QVBoxLayout* mainLayout_;

    /// returns a String representation for LCD Displays of the current time
    QString getTimeString(int);

    void keyPressEvent(QKeyEvent*);

    /// LCDTimeCounter for visual feedback
    QLabel* timeCounter_;

    /// the type of the currently selected item (keyframe or interval)
    QLabel* currentItemLabel_;

    /// the time of the currently selected item (keyframe / interval)
    QLabel* currentItemTime_;

    /// QWidget reimplementation of the QWidget resizeevent
    void resizeEvent(QResizeEvent*);
    /// Scrollarea which encapsulates all Processortimelines
    QScrollArea* scrollArea_;
    /// layout for the Scrollarea
    QVBoxLayout* scrollAreaLayout_;
    QWidget* containerWidget_;
    AnimationOverviewWidget* overviewTimeline_;
    QGroupBox* timeBox_;

    AnimationKeyframeWidget* keyframeWidget_;
    AnimationIntervalWidget* intervalWidget_;
    AnimationPropertyListWidget* listWidget_;

    //builds the union of the key values of all property timelines as graphical items
    void insertAllKeyValuesAsGraphicalItems();

    //puts the graphical key frames as key values into all currently animated property timelines and builds a map, associating the graphics items with the corresponding property key values
    void putKeyframesIntoTimelines();

    // indicates wether there was a change once
    bool changed_;

    //keyframes
    std::vector<KeyframeGraphicsItem*> keyframeGraphicsItems_;

    //the currently animated property timelines
    std::vector<PropertyTimeline*> propertyTimelines_;

    // associates graphical keyframe items with the corresponding key values in all animated property timelines
    std::map<KeyframeGraphicsItem*, std::vector<std::pair<PropertyKeyValueBase*, PropertyTimeline*> > > keyframeAssociations_;

    void removeWidgetsFromLayout(QLayout* layout);

    void resetSelectionOfItems();

    /**
     * If nothing is selected, this sets the list of currently animated property timelines to the scrollable area.
     * If a keyframe is selected, the scroll displays the property widgets for each timeline with the values of the selected key.
     * If an interval is selected, the property widgets of the left and right key as well as the interpolation functions between them are displayed.
     */
    void updateScrollArea();

    bool keyframeSelected_;                     ///< true if a keyframe item is currently selected
    bool intervalSelected_;                     ///< true if an interval between two keyframes is currently selected
    KeyframeGraphicsItem* currentKfgi_;         ///< currently selected keyframe graphics item
    KeyframeGraphicsItem* leftKfgi_;            ///< keyframe to the left of the currently selected interval
    KeyframeGraphicsItem* rightKfgi_;           ///< keyframe to the right of the currently selected interval
    QGraphicsRectItem* intervalGraphicsItem_;              ///< rectangle for the selected interval


protected slots:
    /// renders the current Animation at the given time
    void renderAt(float);
    /// sets the currentFrame
    void currentFrame(int);

    /**
     * Adds a Keyframe at given point in animation core.
     * The values of the keyframe properties is computed by interpolating the left and right neighbor keys (if existing).
     */
    KeyframeGraphicsItem* addKeyframe(QPointF);

    void removeKeyframe(KeyframeGraphicsItem*);

    /**
     * Records the current network settings of the animated properties to an existing keyframe
     */
    void recordToKeyframe(KeyframeGraphicsItem*);

    /**
     * Records the current network settings of the animated properties into a new keyframe
     */
    void recordAt(int);

    ///Invoken when an item is clicked
    void itemClicked(KeyframeGraphicsItem*);
    ///Invoken when an item is released
    void itemReleased(KeyframeGraphicsItem*, bool shift);
    /// Item is moving
    void itemMoving(KeyframeGraphicsItem* kfgi);

    /// removes graphical key items that are not within the current animation duration
    void removeKeysAboveDuration(int);

    void intervalSelectedAt(QPointF);

signals:
    /// orders a zoom to all timelines with zommfactor of the given integer
    void zoomOrder(int);

    void showActiveTimelines();
    /// orders a scene translation via the given matrix
    void sceneOrder(QMatrix);
    /// order to distribute vertical scrollbarpositions
    void scrollBarOrder(int);
    /// signals a framechange
    void currentFrameChanged(int);
    /// used to record the network state which is taken from the voreenve app
    void recordSignal();
    /**
    * signals a resizeEvent. Mainly used to work around some Layout bugs related to QScrollAreas
    * two days of ircing and forum writing couldn't solve this any other way.
    */
    void resizeSignal(int);
    /// The duration of the Animation has changed
    //void durationChanged(int);
    /// The fps of the Animation has changed
    void fpsChanged(int);
    /**
    * signals a resizeEvent which comes from the propertytimelineview. This is used to determine the
    * geometry of the Viewportvisualisation in the OverviewWidget
    */
    void viewResizeSignal(int);

    void setAnimationEditorDuration(int);

    void autoPreview(bool);
    void updatePreviews();

    void playingSignal(bool);

    void deleteTimelineSignal(PropertyTimeline*);
};

} // namespace voreen

#endif
