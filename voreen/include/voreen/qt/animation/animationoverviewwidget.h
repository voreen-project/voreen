/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2019 University of Muenster, Germany,                        *
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

#ifndef ANIMATIONOVERVIEWWIDGET_H
#define ANIMATIONOVERVIEWWIDGET_H

#include "modules/core/processors/output/canvasrenderer.h" //< core module is always available
#include "voreen/core/network/networkevaluator.h"

#include "currentframegraphicsitem.h"
#include "keyframegraphicsitem.h"


#include <QWidget>
#include <QGraphicsScene>
#include <QGraphicsView>

class QGraphicsScene;
class QGraphicsRectItem;
class QPixmap;
class RenderPort;

namespace voreen {

/**
* The View for the Overviewwidget
*/
class AnimationOverviewView : public QGraphicsView {
    Q_OBJECT
public:
   AnimationOverviewView(QGraphicsScene*, QWidget* = 0, int duration = 1200);
   void setBar(QGraphicsRectItem*);

   //void clearFocus();

   void setDuration(int);

public slots:

    /// Sets the current Frame. It's used if currentFrame changes come the AnimationEditor
    void setCurrentFrame(int);
    void offsetCorrection(int);

protected:
   bool barMovement_;
   int relativeBarPosition_;
   /// the currently selected frame
   QGraphicsRectItem* highlightBar_;
   int currentFrame_;
   QGraphicsScene* scene_;
   bool slide_;
   /// GraphicsItem representing the current Frame
   CurrentFrameGraphicsItem* currentFrameGraphicsItem_;
   virtual void contextMenuEvent(QContextMenuEvent*);
   virtual void mousePressEvent(QMouseEvent*);
   virtual void mouseMoveEvent(QMouseEvent*);
   void mouseReleaseEvent(QMouseEvent*);
   int offset_;

   int duration_;

signals:
    void currentFrameChanged(int);

    void contextMenuRequest(QPoint);
    void keyframeContextMenuRequest(KeyframeGraphicsItem*, QPoint);
    //void barMovement(int);

    /// for selecting intervals between keys
    void intervalSelectedAt(QPointF);

    void clearSelection();

};

// ---------------------------------------------------------------------------

/**
* displays a timeline for global animationorientation and intuitiv framechanges
*/
class AnimationOverviewWidget : public QWidget {
    Q_OBJECT
public:
    AnimationOverviewWidget(QWidget* = 0, NetworkEvaluator* = 0);
    // resets the Overview to its start status
    void reset();
    int getCurrentFrame() const;
    int getDuration() const;

    void addKeyframeToScene(KeyframeGraphicsItem* kfgi);
    void removeKeyframeFromScene(KeyframeGraphicsItem* kfgi);

    void addIntervalToScene(QGraphicsItem* item);
    void removeIntervalFromScene(QGraphicsItem* item);

public slots:
    void contextMenuRequest(QPoint);
    void keyframeContextMenuRequest(KeyframeGraphicsItem*, QPoint);
    /// sets the current frame
    void setCurrentFrame(int);
    void sceneOrder(QMatrix);
    void scrollBarOrder(int);
    void setFps(int);
    void setDuration(int);
    void updateViewportRect(int);
    void renderPreviews();
    void setStandardRenderport();
    void autoPreview(bool);
    void updatePreviews();
    void playingStateChanged(bool);

protected:
    QGraphicsScene* overviewScene_;
    AnimationOverviewView* overviewView_;
    int xPosition_;
    QMatrix matrix_;
    /// highlights the area which is shown in the propertytimelinewidgets
    //QGraphicsRectItem* highlight_;
    QGraphicsRectItem* highlightBar_;
    int duration_;
    int zoom_;
    int viewportWidth_;
    int scrollBarPosition_;
    CanvasRenderer* canvasRenderer_;
    NetworkEvaluator* networkEvaluator_;
    int currentFrame_;
    int fps_;
    RenderPort* renderPort_;
    std::vector<QGraphicsPixmapItem*> previews_;
    std::string getTimeString(int sec);
    bool autoPreview_;

    bool currentlyPlaying_;

    /// links GraphicsItem with the corresponding key values of the property timelines
    //std::map<KeyframeGraphicsItem*, std::vector<std::pair<PropertyTimeline*, const PropertyKeyValueBase*> > > keyframes_;

protected slots:
    void selectCanvasRenderer();

signals:
    /// signals a change of the current frame
    void currentFrameChanged(int);
    /// snapshots the current network state at the x position of the mouse
    void snapshotNetwork(int);

    void recordAt(int);
    //void barMovement(int);
    void offsetCorrection(int);

    /// for adding keyframes
    void addKeyframe(QPointF);

    /// for removing keyframes
    void removeKeyframe(KeyframeGraphicsItem*);

    void itemSelected(KeyframeGraphicsItem*);

    void recordToKeyframe(KeyframeGraphicsItem*);

    /// for selecting intervals between keys
    void intervalSelectedAt(QPointF);

    void clearSelection();
};

} // namespace voreen

#endif
