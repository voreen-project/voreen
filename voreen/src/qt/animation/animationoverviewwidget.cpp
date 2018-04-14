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

#include "voreen/qt/animation/animationoverviewwidget.h"
#include "voreen/qt/animation/animationexportwidget.h"
#include "voreen/qt/animation/currentframegraphicsitem.h"
//#include "voreen/core/ports/renderport.h"

#include <QContextMenuEvent>
#include <QDialog>
#include <QGraphicsItem>
#include <QGraphicsRectItem>
#include <QGraphicsScene>
#include <QGraphicsView>
#include <QImage>
#include <QLabel>
#include <QMenu>
#include <QPushButton>
#include <QScrollBar>
#include <QVBoxLayout>
#include <QWheelEvent>
#include <cmath>
#include <iostream>
#include <sstream>

namespace voreen {

#define PREVIEW_SIZE 40

AnimationOverviewView::AnimationOverviewView(QGraphicsScene* qgs, QWidget* parent, int duration)
    : QGraphicsView(qgs, parent)
    , barMovement_(false)
    , scene_(qgs)
    , slide_(false)
    , offset_(0)
    , duration_(duration)
{
    currentFrameGraphicsItem_ = new CurrentFrameGraphicsItem(true, true);
    setRenderHint(QPainter::SmoothPixmapTransform);
    //setCacheMode(CacheBackground);
    //setDragMode(QGraphicsView::ScrollHandDrag);
    setVerticalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
    setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOn);
    horizontalScrollBar()->setContextMenuPolicy(Qt::NoContextMenu);
    setAlignment(Qt::AlignLeft);
    scene_->setBackgroundBrush(Qt::blue);
    QLinearGradient gradient(0, 0, 0, 80);
    gradient.setSpread(QGradient::PadSpread);
    gradient.setColorAt(0, QColor(255,255,255,0));
    gradient.setColorAt(0.5, QColor(150,150,150,200));
    gradient.setColorAt(1, QColor(50,0,50,100));
    scene_->setBackgroundBrush(gradient);
    scene_->addItem(currentFrameGraphicsItem_);
}

void AnimationOverviewView::setBar(QGraphicsRectItem* bar) {
    highlightBar_ = bar;
}

void AnimationOverviewView::setDuration(int duration) {
    duration_ = duration;

    //set duration to keyframe items
    float fDuration = static_cast<float>(duration);
    for (int i = 0; i < scene_->items().size(); ++i)
        if (KeyframeGraphicsItem* item = dynamic_cast<KeyframeGraphicsItem*>(scene_->items().at(i)))
            item->setDuration(fDuration);
}

void AnimationOverviewView::contextMenuEvent(QContextMenuEvent* e) {
    if (KeyframeGraphicsItem* key = dynamic_cast<KeyframeGraphicsItem*>(scene_->itemAt(mapToScene(e->pos()).x(), mapToScene(e->pos()).y(), QTransform())))
        emit keyframeContextMenuRequest(key, e->pos());
    else if (scene_->itemAt(mapToScene(e->pos()).x(), mapToScene(e->pos()).y(), QTransform()) == highlightBar_)
        return;
    else
        emit contextMenuRequest(e->pos());
}

void AnimationOverviewView::setCurrentFrame(int frame) {
    currentFrame_ = frame;
    currentFrameGraphicsItem_->setPos(frame, 0);
}

void AnimationOverviewView::mousePressEvent(QMouseEvent* e) {

    if (scene_->itemAt(mapToScene(e->pos()).x(), mapToScene(e->pos()).y(), QTransform())
        && (scene_->itemAt(mapToScene(e->pos()).x(), mapToScene(e->pos()).y(), QTransform())->boundingRect() == highlightBar_->boundingRect()
            || scene_->itemAt(mapToScene(e->pos()).x(), mapToScene(e->pos()).y(), QTransform())->boundingRect() == currentFrameGraphicsItem_->boundingRect())
        && e->button() == Qt::LeftButton) {
        QScrollBar* scroll = horizontalScrollBar();
        if (e->x()+ scroll->value() >= 0 && e->x()+ scroll->value() <= duration_) {
            currentFrameGraphicsItem_->setPos(e->x() + scroll->value(), 0);
            slide_ = true;
            emit currentFrameChanged(e->x() + scroll->value());
        }
    }
    else if((e->button() == Qt::LeftButton) /*&& !scene_->itemAt(mapToScene(e->pos()).x(), mapToScene(e->pos()).y())*/) {

        if(!(e->modifiers() & Qt::CTRL))
            //interval selected
            emit intervalSelectedAt(mapToScene(e->pos()));
        else {
            clearSelection();
            emit clearSelection();
        }
    }
    QGraphicsView::mousePressEvent(e);
}

void AnimationOverviewView::mouseMoveEvent(QMouseEvent* e) {
    QScrollBar* scroll = horizontalScrollBar();

    if (slide_ && e->x()+ scroll->value() >= 0 && e->x() + scroll->value() <= duration_) {
        currentFrameGraphicsItem_->setPos(e->x() + scroll->value(), 0);
        emit currentFrameChanged(e->x() + scroll->value());
    }
    else if (slide_ && e->x()+ scroll->value() > duration_) {
        currentFrameGraphicsItem_->setPos(duration_, 0);
        emit currentFrameChanged(static_cast<int>(duration_));
    }
    else if (slide_ && e->x()+ scroll->value() < 0) {
        currentFrameGraphicsItem_->setPos(0, 0);
        emit currentFrameChanged(0);
    }

    QGraphicsView::mouseMoveEvent(e);
}

void AnimationOverviewView::mouseReleaseEvent(QMouseEvent* e) {
    slide_ = false;
    barMovement_ = false;
    //offset_ = static_cast<int>(highlightBar_->boundingRect().x());
    QGraphicsView::mouseReleaseEvent(e);
}

void AnimationOverviewView::offsetCorrection(int x) {
    if (!barMovement_)
        offset_ = x;
}

// ---------------------------------------------------------------------------

AnimationOverviewWidget::AnimationOverviewWidget(QWidget* parent, NetworkEvaluator* networkEval)
    : QWidget(parent)
    , duration_(1200)
    , zoom_(1)
    , viewportWidth_(500)
    , scrollBarPosition_(0)
    , canvasRenderer_(0)
    , networkEvaluator_(networkEval)
    , currentFrame_(0)
    , fps_(30)
    , autoPreview_(false)
    , currentlyPlaying_(false)
{
    canvasRenderer_ = 0;
    previews_.clear();
    QHBoxLayout* mainLayout = new QHBoxLayout(this);

    mainLayout->setMargin(1);
    mainLayout->setSpacing(1);
    setFixedHeight(70 + PREVIEW_SIZE);

    overviewScene_ = new QGraphicsScene(this);

    overviewView_ = new AnimationOverviewView(overviewScene_, this, duration_);

    highlightBar_ = new QGraphicsRectItem();
    overviewView_->setBar(highlightBar_);
    highlightBar_->setZValue(2);
    highlightBar_->setPen(QPen(Qt::DotLine));
    QLinearGradient gradient1(0,0,0, PREVIEW_SIZE);
    gradient1.setSpread(QGradient::ReflectSpread);
    gradient1.setColorAt(0.0, QColor(20, 100 ,100, 30));
    gradient1.setColorAt(1.0, QColor(80, 100 ,100, 30));
    QBrush brush(gradient1);
    highlightBar_->setBrush(gradient1);
    overviewScene_->addItem(highlightBar_);

    overviewView_->setStyleSheet("background:transparent");
    connect(overviewView_, SIGNAL(currentFrameChanged(int)), this, SIGNAL(currentFrameChanged(int)));
    connect(overviewView_, SIGNAL(intervalSelectedAt(QPointF)), this, SIGNAL(intervalSelectedAt(QPointF)));
    connect(overviewView_, SIGNAL(clearSelection()), this, SIGNAL(clearSelection()));
    connect(overviewView_, SIGNAL(contextMenuRequest(QPoint)), this, SLOT(contextMenuRequest(QPoint)));
    connect(overviewView_, SIGNAL(keyframeContextMenuRequest(KeyframeGraphicsItem*, QPoint)), this, SLOT(keyframeContextMenuRequest(KeyframeGraphicsItem*, QPoint)));
    //connect(overviewView_, SIGNAL(barMovement(int)), this, SIGNAL(barMovement(int)));
    connect(this, SIGNAL(offsetCorrection(int)), overviewView_, SLOT(offsetCorrection(int)));

    overviewView_->setGeometry(QRect(0, 0, duration_, 50));
    overviewView_->setSceneRect(QRect(0, 0, duration_, 50));
    overviewView_->setCurrentFrame(0);

    mainLayout->addWidget(overviewView_);
    QLabel* nameLabel = new QLabel(this);
    nameLabel->lower();
    nameLabel->move(8, 45);
    nameLabel->setText(QString::fromStdString("Overview"));

    if (canvasRenderer_) {
        const std::vector<Port*> port = canvasRenderer_->getInports();
        if (port.size() == 1) {
            if(dynamic_cast<RenderPort*>(port.at(0)) != 0)
                renderPort_ = dynamic_cast<RenderPort*>(port.at(0));
        }
    }
    setStandardRenderport();
}

int AnimationOverviewWidget::getDuration() const {
    return duration_;
}

int AnimationOverviewWidget::getCurrentFrame() const {
    return currentFrame_;
}

void AnimationOverviewWidget::autoPreview(bool autoPreview) {
    autoPreview_ = autoPreview;

    renderPreviews();
}

void AnimationOverviewWidget::reset() {
    overviewScene_->clearSelection();
    overviewView_->setCurrentFrame(0);
    setFps(30);
    scrollBarOrder(0);

    renderPreviews();
}

void AnimationOverviewWidget::addKeyframeToScene(KeyframeGraphicsItem* kfgi) {
    overviewScene_->addItem(kfgi);
}

void AnimationOverviewWidget::removeKeyframeFromScene(KeyframeGraphicsItem* kfgi) {
    overviewScene_->removeItem(kfgi);
}

void AnimationOverviewWidget::addIntervalToScene(QGraphicsItem* item) {
    overviewScene_->addItem(item);
}

void AnimationOverviewWidget::removeIntervalFromScene(QGraphicsItem* item) {
    overviewScene_->removeItem(item);
}

void AnimationOverviewWidget::contextMenuRequest(QPoint pos) {
    emit clearSelection();
    QMenu* addFrameMenu = new QMenu(this);
    addFrameMenu->setStyleSheet("background:white");
    QAction addFrameAction(tr("Add Keyframe"), this);
    addFrameAction.setDisabled(currentlyPlaying_);
    QAction snapshotAction(tr("Take Snapshot"), this);
    snapshotAction.setDisabled(currentlyPlaying_);
    QAction regeneratePreviews(tr("Update Previews"), this);
    regeneratePreviews.setDisabled(currentlyPlaying_);

    addFrameMenu->addAction(&addFrameAction);
    addFrameMenu->addAction(&snapshotAction);
    addFrameMenu->addSeparator();
    if (autoPreview_)
        addFrameMenu->addAction(&regeneratePreviews);

    QMenu* canvasRendererMenu = new QMenu("Canvasrenderer", this);
    const ProcessorNetwork* network = networkEvaluator_->getProcessorNetwork();
    const std::vector<CanvasRenderer*>& canvasRenderer = network->getProcessorsByType<CanvasRenderer>();
    std::map<QAction*, CanvasRenderer*> menuMap;

    for (size_t i = 0; i < canvasRenderer.size(); ++i) {
        QString canvasRendererName = QString::fromStdString(canvasRenderer.at(i)->getID());
        QAction* action = new QAction(canvasRendererName, canvasRendererMenu);
        action->setDisabled(currentlyPlaying_);
        canvasRendererMenu->addAction(action);
        menuMap[action] = canvasRenderer.at(i);
    }

    addFrameMenu->addMenu(canvasRendererMenu);
    QAction* action = addFrameMenu->exec(QCursor::pos());
    if (action == &addFrameAction) {
        emit addKeyframe(overviewView_->mapToScene(pos));
    }
    else if (action == &snapshotAction) {
        emit recordAt(static_cast<int>(overviewView_->mapToScene(pos).x()));
        if(autoPreview_)
            renderPreviews();
    }
    else if (action == &regeneratePreviews) {
        updatePreviews();
    }
    else if (action != 0){
        canvasRenderer_ = menuMap[action];
        const std::vector<Port*> port = canvasRenderer_->getInports();
        if (port.size() == 1) {
             if(dynamic_cast<RenderPort*>(port.at(0)) != 0)
                 renderPort_ = dynamic_cast<RenderPort*>(port.at(0));
         }
        if(autoPreview_)
            renderPreviews();
    }
}

void AnimationOverviewWidget::keyframeContextMenuRequest(KeyframeGraphicsItem* kfgi, QPoint pos) {
    //deselect previously selected items and select the current
    overviewScene_->clearSelection();
    kfgi->setSelected(true);
    emit itemSelected(kfgi);

    QMenu* keyFrameMenu = new QMenu(this);
    keyFrameMenu->setStyleSheet("background:white");
    QAction removeFrameAction(tr("Remove key frame"), this);
    QAction snapshotAction(tr("Take Snapshot"), this);

    keyFrameMenu->addAction(&removeFrameAction);
    keyFrameMenu->addAction(&snapshotAction);

    QAction* action = keyFrameMenu->exec(QCursor::pos());
    if (action == &removeFrameAction) {
        emit removeKeyframe(kfgi/*overviewView_->mapToScene(pos)*/);
    }
    else if (action == &snapshotAction) {
        emit recordToKeyframe(kfgi);
        //LWARNINGC("voreen.qt.AnimationOverviewWidget", "Snapshot function not supported yet!");
        if(autoPreview_)
            renderPreviews();
    }
}
void AnimationOverviewWidget::setStandardRenderport() {
    const ProcessorNetwork* network = networkEvaluator_->getProcessorNetwork();
    if(network != 0) {
        const std::vector<CanvasRenderer*>& canvasRenderer = network->getProcessorsByType<CanvasRenderer>();
        if (!canvasRenderer.empty()) {
            canvasRenderer_ = *canvasRenderer.begin();
            const std::vector<Port*> port = canvasRenderer_->getInports();
            if (port.size() == 1) {
                if (dynamic_cast<RenderPort*>(port.at(0)) != 0)
                    renderPort_ = dynamic_cast<RenderPort*>(port.at(0));
            }
        }
        else
            renderPort_ = 0;
    }
}

void AnimationOverviewWidget::selectCanvasRenderer() {
   QMenu* menu = new QMenu(this);
   const ProcessorNetwork* network = networkEvaluator_->getProcessorNetwork();
   const std::vector<CanvasRenderer*>& canvasRenderer = network->getProcessorsByType<CanvasRenderer>();
   std::map<QAction*, CanvasRenderer*> menuMap;

   for (size_t i = 0; i < canvasRenderer.size(); ++i) {
       QString canvasRendererName = QString::fromStdString(canvasRenderer.at(i)->getID());
       QAction* action = new QAction(canvasRendererName, menu);
       menu->addAction(action);
       menuMap[action] = canvasRenderer.at(i);
   }
   QAction* result = menu->exec(QCursor::pos());
   if (result != 0) {
       canvasRenderer_ = menuMap[result];
       const std::vector<Port*> port = canvasRenderer_->getInports();
       if (port.size() == 1) {
            if (dynamic_cast<RenderPort*>(port.at(0)) != 0)
                renderPort_ = dynamic_cast<RenderPort*>(port.at(0));
        }
   }
   delete menu;
   if(autoPreview_)
       renderPreviews();
}

void AnimationOverviewWidget::sceneOrder(QMatrix matrix) {
    matrix_ = matrix;
    zoom_ = static_cast<int>(matrix.m11());
}

void AnimationOverviewWidget::playingStateChanged(bool playing) {
    currentlyPlaying_ = playing;
}

void AnimationOverviewWidget::scrollBarOrder(int scrollBarPosition) {
    /*if (zoom_ != 0) {
        scrollBarPosition_ = scrollBarPosition;
        highlight_->setRect(scrollBarPosition_/zoom_, -15, (viewportWidth_) / zoom_ , 85);
        highlightBar_->setRect(scrollBarPosition_/zoom_, -15, (viewportWidth_) / zoom_ , 15);
        emit offsetCorrection(static_cast<int>(highlightBar_->boundingRect().x()));
    }*/
}

void AnimationOverviewWidget::updatePreviews() {
    if (autoPreview_)
        renderPreviews();
}

void AnimationOverviewWidget::renderPreviews() {
    if(autoPreview_) {
        if (renderPort_) {
            int currentFrame = currentFrame_;
            if (!previews_.empty()) {
                for (size_t vs = 0; vs < previews_.size(); vs++) {
                    overviewScene_->removeItem(previews_.at(vs));
                    delete previews_.at(vs);
                }
            }
            previews_.clear();
            int quadSize = PREVIEW_SIZE;
            int stepLength = duration_ / quadSize;
            if (stepLength < quadSize)
                stepLength = quadSize;
            // save the current settings
            tgt::ivec2 size = renderPort_->getSizeOriginProperty()->get();
            renderPort_->requestSize(tgt::ivec2(quadSize,quadSize));
            QPixmap previewPixmap;
            QImage* preview = 0;
            for (int i = 0; i < duration_ - 1; i += stepLength) {
                preview = new QImage(QSize(quadSize, quadSize), QImage::Format_ARGB32);
                emit currentFrameChanged(i);
                networkEvaluator_->process();
                tgt::col4* buffer = renderPort_->readColorBuffer<uint8_t>(); // TODO: catch exceptions
                for (int x = quadSize-1; x >= 0; x-=1) {
                    for (int y = 0; y < quadSize; y+=1) {
                        QColor color(buffer->r, buffer->g, buffer->b, 255);
                        QRgb rgba = color.rgba();
                        preview->setPixel(y,x, rgba);
                        buffer++;
                    }
                }
                const QImage img = preview->copy(0,0, quadSize,quadSize);

                previewPixmap = QPixmap::fromImage(img);
                QGraphicsPixmapItem* pixmapItem = overviewScene_->addPixmap(previewPixmap);
                previews_.push_back(pixmapItem);
                pixmapItem->moveBy(i, 26);
                pixmapItem->setZValue(-1);
                delete preview;
            }
            renderPort_->requestSize(size);
            setCurrentFrame(currentFrame);
            emit currentFrameChanged(currentFrame);
        }
    }
    else {
        if (!previews_.empty()) {
                for (size_t vs = 0; vs < previews_.size(); vs++) {
                    overviewScene_->removeItem(previews_.at(vs));
                    delete previews_.at(vs);
                }
            }
        previews_.clear();
    }
}

void AnimationOverviewWidget::setCurrentFrame(int currentFrame) {
    currentFrame_ = currentFrame;
    overviewView_->setCurrentFrame(currentFrame);
}

void AnimationOverviewWidget::setFps(int fps) {
    QList<QGraphicsItem*> items = overviewScene_->items();
    for (int i = 0; i < items.size(); ++i) {
        if (!(dynamic_cast<CurrentFrameGraphicsItem*>(items.at(i))
                || dynamic_cast<QGraphicsRectItem*>(items.at(i))
                || dynamic_cast<KeyframeGraphicsItem*>(items.at(i))))
            overviewScene_->removeItem(items.at(i));
    }
    fps_ = fps;
    int secs = 0;
    for (int x = 0; x < duration_ + fps; x+=fps) {
        if (x % (2*fps) == 0) {
            overviewScene_->addLine(x, -14/*15*/, x, -2/*35*/);
            QGraphicsTextItem* textItem = overviewScene_->addText(QString::fromStdString(getTimeString(secs)));
            textItem->moveBy(x-14, -5);

        }
        else if (x % fps == 0) {
            overviewScene_->addLine(x, -14/*20*/, x, -7/*30*/);
        }
        secs+=1;
    }
    overviewView_->setCurrentFrame(0);
    //overviewScene_->addRect(0, -15/*25*/, duration_, 1);
    overviewView_->setSceneRect(-2,-2, duration_, 60);

    highlightBar_->setRect(0, -15, duration_, 30/*15*/);
    highlightBar_->show();

    previews_.clear();
}

void AnimationOverviewWidget::updateViewportRect(int viewportWidth) {
    viewportWidth_ = viewportWidth;
    //highlight_->setRect(scrollBarPosition_/zoom_, -15, (viewportWidth_) / zoom_ , 80);
    //highlightBar_->setRect(0, -15, viewportWidth_ , 15);
}

void AnimationOverviewWidget::setDuration(int duration) {
    duration_ = duration;
    overviewView_->setDuration(duration);
    setFps(fps_);
}

std::string AnimationOverviewWidget::getTimeString(int sec) {
    std::stringstream out;
    if (sec < 60) {
        out << "00:";
        if (sec < 10)
            out << "0" << sec;
        else
            out << sec;
    }
    else if (sec < 600) {
        int min = 0;
        while (sec >= 60) {
            sec -= 60;
            min += 1;
        }
        out << "0" << min;
        if (sec < 10)
            out << ":0" << sec;
        else
            out << ":" << sec;
    }
    else {
        int min = 0;
        while (sec > 60) {
            sec -= 60;
            min += 1;
        }
        if (sec < 10)
            out << min << ":0" << sec;
        else
            out << min << ":" << sec;
    }
    return out.str();
}

} //namespace voreen
