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

#include "voreen/qt/widgets/property/transfunc/2d/2dprimitives/transfunc2dprimitivespropertyeditorcanvas.h"

#include "voreen/core/properties/transfunc/2d/2dprimitives/transfunc2dprimitivesproperty.h"
#include "voreen/core/datastructures/transfunc/2d/2dprimitives/transfunc2dprimitives.h"

#include "voreen/core/datastructures/transfunc/1d/1dkeys/utils/transfuncmappingkey.h"


#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/datastructures/volume/volumedisk.h"
#include "voreen/core/datastructures/volume/volumedecorator.h"
#include "voreen/core/datastructures/volume/histogram.h"


#include <QAction>
#include <QApplication>
#include <QColor>
#include <QColorDialog>
#include <QMouseEvent>
#include <QPainter>
#include <QString>
#include <QToolTip>
#include <QImage>

#include <iostream>

namespace  {
    const QString X_AXIS_TEXT = QString("intensity");
    const QString Y_AXIS_TEXT = QString("gradient");

    const int AXIS_OFFSET  = 12;
    const int ARROW_LENGTH = 10;
    const int ARROW_WIDTH  = 3;
    const int KEY_POINT_SIZE = 10;
    const float KEY_SPLIT_FACTOR = 1.5f;
    const int MIN_CELL_SIZE = 8;
}

namespace voreen {
TransFunc2DPrimitivesPropertyEditorCanvas::TransFunc2DPrimitivesPropertyEditorCanvas(QWidget* parent, TransFunc2DPrimitivesProperty* prop)
    : QWidget(parent)
    //interaction
    , selectedKey_(0), selectedPrimitive_(0), isBeingDraged_(false)
    //paint
    , gridSpacing_(0.1f,0.1f), histogram_(0), histogramCache_(0), visibleHistogramRangeX_(0.f,1.f), visibleHistogramRangeY_(0.f,1.f)
    , showHistogram_(true), histogramBrightness_(2.f), showTexture_(true)
    //general
    , tfProp_(prop)
{
    setObjectName("TransFunc2DPrimitivesPropertyEditorCanvas");
    setMouseTracking(true);
    setFocusPolicy(Qt::StrongFocus);
    setFocus();

    // init the context menu
    initializeKeyContextMenu();
    initializePrimitiveContextMenu();
}

TransFunc2DPrimitivesPropertyEditorCanvas::~TransFunc2DPrimitivesPropertyEditorCanvas() {
    delete histogramCache_;
}

void TransFunc2DPrimitivesPropertyEditorCanvas::updateFromProperty() {
    tgtAssert(tfProp_ && tfProp_->get(), "no tf or tf property!");
    if(const VolumeBase* volume = tfProp_->getVolume()) {
        //get histogram
        if(volume->hasDerivedData<VolumeHistogramIntensityGradient>()) {
            Histogram2D* histogram = &volume->getDerivedData<VolumeHistogramIntensityGradient>()->getHistogram(tfProp_->getVolumeChannel());
            if(histogram_ != histogram) {
                histogram_ = histogram;
                delete histogramCache_;
                histogramCache_ = 0;
            }
        }
        else {
            if (tfProp_->getComputeHistogram())
                volume->getDerivedDataThreaded<VolumeHistogramIntensityGradient>();
            histogram_ = 0;
            delete histogramCache_;
            histogramCache_ = 0;
        }
    } else {
        histogram_ = 0;
        delete histogramCache_;
        histogramCache_ = 0;
    }
    if (visibleHistogramRangeX_ != tfProp_->get()->getDomain(0) || visibleHistogramRangeY_ != tfProp_->get()->getDomain(1)) {
        visibleHistogramRangeX_ = tfProp_->get()->getDomain(0); visibleHistogramRangeY_ = tfProp_->get()->getDomain(1);
        delete histogramCache_;
        histogramCache_ = 0;
    }
    update();
}

//-------------------------------------------------------------------------------------------------------------
//      Context menu
//-------------------------------------------------------------------------------------------------------------
void TransFunc2DPrimitivesPropertyEditorCanvas::initializeKeyContextMenu() {
    QAction* keyColorChangeAction = new QAction(tr("Change color of key"), this);
    keyColorChangeAction->setIcon(QIcon(":/qt/icons/colorize.png"));
    keyContextMenu_.addAction(keyColorChangeAction);
    connect(keyColorChangeAction, SIGNAL(triggered()), this, SLOT(keyColorChangeActionSlot()));

    keyContextMenu_.addSeparator();

    QAction* primitiveColorChangeAction = new QAction(tr("Change color of primitive"), this);
    primitiveColorChangeAction->setIcon(QIcon(":/qt/icons/colorize.png"));
    keyContextMenu_.addAction(primitiveColorChangeAction);
    connect(primitiveColorChangeAction, SIGNAL(triggered()), this, SLOT(primitiveColorChangeActionSlot()));

    QAction* deleteAction_ = new QAction(tr("Delete this primitive"), this);
    deleteAction_->setIcon(QIcon(":/qt/icons/eraser.png"));
    keyContextMenu_.addAction(deleteAction_);
    connect(deleteAction_, SIGNAL(triggered()), this, SLOT(deletePrimitiveSlot()));
}

void TransFunc2DPrimitivesPropertyEditorCanvas::initializePrimitiveContextMenu() {
    QAction* primitiveColorChangeAction = new QAction(tr("Change color of primitive"), this);
    primitiveColorChangeAction->setIcon(QIcon(":/qt/icons/colorize.png"));
    primitiveContextMenu_.addAction(primitiveColorChangeAction);
    connect(primitiveColorChangeAction, SIGNAL(triggered()), this, SLOT(primitiveColorChangeActionSlot()));

    QAction* deleteAction_ = new QAction(tr("Delete this primitive"), this);
    deleteAction_->setIcon(QIcon(":/qt/icons/eraser.png"));
    primitiveContextMenu_.addAction(deleteAction_);
    connect(deleteAction_, SIGNAL(triggered()), this, SLOT(deletePrimitiveSlot()));
}

void TransFunc2DPrimitivesPropertyEditorCanvas::keyColorChangeActionSlot() {
    if (!selectedKey_)
        return;

    QColor oldColor = QColor(selectedKey_->color_.r,selectedKey_->color_.g,selectedKey_->color_.b,selectedKey_->color_.a);

    QColor newColor = QColorDialog::getColor(oldColor);
    if (newColor.isValid())
        changeCurrentColorSlot(newColor);
}

void TransFunc2DPrimitivesPropertyEditorCanvas::primitiveColorChangeActionSlot() {
    if (!selectedPrimitive_)
        return;

    //remove selected key to force color update of entire primitive
    selectedKey_ = 0;

    QColor oldColor = QColor(selectedPrimitive_->getControlPoint(0).color_.r,selectedPrimitive_->getControlPoint(0).color_.g,
                             selectedPrimitive_->getControlPoint(0).color_.b,selectedPrimitive_->getControlPoint(0).color_.a);

    QColor newColor = QColorDialog::getColor(oldColor);
    if (newColor.isValid())
        changeCurrentColorSlot(newColor);
}


void TransFunc2DPrimitivesPropertyEditorCanvas::changeCurrentColorSlot(const QColor& c) {
    if ((!selectedKey_ && !selectedPrimitive_) || !c.isValid())
        return;

    tgt::col4 color(c.red(),c.green(),c.blue(),c.alpha());
    bool changedColor = false;
    if(selectedKey_) {
        if(selectedKey_->color_ != color) {
            selectedKey_->color_ = color;
            changedColor = true;
        }
    } else { //case primitive
        selectedPrimitive_->setColor(color);
        changedColor = true;
    }

    if (changedColor) {
        emit colorChangedSignal(c);
        tfProp_->get()->invalidateTexture();
        tfProp_->invalidate();
    }
}

void TransFunc2DPrimitivesPropertyEditorCanvas::deletePrimitiveSlot() {
    if(tfProp_->get()) {
        if(selectedPrimitive_) {
            tfProp_->get()->removePrimitive(selectedPrimitive_);
            selectedKey_ = 0;
            selectedPrimitive_ = 0;
            emit colorChangedSignal(QColor(0,0,0,0));
            emit fuzzinessChangedSignal(-1);
            tfProp_->invalidate();
        }
    }
}

void TransFunc2DPrimitivesPropertyEditorCanvas::updateFuzzinessSlot(int value) {
    if(selectedPrimitive_) {
        selectedPrimitive_->setFuzziness(static_cast<float>(value)/100.f);
        tfProp_->get()->invalidateTexture();
        tfProp_->invalidate();
    }
}

void TransFunc2DPrimitivesPropertyEditorCanvas::newPrimitiveAddedSlot(TransFuncPrimitive* prim) {
    selectedPrimitive_ = 0; selectedKey_ = 0;
    emit colorChangedSignal(QColor(prim->getControlPoint(0).color_.r,prim->getControlPoint(0).color_.g,
                                   prim->getControlPoint(0).color_.b,prim->getControlPoint(0).color_.a));
    emit fuzzinessChangedSignal(prim->getFuzziness()*100);
    //set selection after emit to remove color changes
    selectedPrimitive_ = prim;
    //not needed, since invalidate is called from the editor after emitting the signal.
    //update();
}

//-------------------------------------------------------------------------------------------------------------
//      Mouse events and other
//-------------------------------------------------------------------------------------------------------------
void TransFunc2DPrimitivesPropertyEditorCanvas::mousePressEvent(QMouseEvent* event) {
    //activate interaction mode
    if (event->button() == Qt::LeftButton)
        emit toggleInteractionModeSignal(true);

    event->accept();

    tgt::vec2 sHit = tgt::vec2(event->x(), static_cast<float>(height()) - event->y());
    tgt::vec2 hit = pixelToNormalized(sHit);

    // see if a key was selected
    selectedKey_ = 0; selectedPrimitive_ = 0;
    for(int p = 0; p < tfProp_->get()->getNumPrimitives(); ++p) {
       TransFuncPrimitive* prim = tfProp_->get()->getPrimitive(p);
        for (int c=0; c<prim->getNumControlPoints(); ++c) {
            const TransFuncPrimitiveControlPoint& key = prim->getControlPoint(c);
            tgt::vec2 sp = normalizedToPixel(key.position_);
            if (sHit.x > sp.x - KEY_POINT_SIZE && sHit.x < sp.x + KEY_POINT_SIZE &&
                sHit.y > sp.y - KEY_POINT_SIZE && sHit.y < sp.y + KEY_POINT_SIZE)
            {
                 emit colorChangedSignal(QColor(key.color_.r,key.color_.g,key.color_.b,key.color_.a));
                 emit fuzzinessChangedSignal(prim->getFuzziness()*100);
                 selectedKey_ = const_cast<TransFuncPrimitiveControlPoint*>(&key);
                 selectedPrimitive_ = prim;
                 break;
            }
        }
        //end loop, if key has been found
        if(selectedKey_)
            break;
    }

    //check, if a primitive is been hit
    if(selectedPrimitive_ == 0) {
        TransFuncPrimitive* selectedPrimitiveTmp = tfProp_->get()->getPrimitive(hit);
        if(selectedPrimitiveTmp) {
            emit fuzzinessChangedSignal(selectedPrimitiveTmp->getFuzziness()*100);
            TransFuncPrimitiveControlPoint key = selectedPrimitiveTmp->getControlPoint(0);
            emit colorChangedSignal(QColor(key.color_.r,key.color_.g,key.color_.b,key.color_.a));
            selectedPrimitive_ = selectedPrimitiveTmp;
        }
    }

    //handle right click
    if (event->button() == Qt::RightButton) {
        if (selectedKey_) {
            keyContextMenu_.popup(event->globalPos());
        } else if(selectedPrimitive_) {
            primitiveContextMenu_.popup(event->globalPos());
        } else {
            //reset color
            emit colorChangedSignal(QColor(0,0,0,0));
            emit fuzzinessChangedSignal(-1);
        }
        update();
        return;
    }

    if (selectedPrimitive_ != 0 && event->button() == Qt::LeftButton) {
        isBeingDraged_ = true;
        //keep values within valid range
        lastMousePos_ = event->pos();
        hit = tgt::clamp(hit, 0.f, 1.f);
        updateToolTipCoordinates(event->pos(), hit);
        update();
        return;
    }

    //nothing hit or unknown button
    emit colorChangedSignal(QColor(0,0,0,0));
    emit fuzzinessChangedSignal(-1);
    update();
}

void TransFunc2DPrimitivesPropertyEditorCanvas::mouseMoveEvent(QMouseEvent* event) {
    event->accept();

    tgt::vec2 pixelHit = tgt::vec2(event->x(), static_cast<float>(height()) - event->y());
    tgt::vec2 normalizedHit = pixelToNormalized(pixelHit);

    if (isBeingDraged_) {
        setCursor(Qt::ClosedHandCursor);
        if(selectedKey_) {
            selectedKey_->position_ = normalizedHit;
        } else if(selectedPrimitive_) {
            selectedPrimitive_->move(normalizedHit - pixelToNormalized(tgt::vec2(lastMousePos_.x(), static_cast<float>(height()) - lastMousePos_.y())));
        } else {
            tgtAssert(false, "should not get here!");
        }
        updateToolTipCoordinates(event->pos(), tgt::clamp(normalizedHit,0.f,1.f));
        tfProp_->get()->invalidateTexture();
        tfProp_->invalidate(); //update called from tf invalidate
        lastMousePos_ = event->pos();
    } else {
        //Nothing to do here
        unsetCursor();
    }
}

void TransFunc2DPrimitivesPropertyEditorCanvas::mouseReleaseEvent(QMouseEvent* event) {
    event->accept();
    if (event->button() == Qt::LeftButton) {
        isBeingDraged_ = false;
        QToolTip::hideText();
        update();
        emit toggleInteractionModeSignal(false);
    }
}

void TransFunc2DPrimitivesPropertyEditorCanvas::mouseDoubleClickEvent(QMouseEvent *event) {
    event->accept();
    if (event->button() == Qt::LeftButton) {
        if(selectedKey_)
            keyColorChangeActionSlot();
        else if (selectedPrimitive_)
            primitiveColorChangeActionSlot();
    }
}

/*void TransFunc2DPrimitivesPropertyEditorCanvas::keyPressEvent(QKeyEvent* event) {
    if (event->key() == Qt::Key_Shift                    && underMouse() &&
        isLineHit(tgt::vec2(mousePos_.x(), mousePos_.y())) >= 0 && !isKeyBeingDraged_)
    {
        setCursor(Qt::SizeVerCursor);
    }
}*/

void TransFunc2DPrimitivesPropertyEditorCanvas::keyReleaseEvent(QKeyEvent* event) {
    unsetCursor(); //probably not needed
    if (event->key() == Qt::Key_Delete && (selectedKey_ != 0 || selectedPrimitive_ != 0)) {
        event->accept();
        deletePrimitiveSlot();
    }
}

void TransFunc2DPrimitivesPropertyEditorCanvas::updateToolTipCoordinates(QPoint pos, tgt::vec2 values) {
    std::ostringstream os;
    os.precision(2);
    os.setf(std::ios::fixed, std::ios::floatfield);

    float intensity = values.x, gradient = values.y;
    if(tfProp_->get()) {
        intensity  = tfProp_->get()->getDomain(0).x + (tfProp_->get()->getDomain(0).y - tfProp_->get()->getDomain(0).x) * intensity;
        gradient   = tfProp_->get()->getDomain(1).x + (tfProp_->get()->getDomain(1).y - tfProp_->get()->getDomain(1).x) * gradient;
    }
    os << intensity << " / " << gradient; // intensity / gradient
    QToolTip::showText(mapToGlobal(pos), QString(os.str().c_str()));
}

//-------------------------------------------------------------------------------------------------------------
//      Paint Functions
//-------------------------------------------------------------------------------------------------------------
void TransFunc2DPrimitivesPropertyEditorCanvas::paintEvent(QPaintEvent* event) {
    tgtAssert(tfProp_ && tfProp_->get(), "no tf or property!");

    //the histogram is automatically painted onto this widget
    //we do not need to call the paintevent for the Histogrampainter directly
    event->accept();

    QPainter paint(this);

    // put origin in lower lefthand corner
    QMatrix m;
    m.translate(0.0, static_cast<float>(height())-1);
    m.scale(1.f, -1.f);
    paint.setMatrix(m);
    //draw white field
    paint.setMatrixEnabled(true);
    paint.setRenderHint(QPainter::Antialiasing, false);
    paint.setPen(Qt::NoPen);
    paint.setBrush(Qt::white);
    paint.drawRect(0, 0, width() - 1, height() - 1);

    //draw the grid
    drawCheckBoard(&paint);
    //draw tf
    if(showTexture_)
        drawTransFunc(&paint);
    //draw histogram
    if(showHistogram_)
        drawHistogram(&paint);
    //draw axis
    drawAxes(&paint);
    //draw threshold
    drawThreshold(&paint);
    //draw mapping keys
    drawMappingKeys(&paint);

    paint.setRenderHint(QPainter::Antialiasing, false);
    paint.setPen(Qt::lightGray);
    paint.setBrush(Qt::NoBrush);
    paint.drawRect(0, 0, width() - 1, height() - 1);
}

void TransFunc2DPrimitivesPropertyEditorCanvas::drawCheckBoard(QPainter* painter) {
    // draw grid
    painter->setPen(QColor(220, 220, 220));
    painter->setBrush(QColor(220, 220, 220));
    painter->setRenderHint(QPainter::Antialiasing, false);

    tgt::vec2 minCoord = normalizedToPixel(tgt::vec2(0.f,0.f));
    tgt::vec2 maxCoord = normalizedToPixel(tgt::vec2(1.f,1.f));

    if(showTexture_) {
        int cX = 0, cY = 0;
        tgt::vec2 size = normalizedToPixel(tgt::vec2(gridSpacing_.x, gridSpacing_.y)) - minCoord;
        for (float x = 0.f; x < 1.f - gridSpacing_.x*0.5; x += gridSpacing_.x) {
            for (float y = 0.f; y < 1.f - gridSpacing_.y*0.5; y +=gridSpacing_.y) {
                if((cX+cY) % 2 == 0) {
                    tgt::vec2 start = normalizedToPixel(tgt::vec2(x, y));
                    tgt::vec2 tmp = normalizedToPixel(tgt::vec2(x+gridSpacing_.x, y+gridSpacing_.y));
                    painter->drawRect(start.x, start.y,size.x,size.y);
                }
                cY++;
            }
            cY = 0;
            cX++;
        }
        painter->drawLine(maxCoord.x,minCoord.y,maxCoord.x,maxCoord.y);
        painter->drawLine(minCoord.x,maxCoord.y,maxCoord.x,maxCoord.y);
    } else {
        for (float x = gridSpacing_.x; x < 1.f + gridSpacing_.x*0.5; x += gridSpacing_.x) {
            tgt::vec2 tmp = normalizedToPixel(tgt::vec2(x,0.f));
            painter->drawLine(tmp.x,minCoord.y,tmp.x,maxCoord.y);
        }
        for (float y = gridSpacing_.y; y < 1.f + gridSpacing_.y*0.5; y += gridSpacing_.y) {
            tgt::vec2 tmp = normalizedToPixel(tgt::vec2(0.f,y));
            painter->drawLine(minCoord.x,tmp.y,maxCoord.x,tmp.y);
        }
    }
}

void TransFunc2DPrimitivesPropertyEditorCanvas::drawTransFunc(QPainter* painter) {
    QPixmap icon(256,256);
    icon.fill(Qt::transparent);
    //icon.set
    QPainter pmPainter(&icon);
    for(int x = 0; x < 256; x++) {
        for(int y = 0; y < 256; y++) {
            tgt::Color color = tfProp_->get()->getTexture()->texelAsFloat(tgt::svec2(static_cast<size_t>((static_cast<float>(x)/255.f)*(tfProp_->get()->getDimensions().x-1)),
                                                                                     static_cast<size_t>((static_cast<float>(y)/255.f)*(tfProp_->get()->getDimensions().y-1))));
            pmPainter.setPen(QColor(color.r*255,color.g*255,color.b*255,color.a*255));
            pmPainter.drawPoint(x,y);
        }
    }

    tgt::vec2 start = normalizedToPixel(tgt::vec2::zero);
    tgt::vec2 end = normalizedToPixel(tgt::vec2::one);
    painter->drawPixmap(start.x,start.y,end.x-start.x,end.y-start.y, icon);
}

void TransFunc2DPrimitivesPropertyEditorCanvas::drawAxes(QPainter* painter) {
    // draw x and y axes
    painter->setRenderHint(QPainter::Antialiasing, true);
    painter->setPen(Qt::gray);
    painter->setBrush(Qt::gray);

    // draw axes independently from visible range
    tgt::vec2 origin = normalizedToPixel(tgt::vec2(0.f, 0.f));
    origin.x = floor(origin.x) + 0.5f;
    origin.y = floor(origin.y) - 0.5f;

    painter->setRenderHint(QPainter::Antialiasing, true);

    painter->drawLine(QPointF(AXIS_OFFSET, origin.y),
                   QPointF(width() - AXIS_OFFSET, origin.y));

    painter->drawLine(QPointF(origin.x, AXIS_OFFSET),
                   QPointF(origin.x, height() - AXIS_OFFSET));

    QPointF arrow[3];
    arrow[0] = QPointF(origin.x, height() - AXIS_OFFSET);
    arrow[1] = QPointF(origin.x + ARROW_WIDTH, height() - AXIS_OFFSET - ARROW_LENGTH);
    arrow[2] = QPointF(origin.x - ARROW_WIDTH, height() - AXIS_OFFSET - ARROW_LENGTH);

    painter->drawConvexPolygon(arrow, 3);

    arrow[0] = QPointF(width() - AXIS_OFFSET, origin.y);
    arrow[1] = QPointF(width() - AXIS_OFFSET - ARROW_LENGTH, origin.y - ARROW_WIDTH);
    arrow[2] = QPointF(width() - AXIS_OFFSET - ARROW_LENGTH, origin.y + ARROW_WIDTH);

    painter->drawConvexPolygon(arrow, 3);

    painter->scale(-1.f, 1.f);
    painter->rotate(180.f);
    painter->drawText(static_cast<int>(width() - painter->fontMetrics().width(X_AXIS_TEXT) - 2.78f * AXIS_OFFSET), static_cast<int>(-1 * (origin.y+1.f - 0.8f * AXIS_OFFSET)), X_AXIS_TEXT);
    painter->drawText(static_cast<int>(1.6f * AXIS_OFFSET), static_cast<int>(-1 * (height() - 1.85f * AXIS_OFFSET)), Y_AXIS_TEXT);
    painter->rotate(180.f);
    painter->scale(-1.f, 1.f);
}

void TransFunc2DPrimitivesPropertyEditorCanvas::drawHistogram(QPainter* painter) {
    if (!histogram_) {
        painter->setMatrixEnabled(false);
        painter->setPen(Qt::red);
        painter->drawText(QRectF(0, 7, width() - 1, height() - 8), tr("No volume or calculating histogram"), QTextOption(Qt::AlignHCenter));
        painter->setMatrixEnabled(true);
    } else {
        int numBucketsX = histogram_->getNumBuckets(0);
        int numBucketsY = histogram_->getNumBuckets(1);

        if (histogramCache_ == 0 || histogramCache_->rect() != rect()) {
            int maxBucket = histogram_->getMaxBucket();
            delete histogramCache_;
            histogramCache_ = new QPixmap(numBucketsX,numBucketsY);
            histogramCache_->fill(Qt::transparent);

            QPainter pmPainter(histogramCache_);
            for(int x = 0; x < 256; x++) {
                for(int y = 0; y < 256; y++) {
                    float value = logf(static_cast<float>(1+histogram_->getBucket(x,y))) / (logf(1+maxBucket)+0.0001f); // to not divide by zero
                    pmPainter.setPen(QColor(255, 135, 135,tgt::clamp(static_cast<int>(value*255.f*histogramBrightness_),0,255)));
                    pmPainter.drawPoint(x,y);
                }
            }

        }
        tgt::vec2 start = normalizedToPixel(tgt::vec2::zero), end = normalizedToPixel(tgt::vec2::one);
        float widthInPixels = end.x-start.x, heightInPixels = end.y-start.y;

        float bucketStepSizeX = (histogram_->getMaxValue(0) - histogram_->getMinValue(0)) / numBucketsX;
        float pixelStepSizeX =  (tfProp_->get()->getDomain(0).y - tfProp_->get()->getDomain(0).x) /widthInPixels;
        float bucketStepSizeY = (histogram_->getMaxValue(1) - histogram_->getMinValue(1)) / numBucketsY;
        float pixelStepSizeY =  (tfProp_->get()->getDomain(1).y - tfProp_->get()->getDomain(1).x) /heightInPixels;

        float mapStartX, mapSizeX, histoStartX, histoSizeX;
        if(tfProp_->get()->getDomain(0).x >= histogram_->getMinValue(0)) {
            mapStartX = start.x;
            histoStartX = (tfProp_->get()->getDomain(0).x - histogram_->getMinValue(0)) / bucketStepSizeX;
        } else {
            mapStartX = (histogram_->getMinValue(0) - tfProp_->get()->getDomain(0).x) / pixelStepSizeX;
            histoStartX = 0;
        }
        if(tfProp_->get()->getDomain(0).y <= histogram_->getMaxValue(0)) {
            mapSizeX = end.x - mapStartX;
            histoSizeX = (tfProp_->get()->getDomain(0).y - histogram_->getMinValue(0)) / bucketStepSizeX;
        } else {
            mapSizeX = (histogram_->getMaxValue(0) - tfProp_->get()->getDomain(0).x) / pixelStepSizeX - mapStartX;
            histoSizeX = numBucketsX - histoStartX;
        }

        float mapStartY, mapSizeY, histoStartY, histoSizeY;
        if(tfProp_->get()->getDomain(1).x >= histogram_->getMinValue(1)) {
            mapStartY = start.y;
            histoStartY = (tfProp_->get()->getDomain(1).x - histogram_->getMinValue(1)) / bucketStepSizeY;
        } else {
            mapStartY = (histogram_->getMinValue(1) - tfProp_->get()->getDomain(1).x) / pixelStepSizeY;
            histoStartY = 0;
        }
        if(tfProp_->get()->getDomain(1).y <= histogram_->getMaxValue(1)) {
            mapSizeY = end.y - mapStartY;
            histoSizeY = (tfProp_->get()->getDomain(1).y - histogram_->getMinValue(1)) / bucketStepSizeY;
        } else {
            mapSizeY = (histogram_->getMaxValue(1) - tfProp_->get()->getDomain(1).x) / pixelStepSizeY - mapStartY;
            histoSizeY = numBucketsY - histoStartY;
        }

        painter->drawPixmap(static_cast<int>(mapStartX),static_cast<int>(mapStartY),static_cast<int>(mapSizeX),static_cast<int>(mapSizeY),*histogramCache_,
                 static_cast<int>(histoStartX),static_cast<int>(histoStartY),static_cast<int>(histoSizeX),static_cast<int>(histoSizeY));
    }
}

void TransFunc2DPrimitivesPropertyEditorCanvas::drawMappingKeys(QPainter* painter) {
    //draw lines
        //setLine Color
        QPen pen = QPen(Qt::darkRed);
        pen.setWidthF(1.5f);
        painter->setPen(pen);
        painter->setBrush(Qt::NoBrush);

    for (int i = 0; i < tfProp_->get()->getNumPrimitives(); ++i) {
        const TransFuncPrimitive* prim = tfProp_->get()->getPrimitive(i);
        pen.setWidthF((prim == selectedPrimitive_) && (selectedKey_ == 0) ? 3.f : 1.5f);
        painter->setPen(pen);
        if(dynamic_cast<const TransFuncQuad*>(prim) || dynamic_cast<const TransFuncTriangle*>(prim)) { //quad or triangle
            tgt::vec2 oldPos = normalizedToPixel(prim->getControlPoint(prim->getNumControlPoints()-1).position_);
            for(int j = 0; j < prim->getNumControlPoints(); j++) {
                tgt::vec2 newPos = normalizedToPixel(prim->getControlPoint(j).position_);
                painter->drawLine(QPointF(oldPos.x /*+ 1.f*/, oldPos.y),
                           QPointF(newPos.x /*- 1.f*/, newPos.y));
                oldPos = newPos;
            }
        } else if(dynamic_cast<const TransFuncBanana*>(prim)) { //banana
            tgt::vec2 p0 = normalizedToPixel(prim->getControlPoint(0).position_);
            tgt::vec2 p1 = normalizedToPixel((2.f * prim->getControlPoint(1).position_) - (0.5f * prim->getControlPoint(0).position_) - (0.5f * prim->getControlPoint(3).position_));
            tgt::vec2 p2 = normalizedToPixel( (2.f * prim->getControlPoint(2).position_) - (0.5f * prim->getControlPoint(0).position_) - (0.5f * prim->getControlPoint(3).position_));
            tgt::vec2 p3 = normalizedToPixel(prim->getControlPoint(3).position_);
            QPainterPath path(QPointF(p0.x,p0.y));
            path.quadTo(p1.x,p1.y,p3.x,p3.y);
            path.quadTo(p2.x,p2.y,p0.x,p0.y);
            painter->drawPath(path);
        } else {
            tgtAssert(false,"Unexpected primitive!");
            LWARNINGC("TransFunc2DPrimitivesPropertyEditorCanvas","Unexpected Primitive!");
        }
    }

    //draw keys
    bool highlight = false;
    for (int i = 0; i < tfProp_->get()->getNumPrimitives(); ++i) {
        const TransFuncPrimitive* prim = tfProp_->get()->getPrimitive(i);
        highlight = (prim == selectedPrimitive_) && (selectedKey_ == 0);
        for(int j = 0; j < prim->getNumControlPoints(); j++) {
            const TransFuncPrimitiveControlPoint& point = prim->getControlPoint(j);
            drawControlPoint(*painter, point.color_, normalizedToPixel(point.position_), highlight || (&point == selectedKey_));
        }
    }
}

void TransFunc2DPrimitivesPropertyEditorCanvas::drawControlPoint(QPainter& paint, const tgt::col4& tgtcolor, const tgt::vec2& p, bool selected) {
    paint.setBrush(QColor(tgtcolor.r,tgtcolor.g,tgtcolor.b,tgtcolor.a));
    QPen pen(QBrush(Qt::darkGray), Qt::SolidLine);
    if (selected)
        pen.setWidth(3); //default 1
    paint.setPen(pen);
    paint.drawEllipse(QRectF(p.x - KEY_POINT_SIZE/2, p.y - KEY_POINT_SIZE/2,
                                 KEY_POINT_SIZE, KEY_POINT_SIZE));
}

void TransFunc2DPrimitivesPropertyEditorCanvas::drawThreshold(QPainter* painter) {
    // grey out threshold area
    tgt::vec2 origin = normalizedToPixel(tgt::vec2(0.f));
    painter->setBrush(QBrush(QColor(192, 192, 192, 230), Qt::SolidPattern));
    painter->setPen(Qt::NoPen);
    tgt::vec2 upperRight = normalizedToPixel(tgt::vec2(1.f));
    tgt::vec2 lowerLeft = normalizedToPixel(tgt::vec2(0.f));
    int w = static_cast<int>(upperRight.x - lowerLeft.x);
    int h = static_cast<int>(upperRight.y - lowerLeft.y);
    tgt::ivec2 thresholdXVisible(tfProp_->get()->getThreshold(0).x * w, tfProp_->get()->getThreshold(0).y * w);
    tgt::ivec2 thresholdYVisible(tfProp_->get()->getThreshold(1).x * h, tfProp_->get()->getThreshold(1).y * h);

    //draw x
    if (tfProp_->get()->getThreshold(0).x > 0.f) {
        painter->drawRect(static_cast<int>(origin.x), static_cast<int>(origin.y),
                       thresholdXVisible.x, h);
    }
    if (tfProp_->get()->getThreshold(0).y < 1.f) {
        painter->drawRect(static_cast<int>(origin.x + thresholdXVisible.y), static_cast<int>(origin.y),
            w - thresholdXVisible.y, h);
    }
    //draw y
    if (tfProp_->get()->getThreshold(1).x > 0.f) {
        painter->drawRect(static_cast<int>(origin.x + thresholdXVisible.x), static_cast<int>(origin.y),
                          thresholdXVisible.y - thresholdXVisible.x, thresholdYVisible.x);
    }
    if (tfProp_->get()->getThreshold(1).y < 1.f) {
        painter->drawRect(static_cast<int>(origin.x + thresholdXVisible.x),static_cast<int>(origin.y + thresholdYVisible.y),
                           thresholdXVisible.y - thresholdXVisible.x, h - thresholdYVisible.y);
    }

}

void TransFunc2DPrimitivesPropertyEditorCanvas::toggleHistogram(bool state) {
    showHistogram_ = state;
    update();
}

void TransFunc2DPrimitivesPropertyEditorCanvas::setHistogramBrightness(int value) {
    float brightness = static_cast<float>(value)/10.f;

    if(histogramBrightness_ != brightness) {
        histogramBrightness_ = brightness;
        delete histogramCache_; histogramCache_ = 0;
        update();
    }
}

void TransFunc2DPrimitivesPropertyEditorCanvas::toggleTexture(bool state) {
    showTexture_ = state;
    update();
}


//-------------------------------------------------------------------------------------------------------------
//      Helper Functions
//-------------------------------------------------------------------------------------------------------------
tgt::vec2 TransFunc2DPrimitivesPropertyEditorCanvas::normalizedToPixel(tgt::vec2 p) {
    float sx = p.x * (static_cast<float>(width())  - 2.f * AXIS_OFFSET - 1.5f * ARROW_LENGTH) + AXIS_OFFSET;
    float sy = p.y * (static_cast<float>(height()) - 2.f * AXIS_OFFSET - 1.5f * ARROW_LENGTH) + AXIS_OFFSET;
    return tgt::vec2(sx, sy);
}

tgt::vec2 TransFunc2DPrimitivesPropertyEditorCanvas::pixelToNormalized(tgt::vec2 p) {
    float wx = (p.x - AXIS_OFFSET) / (static_cast<float>(width())  - 2.f * AXIS_OFFSET - 1.5f * ARROW_LENGTH);
    float wy = (p.y - AXIS_OFFSET) / (static_cast<float>(height()) - 2.f * AXIS_OFFSET - 1.5f * ARROW_LENGTH);
    return tgt::vec2(wx, wy);
}

} // namespace voreen
