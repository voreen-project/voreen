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

#include "voreen/qt/widgets/property/transfunc/1d/1dkeys/transfunc1dkeyspropertyeditorcanvas.h"

#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/datastructures/volume/volumedisk.h"
#include "voreen/core/datastructures/volume/volumedecorator.h"
#include "voreen/core/datastructures/volume/histogram.h"
#include "voreen/core/datastructures/transfunc/1d/1dkeys/transfunc1dkeys.h"
#include "voreen/core/datastructures/transfunc/1d/1dkeys/utils/transfuncmappingkey.h"

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
    const QString Y_AXIS_TEXT = QString("opacity");

    const int AXIS_OFFSET  = 12;
    const int ARROW_LENGTH = 10;
    const int ARROW_WIDTH  = 3;
    const int KEY_POINT_SIZE = 10;
    const float KEY_SPLIT_FACTOR = 1.5f;
    const int MIN_CELL_SIZE = 8;
}

namespace voreen {
TransFunc1DKeysPropertyEditorCanvas::TransFunc1DKeysPropertyEditorCanvas(QWidget* parent, TransFunc1DKeysProperty* prop)
    : QWidget(parent)
    //context
    , splitMergeAction_(0), zeroAction_(0), deleteAction_(0)
    //interaction
    , selectedKey_(0), isLeftPartSelected_(true), isKeyBeingDraged_(false)
    , dragLine_(-1), dragLineAlphaLeft_(-1.f), dragLineAlphaRight_(-1.f)
    //paint
    , gridSpacing_(0.1f,0.1f), histogram_(0), histogramCache_(0), visibleHistogramRange_(0.f,1.f)
    , showHistogram_(true), showTexture_(true)
    //general
    , tfProp_(prop)
{
    setObjectName("TransFunc1DKeysPropertyEditorCanvas");
    setMouseTracking(true);
    setFocusPolicy(Qt::StrongFocus);
    setFocus();

    // init the context menu
    initializeKeyContextMenu();
}

TransFunc1DKeysPropertyEditorCanvas::~TransFunc1DKeysPropertyEditorCanvas() {
    delete histogramCache_;
}

void TransFunc1DKeysPropertyEditorCanvas::updateFromProperty() {
    tgtAssert(tfProp_ && tfProp_->get(), "no tf or tf property!");
    if(const VolumeBase* volume = tfProp_->getVolume()) {
        //get histogram
        if(volume->hasDerivedData<VolumeHistogramIntensity>()) {
            Histogram1D* histogram = &volume->getDerivedData<VolumeHistogramIntensity>()->getHistogram(tfProp_->getVolumeChannel());
            if(histogram_ != histogram) {
                histogram_ = histogram;
                delete histogramCache_;
                histogramCache_ = 0;
            }
        }
        else {
            if (tfProp_->getComputeHistogram())
                volume->getDerivedDataThreaded<VolumeHistogramIntensity>();
            histogram_ = 0;
            delete histogramCache_;
            histogramCache_ = 0;
        }
    } else {
        histogram_ = 0;
        delete histogramCache_;
        histogramCache_ = 0;
    }
    if (visibleHistogramRange_ != tfProp_->get()->getDomain()) {
        visibleHistogramRange_ = tfProp_->get()->getDomain();
        delete histogramCache_;
        histogramCache_ = 0;
    }
    update();
}

//-------------------------------------------------------------------------------------------------------------
//      Context menu
//-------------------------------------------------------------------------------------------------------------
void TransFunc1DKeysPropertyEditorCanvas::initializeKeyContextMenu() {
    QAction* colorChangeAction = new QAction(tr("Change color of key"), this);
    colorChangeAction->setIcon(QIcon(":/qt/icons/colorize.png"));
    keyContextMenu_.addAction(colorChangeAction);
    connect(colorChangeAction, SIGNAL(triggered()), this, SLOT(colorChangeActionSlot()));

    keyContextMenu_.addSeparator();

    splitMergeAction_ = new QAction(tr(""), this); // Text will be set later
    keyContextMenu_.addAction(splitMergeAction_);
    connect(splitMergeAction_, SIGNAL(triggered()), this, SLOT(splitMergeKeysSlot()));

    zeroAction_ = new QAction("", this); // Text will be set later
    keyContextMenu_.addAction(zeroAction_);
    connect(zeroAction_, SIGNAL(triggered()), this, SLOT(zeroKeySlot()));

    keyContextMenu_.addSeparator();

    deleteAction_ = new QAction(tr("Delete this key"), this);
    deleteAction_->setIcon(QIcon(":/qt/icons/eraser.png"));
    keyContextMenu_.addAction(deleteAction_);
    connect(deleteAction_, SIGNAL(triggered()), this, SLOT(deleteKeySlot()));
}

void TransFunc1DKeysPropertyEditorCanvas::showKeyContextMenu(QMouseEvent* event) {
    tgtAssert(tfProp_ && tfProp_->get(), "no tf or tf property!");
    // Set context-dependent text for menu items

    // Split/merge
    QString splitMergeText;
    if (selectedKey_->isSplit())
        splitMergeText = tr("Merge this key");
    else
        splitMergeText = tr("Split this key");
    splitMergeAction_->setText(splitMergeText);

    // Zero/unzero
    QString zeroText;
    if (isLeftPartSelected_)
        zeroText = tr("Zero to the left");
    else
        zeroText = tr("Zero to the right");
    zeroAction_->setText(zeroText);

    // allow deletion of keys only if there are at least two of them
    deleteAction_->setEnabled(tfProp_->get()->getNumKeys() > 1);

    keyContextMenu_.popup(event->globalPos());
}

void TransFunc1DKeysPropertyEditorCanvas::colorChangeActionSlot() {
    if (!selectedKey_)
        return;

    QColor oldColor;
    if (selectedKey_->isSplit() && !isLeftPartSelected_)
        oldColor = Col2QColor( selectedKey_->getColorR() );
    else
        oldColor = Col2QColor( selectedKey_->getColorL() );

    QColor newColor = QColorDialog::getColor(oldColor, 0);
    if (newColor.isValid())
        changeCurrentColorSlot(newColor);
}

void TransFunc1DKeysPropertyEditorCanvas::changeCurrentColorSlot(const QColor& c) {
    if (!selectedKey_ || !c.isValid())
        return;

    tgt::col4 tgtcolor = QColor2Col(c);
    bool changedColor = false;
    if (selectedKey_->isSplit() && !isLeftPartSelected_) {
        tgtcolor.a = selectedKey_->getColorR().a;
        if (selectedKey_->getColorR() != tgtcolor) {
            selectedKey_->setColorR(tgtcolor);
            changedColor = true;
        }
    }
    else {
        tgtcolor.a = selectedKey_->getColorL().a;
        if (selectedKey_->getColorL() != tgtcolor) {
            selectedKey_->setColorL(tgtcolor);
            changedColor = true;
        }
    }

    if (changedColor) {
        tfProp_->get()->updateKey(selectedKey_);
        tfProp_->invalidate();
        emit colorChangedSignal(c);
    }
}

void TransFunc1DKeysPropertyEditorCanvas::splitMergeKeysSlot() {
    if (!selectedKey_)
        return;

    selectedKey_->setSplit(!selectedKey_->isSplit(), isLeftPartSelected_);
    tfProp_->invalidate();
}

void TransFunc1DKeysPropertyEditorCanvas::zeroKeySlot() {
    if (!selectedKey_)
        return;

    TransFuncMappingKey* otherKey = getOtherKey(selectedKey_, isLeftPartSelected_);
    if (otherKey) {
        if (!otherKey->isSplit())
            otherKey->setSplit(true);
        if (isLeftPartSelected_)
            otherKey->setAlphaR(0.0);
        else
            otherKey->setAlphaL(0.0);
    }

    if (!selectedKey_->isSplit())
        selectedKey_->setSplit(true);

    if (isLeftPartSelected_)
        selectedKey_->setAlphaL(0.0);
    else
        selectedKey_->setAlphaR(0.0);

    tfProp_->invalidate();
}

void TransFunc1DKeysPropertyEditorCanvas::deleteKeySlot() {
    tgtAssert(tfProp_ && tfProp_->get(), "no tf or tf property!");

    // allow deletion of keys only if there are at least two of them
    if (!selectedKey_ || tfProp_->get()->getNumKeys() < 2)
        return;

    tfProp_->get()->removeKey(selectedKey_);
    selectedKey_ = 0;
    emit colorChangedSignal(QColor(0,0,0,0));

    tfProp_->invalidate();
}

//-------------------------------------------------------------------------------------------------------------
//      Mouse events and other
//-------------------------------------------------------------------------------------------------------------
void TransFunc1DKeysPropertyEditorCanvas::mousePressEvent(QMouseEvent* event) {
    //activate interaction mode
    if (event->button() == Qt::LeftButton)
        emit toggleInteractionModeSignal(true);

    event->accept();

    //check, if a line is been hit
    dragLine_ = isLineHit(tgt::vec2(event->x(), event->y()));
    if (dragLine_ >= 0 && event->modifiers() == Qt::ShiftModifier) {
        dragLineStartY_ = event->y();
        return;
    }

    tgt::vec2 sHit = tgt::vec2(event->x(), static_cast<float>(height()) - event->y());
    tgt::vec2 hit = stow(sHit);

    // see if a key was selected
    selectedKey_ = 0;
    for (int i=0; i<tfProp_->get()->getNumKeys(); ++i) {
        TransFuncMappingKey* key = tfProp_->get()->getKey(i);
        tgt::vec2 sp = wtos(tgt::vec2(key->getIntensity(), key->getColorL().a / 255.0));
        tgt::vec2 spr = wtos(tgt::vec2(key->getIntensity(), key->getColorR().a / 255.0));
        if (key->isSplit()) {
            if (sHit.x > sp.x - KEY_SPLIT_FACTOR * KEY_POINT_SIZE && sHit.x <= sp.x &&
                sHit.y > sp.y - KEY_POINT_SIZE && sHit.y < sp.y + KEY_POINT_SIZE)
            {
                selectedKey_ = key;
                isLeftPartSelected_ = true;
            }
            if (sHit.x >= spr.x && sHit.x < spr.x + KEY_SPLIT_FACTOR * KEY_POINT_SIZE &&
                sHit.y > spr.y - KEY_POINT_SIZE && sHit.y < spr.y + KEY_POINT_SIZE)
            {
                selectedKey_ = key;
                isLeftPartSelected_ = false;
            }
        }
        else {
            if (sHit.x > sp.x - KEY_POINT_SIZE && sHit.x < sp.x + KEY_POINT_SIZE &&
                sHit.y > sp.y - KEY_POINT_SIZE && sHit.y < sp.y + KEY_POINT_SIZE)
            {
                selectedKey_ = key;
                isLeftPartSelected_ = false;
            }
        }
    }


    if (event->button() == Qt::RightButton) {
        if (selectedKey_) {
            showKeyContextMenu(event);
            emit colorChangedSignal(Col2QColor(selectedKey_->getColorL()));
        } else {
            emit colorChangedSignal(QColor(0,0,0,0));
        }
        update();
        return;
    }

    if (selectedKey_ != 0 && event->button() == Qt::LeftButton) {
        isKeyBeingDraged_ = true;
        //keep values within valid range
        hit = tgt::clamp(hit, 0.f, 1.f);
        updateToolTipCoordinates(event->pos(), hit);
        if (selectedKey_->isSplit() && !isLeftPartSelected_)
            emit colorChangedSignal(Col2QColor(selectedKey_->getColorR()));
        else
            emit colorChangedSignal(Col2QColor(selectedKey_->getColorL()));
        update();
        return;
    }

    // no key was selected -> insert new key
    if (hit.x >= 0.f && hit.x <= 1.f &&
        hit.y >= 0.f && hit.y <= 1.f &&
        event->button() == Qt::LeftButton)
    {
        insertNewKey(hit); //calls invalidate
        isKeyBeingDraged_ = true;
        dragLine_ = -1;
        updateToolTipCoordinates(event->pos(), hit);
        update();
        emit colorChangedSignal(Col2QColor(selectedKey_->getColorL()));

    }
}

void TransFunc1DKeysPropertyEditorCanvas::mouseMoveEvent(QMouseEvent* event) {
    event->accept();
    mousePos_ = event->pos();

    tgt::vec2 sHit = tgt::vec2(event->x(), static_cast<float>(height()) - event->y());
    tgt::vec2 hit = stow(sHit);


    if (!isKeyBeingDraged_ && isLineHit(tgt::vec2(event->x(), event->y())) >= 0 && event->modifiers() == Qt::ShiftModifier)
        setCursor(Qt::SizeVerCursor);
    else
        unsetCursor();

    if (dragLine_ >= 0) {
        // a line between 2 keys is moved (shift modifier was used)
        float delta = dragLineStartY_ - event->y();
        dragLineStartY_ = event->y();
        //left key
        TransFuncMappingKey* key = tfProp_->get()->getKey(dragLine_);
        if (dragLineAlphaLeft_ == -1.f)
            dragLineAlphaLeft_ = key->isSplit() ? key->getAlphaR() : key->getAlphaL();
        dragLineAlphaLeft_ = wtos(tgt::vec2(dragLineAlphaLeft_)).y;
        dragLineAlphaLeft_ += delta;
        dragLineAlphaLeft_ = stow(tgt::vec2(dragLineAlphaLeft_)).y;
        if (dragLineAlphaLeft_ < 0.f)
            dragLineAlphaLeft_ = 0.f;
        if (dragLineAlphaLeft_ > 1.f)
            dragLineAlphaLeft_ = 1.f;
        key->setAlphaR(dragLineAlphaLeft_);
        tfProp_->get()->updateKey(key);
        if (tfProp_->get()->getNumKeys() >= dragLine_+1) {
            //right key - when existing
            key = tfProp_->get()->getKey(dragLine_+1);
            if (dragLineAlphaRight_ == -1.f)
                dragLineAlphaRight_ = key->getAlphaL();
            dragLineAlphaRight_ = wtos(tgt::vec2(dragLineAlphaRight_)).y;
            dragLineAlphaRight_ += delta;
            dragLineAlphaRight_ = stow(tgt::vec2(dragLineAlphaRight_)).y;
            if (dragLineAlphaRight_ < 0.f)
                dragLineAlphaRight_ = 0.f;
            if (dragLineAlphaRight_ > 1.f)
                dragLineAlphaRight_ = 1.f;
            key->setAlphaL(dragLineAlphaRight_);
            tfProp_->get()->updateKey(key);
        }
        update();
        tfProp_->invalidate();
        return;
    }

    // return when no key was inserted or selected
    if (!isKeyBeingDraged_)
        return;

    // keep location within valid texture coord range
    hit = tgt::clamp(hit, 0.f, 1.f);

    if (selectedKey_ != 0) {
        updateToolTipCoordinates(event->pos(), hit);
        if (event->modifiers() != Qt::ShiftModifier) {
            selectedKey_->setIntensity(hit.x);
        }
        if (event->modifiers() != Qt::ControlModifier) {
            if (selectedKey_->isSplit()) {
                if (isLeftPartSelected_)
                    selectedKey_->setAlphaL(hit.y);
                else
                    selectedKey_->setAlphaR(hit.y);
            }
            else
                selectedKey_->setAlphaL(hit.y);
        }
        bool selectedFound = false;
        for (size_t i = 0; i < tfProp_->get()->getKeys().size(); ++i) {
            TransFuncMappingKey* key = tfProp_->get()->getKey(static_cast<int>(i));
            //is the tf key the selected one?
            if (key == selectedKey_) {
                selectedFound = true;
                continue;
            }
            if (selectedFound) {
                //change intensity of key if its lower than the intensity of selectedKey_
                if (key->getIntensity() < selectedKey_->getIntensity())
                    key->setIntensity(selectedKey_->getIntensity());
            }
            else {
                //change intensity of key if its higher than the intensity of selectedKey_
                if (key->getIntensity() > selectedKey_->getIntensity())
                    key->setIntensity(selectedKey_->getIntensity());
            }
        }
        tfProp_->get()->updateKey(selectedKey_);

        update();
        tfProp_->invalidate();
    }
}

void TransFunc1DKeysPropertyEditorCanvas::mouseReleaseEvent(QMouseEvent* event) {
    event->accept();
    if (event->button() == Qt::LeftButton) {
        isKeyBeingDraged_ = false;
        dragLine_ = -1;
        dragLineAlphaLeft_ = -1.f;
        dragLineAlphaRight_ = -1.f;
        QToolTip::hideText();
        update();
        emit toggleInteractionModeSignal(false);
    }
}

void TransFunc1DKeysPropertyEditorCanvas::mouseDoubleClickEvent(QMouseEvent *event) {
    event->accept();
    if (event->button() == Qt::LeftButton)
        colorChangeActionSlot();
}

void TransFunc1DKeysPropertyEditorCanvas::keyPressEvent(QKeyEvent* event) {
    if (event->key() == Qt::Key_Shift                    && underMouse() &&
        isLineHit(tgt::vec2(mousePos_.x(), mousePos_.y())) >= 0 && !isKeyBeingDraged_)
    {
        setCursor(Qt::SizeVerCursor);
    }
}

void TransFunc1DKeysPropertyEditorCanvas::keyReleaseEvent(QKeyEvent* event) {
    unsetCursor();
    if (event->key() == Qt::Key_Delete && selectedKey_ != 0) {
        event->accept();
        deleteKeySlot();
    }
}

int TransFunc1DKeysPropertyEditorCanvas::isLineHit(const tgt::vec2& pixelCoordinates) {
    tgtAssert(tfProp_ && tfProp_->get(), "no tf or property!");

    int hit = -1;
    tgt::vec2 sHit = tgt::vec2(pixelCoordinates.x, static_cast<float>(height()) - pixelCoordinates.y);
    tgt::vec2 old;
    for (int i=0; i < tfProp_->get()->getNumKeys(); ++i) {
        TransFuncMappingKey* key = tfProp_->get()->getKey(i);
        if (i > 0) {
            tgt::vec2 p1 = tgt::vec2(old.x + 1.f, old.y);
            tgt::vec2 p2 = tgt::vec2(pixelCoordinates.x - 1.f, pixelCoordinates.y);
            float s = (p2.y - p1.y) / (p2.x - p1.x);
            int a = static_cast<int>(p1.y + (sHit.x - p1.x) * s);
            if ((sHit.x >= p1.x+10) && (sHit.x <= p2.x-10) && (abs(static_cast<int>(sHit.y) - a) < 5)) {
                hit = i - 1;
            }
        }

        old = pixelCoordinates;
        if (key->isSplit())
            old = wtos(tgt::vec2(key->getIntensity(), key->getColorR().a / 255.f));
    }
    return hit;
}

void TransFunc1DKeysPropertyEditorCanvas::updateToolTipCoordinates(QPoint pos, tgt::vec2 values) {
    std::ostringstream os;
    os.precision(2);
    os.setf(std::ios::fixed, std::ios::floatfield);

    float intensity = values.x;
    if(tfProp_->get()) {
        intensity  = tfProp_->get()->getDomain().x + (tfProp_->get()->getDomain().y - tfProp_->get()->getDomain().x) * intensity;
    }
    os << intensity << " / " << values.y; // intensity / alpha
    QToolTip::showText(mapToGlobal(pos), QString(os.str().c_str()));
}

//-------------------------------------------------------------------------------------------------------------
//      Paint Functions
//-------------------------------------------------------------------------------------------------------------
void TransFunc1DKeysPropertyEditorCanvas::paintEvent(QPaintEvent* event) {
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

void TransFunc1DKeysPropertyEditorCanvas::drawCheckBoard(QPainter* painter) {
    // draw grid
    painter->setPen(QColor(220, 220, 220));
    painter->setBrush(QColor(220, 220, 220));
    painter->setRenderHint(QPainter::Antialiasing, false);

    tgt::vec2 minCoord = wtos(tgt::vec2(0.f,0.f));
    tgt::vec2 maxCoord = wtos(tgt::vec2(1.f,1.f));

    if(showTexture_) {
        int cX = 0, cY = 0;
        tgt::vec2 size = wtos(tgt::vec2(gridSpacing_.x, gridSpacing_.y)) - minCoord;
        for (float x = 0.f; x < 1.f - gridSpacing_.x*0.5; x += gridSpacing_.x) {
            for (float y = 0.f; y < 1.f - gridSpacing_.y*0.5; y +=gridSpacing_.y) {
                if((cX+cY) % 2 == 0) {
                    tgt::vec2 start = wtos(tgt::vec2(x, y));
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
            tgt::vec2 tmp = wtos(tgt::vec2(x,0.f));
            painter->drawLine(tmp.x,minCoord.y,tmp.x,maxCoord.y);
        }
        for (float y = gridSpacing_.y; y < 1.f + gridSpacing_.y*0.5; y += gridSpacing_.y) {
            tgt::vec2 tmp = wtos(tgt::vec2(0.f,y));
            painter->drawLine(minCoord.x,tmp.y,maxCoord.x,tmp.y);
        }
    }
}

void TransFunc1DKeysPropertyEditorCanvas::drawTransFunc(QPainter* painter) {
    QPixmap icon(256,2);
    icon.fill(Qt::transparent);
    //icon.set
    QPainter pmPainter(&icon);
    for(int i = 0; i < 256; i++) {
        tgt::Color color = tfProp_->get()->getTexture()->texelAsFloat(static_cast<size_t>((static_cast<float>(i)/255.f)*(tfProp_->get()->getDimensions().x-1)));
        pmPainter.setPen(QColor(color.r*255,color.g*255,color.b*255,color.a*255));
        pmPainter.drawLine(i,0,i,1);
    }
    tgt::vec2 start = wtos(tgt::vec2::zero);
    tgt::vec2 end = wtos(tgt::vec2::one);
    painter->drawPixmap(start.x,start.y,end.x-start.x,end.y-start.y, icon);
}

void TransFunc1DKeysPropertyEditorCanvas::drawAxes(QPainter* painter) {
    // draw x and y axes
    painter->setRenderHint(QPainter::Antialiasing, true);
    painter->setPen(Qt::gray);
    painter->setBrush(Qt::gray);

    // draw axes independently from visible range
    tgt::vec2 origin = wtos(tgt::vec2(0.f, 0.f));
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
    painter->drawText(static_cast<int>(width() - painter->fontMetrics().width(X_AXIS_TEXT) - 2.78f * AXIS_OFFSET), static_cast<int>(-1 * (origin.y+1.0f - 0.8f * AXIS_OFFSET)), X_AXIS_TEXT);
    painter->drawText(static_cast<int>(1.6f * AXIS_OFFSET), static_cast<int>(-1 * (height() - 1.85f * AXIS_OFFSET)), Y_AXIS_TEXT);
    painter->rotate(180.f);
    painter->scale(-1.f, 1.f);
}

void TransFunc1DKeysPropertyEditorCanvas::drawHistogram(QPainter* painter) {
    if (!histogram_) {
        painter->setMatrixEnabled(false);
        painter->setPen(Qt::red);
        painter->drawText(QRectF(0, 7, width() - 1, height() - 8), tr("No volume or calculating histogram"), QTextOption(Qt::AlignHCenter));
        painter->setMatrixEnabled(true);
    } else {
        if (histogramCache_ == 0 || histogramCache_->rect() != rect()) {
            delete histogramCache_;
            histogramCache_ = new QPixmap(rect().size());
            histogramCache_->fill(Qt::transparent);

            QPainter paint(histogramCache_);

            if (histogram_) {
                // draw histogram
                paint.setPen(Qt::NoPen);
                paint.setBrush(QColor(255, 135, 135, 255)); //200 0 0 120
                paint.setRenderHint(QPainter::Antialiasing, true);

                int histogramWidth = static_cast<int>(histogram_->getNumBuckets());
                tgt::vec2 p;

                // Qt can't handle polygons that have more than 65536 points
                // so we have to split the polygon
                int maxSize = 65536; //max size of polygon
                std::vector<QPointF*> points;
                int vi = 0; //iterator in the points vector
                points.push_back(new QPointF[maxSize]);
                int count = 0;

                for (int x=0; x < histogramWidth; ++x) {
                    float xpos = static_cast<float>(x) / histogramWidth;
                    xpos = histogram_->getMinValue() + (histogram_->getMaxValue() - histogram_->getMinValue()) * xpos;
                    // Do some simple clipping here, as the automatic clipping of drawPolygon()
                    // gets very slow if lots of polygons have to be clipped away, e.g. when
                    // zooming to small part of the histogram.
                    if (xpos >= visibleHistogramRange_.x && xpos <= visibleHistogramRange_.y) {
                        //open new list, if old one is full
                        if( count == maxSize-2 ){
                            count = 0;
                            points.push_back(new QPointF[maxSize]);
                            vi++;
                            //copy last point to connect two polygons
                            points[vi][count].rx() = p.x;
                            points[vi][count].ry() = p.y;
                            count++;
                        }
                        float value = (true ? histogram_->getBucketLogNormalized(x) : histogram_->getBucketNormalized(x));
                        p = histoWtos(tgt::vec2(xpos, value));

                        // optimization: if the y-coord has not changed from the two last points
                        // then just update the last point's x-coord to the current one
                        if( (count >= 2 ) && (points[vi][count - 2].ry() == p.y) && (points[vi][count - 1].ry() == p.y) && (count >= 2) ){
                            points[vi][count - 1].rx() = p.x;
                        } else {
                            points[vi][count].rx() = p.x;
                            points[vi][count].ry() = p.y;
                            count++;
                        }
                    }
                }

                for(size_t i = 0; i < points.size(); ++i){
                    if (count > 0) {
                        if (static_cast<int>(i) == vi){
                            // needed for a closed polygon
                            p = histoWtos(tgt::vec2(0.f, 0.f));
                            points[i][count].rx() = points[i][count -1].rx();
                            points[i][count].ry() = p.y;
                            count++;
                            points[i][count].rx() = points[i][0].rx();
                            points[i][count].ry() = p.y;
                            count++;

                            paint.drawPolygon(points[i], count);
                        } else {
                            // needed for a closed polygon
                            p = histoWtos(tgt::vec2(0.f, 0.f));
                            points[i][maxSize - 2].rx() = points[i][maxSize - 3].rx();
                            points[i][maxSize - 2].ry() = p.y;
                            points[i][maxSize - 1].rx() = points[i][0].rx();
                            points[i][maxSize - 1].ry() = p.y;

                            paint.drawPolygon(points[i], maxSize);
                        }
                    }
                }

                for (size_t i = 0; i < points.size(); ++i)
                    delete[] points[i];
            }
        }
        painter->drawPixmap(0,0,*histogramCache_);
    }
}

void TransFunc1DKeysPropertyEditorCanvas::drawMappingKeys(QPainter* painter) {
    QPen pen = QPen(Qt::darkRed);
    pen.setWidthF(1.5f);
    painter->setPen(pen);

    tgt::vec2 old(0.0f);
    for (int i=0; i<tfProp_->get()->getNumKeys(); ++i) {
        TransFuncMappingKey *key = tfProp_->get()->getKey(i);
        tgt::vec2 p = wtos(tgt::vec2(key->getIntensity(), key->getColorL().a / 255.f));
        if (i == 0)  {
            if (tfProp_->get()->getKey(0)->getIntensity() > 0.f)
                painter->drawLine(QPointF(wtos(tgt::vec2(0.f, 0.f)).x, p.y),
                               QPointF(p.x - 1.f, p.y));
        }
        else {
            painter->drawLine(QPointF(old.x + 1.f, old.y),
                           QPointF(p.x - 1.f, p.y));
        }
        old = p;
        if (key->isSplit())
            old = wtos(tgt::vec2(key->getIntensity(), key->getColorR().a / 255.f));
    }
    if (tfProp_->get()->getNumKeys() > 0 && (tfProp_->get()->getKey(tfProp_->get()->getNumKeys()-1)->getIntensity() < 1.f)) {
        painter->drawLine(QPointF(old.x + 1.f, old.y),
                       QPointF(wtos(tgt::vec2(1.f, 0.f)).x, old.y));
    }

    if (1.f != 0.f)
        paintKeys(*painter);
}

void TransFunc1DKeysPropertyEditorCanvas::drawThreshold(QPainter* painter) {
    // grey out threshold area
    tgt::vec2 origin = wtos(tgt::vec2(0.f));
    painter->setBrush(QBrush(QColor(192, 192, 192, 230), Qt::SolidPattern));
    painter->setPen(Qt::NoPen);
    tgt::vec2 upperRight = wtos(tgt::vec2(1.f));
    tgt::vec2 lowerLeft = wtos(tgt::vec2(0.f));
    int w = static_cast<int>(upperRight.x - lowerLeft.x);
    int h = static_cast<int>(upperRight.y - lowerLeft.y);

    if (tfProp_->get()->getThreshold().x > 0.f) {
        painter->drawRect(static_cast<int>(origin.x), static_cast<int>(origin.y),
                       static_cast<int>(tfProp_->get()->getThreshold().x * w + 1), h);
    }
    if (tfProp_->get()->getThreshold().y < 1.f) {
        painter->drawRect(static_cast<int>(origin.x + floor(tfProp_->get()->getThreshold().y * w)),
                       static_cast<int>(origin.y), static_cast<int>((1 - tfProp_->get()->getThreshold().y) * w + 1), h);
    }
}

void TransFunc1DKeysPropertyEditorCanvas::toggleHistogram(bool state) {
    showHistogram_ = state;
    update();
}

void TransFunc1DKeysPropertyEditorCanvas::toggleTexture(bool state) {
    showTexture_ = state;
    update();
}

//-------------------------------------------------------------------------------------------------------------
//      Mapping Key Functions
//-------------------------------------------------------------------------------------------------------------
void TransFunc1DKeysPropertyEditorCanvas::insertNewKey(tgt::vec2& hit) {

    if (!tfProp_->get())
        return;

    hit = tgt::clamp(hit, 0.f, 1.f);

    TransFuncMappingKey* key = new TransFuncMappingKey(hit.x, QColor2Col(Qt::lightGray));

    tfProp_->get()->addKey(key);
    TransFuncMappingKey* leftKey = getOtherKey(key, true);
    TransFuncMappingKey* rightKey = getOtherKey(key, false);

    // interpolate color of inserted key from neighbouring keys
    // (weighted by distance)
    // the alpha value is determined by hit.y
    tgt::col4 keyColor;
    if (!leftKey && !rightKey)
        keyColor = tgt::vec4(0.f);
    else if (!leftKey)
        keyColor = rightKey->getColorL();
    else if (!rightKey)
        keyColor = leftKey->getColorR();
    else {
        float leftSource = leftKey->getIntensity();
        float rightSource = rightKey->getIntensity();
        float distSource = rightSource - leftSource;
        tgt::vec4 leftColor = static_cast<tgt::vec4>(leftKey->getColorR());
        tgt::vec4 rightColor = static_cast<tgt::vec4>(rightKey->getColorL());

        keyColor = static_cast<tgt::col4>(
            leftColor* ( (distSource-(hit.x-leftSource))/distSource  ) +
            rightColor*( (distSource-(rightSource-hit.x))/distSource ) );
    }
    key->setColorL(keyColor);
    //overwrite alpha value with clicked position
    key->setAlphaL(hit.y);

    selectedKey_ = key;

    tfProp_->invalidate();
}

TransFuncMappingKey* TransFunc1DKeysPropertyEditorCanvas::getOtherKey(TransFuncMappingKey* selectedKey, bool selectedLeftPart) {

    if (!tfProp_->get())
        return 0;

    TransFuncMappingKey* otherKey = 0;
    for (int i=0; i < tfProp_->get()->getNumKeys(); ++i) {
        if ((selectedLeftPart && i < tfProp_->get()->getNumKeys() - 1 && tfProp_->get()->getKey(i + 1) == selectedKey) ||
            (!selectedLeftPart && i > 0 && tfProp_->get()->getKey(i - 1) == selectedKey))
        {
            otherKey = tfProp_->get()->getKey(i);
        }
    }
    return otherKey;
}

void TransFunc1DKeysPropertyEditorCanvas::paintKeys(QPainter& paint) {

    if (!tfProp_->get())
        return;

    for (int i=0; i<tfProp_->get()->getNumKeys(); ++i) {
        TransFuncMappingKey *key = tfProp_->get()->getKey(i);
        tgt::vec2 p = wtos(tgt::vec2(key->getIntensity(), key->getColorL().a / 255.0));
        int props;
        if (key->isSplit()) {
            props = MARKER_LEFT;
            if (key == selectedKey_ && isLeftPartSelected_)
                props |= MARKER_SELECTED;

            drawMarker(paint, key->getColorL(), p, props);

            p = wtos(tgt::vec2(key->getIntensity(), key->getColorR().a / 255.0));
            props = MARKER_RIGHT;
            if (key == selectedKey_ && !isLeftPartSelected_)
                props |= MARKER_SELECTED;

            drawMarker(paint, key->getColorR(), p, props);
        }
        else {
            props = MARKER_NORMAL;
            if (key == selectedKey_)
                props |= MARKER_SELECTED;
            drawMarker(paint, key->getColorL(), p, props);
        }
    }
}

void TransFunc1DKeysPropertyEditorCanvas::drawMarker(QPainter& paint, const tgt::col4& tgtcolor, const tgt::vec2& p, int props) {
    paint.setBrush(Col2QColor(tgtcolor));

    QPen pen(QBrush(Qt::darkGray), Qt::SolidLine);
    if (props & MARKER_SELECTED)
        pen.setWidth(3);
    paint.setPen(pen);

    if (props & MARKER_LEFT) {
        paint.drawPie(QRectF(p.x - KEY_SPLIT_FACTOR * KEY_POINT_SIZE/2, p.y - KEY_POINT_SIZE/2,
                             KEY_SPLIT_FACTOR * KEY_POINT_SIZE, KEY_POINT_SIZE),
                      90 * 16, 180 * 16);
    }
    else if (props & MARKER_RIGHT) {
        paint.drawPie(QRectF(p.x - KEY_SPLIT_FACTOR * KEY_POINT_SIZE/2, p.y - KEY_POINT_SIZE/2,
                             KEY_SPLIT_FACTOR * KEY_POINT_SIZE, KEY_POINT_SIZE),
                      270 * 16, 180 * 16);
    }
    else {
        paint.drawEllipse(QRectF(p.x - KEY_POINT_SIZE/2, p.y - KEY_POINT_SIZE/2,
                                 KEY_POINT_SIZE, KEY_POINT_SIZE));
    }
}

//-------------------------------------------------------------------------------------------------------------
//      Helper Functions
//-------------------------------------------------------------------------------------------------------------
tgt::vec2 TransFunc1DKeysPropertyEditorCanvas::wtos(tgt::vec2 p) {
    float sx = p.x * (static_cast<float>(width())  - 2.f * AXIS_OFFSET - 1.5f * ARROW_LENGTH) + AXIS_OFFSET;
    float sy = p.y * (static_cast<float>(height()) - 2.f * AXIS_OFFSET - 1.5f * ARROW_LENGTH) + AXIS_OFFSET;
    return tgt::vec2(sx, sy);
}

tgt::vec2 TransFunc1DKeysPropertyEditorCanvas::histoWtos(tgt::vec2 p) {
    float sx = (p.x - visibleHistogramRange_.x) / (visibleHistogramRange_.y - visibleHistogramRange_.x) * (static_cast<float>(width())  - 2.f * AXIS_OFFSET - 1.5f * ARROW_LENGTH) + AXIS_OFFSET;
    float sy = p.y * (static_cast<float>(height()) - 2.f * AXIS_OFFSET - 1.5f * ARROW_LENGTH) + AXIS_OFFSET;
    return tgt::vec2(sx, sy);
}

tgt::vec2 TransFunc1DKeysPropertyEditorCanvas::stow(tgt::vec2 p) {
    float wx = (p.x - AXIS_OFFSET) / (static_cast<float>(width())  - 2.f * AXIS_OFFSET - 1.5f * ARROW_LENGTH);
    float wy = (p.y - AXIS_OFFSET) / (static_cast<float>(height()) - 2.f * AXIS_OFFSET - 1.5f * ARROW_LENGTH);
    return tgt::vec2(wx, wy);
}

QColor TransFunc1DKeysPropertyEditorCanvas::Col2QColor(const tgt::col4& color) {
    return QColor(color.r, color.g, color.b); // ignore alpha
}

tgt::col4 TransFunc1DKeysPropertyEditorCanvas::QColor2Col(const QColor& color) {
    return tgt::col4(color.red(), color.green(), color.blue(), 255); // ignore alpha
}

} // namespace voreen
