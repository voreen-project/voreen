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

#include "voreen/qt/widgets/property/transfunc/2d/transfunc2dpropertywidget.h"

#include "voreen/core/datastructures/transfunc/2d/transfunc2d.h"
#include "voreen/core/datastructures/volume/histogram.h"
#include "voreen/core/utils/stringutils.h"

#include "voreen/core/voreenapplication.h"

#include "tgt/shadermanager.h"
#include "tgt/gpucapabilities.h"
#include "tgt/glcanvas.h"
#include "tgt/matrixstack.h"

#include <sstream>

#include <QMouseEvent>
#include <QApplication>
#include <QFont>
#include <QFontMetrics>
#include <QToolTip>
namespace voreen {

//Magic numbers
     const int CANVAS_BORDER_LEFT_SIZE   = 12;
     const int CANVAS_BORDER_RIGHT_SIZE  =  8;
     const int CANVAS_BORDER_TOP_SIZE    =  8;
     const int CANVAS_BORDER_BOTTOM_SIZE = 12;

     const int SLIDER_GAMMA_X_SIZE  = CANVAS_BORDER_TOP_SIZE - 1;
     const int SLIDER_GAMMA_Y_SIZE  = CANVAS_BORDER_RIGHT_SIZE - 1;
     const int SLIDER_DOMAIN_X_SIZE = CANVAS_BORDER_BOTTOM_SIZE -1;
     const int SLIDER_DOMAIN_Y_SIZE = CANVAS_BORDER_LEFT_SIZE - 1;
//--------------------------------
//--------------------------------
// base and tgt::Paint functions
//--------------------------------
//--------------------------------
TransFunc2DPropertyWidgetPainter::TransFunc2DPropertyWidgetPainter(tgt::GLCanvas* canvas, QColor backgroundColor)
    : TransFuncPropertyWidgetPainterBase(canvas)
    //base member
    , tf_(0), histogram_(0)
    //slider
    , gammaSliderX_(0), gammaSliderY_(0)
    , domainSliderX_(0), domainSliderY_(0)
    , pressedSlider_(NO_SLIDER)
    , mousePressedPosition_(-1) , domainLowerPressedPosition_(-1), domainUpperPressedPosition_(-1)
    //helper
    , isInitialized_(false), renderHistogram_(true), renderGammaSlider_(true), renderDomainSlider_(true)
    , logarithmic_(true)
    , backgroundColor_(backgroundColor)
{
    visibleDomainValues_[0] = tgt::vec2(0.f,1.f); visibleDomainValues_[1] = tgt::vec2(0.f,1.f);
    //create slider and connect slots
    gammaSliderX_ = new TransFuncPropertyWidgetPainterGammaSlider(0,TransFuncPropertyWidgetPainterGammaSlider::SO_TOP, SLIDER_GAMMA_X_SIZE);
    gammaSliderY_ = new TransFuncPropertyWidgetPainterGammaSlider(1,TransFuncPropertyWidgetPainterGammaSlider::SO_RIGHT, SLIDER_GAMMA_Y_SIZE);

    domainSliderX_ = new TransFuncPropertyWidgetPainterDomainSlider(0,TransFuncPropertyWidgetPainterGammaSlider::SO_BOTTOM,SLIDER_DOMAIN_X_SIZE);
    domainSliderY_ = new TransFuncPropertyWidgetPainterDomainSlider(1,TransFuncPropertyWidgetPainterGammaSlider::SO_LEFT,SLIDER_DOMAIN_Y_SIZE);

    connect(gammaSliderX_,SIGNAL(gammaChanged(float, size_t)),this,SLOT(gammaSlot(float,size_t)));
    connect(gammaSliderY_,SIGNAL(gammaChanged(float, size_t)),this,SLOT(gammaSlot(float,size_t)));
    connect(domainSliderX_,SIGNAL(domainChanged(tgt::vec2, size_t)),this,SLOT(domainSlot(tgt::vec2,size_t)));
    connect(domainSliderY_,SIGNAL(domainChanged(tgt::vec2, size_t)),this,SLOT(domainSlot(tgt::vec2,size_t)));
}

TransFunc2DPropertyWidgetPainter::~TransFunc2DPropertyWidgetPainter() {
    delete gammaSliderX_;  delete gammaSliderY_;
    delete domainSliderX_; delete domainSliderY_;
}

void TransFunc2DPropertyWidgetPainter::initialize() {
    tgt::Painter::initialize();
}

void TransFunc2DPropertyWidgetPainter::sizeChanged(const tgt::ivec2& _size) {
    tgt::ivec2 size = canvas_->getPhysicalSize();

    domainSliderX_->setCanvasSize(size); domainSliderY_->setCanvasSize(size);
    gammaSliderX_->setCanvasSize(size);  gammaSliderY_->setCanvasSize(size);

    domainSliderX_->setSubCanvasSize(size - tgt::ivec2(CANVAS_BORDER_LEFT_SIZE+CANVAS_BORDER_RIGHT_SIZE,
                                                       CANVAS_BORDER_TOP_SIZE+CANVAS_BORDER_BOTTOM_SIZE));
    domainSliderY_->setSubCanvasSize(size - tgt::ivec2(CANVAS_BORDER_LEFT_SIZE+CANVAS_BORDER_RIGHT_SIZE,
                                                       CANVAS_BORDER_TOP_SIZE+CANVAS_BORDER_BOTTOM_SIZE));
    domainSliderX_->setSubCanvasOffset(tgt::ivec2(CANVAS_BORDER_LEFT_SIZE, CANVAS_BORDER_BOTTOM_SIZE));
    domainSliderY_->setSubCanvasOffset(tgt::ivec2(CANVAS_BORDER_LEFT_SIZE, CANVAS_BORDER_BOTTOM_SIZE));

    gammaSliderX_->setSubCanvasSize(tgt::ivec2(domainSliderX_->getPosition().y-domainSliderX_->getPosition().x+1));
    gammaSliderX_->setSubCanvasOffset(tgt::ivec2(domainSliderX_->getPosition().x));
    gammaSliderY_->setSubCanvasSize(tgt::ivec2(domainSliderY_->getPosition().y-domainSliderY_->getPosition().x+1));
    gammaSliderY_->setSubCanvasOffset(tgt::ivec2(domainSliderY_->getPosition().x));
}

//--------------------------------
//--------------------------------
//          paint
//--------------------------------
//--------------------------------
void TransFunc2DPropertyWidgetPainter::paint() {
    tgtAssert(initialized_, "painter not initialized!");

    tgt::ivec2 size = canvas_->getPhysicalSize();

    //fill canvas with background color
    glViewport(0,0,size.x,size.y);
    glClearColor(backgroundColor_.redF(),backgroundColor_.greenF(),backgroundColor_.blueF(),1.f);
    glClear(GL_COLOR_BUFFER_BIT);

    glViewport(CANVAS_BORDER_LEFT_SIZE, CANVAS_BORDER_BOTTOM_SIZE, size.x - (CANVAS_BORDER_LEFT_SIZE+CANVAS_BORDER_RIGHT_SIZE),
                                                                   size.y - (CANVAS_BORDER_TOP_SIZE+CANVAS_BORDER_BOTTOM_SIZE));
    MatStack.matrixMode(tgt::MatrixStack::PROJECTION);
    MatStack.loadIdentity();
    MatStack.multMatrix(tgt::mat4::createOrtho(0.f, 1.f, 0.f, 1.f, -2.f, 1.f));
    MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    //draw check board
    drawCheckBoard();

    if(tf_) {
        MatStack.scale(1.0f / (visibleDomainValues_[0].y - visibleDomainValues_[0].x), 1.0f / (visibleDomainValues_[1].y - visibleDomainValues_[1].x), 1.0f);
        MatStack.translate(-visibleDomainValues_[0].x, -visibleDomainValues_[1].x, 0.0f);

        //draw transfer function
        drawTransferFunction();

        //draw histogram
        if(histogram_ && renderHistogram_) {
            drawHistogram();
        }

        //draw threshold
        drawThreshold();

        //the sliders need a ortho projection in pixels
        glViewport(0,0,size.x,size.y);
        MatStack.loadIdentity();
        MatStack.matrixMode(tgt::MatrixStack::PROJECTION);
        MatStack.loadIdentity();
        MatStack.multMatrix(tgt::mat4::createOrtho(0.f, static_cast<float>(size.x-1), 0.f, static_cast<float>(size.y-1), -2.f, 1.f));
        MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
        //render slider
        glDisable(GL_DEPTH_TEST);
        if(renderGammaSlider_) {
            gammaSliderX_->paint();
            gammaSliderY_->paint();
        }
        if(renderDomainSlider_) {
            domainSliderX_->paint();
            domainSliderY_->paint();
        }
        glColor4f(1.f,1.f,1.f,1.f);
        glEnable(GL_DEPTH_TEST);
    }
    MatStack.matrixMode(tgt::MatrixStack::PROJECTION);
    MatStack.loadIdentity();
    MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
}

void TransFunc2DPropertyWidgetPainter::drawCheckBoard() {
    int xDim = 10, yDim = 10;
    float xInc = 1.f / xDim, yInc = 1.f / yDim;
    // paint checkerboard
    glBegin(GL_QUADS);
        for (int x = 0 ; x < xDim ; ++x) {
            for (int y = 0 ; y < yDim ; ++y) {
                if ((x+y) % 2)
                    glColor3f(0.6f, 0.6f, 0.6f);
                else
                    glColor3f(1.f, 1.f, 1.f);
                glVertex3f( x      * xInc,  y      * yInc,  -0.5f);  // Bottom Left
                glVertex3f((x + 1) * xInc,  y      * yInc,  -0.5f);  // Bottom Right
                glVertex3f((x + 1) * xInc, (y + 1) * yInc,  -0.5f);  // Top Right
                glVertex3f( x      * xInc, (y + 1) * yInc,  -0.5f);  // Top Left
            }
        }
    glEnd();

}

void TransFunc2DPropertyWidgetPainter::drawTransferFunction() {
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glDisable(GL_DEPTH_TEST);

    glActiveTexture(GL_TEXTURE0);

    tf_->getTexture()->bind();
    tf_->getTexture()->enable();

    MatStack.matrixMode(tgt::MatrixStack::TEXTURE);
    MatStack.loadIdentity();
    MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);

    tgt::vec2 domX = tf_->getDomain(0);
    tgt::vec2 domY = tf_->getDomain(1);

    tgt::vec2 texCoordsX((visibleDomainValues_[0].x - domX.x) /(domX.y-domX.x),(visibleDomainValues_[0].y - domX.x) /(domX.y-domX.x));
    tgt::vec2 texCoordsY((visibleDomainValues_[1].x - domY.x) /(domY.y-domY.x),(visibleDomainValues_[1].y - domY.x) /(domY.y-domY.x));
    //draw tf with alpha
    glBegin(GL_QUADS);
        glTexCoord2f(texCoordsX.x,texCoordsY.x);
        glVertex3f(visibleDomainValues_[0].x, visibleDomainValues_[1].x, -0.5f);
        glTexCoord2f(texCoordsX.y,texCoordsY.x);
        glVertex3f(visibleDomainValues_[0].y, visibleDomainValues_[1].x, -0.5f);
        glTexCoord2f(texCoordsX.y,texCoordsY.y);
        glVertex3f(visibleDomainValues_[0].y, visibleDomainValues_[1].y, -0.5f);
        glTexCoord2f(texCoordsX.x,texCoordsY.y);
        glVertex3f(visibleDomainValues_[0].x, visibleDomainValues_[1].y, -0.5f);
    glEnd();

    glBlendFunc(GL_ONE, GL_ZERO);
    glDisable(GL_BLEND);
    tf_->getTexture()->disable();
    glEnable(GL_DEPTH_TEST);
}

void TransFunc2DPropertyWidgetPainter::drawHistogram() {
    float max = histogram_->getMaxBucket();

    //draw histogram, if it is not empty
    if(max > 0) {

        glDisable(GL_DEPTH_TEST);
        glColor4f(1.f, 0.53f, 0.53f,1.f);

        int nBucketsX = (int)histogram_->getNumBuckets(0);
        int nBucketsY = (int)histogram_->getNumBuckets(1);
        float bucketSizeX = (histogram_->getMaxValue(0) - histogram_->getMinValue(0)) / nBucketsX;
        float bucketSizeY = (histogram_->getMaxValue(1) - histogram_->getMinValue(1)) / nBucketsY;
        float offsetX = histogram_->getMinValue(0);
        float offsetY = histogram_->getMinValue(1);

        glBegin(GL_QUADS);
            for (int x = 0 ; x < nBucketsX ; ++x) {
                for (int y = 0 ; y < nBucketsY ; ++y) {

                    uint64_t bucket = histogram_->getBucket(x,y);
                    if(bucket > 0) {
                        glVertex3f(offsetX + ( x    * bucketSizeX), offsetY + ( y    * bucketSizeY),  -0.5f);  // Bottom Left
                        glVertex3f(offsetX + ((x+1) * bucketSizeX), offsetY + ( y    * bucketSizeY),  -0.5f);  // Bottom Right
                        glVertex3f(offsetX + ((x+1) * bucketSizeX), offsetY + ((y+1) * bucketSizeY),  -0.5f);  // Top Right
                        glVertex3f(offsetX + ( x    * bucketSizeX), offsetY + ((y+1) * bucketSizeY),  -0.5f);  // Top Left
                    }
                }
            }
        glEnd();
        glEnable(GL_DEPTH_TEST);
    }
    //logarythmic
    /*max = logf(max)
    int bucket = histogram_->getBucket(x,y);
    alpha = logf(static_cast<float>(bucket)) / max;
    */
}

void TransFunc2DPropertyWidgetPainter::drawThreshold() {
    //paint threshold quads
    tgt::vec2 threshX = tf_->getThreshold(0), threshY = tf_->getThreshold(1);
    tgt::vec2 domX = tf_->getDomain(0),       domY = tf_->getDomain(1);
    tgt::vec2 visibleXThresholdInterval(threshX.x*(domX.y-domX.x)+domX.x,threshX.y*(domX.y-domX.x)+domX.x);
    tgt::vec2 visibleYThresholdInterval(threshY.x*(domY.y-domY.x)+domY.x,threshY.y*(domY.y-domY.x)+domY.x);
    glDisable(GL_DEPTH_TEST);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glBegin(GL_QUADS);
        glColor4f(0.75f, 0.75f, 0.75f, 0.9f);
        // Quad left of threshold:
        if(threshX.x > 0.f) {
            glVertex3f(visibleDomainValues_[0].x,   visibleDomainValues_[1].x, -0.5f);
            glVertex3f(visibleXThresholdInterval.x, visibleDomainValues_[1].x, -0.5f);
            glVertex3f(visibleXThresholdInterval.x, visibleDomainValues_[1].y, -0.5f);
            glVertex3f(visibleDomainValues_[0].x,   visibleDomainValues_[1].y, -0.5f);
        }
        // Quad right of threshold:
        if(threshX.y < 1.f) {
            glVertex3f(visibleXThresholdInterval.y, visibleDomainValues_[1].x, -0.5f);
            glVertex3f(visibleDomainValues_[0].y,   visibleDomainValues_[1].x, -0.5f);
            glVertex3f(visibleDomainValues_[0].y,   visibleDomainValues_[1].y, -0.5f);
            glVertex3f(visibleXThresholdInterval.y, visibleDomainValues_[1].y, -0.5f);
        }
        // Quad bottom of threshold:
        if(threshY.x > 0.f) {
            glVertex3f(visibleXThresholdInterval.x, visibleDomainValues_[1].x, -0.5f);
            glVertex3f(visibleXThresholdInterval.y, visibleDomainValues_[1].x, -0.5f);
            glVertex3f(visibleXThresholdInterval.y, visibleYThresholdInterval.x, -0.5f);
            glVertex3f(visibleXThresholdInterval.x, visibleYThresholdInterval.x, -0.5f);
        }
        // Quad top of threshold:
        if(threshY.y < 1.f) {
            glVertex3f(visibleXThresholdInterval.x, visibleYThresholdInterval.y, -0.5f);
            glVertex3f(visibleXThresholdInterval.y, visibleYThresholdInterval.y, -0.5f);
            glVertex3f(visibleXThresholdInterval.y, visibleDomainValues_[1].y, -0.5f);
            glVertex3f(visibleXThresholdInterval.x, visibleDomainValues_[1].y, -0.5f);
        }
        glColor4f(1.f,1.f,1.f,1.f);
    glEnd();

    glEnable(GL_DEPTH_TEST);
    glBlendFunc(GL_ONE, GL_ZERO);
    glDisable(GL_BLEND);
}

//--------------------------------
//--------------------------------
//      getter and setter
//--------------------------------
//--------------------------------
void TransFunc2DPropertyWidgetPainter::setTransFunc(TransFuncBase* tf) {
    if(TransFunc2D* tf2d = dynamic_cast<TransFunc2D*>(tf)) {
        if (tf_ != tf) {
            tf_ = tf2d;
            resetZoom();
        }
        if(tf_ != 0) {
            bool visibleDomainChanged = false;
            if(tf2d->getDomain(0).x < visibleDomainValues_[0].x) {
                visibleDomainValues_[0].x = tf2d->getDomain(0).x;
                visibleDomainChanged = true;
            }
            if(tf2d->getDomain(0).y > visibleDomainValues_[0].y) {
                visibleDomainValues_[0].y = tf2d->getDomain(0).y;
                visibleDomainChanged = true;
            }
            if(tf2d->getDomain(1).x < visibleDomainValues_[1].x) {
                visibleDomainValues_[1].x = tf2d->getDomain(1).x;
                visibleDomainChanged = true;
            }
            if(tf2d->getDomain(1).y > visibleDomainValues_[1].y) {
                visibleDomainValues_[1].y = tf2d->getDomain(1).y;
                visibleDomainChanged = true;
            }
            //update meta data
            if(visibleDomainChanged)
                emit storeZoomMetaDataSignal();

            //update Slider
            domainSliderX_->setValueRange(tgt::vec2(visibleDomainValues_[0].x,visibleDomainValues_[0].y));
            domainSliderX_->setDomainValue(tf2d->getDomain(0));
            domainSliderY_->setValueRange(tgt::vec2(visibleDomainValues_[1].x,visibleDomainValues_[1].y));
            domainSliderY_->setDomainValue(tf2d->getDomain(1));

            gammaSliderX_->setSubCanvasSize(tgt::ivec2(domainSliderX_->getPosition().y-domainSliderX_->getPosition().x+1));
            gammaSliderX_->setSubCanvasOffset(tgt::ivec2(domainSliderX_->getPosition().x));
            gammaSliderX_->setGammaValue(tf2d->getGammaValue().x);
            gammaSliderY_->setSubCanvasSize(tgt::ivec2(domainSliderY_->getPosition().y-domainSliderY_->getPosition().x+1));
            gammaSliderY_->setSubCanvasOffset(tgt::ivec2(domainSliderY_->getPosition().x));
            gammaSliderY_->setGammaValue(tf2d->getGammaValue().y);
        }
    }
}

void TransFunc2DPropertyWidgetPainter::setHistogram(const Histogram* histogram) {
    histogram_ = static_cast<const Histogram2D*>(histogram);
}

//--------------------------------------------
//--------------------------------------------
//          Mouse Events
//--------------------------------------------
//--------------------------------------------
void TransFunc2DPropertyWidgetPainter::mousePressEvent(tgt::MouseEvent* e) {
    //TODO: does not work
    //emit showToolTip(e-> coord(),QString("test"));
    if(renderGammaSlider_ && gammaSliderX_->isHit(e->coord())) {
        pressedSlider_ = GAMMA_SLIDER_X;
    } else
    if(renderGammaSlider_ && gammaSliderY_->isHit(e->coord())) {
        pressedSlider_ = GAMMA_SLIDER_Y;
    } else
    if(renderDomainSlider_ && domainSliderX_->isHit(e->coord())) {
        if(domainSliderX_->isLeftHit(e->coord()))
            pressedSlider_ = DOMAIN_LEFT_SLIDER_X;
        else
            pressedSlider_ = DOMAIN_RIGHT_SLIDER_X;
    } else
    if(renderDomainSlider_ && domainSliderY_->isHit(e->coord())) {
        if(domainSliderY_->isLeftHit(e->coord()))
            pressedSlider_ = DOMAIN_LEFT_SLIDER_Y;
        else
            pressedSlider_ = DOMAIN_RIGHT_SLIDER_Y;
    } else
    if(renderDomainSlider_ && domainSliderX_->isIntervalHit(e->coord())) {
        pressedSlider_ = DOMAIN_BOTH_SLIDER_X;
        mousePressedPosition_ = e->coord();
        domainLowerPressedPosition_ = domainSliderX_->getPosition().x;
        domainUpperPressedPosition_ = domainSliderX_->getPosition().y;
    } else
    if(renderDomainSlider_ && domainSliderY_->isIntervalHit(e->coord())) {
        pressedSlider_ = DOMAIN_BOTH_SLIDER_Y;
        mousePressedPosition_ = e->coord();
        domainLowerPressedPosition_ = domainSliderY_->getPosition().x;
        domainUpperPressedPosition_ = domainSliderY_->getPosition().y;
    } else
        pressedSlider_ = NO_SLIDER;

    if(pressedSlider_ != NO_SLIDER)
        emit interaction(true);
}

void TransFunc2DPropertyWidgetPainter::mouseMoveEvent(tgt::MouseEvent* e) {
    switch(pressedSlider_) {
    case GAMMA_SLIDER_X:
        gammaSliderX_->setPosition(e->coord());
        canvas_->update();
    break;
    case GAMMA_SLIDER_Y:
        gammaSliderY_->setPosition(e->coord());
        canvas_->update();
    break;
    case DOMAIN_LEFT_SLIDER_X:
        domainSliderX_->setLeftPosition(e->coord());
        canvas_->update();
    break;
    case DOMAIN_LEFT_SLIDER_Y:
        domainSliderY_->setLeftPosition(e->coord());
        canvas_->update();
    break;
    case DOMAIN_RIGHT_SLIDER_X:
        domainSliderX_->setRightPosition(e->coord());
        canvas_->update();
    break;
    case DOMAIN_RIGHT_SLIDER_Y:
        domainSliderY_->setRightPosition(e->coord());
        canvas_->update();
    break;
    case DOMAIN_BOTH_SLIDER_X:
        tgtAssert(mousePressedPosition_ != tgt::ivec2(-1), "mousePressedPosition has not been set!");
        tgtAssert(domainLowerPressedPosition_ != -1, "domainLeftPressedPosition has not been set!");
        tgtAssert(domainUpperPressedPosition_ != -1, "domainRightPressedPosition has not been set!");
        domainSliderX_->setLeftPosition( tgt::ivec2(domainLowerPressedPosition_, canvas_->getPhysicalSize().y -1 - domainLowerPressedPosition_) - (mousePressedPosition_ - e->coord()));
        domainSliderX_->setRightPosition(tgt::ivec2(domainUpperPressedPosition_, canvas_->getPhysicalSize().y -1 - domainUpperPressedPosition_) - (mousePressedPosition_ -e->coord()));
        canvas_->update();
    break;
    case DOMAIN_BOTH_SLIDER_Y:
        tgtAssert(mousePressedPosition_ != tgt::ivec2(-1), "mousePressedPosition has not been set!");
        tgtAssert(domainLowerPressedPosition_ != -1, "domainLeftPressedPosition has not been set!");
        tgtAssert(domainUpperPressedPosition_ != -1, "domainRightPressedPosition has not been set!");
        domainSliderY_->setLeftPosition( tgt::ivec2(domainLowerPressedPosition_, canvas_->getPhysicalSize().y -1 - domainLowerPressedPosition_) - (mousePressedPosition_ -e->coord()));
        domainSliderY_->setRightPosition(tgt::ivec2(domainUpperPressedPosition_, canvas_->getPhysicalSize().y -1 - domainUpperPressedPosition_) - (mousePressedPosition_ -e->coord()));
        canvas_->update();
    break;
    default:
        // do nothing
    break;
    }
    //creates and shows tooltip
    createInfoToolTip(QPoint(e->coord().x,e->coord().y));
}

void TransFunc2DPropertyWidgetPainter::mouseReleaseEvent(tgt::MouseEvent* e) {
    pressedSlider_ = NO_SLIDER;
    mousePressedPosition_ = tgt::ivec2(-1);
    domainLowerPressedPosition_ = -1;
    domainUpperPressedPosition_ = -1;
    emit hideInfoToolTip();
    emit interaction(false);
}

void TransFunc2DPropertyWidgetPainter::mouseDoubleClickEvent(tgt::MouseEvent* e) {
    if(renderGammaSlider_) {
        if(gammaSliderX_->isHit(e->coord())) {
             gammaSliderX_->setGammaValue(1.f);
             gammaSlot(1.f,0);
             canvas_->update();
        } else
        if (gammaSliderY_->isHit(e->coord())) {
             gammaSliderY_->setGammaValue(1.f);
             gammaSlot(1.f,1);
             canvas_->update();
        }
    }
}

void TransFunc2DPropertyWidgetPainter::createInfoToolTip(QPoint mousePos) {
    //info tool tip output
    QPoint ittPos;    //x: center of tip, y: lower border
    std::stringstream ittText;
    QString qText;
    //helper to determine tool tip position
    int ittLeftBorder = 6; // empty space on left side (estimated)
    int ittWidth = 0;
    QFont ttFont = QToolTip::font();
    QFont bttFont = ttFont;
    bttFont.setBold(true);
    QFontMetrics ttMetr(ttFont);
    QFontMetrics bttMetr(bttFont);
    //percision
    /*float interval = domainSlider_->getDomainValue().y - domainSlider_->getDomainValue().x;
    int percision = 3;
    if(interval > 1000.f)
        percision = 0;
    else if(interval > 100.f)
        percision = 1;
    else if(interval > 10.f)
        percision = 2;*/

    switch(pressedSlider_){
    case GAMMA_SLIDER_X:
        //position
        ittPos.setX(mousePos.x());
        ittPos.setY(-bttMetr.height()+3);
        //gamma
        qText = QString(1,QChar(0x03b3)) + QString("=");
        ittWidth += bttMetr.width(qText);
        //value
        ittText << ftos(gammaSliderX_->getGammaValue(),2) << "</b>";
        ittWidth += bttMetr.width(ftos(gammaSliderX_->getGammaValue(),2).c_str());
        qText = QString("<b>") + qText + QString(ittText.str().c_str());
        break;
    case GAMMA_SLIDER_Y:
        //gamma
        qText = QString(1,QChar(0x03b3)) + QString("=");
        ittWidth += bttMetr.width(qText);
        //value
        ittText << ftos(gammaSliderY_->getGammaValue(),2) << "</b>";
        ittWidth += bttMetr.width(ftos(gammaSliderY_->getGammaValue(),2).c_str());
        qText = QString("<b>") + qText + QString(ittText.str().c_str());
        //position
        ittPos.setX(canvas_->getPhysicalWidth() - ittWidth/2 - 3);
        ittPos.setY(mousePos.y());
        break;

        /*case DOMAIN_LEFT_SLIDER:
        ittPos.setX(canvas_->getWidth() / 2);
        ittPos.setY(canvas_->getHeight()+10);
        ittText << "<nobr>Domain: [ <b>" << ftos(domainSlider_->getDomainValue().x,percision) << "</b> , " <<
                                    ftos(domainSlider_->getDomainValue().y,percision) << " ]</nobr>";
        ittWidth += ttMetr.width("Domain: [  ,  ]");
        ittWidth += bttMetr.width(ftos(domainSlider_->getDomainValue().x,percision).c_str());
        ittWidth += ttMetr.width(ftos(domainSlider_->getDomainValue().y,percision).c_str());
        qText = QString(ittText.str().c_str());
        break;
    case DOMAIN_RIGHT_SLIDER:
        ittPos.setX(canvas_->getWidth() / 2);
        ittPos.setY(canvas_->getHeight()+10);
        ittText << "<nobr>Domain: [ " + ftos(domainSlider_->getDomainValue().x,percision) << " , <b>" <<
                                ftos(domainSlider_->getDomainValue().y,percision) << "</b> ]</nobr>";
        ittWidth += ttMetr.width("Domain: [  ,  ]");
        ittWidth += ttMetr.width(ftos(domainSlider_->getDomainValue().x,percision).c_str());
        ittWidth += bttMetr.width(ftos(domainSlider_->getDomainValue().y,percision).c_str());
        qText = QString(ittText.str().c_str());
        break;
    case DOMAIN_BOTH_SLIDER:
        ittPos.setX(canvas_->getWidth() / 2);
        ittPos.setY(canvas_->getHeight()+10);
        ittText << "<nobr>Domain: [ <b>" << ftos(domainSlider_->getDomainValue().x,percision) << "</b> , <b>" +
                                    ftos(domainSlider_->getDomainValue().y,percision) << "</b> ]</nobr>";
        ittWidth += ttMetr.width("Domain: [  ,  ]");
        ittWidth += bttMetr.width(ftos(domainSlider_->getDomainValue().x,percision).c_str());
        ittWidth += bttMetr.width(ftos(domainSlider_->getDomainValue().y,percision).c_str());
        qText = QString(ittText.str().c_str());
        break;*/
    default:
        return;
        break;
    }
    //shift info tool tip
    ittPos.setX(ittPos.x()-ittWidth/2-ittLeftBorder);
    //show tip
    emit showInfoToolTip(ittPos,qText);
}

void TransFunc2DPropertyWidgetPainter::gammaSlot(float gamma, size_t dimension) {
    tgt::vec2 tmpGamma = tf_->getGammaValue();
    if(tmpGamma[dimension] != gamma) {
        tmpGamma[dimension] = gamma;
        tf_->setGammaValue(tmpGamma);
        emit changedGamma();
    }
}
void TransFunc2DPropertyWidgetPainter::domainSlot(tgt::vec2 domain, size_t dimension) {
    tgt::vec2 tmpDomain = tf_->getDomain(dimension);
    if(tmpDomain != domain) {
        tf_->setDomain(domain,dimension);
        //move gamma slider
        if(dimension == 0) {
            gammaSliderX_->setSubCanvasSize(tgt::ivec2(domainSliderX_->getPosition().y-domainSliderX_->getPosition().x+1));
            gammaSliderX_->setSubCanvasOffset(tgt::ivec2(domainSliderX_->getPosition().x));
        }else {
            gammaSliderY_->setSubCanvasSize(tgt::ivec2(domainSliderY_->getPosition().y-domainSliderY_->getPosition().x+1));
            gammaSliderY_->setSubCanvasOffset(tgt::ivec2(domainSliderY_->getPosition().x));
        }
        emit changedDomain();
    }
}

//--------------------------------------------
//--------------------------------------------
//  zoom functions
//--------------------------------------------
//--------------------------------------------
void TransFunc2DPropertyWidgetPainter::zoomIn() {
    if(tf_) {
        float domCenterX = (tf_->getDomain(0).x + tf_->getDomain(0).y) * 0.5f;
        float domCenterY = (tf_->getDomain(1).x + tf_->getDomain(1).y) * 0.5f;
        float domSizeX = tf_->getDomain(0).y - tf_->getDomain(0).x;
        float domSizeY = tf_->getDomain(1).y - tf_->getDomain(1).x;
        float viewSizeX = visibleDomainValues_[0].y - visibleDomainValues_[0].x;
        float viewSizeY = visibleDomainValues_[1].y - visibleDomainValues_[1].x;

        viewSizeX *= 0.5f; viewSizeY *= 0.5f;
        if(viewSizeX < domSizeX) viewSizeX = domSizeX;
        if(viewSizeY < domSizeY) viewSizeY = domSizeY;

        visibleDomainValues_[0].x = domCenterX - (viewSizeX * 0.5f); visibleDomainValues_[0].y = domCenterX + (viewSizeX * 0.5f);
        visibleDomainValues_[1].x = domCenterY - (viewSizeY * 0.5f); visibleDomainValues_[1].y = domCenterY + (viewSizeY * 0.5f);

        domainSliderX_->setValueRange(tgt::vec2(visibleDomainValues_[0].x,visibleDomainValues_[0].y));
        gammaSliderX_->setSubCanvasSize(tgt::ivec2(domainSliderX_->getPosition().y-domainSliderX_->getPosition().x+1));
        gammaSliderX_->setSubCanvasOffset(tgt::ivec2(domainSliderX_->getPosition().x));
        domainSliderY_->setValueRange(tgt::vec2(visibleDomainValues_[1].x,visibleDomainValues_[1].y));
        gammaSliderY_->setSubCanvasSize(tgt::ivec2(domainSliderY_->getPosition().y-domainSliderY_->getPosition().x+1));
        gammaSliderY_->setSubCanvasOffset(tgt::ivec2(domainSliderY_->getPosition().x));

        //repaint
        canvas_->update();
    }
}

void TransFunc2DPropertyWidgetPainter::zoomOut() {
    if(tf_) {
        float domCenterX = (tf_->getDomain(0).x + tf_->getDomain(0).y) * 0.5f;
        float domCenterY = (tf_->getDomain(1).x + tf_->getDomain(1).y) * 0.5f;
        float viewSizeX = visibleDomainValues_[0].y - visibleDomainValues_[0].x;
        float viewSizeY = visibleDomainValues_[1].y - visibleDomainValues_[1].x;
        float viewCenterX = domCenterX;
        float viewCenterY = domCenterY;

        viewSizeX *= 2.0f; viewSizeY *= 2.0f;

        visibleDomainValues_[0].x = viewCenterX - (viewSizeX * 0.5f); visibleDomainValues_[0].y = viewCenterX + (viewSizeX * 0.5f);
        visibleDomainValues_[1].x = viewCenterY - (viewSizeY * 0.5f); visibleDomainValues_[1].y = viewCenterY + (viewSizeY * 0.5f);

        domainSliderX_->setValueRange(tgt::vec2(visibleDomainValues_[0].x,visibleDomainValues_[0].y));
        gammaSliderX_->setSubCanvasSize(tgt::ivec2(domainSliderX_->getPosition().y-domainSliderX_->getPosition().x+1));
        gammaSliderX_->setSubCanvasOffset(tgt::ivec2(domainSliderX_->getPosition().x));
        domainSliderY_->setValueRange(tgt::vec2(visibleDomainValues_[1].x,visibleDomainValues_[1].y));
        gammaSliderY_->setSubCanvasSize(tgt::ivec2(domainSliderY_->getPosition().y-domainSliderY_->getPosition().x+1));
        gammaSliderY_->setSubCanvasOffset(tgt::ivec2(domainSliderY_->getPosition().x));

        //repaint
        canvas_->update();
    }
}

void TransFunc2DPropertyWidgetPainter::resetZoom() {
    if(tf_) {
        visibleDomainValues_[0] = tf_->getDomain(0);
        visibleDomainValues_[1] = tf_->getDomain(1);

        domainSliderX_->setValueRange(tgt::vec2(visibleDomainValues_[0].x,visibleDomainValues_[0].y));
        gammaSliderX_->setSubCanvasSize(tgt::ivec2(domainSliderX_->getPosition().y-domainSliderX_->getPosition().x+1));
        gammaSliderX_->setSubCanvasOffset(tgt::ivec2(domainSliderX_->getPosition().x));
        domainSliderY_->setValueRange(tgt::vec2(visibleDomainValues_[1].x,visibleDomainValues_[1].y));
        gammaSliderY_->setSubCanvasSize(tgt::ivec2(domainSliderY_->getPosition().y-domainSliderY_->getPosition().x+1));
        gammaSliderY_->setSubCanvasOffset(tgt::ivec2(domainSliderY_->getPosition().x));

        //repaint
        canvas_->update();
    }
}

float TransFunc2DPropertyWidgetPainter::getMinVisibleDomainValue(size_t dimension) const {
    tgtAssert(dimension <= 1, "wrong dimension");
    return visibleDomainValues_[dimension].x;
}

float TransFunc2DPropertyWidgetPainter::getMaxVisibleDomainValue(size_t dimension) const {
    tgtAssert(dimension <= 1, "wrong dimension");
    return visibleDomainValues_[dimension].y;
}

void TransFunc2DPropertyWidgetPainter::setVisibleDomainValues(tgt::vec2 visible1Domain, tgt::vec2 visible2Domain) {
    visibleDomainValues_[0] = visible1Domain;
    visibleDomainValues_[1] = visible2Domain;

    domainSliderX_->setValueRange(tgt::vec2(visibleDomainValues_[0].x,visibleDomainValues_[0].y));
    gammaSliderX_->setSubCanvasSize(tgt::ivec2(domainSliderX_->getPosition().y-domainSliderX_->getPosition().x+1));
    gammaSliderX_->setSubCanvasOffset(tgt::ivec2(domainSliderX_->getPosition().x));
    domainSliderY_->setValueRange(tgt::vec2(visibleDomainValues_[1].x,visibleDomainValues_[1].y));
    gammaSliderY_->setSubCanvasSize(tgt::ivec2(domainSliderY_->getPosition().y-domainSliderY_->getPosition().x+1));
    gammaSliderY_->setSubCanvasOffset(tgt::ivec2(domainSliderY_->getPosition().x));

    //repaint
    canvas_->update();
}

} // namespace voreen
