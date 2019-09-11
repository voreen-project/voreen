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

#include "voreen/qt/widgets/property/transfunc/1d/transfunc1dpropertywidget.h"

#include "voreen/core/datastructures/transfunc/1d/transfunc1d.h"
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

//--------------------------------
//--------------------------------
// base and tgt::Paint functions
//--------------------------------
//--------------------------------
TransFunc1DPropertyWidgetPainter::TransFunc1DPropertyWidgetPainter(tgt::GLCanvas* canvas)
    : TransFuncPropertyWidgetPainterBase(canvas)
    //base member
    , tf_(0), histogram_(0)
    //slider
    , gammaSlider_(0)
    , domainSlider_(0)
    , pressedSlider_(NO_SLIDER)
    , mousePressedPosition_(-1) , domainLeftPressedPosition_(-1), domainRightPressedPosition_(-1)
    //helper
    , isInitialized_(false), renderHistogram_(true), renderGammaSlider_(true), renderDomainSlider_(true)
    , logarithmic_(true), visibleDomainValues_(0.f,1.f)
{
    //create slider and connect slots
    gammaSlider_ = new TransFuncPropertyWidgetPainterGammaSlider(0,TransFuncPropertyWidgetPainterGammaSlider::SO_TOP, 5);
    domainSlider_ = new TransFuncPropertyWidgetPainterDomainSlider(0,TransFuncPropertyWidgetPainterGammaSlider::SO_BOTTOM, 10);
    connect(gammaSlider_,SIGNAL(gammaChanged(float, size_t)),this,SLOT(gammaSlot(float)));
    connect(domainSlider_,SIGNAL(domainChanged(tgt::vec2, size_t)),this,SLOT(domainSlot(tgt::vec2)));
}

TransFunc1DPropertyWidgetPainter::~TransFunc1DPropertyWidgetPainter() {
    delete gammaSlider_;
    delete domainSlider_;
}

void TransFunc1DPropertyWidgetPainter::initialize() {
    tgt::Painter::initialize();
}

void TransFunc1DPropertyWidgetPainter::sizeChanged(const tgt::ivec2& _size) {
    tgt::ivec2 size = canvas_->getPhysicalSize();
    domainSlider_->setCanvasSize(size);
    gammaSlider_->setCanvasSize(size);

    domainSlider_->setSubCanvasSize(size);
    gammaSlider_->setSubCanvasSize(tgt::ivec2(domainSlider_->getPosition().y-domainSlider_->getPosition().x+1));
    gammaSlider_->setSubCanvasOffset(tgt::ivec2(domainSlider_->getPosition().x));
}

//--------------------------------
//--------------------------------
//          paint
//--------------------------------
//--------------------------------
void TransFunc1DPropertyWidgetPainter::paint() {
    tgtAssert(initialized_, "painter not initialized!");

    tgt::ivec2 size = canvas_->getPhysicalSize();

    glViewport(0, 0, size.x, size.y);
    MatStack.matrixMode(tgt::MatrixStack::PROJECTION);
    MatStack.loadIdentity();
    MatStack.multMatrix(tgt::mat4::createOrtho(0.f, 1.f, 0.f, 1.f, -2.f, 1.f));
    MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    //draw check board
    drawCheckBoard();

    if(tf_) {
        MatStack.scale(1.0f / (visibleDomainValues_.y - visibleDomainValues_.x), 1.0f, 1.0f);
        MatStack.translate(-visibleDomainValues_.x, 0.0f, 0.0f);

        //draw transfer function
        drawTransferFunction();

        //draw histogram
        if(histogram_ && renderHistogram_) {
            drawHistogram();
        }

        //draw threshold
        drawThreshold();

        //the sliders need a ortho projection in pixels
        MatStack.loadIdentity();
        MatStack.matrixMode(tgt::MatrixStack::PROJECTION);
        MatStack.loadIdentity();
        MatStack.multMatrix(tgt::mat4::createOrtho(0.f, static_cast<float>(size.x-1), 0.f, static_cast<float>(size.y-1), -2.f, 1.f));
        MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
        //render slider
        glDisable(GL_DEPTH_TEST);
        if(renderGammaSlider_)
            gammaSlider_->paint();
        if(renderDomainSlider_)
            domainSlider_->paint();
        glColor4f(1.f,1.f,1.f,1.f);
        glEnable(GL_DEPTH_TEST);
    }
    MatStack.matrixMode(tgt::MatrixStack::PROJECTION);
    MatStack.loadIdentity();
    MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
}

void TransFunc1DPropertyWidgetPainter::drawCheckBoard() {
    float inc = 0.1f;
    // paint checkerboard
    for (int i = 0 ; i < 10 ; ++i) {
        glBegin(GL_QUADS);
            // Front Face
            if (i % 2)
                glColor3f(0.6f, 0.6f, 0.6f);
            else
                glColor3f(1.f, 1.f, 1.f);
            glVertex3f( i      * inc, 0.0f,  -0.5f);  // Bottom Left
            glVertex3f((i + 1) * inc, 0.0f,  -0.5f);  // Bottom Right
            glVertex3f((i + 1) * inc, 0.5f,  -0.5f);  // Top Right
            glVertex3f( i      * inc, 0.5f,  -0.5f);  // Top Left
        glEnd();
        glBegin(GL_QUADS);
            // Front Face
            if (i % 2)
                glColor3f(1.f, 1.f, 1.f);
            else
                glColor3f(0.6f, 0.6f, 0.6f);
            glVertex3f( i      * inc, 0.5f,  -0.5f);  // Bottom Left
            glVertex3f((i + 1) * inc, 0.5f,  -0.5f);  // Bottom Right
            glVertex3f((i + 1) * inc, 1.0f,  -0.5f);  // Top Right
            glVertex3f( i      * inc, 1.0f,  -0.5f);  // Top Left
        glEnd();
    }
}

void TransFunc1DPropertyWidgetPainter::drawTransferFunction() {
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glDisable(GL_DEPTH_TEST);

    glActiveTexture(GL_TEXTURE0);

    tf_->getTexture()->bind();
    tf_->getTexture()->enable();

    MatStack.matrixMode(tgt::MatrixStack::TEXTURE);
    MatStack.loadIdentity();
    MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);

    tgt::vec2 dom = tf_->getDomain();

    //draw tf with alpha
    glBegin(GL_QUADS);
        // Quad left of domain:
        glTexCoord1f(0.f);
        glVertex3f(visibleDomainValues_.x, 0.f, -0.5f);
        glVertex3f(dom.x, 0.f, -0.5f);
        glVertex3f(dom.x, 1.f, -0.5f);
        glVertex3f(visibleDomainValues_.x, 1.f, -0.5f);
        // Inside domain:
        glTexCoord1f(0.f);
        glVertex3f(dom.x, 0.f, -0.5f);
        glTexCoord1f(1.f);
        glVertex3f(dom.y, 0.f, -0.5f);
        glVertex3f(dom.y, 1.f, -0.5f);
        glTexCoord1f(0.f);
        glVertex3f(dom.x, 1.f, -0.5f);
        // Quad right of domain:
        glTexCoord1f(1.f);
        glVertex3f(dom.y, 0.f, -0.5f);
        glVertex3f(visibleDomainValues_.y, 0.f, -0.5f);
        glVertex3f(visibleDomainValues_.y, 1.f, -0.5f);
        glVertex3f(dom.y, 1.f, -0.5f);
    glEnd();

    //draw tf without alpha
    glBlendFunc(GL_ONE, GL_ZERO);
    glDisable(GL_BLEND);
    float height = 22.f/30.f;
    glBegin(GL_QUADS);
        glColor4f(1.f, 1.f, 1.f, 1.f);
        // Quad left of domain:
        glTexCoord1f(0.f);
        glVertex3f(visibleDomainValues_.x, height, -0.5f);
        glVertex3f(dom.x, height, -0.5f);
        glVertex3f(dom.x, 1.f, -0.5f);
        glVertex3f(visibleDomainValues_.x, 1.f, -0.5f);
        // Inside domain:
        glTexCoord1f(0.f);
        glVertex3f(dom.x, height, -0.5f);
        glTexCoord1f(1.f);
        glVertex3f(dom.y, height, -0.5f);
        glVertex3f(dom.y, 1.f, -0.5f);
        glTexCoord1f(0.f);
        glVertex3f(dom.x, 1.f, -0.5f);
        // Quad right of domain:
        glTexCoord1f(1.f);
        glVertex3f(dom.y, height, -0.5f);
        glVertex3f(visibleDomainValues_.y, height, -0.5f);
        glVertex3f(visibleDomainValues_.y, 1.f, -0.5f);
        glVertex3f(dom.y, 1.f, -0.5f);
    glEnd();

    tf_->getTexture()->disable();

    glBegin(GL_LINES);
        glColor3f(0.2f,0.2f,0.2f);
        glVertex3f(visibleDomainValues_.x,height-(1.f/30.f)+0.02,-0.5f);
        glVertex3f(visibleDomainValues_.y,height-(1.f/30.f)+0.02,-0.5f);
    glEnd();
    glEnable(GL_DEPTH_TEST);
}

void TransFunc1DPropertyWidgetPainter::drawHistogram() {
    // paint histogram
        glDisable(GL_DEPTH_TEST);
        int nBuckets = (int)histogram_->getNumBuckets();

        glColor3f(1.f, 0.53f, 0.53f);
        float bucketSize = (histogram_->getMaxValue() - histogram_->getMinValue()) / nBuckets;
        float offset = histogram_->getMinValue();
        float max = histogram_->getMaxBucket();
        if(max > 0) {
            if(logarithmic_)
                max = logf(max);
            glBegin(GL_QUADS);
            for (int i = 0 ; i < nBuckets ; ++i) {
                float y = 0.0f;
                int bucket = histogram_->getBucket(i);
                if(bucket > 0) {
                    if(logarithmic_) {
                        y = logf(static_cast<float>(bucket)) / max;
                    }
                    else {
                        y = bucket / max;
                    }
                }

                //only 2/3 should be histogram
                y *= 2.f/3.f;

                glVertex3f(offset + (i * bucketSize), 0.0f,  -0.5f);  // Bottom Left
                glVertex3f(offset + ((i+1) * bucketSize), 0.0f,  -0.5f);  // Bottom Right
                glVertex3f(offset + ((i+1) * bucketSize), y,  -0.5f);  // Top Right
                glVertex3f(offset + (i * bucketSize), y,  -0.5f);  // Top Left
            }
            glEnd();
        }
        glEnable(GL_DEPTH_TEST);
}

void TransFunc1DPropertyWidgetPainter::drawThreshold() {
    //paint threshold quads
    tgt::vec2 thresh = tf_->getThreshold();
    tgt::vec2 dom = tf_->getDomain();
    glDisable(GL_DEPTH_TEST);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glBegin(GL_QUADS);
    glColor4f(0.75f, 0.75f, 0.75f, 0.9f);

    if(thresh.x > 0.f) {
    // Quad left of threshold:
        glVertex3f(visibleDomainValues_.x, 0.f, -0.5f);
        glVertex3f(thresh.x*(dom.y-dom.x)+dom.x, 0.f, -0.5f);
        glVertex3f(thresh.x*(dom.y-dom.x)+dom.x, 1.f, -0.5f);
        glVertex3f(visibleDomainValues_.x, 1.f, -0.5f);
    }

    if(thresh.y < 1.f) {
    // Quad right of threshold:
    glVertex3f(thresh.y*(dom.y-dom.x)+dom.x, 0.f, -0.5f);
    glVertex3f(visibleDomainValues_.y, 0.f, -0.5f);
    glVertex3f(visibleDomainValues_.y, 1.f, -0.5f);
    glVertex3f(thresh.y*(dom.y-dom.x)+dom.x, 1.f, -0.5f);
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
void TransFunc1DPropertyWidgetPainter::setTransFunc(TransFuncBase* tf) {
    if(TransFunc1D* tf1d = dynamic_cast<TransFunc1D*>(tf)) {
        if (tf_ != tf) {
            tf_ = tf1d;
            resetZoom();
        }
        bool visibleDomainChanged = false;
        if(tf1d->getDomain().x < visibleDomainValues_.x) {
            visibleDomainValues_.x = tf1d->getDomain().x;
            visibleDomainChanged = true;
        }
        if(tf1d->getDomain().y > visibleDomainValues_.y) {
            visibleDomainValues_.y = tf1d->getDomain().y;
            visibleDomainChanged = true;
        }
        //update meta data
        if(visibleDomainChanged)
            emit storeZoomMetaDataSignal();

        //update Slider
        domainSlider_->setValueRange(tgt::vec2(visibleDomainValues_.x,visibleDomainValues_.y));
        domainSlider_->setDomainValue(tf1d->getDomain());

        gammaSlider_->setSubCanvasSize(tgt::ivec2(domainSlider_->getPosition().y-domainSlider_->getPosition().x+1,0));
        gammaSlider_->setSubCanvasOffset(tgt::ivec2(domainSlider_->getPosition().x,0));
        gammaSlider_->setGammaValue(tf1d->getGammaValue());
    }
}

void TransFunc1DPropertyWidgetPainter::setHistogram(const Histogram* histogram) {
    if(histogram_ != histogram) {
        histogram_ = static_cast<const Histogram1D*>(histogram);
    }
}

//--------------------------------------------
//--------------------------------------------
//          Mouse Events
//--------------------------------------------
//--------------------------------------------
void TransFunc1DPropertyWidgetPainter::mousePressEvent(tgt::MouseEvent* e) {
    //TODO: does not work
    //emit showToolTip(e-> coord(),QString("test"));
    if(renderGammaSlider_ && gammaSlider_->isHit(e->coord())) {
        pressedSlider_ = GAMMA_SLIDER;
    } else if(renderDomainSlider_ && domainSlider_->isHit(e->coord())) {
        if(domainSlider_->isLeftHit(e->coord()))
            pressedSlider_ = DOMAIN_LEFT_SLIDER;
        else
            pressedSlider_ = DOMAIN_RIGHT_SLIDER;
    } else if(renderDomainSlider_ && (e->coord().x <= domainSlider_->getPosition().y) && (e->coord().x >= domainSlider_->getPosition().x)) {
        pressedSlider_ = DOMAIN_BOTH_SLIDER;
        mousePressedPosition_ = e->coord().x;
        domainLeftPressedPosition_ = domainSlider_->getPosition().x;
        domainRightPressedPosition_ = domainSlider_->getPosition().y;
    } else
        pressedSlider_ = NO_SLIDER;

    if(pressedSlider_ != NO_SLIDER)
        emit interaction(true);
}

void TransFunc1DPropertyWidgetPainter::mouseMoveEvent(tgt::MouseEvent* e) {
    switch(pressedSlider_) {
    case GAMMA_SLIDER:
        gammaSlider_->setPosition(e->coord());
        canvas_->update();
    break;
    case DOMAIN_LEFT_SLIDER:
        domainSlider_->setLeftPosition(e->coord());
        canvas_->update();
    break;
    case DOMAIN_RIGHT_SLIDER:
        domainSlider_->setRightPosition(e->coord());
        canvas_->update();
    break;
    case DOMAIN_BOTH_SLIDER:
        tgtAssert(mousePressedPosition_ != -1, "mousePressedPosition has not been set!");
        tgtAssert(domainLeftPressedPosition_ != -1, "domainLeftPressedPosition has not been set!");
        tgtAssert(domainRightPressedPosition_ != -1, "domainRightPressedPosition has not been set!");
        domainSlider_->setLeftPosition( tgt::ivec2(domainLeftPressedPosition_ - (mousePressedPosition_ -e->coord().x)));
        domainSlider_->setRightPosition(tgt::ivec2(domainRightPressedPosition_ - (mousePressedPosition_ -e->coord().x)));
        canvas_->update();
    break;
    default:
        // do nothing
    break;
    }
    //creates and shows tooltip
    createInfoToolTip(QPoint(e->coord().x,e->coord().y));
}

void TransFunc1DPropertyWidgetPainter::mouseReleaseEvent(tgt::MouseEvent* e) {
    pressedSlider_ = NO_SLIDER;
    mousePressedPosition_ = -1;
    domainLeftPressedPosition_ = -1;
    domainRightPressedPosition_ = -1;
    emit hideInfoToolTip();
    emit interaction(false);
}

void TransFunc1DPropertyWidgetPainter::mouseDoubleClickEvent(tgt::MouseEvent* e) {
     if(renderGammaSlider_ && gammaSlider_->isHit(e->coord())) {
         gammaSlider_->setGammaValue(1.f);
         gammaSlot(1.f);
         canvas_->update();
     }
}

void TransFunc1DPropertyWidgetPainter::createInfoToolTip(QPoint mousePos) {
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
    float interval = domainSlider_->getDomainValue().y - domainSlider_->getDomainValue().x;
    int percision = 3;
    if(interval > 1000.f)
        percision = 0;
    else if(interval > 100.f)
        percision = 1;
    else if(interval > 10.f)
        percision = 2;

    switch(pressedSlider_){
    case GAMMA_SLIDER:
        //position
        ittPos.setX(mousePos.x());
        ittPos.setY(-bttMetr.height()+3);
        //gamma
        qText = QString(1,QChar(0x03b3)) + QString("=");
        ittWidth += bttMetr.width(qText);
        //value
        ittText << ftos(gammaSlider_->getGammaValue(),2) << "</b>";
        ittWidth += bttMetr.width(ftos(gammaSlider_->getGammaValue(),2).c_str());
        qText = QString("<b>") + qText + QString(ittText.str().c_str());
        break;
    case DOMAIN_LEFT_SLIDER:
        ittPos.setX(canvas_->getPhysicalWidth() / 2);
        ittPos.setY(canvas_->getPhysicalHeight()+10);
        ittText << "<nobr>Domain: [ <b>" << ftos(domainSlider_->getDomainValue().x,percision) << "</b> , " <<
                                    ftos(domainSlider_->getDomainValue().y,percision) << " ]</nobr>";
        ittWidth += ttMetr.width("Domain: [  ,  ]");
        ittWidth += bttMetr.width(ftos(domainSlider_->getDomainValue().x,percision).c_str());
        ittWidth += ttMetr.width(ftos(domainSlider_->getDomainValue().y,percision).c_str());
        qText = QString(ittText.str().c_str());
        break;
    case DOMAIN_RIGHT_SLIDER:
        ittPos.setX(canvas_->getPhysicalWidth() / 2);
        ittPos.setY(canvas_->getPhysicalHeight()+10);
        ittText << "<nobr>Domain: [ " + ftos(domainSlider_->getDomainValue().x,percision) << " , <b>" <<
                                ftos(domainSlider_->getDomainValue().y,percision) << "</b> ]</nobr>";
        ittWidth += ttMetr.width("Domain: [  ,  ]");
        ittWidth += ttMetr.width(ftos(domainSlider_->getDomainValue().x,percision).c_str());
        ittWidth += bttMetr.width(ftos(domainSlider_->getDomainValue().y,percision).c_str());
        qText = QString(ittText.str().c_str());
        break;
    case DOMAIN_BOTH_SLIDER:
        ittPos.setX(canvas_->getPhysicalWidth() / 2);
        ittPos.setY(canvas_->getPhysicalHeight()+10);
        ittText << "<nobr>Domain: [ <b>" << ftos(domainSlider_->getDomainValue().x,percision) << "</b> , <b>" +
                                    ftos(domainSlider_->getDomainValue().y,percision) << "</b> ]</nobr>";
        ittWidth += ttMetr.width("Domain: [  ,  ]");
        ittWidth += bttMetr.width(ftos(domainSlider_->getDomainValue().x,percision).c_str());
        ittWidth += bttMetr.width(ftos(domainSlider_->getDomainValue().y,percision).c_str());
        qText = QString(ittText.str().c_str());
        break;
    default:
        return;
        break;
    }
    //shift info tool tip
    ittPos.setX(ittPos.x()-ittWidth/2-ittLeftBorder);
    //show tip
    emit showInfoToolTip(ittPos,qText);
}

void TransFunc1DPropertyWidgetPainter::gammaSlot(float gamma) {
    tf_->setGammaValue(gamma);
    emit changedGamma();
}
void TransFunc1DPropertyWidgetPainter::domainSlot(tgt::vec2 domain) {
    tf_->setDomain(domain);
    //move gamma slider
    gammaSlider_->setSubCanvasSize(tgt::ivec2(domainSlider_->getPosition().y-domainSlider_->getPosition().x+1));
    gammaSlider_->setSubCanvasOffset(tgt::ivec2(domainSlider_->getPosition().x));
    emit changedDomain();
}

//--------------------------------------------
//--------------------------------------------
//  zoom functions
//--------------------------------------------
//--------------------------------------------
void TransFunc1DPropertyWidgetPainter::zoomIn() {
    if(tf_) {
        float domCenter = (tf_->getDomain().x + tf_->getDomain().y) * 0.5f;
        float domSize = (tf_->getDomain().y - tf_->getDomain().x);
        float viewSize = visibleDomainValues_.y - visibleDomainValues_.x;
        viewSize *= 0.5f;
        if(viewSize < domSize)
            viewSize = domSize;

        visibleDomainValues_.x = domCenter - (viewSize * 0.5f);
        visibleDomainValues_.y = domCenter + (viewSize * 0.5f);

        domainSlider_->setValueRange(tgt::vec2(visibleDomainValues_.x,visibleDomainValues_.y));
        gammaSlider_->setSubCanvasSize(tgt::ivec2(domainSlider_->getPosition().y-domainSlider_->getPosition().x+1));
        gammaSlider_->setSubCanvasOffset(tgt::ivec2(domainSlider_->getPosition().x));
        //repaint
        canvas_->update();
    }
}

void TransFunc1DPropertyWidgetPainter::zoomOut() {
    if(tf_) {
        float domCenter = (tf_->getDomain().x + tf_->getDomain().y) * 0.5f;
        float domSize = (tf_->getDomain().y - tf_->getDomain().x);
        float viewSize = visibleDomainValues_.y - visibleDomainValues_.x;
        float viewCenter = domCenter;

        viewSize *= 2.0f;

        visibleDomainValues_.x = viewCenter - (viewSize * 0.5f);
        visibleDomainValues_.y = viewCenter + (viewSize * 0.5f);

        domainSlider_->setValueRange(tgt::vec2(visibleDomainValues_.x,visibleDomainValues_.y));
        gammaSlider_->setSubCanvasSize(tgt::ivec2(domainSlider_->getPosition().y-domainSlider_->getPosition().x+1));
        gammaSlider_->setSubCanvasOffset(tgt::ivec2(domainSlider_->getPosition().x));
        //repaint
        canvas_->update();
    }
}

void TransFunc1DPropertyWidgetPainter::resetZoom() {
    if(tf_) {
        visibleDomainValues_.x = tf_->getDomain().x;
        visibleDomainValues_.y = tf_->getDomain().y;

        domainSlider_->setValueRange(tgt::vec2(visibleDomainValues_.x,visibleDomainValues_.y));
        gammaSlider_->setSubCanvasSize(tgt::ivec2(domainSlider_->getPosition().y-domainSlider_->getPosition().x+1));
        gammaSlider_->setSubCanvasOffset(tgt::ivec2(domainSlider_->getPosition().x));
        //repaint
        canvas_->update();
    }
}

float TransFunc1DPropertyWidgetPainter::getMinVisibleDomainValue() const {
    return visibleDomainValues_.x;
}

float TransFunc1DPropertyWidgetPainter::getMaxVisibleDomainValue() const {
    return visibleDomainValues_.y;
}

void TransFunc1DPropertyWidgetPainter::setVisibleDomainValues(tgt::vec2 visibleDomain) {
    visibleDomainValues_ = visibleDomain;

    domainSlider_->setValueRange(tgt::vec2(visibleDomainValues_.x,visibleDomainValues_.y));
    gammaSlider_->setSubCanvasSize(tgt::ivec2(domainSlider_->getPosition().y-domainSlider_->getPosition().x+1));
    gammaSlider_->setSubCanvasOffset(tgt::ivec2(domainSlider_->getPosition().x));
    //repaint
    canvas_->update();
}

} // namespace voreen
