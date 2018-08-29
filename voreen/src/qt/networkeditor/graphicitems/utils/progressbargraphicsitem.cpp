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

#include "voreen/qt/networkeditor/graphicitems/utils/progressbargraphicsitem.h"

#include "voreen/qt/networkeditor/graphicitems/core/portownergraphicsitem.h"
#include "voreen/qt/networkeditor/styles/nwestyle_base.h"

#include "voreen/core/voreenapplication.h"
#include "voreen/core/datastructures/callback/lambdacallback.h"
#include "voreen/core/utils/commandqueue.h"

#include <QApplication>
#include <QGraphicsScene>
#include <QLinearGradient>
#include <QPainter>
#include <QThread>

namespace voreen {

ProgressBarGraphicsItem::ProgressBarGraphicsItem(PortOwnerGraphicsItem* parent, const QPointF& center, qreal width, qreal height)
    : NWEBaseGraphicsItem(parent->getNetworkEditor())
    , backgroundColor1_(3, 3, 3)
    , backgroundColor2_(21, 21, 21)
{
    setParentItem(parent);
    resize(center, width, height);
    time_.start();

    upperForegroundColor1_ = QColor(176, 191, 217);
    upperForegroundColor2_ = QColor(60, 109, 194);
    lowerForegroundColor_ = QColor(22, 87, 199);
}

ProgressBarGraphicsItem::~ProgressBarGraphicsItem() {
    VoreenApplication::app()->getCommandQueue()->removeAll(this);
}

//---------------------------------------------------------------------------------------------------------------
//                  getter and setter
//---------------------------------------------------------------------------------------------------------------
qreal ProgressBarGraphicsItem::getWidth() const {
    return width_;
}

qreal ProgressBarGraphicsItem::getHeight() const {
    return height_;
}

//---------------------------------------------------------------------------------------------------------------
//                  nwebasegraphicsitem functions
//---------------------------------------------------------------------------------------------------------------
QRectF ProgressBarGraphicsItem::boundingRect() const{
    return currentStyle()->ProgressBarGI_boundingRect(this);
}

QPainterPath ProgressBarGraphicsItem::shape() const {
    return currentStyle()->ProgressBarGI_shape(this);
}

void ProgressBarGraphicsItem::initializePaintSettings() {
    currentStyle()->ProgressBarGI_initializePaintSettings(this);
}

void ProgressBarGraphicsItem::mainPaint(QPainter* painter, const QStyleOptionGraphicsItem* option, QWidget* widget, NWEItemSettings& setting) {
    currentStyle()->ProgressBarGI_paint(this,painter,option,widget,setting);
}

//---------------------------------------------------------------------------------------------------------------
//                  progressbar functions
//---------------------------------------------------------------------------------------------------------------


void ProgressBarGraphicsItem::resize(const QPointF& center, qreal width, qreal height) {
    width_ = width;
    height_ = height;
    qreal x = center.x() - width/2.0;
    qreal y = center.y() - height/2.0;

    setPos(x, y);
}

void ProgressBarGraphicsItem::show() {
    QGraphicsItem::show();
}

void ProgressBarGraphicsItem::hide() {
    QGraphicsItem::hide();
}

void ProgressBarGraphicsItem::forceUpdate() {
    update();
}

void ProgressBarGraphicsItem::update() {

    // Remove old commands.
    VoreenApplication::app()->getCommandQueue()->removeAll(this);

    // Setup command for the update
    LambdaFunctionCallback command([this] {
        if (time_.elapsed() > 500 || progress_ == 1.f || progress_ == 0.f) {
            NWEBaseGraphicsItem::update();
            scene()->update();
            //qApp->processEvents(); // responsible for crashes!
            time_.restart();
        }
    });

    bool isGuiThread = (QThread::currentThread() == QCoreApplication::instance()->thread());
    if (!isGuiThread) {
        VoreenApplication::app()->getCommandQueue()->enqueue(this, command);
    }
    else {
        command.exec();
    }
}


} // namespace
