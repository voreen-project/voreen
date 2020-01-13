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

#include "voreen/qt/mainwindow/voreenqtsplashscreen.h"

#include <QPainter>
#include <QPen>

namespace voreen {

VoreenQtSplashScreen::VoreenQtSplashScreen(std::string filepath, bool usesProgressbar, QColor progressbarMessageColor,
                                          QColor progressbarBoarderColor, QColor progressbarForgroundColor,
                                          QPoint progressbarLLF, QSize  progressbarSize)
    : QSplashScreen(QPixmap(filepath.c_str()))
    , usesProgressbar_(usesProgressbar)
    , progressbarMessage_("")
    , progressbarValue_(0.f)
    , progressbarMessageColor_(progressbarMessageColor)
    , progressbarBoarderColor_(progressbarBoarderColor)
    , progressbarForgroundColor_(progressbarForgroundColor)
    , progressbarLLF_(progressbarLLF)
    , progressbarSize_(progressbarSize)
{
    setWindowFlags(windowFlags() | Qt::WindowStaysOnTopHint);
}

void VoreenQtSplashScreen::drawContents(QPainter* painter) {
    // progressbar
    if(usesProgressbar_) {
        //define variabels
        QLinearGradient gradient(progressbarLLF_.x(), progressbarLLF_.y(),
                                 progressbarLLF_.x()+progressbarSize_.width()*progressbarValue_, progressbarLLF_.x()+progressbarSize_.height());
        progressbarForgroundColor_.setAlpha(75);
        gradient.setColorAt(0, progressbarForgroundColor_);
        progressbarForgroundColor_.setAlpha(175);
        gradient.setColorAt(1, progressbarForgroundColor_);
        QRect r = rect();
        QFont font = painter->font();
        r.setRect(progressbarLLF_.x(), progressbarLLF_.y(), progressbarSize_.width(), progressbarSize_.height());
        //draw boarder
        QPen pen(progressbarBoarderColor_);
        painter->setPen(pen);
        painter->drawRect(r);
        //draw progress
        painter->setBrush(gradient);
        pen.setColor(QColor(0,0,0,0));
        painter->setPen(pen);
        r.setRect(progressbarLLF_.x()+1, progressbarLLF_.y()+1,
                 (progressbarSize_.width()-1)*progressbarValue_, progressbarSize_.height()-1);
        painter->drawRect(r);
        //draw message
        r.setRect(progressbarLLF_.x()+5, progressbarLLF_.y(), progressbarSize_.width()-5, progressbarSize_.height());
        font.setPointSize(progressbarSize_.height()/2);
        painter->setFont(font);
        pen.setColor(progressbarMessageColor_);
        painter->setPen(pen);
        painter->drawText(r, Qt::AlignLeft | Qt::AlignVCenter, progressbarMessage_);
    }
}

void VoreenQtSplashScreen::updateProgressMessage(const QString& message, qreal progress) {
    progressbarMessage_ = message;
    if(0.f <= progress && progress <= 1.f)
        progressbarValue_ = progress;
    QSplashScreen::repaint();
}

void VoreenQtSplashScreen::updateProgressValue(qreal progress) {
    if(0.f <= progress && progress <= 1.f)
        progressbarValue_ = progress;
    QSplashScreen::repaint();
}

 } //end namespace
