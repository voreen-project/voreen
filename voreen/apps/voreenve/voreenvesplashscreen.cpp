/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2021 University of Muenster, Germany,                        *
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

#include "voreenvesplashscreen.h"

#include "voreen/core/version.h"

#include <QPainter>
#include <QPen>
#include <QPixmap>

namespace voreen {

VoreenVESplashScreen::VoreenVESplashScreen(std::string filepath)
    : VoreenQtSplashScreen(filepath, true, Qt::white, Qt::darkGray,
                           QColor(255, 192, 0), QPoint(219, 253), QSize(255, 19))
{}

void VoreenVESplashScreen::drawContents(QPainter* painter) {
    //draws progressbar
    VoreenQtSplashScreen::drawContents(painter);
    // version
    QPen pen( Qt::white );
    QRect r = rect();
    QFont font = painter->font();
    font.setPointSize(11);
    painter->setFont(font);
    r.setRect(20, 94, r.width() - 10, r.height() - 10);
    std::string version = "Version " + VoreenVersion::getVersion();
    painter->setPen(pen);
    painter->drawText(r, Qt::AlignLeft, version.c_str());

    // url
    font.setPointSize(9);
    painter->setFont(font);
    pen.setColor(Qt::lightGray);
    painter->setPen(pen);
    r = rect();
    r.setRect(3, 0, 300, r.height() - 2);
    painter->drawText(r, Qt::AlignLeft | Qt::AlignBottom, QString("voreen.uni-muenster.de"));
}

} //end namespace
