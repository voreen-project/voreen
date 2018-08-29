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

#ifndef VRN_VOREENQTSPLASHSCREEN_H
#define VRN_VOREENQTSPLASHSCREEN_H

#include "voreen/qt/voreenqtapi.h"

#include <QSplashScreen>

class QPainter;
class QPixmap;

namespace voreen {

/**
 * A splash screen shown during during start up.
 */
class VRN_QT_API VoreenQtSplashScreen : public QSplashScreen {
public:
    /**
    * Constructor
    * @param filepath path to the screen image
    */
    VoreenQtSplashScreen(std::string filepath, bool usesProgressbar = false, QColor progressbarMessageColor = Qt::white,
                        QColor progressbarBoarderColor = Qt::black, QColor progressbarForgroundColor = QColor(175, 175, 175),
                        QPoint progressbarLLF = QPoint(0, 0), QSize  progressbarSize = QSize(100,20));

protected:
    /**
     * Draws the progress bar. Can be overritten to add more informations.
     */
    virtual void drawContents(QPainter* painter);
public:
    void updateProgressMessage(const QString& message, qreal progress = -1.0f);
    void updateProgressValue(qreal progress);

protected:
    ///progressbar
    bool usesProgressbar_;       ///< screen is using a progressbar?
    QString progressbarMessage_; ///< massage inside the progressbar
    qreal progressbarValue_;     ///< value of progres inside progressbar

    QColor progressbarMessageColor_;    ///< spashscreen specific values
    QColor progressbarBoarderColor_;    ///< spashscreen specific values
    QColor progressbarForgroundColor_;  ///< spashscreen specific values
    QPoint progressbarLLF_;             ///< spashscreen specific values
    QSize  progressbarSize_;            ///< spashscreen specific values
};

} //end namespace

#endif
