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

#ifndef VRN_VOREENSETTINGSDIALOG_H
#define VRN_VOREENSETTINGSDIALOG_H

#include <QWidget>

#include "voreen/qt/voreenqtapi.h"

#include "voreen/qt/mainwindow/voreenqtmainwindow.h"

class QPushButton;
class QVBoxLayout;

namespace voreen {

class VRN_QT_API VoreenSettingsDialog : public QWidget {
Q_OBJECT

public:
    VoreenSettingsDialog(QWidget* parent, VoreenQtMainWindow* mainwindow = 0);

signals:
    void closeSettings();

private:
    QVBoxLayout* widgetLayout_;

    QPushButton* closeButton_;
    QPushButton* resetButton_;

    VoreenQtMainWindow* mainWindow_;

public slots:
    /**
     * Updates the visibility of all widgets according to their visibility in application or network mode.
     */
    void updateWidgetVisibilities(VoreenQtMainWindow::GuiMode mode);


private slots:
    void resetSettings();
};

} //namespace voreen

#endif // VRN_VOREENSETTINGSDIALOG_H

