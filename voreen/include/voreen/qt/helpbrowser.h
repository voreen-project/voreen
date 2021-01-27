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

#ifndef VRN_HELPBROWSER_H
#define VRN_HELPBROWSER_H

#include <QMainWindow>
#include <QUrl>
#include "voreen/qt/voreenqtapi.h"

class QTextBrowser;

namespace voreen {

/**
 * A simple HTML browser
 */
class VRN_QT_API HelpBrowser : public QMainWindow {
Q_OBJECT
public:
    HelpBrowser(const QUrl& startPage = QUrl(), const QString& title ="", QWidget* parent = 0);

    void createAndConnectToolbar();

private:
    QTextBrowser* browser_;
};

} // namespace voreen

#endif // VRN_HELPBROWSER_H
