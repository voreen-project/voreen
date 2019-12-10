/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2019 University of Muenster, Germany,                        *
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

#ifndef VRN_ENTEREXITPUSHBUTTON_H
#define VRN_ENTEREXITPUSHBUTTON_H

#include <QPushButton>

#include "voreen/qt/voreenqtapi.h"

namespace voreen {

class VRN_QT_API EnterExitPushButton : public QPushButton {
Q_OBJECT
public:
    EnterExitPushButton(QWidget* parent = 0);
    EnterExitPushButton(const QString& text, QWidget* parent = 0);
    EnterExitPushButton(const QIcon& icon, const QString& text, QWidget* parent = 0);

protected:
    void enterEvent(QEvent* event);
    void leaveEvent(QEvent* event);

signals:
    void enterEventSignal();
    void leaveEventSignal();
};

} // namespace

#endif // VRN_ENTEREXITPUSHBUTTON_H
