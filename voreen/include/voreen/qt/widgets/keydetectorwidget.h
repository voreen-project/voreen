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

#ifndef VRN_KEYDETECTORWIDGET_H
#define VRN_KEYDETECTORWIDGET_H

#include <QLineEdit>
#include <QKeyEvent>

#include "tgt/event/keyevent.h"

namespace voreen {

class ModifierDetectorWidget : public QLineEdit {
Q_OBJECT
public:
    ModifierDetectorWidget(QWidget* parent = 0);

    static QString getStringForModifier(Qt::KeyboardModifiers modifier);

    static tgt::Event::Modifier getTGTModifierFromQt(Qt::KeyboardModifiers qtKey);
    static Qt::KeyboardModifiers getQtModifierFromTGT(tgt::Event::Modifier tgtKey);

signals:
    void modifierChanged(Qt::KeyboardModifiers);

protected:
    void keyPressEvent(QKeyEvent* event);
};

class KeyDetectorWidget : public QLineEdit {
Q_OBJECT
public:
    KeyDetectorWidget(bool ignoreModifierKeys = true, QWidget* parent = 0);

    static QString getStringForKey(int key, bool ignoreModifierKeys = true);

    static tgt::KeyEvent::KeyCode getTGTKeyFromQt(int qtKey);
    static int getQtKeyFromTGT(tgt::KeyEvent::KeyCode tgtKey);

signals:
    void keyChanged(int);

protected:
    void keyPressEvent(QKeyEvent* event);

    bool ignoreModifierKeys_;
};

} // namespace

#endif // VRN_KEYDETECTORWIDGET_H
