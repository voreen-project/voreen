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

#include "voreen/qt/widgets/clickablelabel.h"

#include <QPalette>

namespace voreen {

ClickableLabel::ClickableLabel(const char* text, QWidget* parent, Qt::WindowFlags f, Qt::CursorShape hoverCursor)
    : QLabel (text, parent, f),
      hoverCursor_(hoverCursor)
{
    state = false;
    text_ = (std::string)text;

    QPalette palette = QPalette();
    palette.setColor(QPalette::Dark, QColor(140, 140, 140));
    setPalette(palette);

    setForegroundRole(QPalette::Dark);
    expanded_ = false;
}

void ClickableLabel::mousePressEvent(QMouseEvent* e) {
    state = true;
    emit pressed();
    e->ignore();
}

void ClickableLabel::mouseReleaseEvent(QMouseEvent* e) {
    state = false;
    emit released();
    emit clicked();
    e->ignore();
}

void ClickableLabel::mouseDoubleClickEvent(QMouseEvent* e) {
    state = false;
    emit doubleClicked();
    e->ignore();
}

void ClickableLabel::enterEvent(QEvent* e) {
    setCursor(hoverCursor_);
    setForegroundRole(QPalette::Text);
    e->ignore();
}

void ClickableLabel::leaveEvent(QEvent* e) {
    unsetCursor();
    if (!expanded_)
        setForegroundRole(QPalette::Dark);
    else
        setForegroundRole(QPalette::Text);
    e->ignore();
}

void ClickableLabel::toggleExpanded() {
    expanded_ = !expanded_;
    if (expanded_)
        setForegroundRole(QPalette::Text);
    else
        setForegroundRole(QPalette::Dark);
}

} // namespace voreen
