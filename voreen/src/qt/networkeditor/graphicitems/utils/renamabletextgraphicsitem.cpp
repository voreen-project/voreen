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

#include "voreen/qt/networkeditor/graphicitems/utils/renamabletextgraphicsitem.h"

#include <QKeyEvent>
#include <QFont>

namespace voreen {

RenamableTextGraphicsItem::RenamableTextGraphicsItem(const QString& text, QGraphicsItem* parent)
    : QGraphicsTextItem(text, parent)
{
    setFlag(ItemIsSelectable, false);
    setDefaultTextColor(Qt::white);
    setFont(QFont("Helvetica", 10));
}

int RenamableTextGraphicsItem::type() const {
    return UserTypesRenamableTextGraphicsItem;
}

void RenamableTextGraphicsItem::setFocus(Qt::FocusReason focusReason) {
    // save old text
    previousText_ = toPlainText();
    QGraphicsTextItem::setFocus(focusReason);
}

void RenamableTextGraphicsItem::setPlainText(const QString& text) {
    previousText_ = text;
    QGraphicsTextItem::setPlainText(text);
}

void RenamableTextGraphicsItem::keyPressEvent(QKeyEvent* event) {
    event->accept();
    if (event->key() == Qt::Key_Escape) {
        // restore saved text
        setPlainText(previousText_);
        emit textChanged();
        emit renameFinished();
    }
    else if ((event->key() == Qt::Key_Return) && (event->modifiers() == Qt::NoModifier)) {
        previousText_ = toPlainText();
        setPlainText(previousText_); // clears the selection as textCursor().clearSelection() should
        emit textChanged();
        emit renameFinished();
    }
    else {
        QGraphicsTextItem::keyPressEvent(event);
        emit textChanged();
    }
}

void RenamableTextGraphicsItem::focusOutEvent(QFocusEvent*) {
    previousText_ = toPlainText();
    setPlainText(previousText_); // clears the selection as textCursor().clearSelection() should
    emit textChanged();
    emit renameFinished();
}

} // namespace
