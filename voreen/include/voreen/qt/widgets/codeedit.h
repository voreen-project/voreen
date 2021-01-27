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

#ifndef VRN_CODEEDIT_H
#define VRN_CODEEDIT_H

#include "voreen/qt/voreenqtapi.h"

#include <QTextBrowser>
#include <QObject>

class QPaintEvent;
class QResizeEvent;
class QSize;
class QWidget;

class StatusArea;

class VRN_QT_API CodeEdit : public QTextBrowser {
Q_OBJECT

public:
    CodeEdit(bool readonly, int fontSize, QWidget* parent = 0);
    virtual ~CodeEdit();

    void statusAreaPaintEvent(QPaintEvent *event);
    int statusAreaWidth();
    void updateFontSize(unsigned char s);
    void setShaderSource(const std::string& source);

public slots:
    void moveCursorToPosition(int line, int col = -1);
    void updateHighlight();

protected:
    void resizeEvent(QResizeEvent* event);
    void paintEvent(QPaintEvent* event);
    void wheelEvent(QWheelEvent* event);

private slots:
    void updateStatusAreaWidth();
    void highlightCurrentLine();
    void processCurrentLine();

private:
    bool textEditLock_;
    StatusArea* statusArea_;
};


class StatusArea : public QWidget {
public:
    StatusArea(CodeEdit* edit) : QWidget(edit), edit_(edit) {}

    QSize sizeHint() const {
        return QSize(edit_->statusAreaWidth(), 0);
    }

protected:
    void paintEvent(QPaintEvent *event) {
        edit_->statusAreaPaintEvent(event);
    }

private:
    CodeEdit* edit_;
};

#endif // VRN_CODEEDIT_H
