/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany.                        *
 * Visualization and Computer Graphics Group <http://viscg.uni-muenster.de>        *
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

#ifndef CLICKABLELABEL_H
#define CLICKABLELABEL_H

#include <QLabel>
#include <QMouseEvent>

#include <string>

namespace voreen {

class ClickableLabel : public QLabel {
    Q_OBJECT
public:
    ClickableLabel(const char * text, QWidget * parent, Qt::WindowFlags f = 0,
                    Qt::CursorShape hoverCursor = Qt::PointingHandCursor);

    virtual void mousePressEvent(QMouseEvent* e);
    virtual void mouseReleaseEvent(QMouseEvent* e);
    virtual void mouseDoubleClickEvent(QMouseEvent* e);
    virtual void enterEvent(QEvent *e);
    virtual void leaveEvent(QEvent *e);

public slots:
    virtual void toggleExpanded();

signals:
    void clicked();
    void doubleClicked();
    void pressed();
    void released();

protected:

    bool state; ///< If true, we've had a pressed event and are waiting for a released event
    bool expanded_;
    std::string text_;
    Qt::CursorShape hoverCursor_;

private:
    ClickableLabel(const ClickableLabel& rhs);
    ClickableLabel& operator=(const ClickableLabel& rhs);
};

} // namespace voreen

#endif // CLICKABLELABEL_H
