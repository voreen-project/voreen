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

#ifndef VRN_VOREENTOOLWINDOW_H
#define VRN_VOREENTOOLWINDOW_H

#include <QWidget>
#include <QDockWidget>
#include "voreen/qt/voreenqtapi.h"

namespace voreen {

class VRN_QT_API VoreenToolWindowTitle : public QWidget {
Q_OBJECT
public:
    VoreenToolWindowTitle(QDockWidget* parent, bool dockable=true);
    virtual QSize sizeHint() const;
    virtual QSize minimumSizeHint() const;
protected:
    virtual void paintEvent(QPaintEvent *event);
    void mousePressEvent(QMouseEvent* event);
    void mouseDoubleClickEvent(QMouseEvent* event);
private:
    bool dockable_;
    QPixmap closeButton_;
    QPixmap undockButton_;
    QPixmap maximizeButton_;
};

//----------------------------------------------------------------------------------------------------------------

class VRN_QT_API VoreenToolWindow : public QDockWidget {
Q_OBJECT
public:
    /**
     * Creates a tool window used in VoreenVE. The window title is extracted from the action.
     * @note name is not the title. Its the Qt object name. Use updateTitle to set the title
     */
    VoreenToolWindow(QAction* action, QWidget* parent, QWidget* child, const QString& name = "", bool dockable=true);
    /**
     * Sets the window title.
     */
    void updateTitle(QString title, QIcon icon);

    /// Returns the action associated with this tool.
    QAction* action() const;
    /// Returns the widget that is wrapped by the tool windows.
    QWidget* child() const;

private:
    QAction* action_;
    QWidget* child_;
};

} // namespace

#endif // VRN_VOREENTOOLWINDOW_H
