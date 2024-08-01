/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2024 University of Muenster, Germany,                        *
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

#ifndef VRN_QPROPERTYWIDGETWITHTOOLWINDOW_H
#define VRN_QPROPERTYWIDGETWITHTOOLWINDOW_H

#include "voreen/qt/widgets/property/qpropertywidget.h"

namespace voreen {

class VoreenToolWindow;

/**
 * This class implements the QPropertyWidget interface and adds the functionality of a
 * VoreenToolWindow to it.
 * The class handles the construction of the ToolWindow.
 * @note the construction of the widget inside the tool window is pure virtual.
 */
class VRN_QT_API QPropertyWidgetWithToolWindow : public QPropertyWidget {
public:
    QPropertyWidgetWithToolWindow(Property* prop, QWidget* parent = 0, bool showNameLabel = true, bool isToolWindowResizable = true);
    virtual ~QPropertyWidgetWithToolWindow();

    /**
     * Stores the state of the tool window.
     * @override PropertyWidget
     */
    virtual MetaDataBase* getWidgetMetaData() const;

protected:
    /** Creates the widget for the tool window. */
    virtual QWidget* createToolWindowWidget() = 0;
    /** Sets custom tool window properties. */
    virtual void customizeToolWindow() = 0;

    /** Creates the tool window. */
    void createToolWindow(Qt::DockWidgetArea area, const QString& titlePostfix = "", const int& initialWidth = -1, const int& initialHeight = -1);
    /** Queries the property's meta data for the tool window's visibility state. */
    bool isToolWindowVisibleOnStartup() const;

    /** Toggles the tool window visibility */
    void toggleToolWindow();

    //--------------------
    //  Member
    //--------------------
    VoreenToolWindow* toolWindow_;      ///< the associated tool window
private:
    bool isToolWindowResizable_;        ///< is the window resizeable?
};

} // namespace

#endif // VRN_QPROPERTYWIDGETWITHTOOLWINDOW_H
