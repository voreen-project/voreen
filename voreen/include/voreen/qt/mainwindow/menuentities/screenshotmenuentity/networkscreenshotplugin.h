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

#ifndef VRN_NETWORKSCREENSHOTPLUGIN_H
#define VRN_NETWORKSCREENSHOTPLUGIN_H

#include "voreen/qt/mainwindow/menuentities/screenshotmenuentity/screenshotpluginbase.h"

namespace voreen {

    class NetworkEditor;

/**
 * This plugin is used by the ScreenshotMenuEntity to save images of the network.
 */
class VRN_QT_API NetworkScreenshotPlugin : public ScreenshotPluginBase {
    Q_OBJECT
public:
    /** Constructor */
    NetworkScreenshotPlugin(QWidget* parent, Workspace* workspace, NetworkEditor* networkEditorWidget);

protected:
    /** @override */
    virtual void saveScreenshot(const QString& filename, int width, int height);
    /** Shortcut to save the network image in native size. Unused yet. */
    void saveScreenshot(const QString& filename);

    //-----------------
    //  Callbacks
    //-----------------
public slots:
    virtual void nativeSizeHasChangedSlot();

    //-----------------
    //  Members
    //-----------------
protected:
    NetworkEditor* networkEditorWidget_;    ///< pointer to the used editor
};

} // namespace voreen

#endif // VRN_NETWORKSCREENSHOTPLUGIN_H
