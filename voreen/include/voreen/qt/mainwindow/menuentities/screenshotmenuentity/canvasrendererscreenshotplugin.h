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

#ifndef VRN_CANVASRENDERERSCREENSHOTPLUGIN_H
#define VRN_CANVASRENDERERSCREENSHOTPLUGIN_H

#include "voreen/qt/mainwindow/menuentities/screenshotmenuentity/screenshotpluginbase.h"

namespace voreen {

    class CanvasRenderer;

/**
 * This plugin is used by the ScreenshotMenuEntity to save images of a CanvasRenderer.
 */
class VRN_QT_API CanvasRendererScreenshotPlugin : public ScreenshotPluginBase {
    Q_OBJECT
    friend class CanvasRendererWidget;
    friend class ScreenshotMenuEntity;
public:
    /** Constructor */
    CanvasRendererScreenshotPlugin(QWidget* parent, Workspace* workspace, CanvasRenderer* canvasRenderer);

protected:
    /** @override */
    virtual void saveScreenshot(const QString& filename, int width, int height);
    /** Called from ScreenshotMenuEntity. */
    void canvasRendererHasBeenRenamed();
    //-----------------
    //  Callbacks
    //-----------------
public slots:
    virtual void nativeSizeHasChangedSlot();

    //-----------------
    //  Members
    //-----------------
protected:
    CanvasRenderer* canvasRenderer_;    ///< pointer to the used editor
};

} // namespace voreen

#endif // VRN_CANVASRENDERERSCREENSHOTPLUGIN_H
