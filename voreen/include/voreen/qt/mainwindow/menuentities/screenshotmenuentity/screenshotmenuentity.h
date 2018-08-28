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

#ifndef VRN_SCREENSHOTMENUENTITY_H
#define VRN_SCREENSHOTMENUENTITY_H

#include "../voreenqtmenuentity.h"
#include "voreen/core/network/workspace.h"
#include "voreen/qt/mainwindow/voreenqtmainwindow.h"

#include <QObject>
#include <map>

namespace voreen {

    class CanvasRenderer;
    class CanvasRendererScreenshotPlugin;
    class NetworkScreenshotPlugin;

/**
 * Menu entity for the screenshot tools
 * It can be used to take a screenshot of the entire network (@see NetworkScreenshotPlugin) or of each
 * existing canvas of the network (@see CanvasRendererScreenshotPlugin).
 */
class ScreenshotMenuEntity : public QObject, public VoreenQtMenuEntity, public UsesWorkspace, public ProcessorNetworkObserver {
    Q_OBJECT
    //friend class VoreenVEMainWindow;

public:
    /** Constructor */
    ScreenshotMenuEntity();
    /** Destructor */
    virtual ~ScreenshotMenuEntity();

    //------------------------------
    //  VoreenQtMenuEntity
    //------------------------------
    virtual QIcon getIcon() const {return QIcon(":/qt/icons/screenshot.png");}
    virtual std::string getName() const {return "Screenshot";}
    virtual std::string getShortCut() const {return "";}
    virtual MenuCategory getMenuCategory() const {return MC_BASIC_TOOL;}
    virtual bool getDefaultVisibilityInApplicationMode() const {return false;}
    virtual void setMainWindow(VoreenQtMainWindow* mainWindow);

    //------------------------------
    //  UsesWorkspace
    //------------------------------
    virtual void setWorkspace(Workspace* workspace);

    //------------------------------
    //  ProcessorNetworkObserver
    //------------------------------
    /** Nothing to do by default */
    virtual void networkChanged() {}
    /** Adjust menu, if a Canvasrenderer has been added. */
    virtual void processorAdded(const Processor* processor);
    /** Adjust menu and delete plugin, if a Canvasrenderer is about to be removed. */
    virtual void processorRemoved(const Processor* processor);
    /** Update plugin name and meta data, if a CanvasRenderer has been renamed. */
    virtual void processorRenamed(const Processor* processor, const std::string& prevName);

    //------------------------------
    //  Menu and Main-Functions
    //------------------------------
protected:
    virtual QAction* createMenuAction();
protected slots:
    void adjustScreenshotMenuSlot(VoreenQtMainWindow::GuiMode mode);
    void screenshotActionTriggeredSlot();

    //------------------------------
    //  Members
    //------------------------------
private:
     NetworkScreenshotPlugin* networkScreenshotPlugin_; ///< Network Plugin used for network screenshots
     std::map<const CanvasRenderer*,CanvasRendererScreenshotPlugin*> canvasRendererScreenshotPluginMap_;     ///< map between renderer and plugins
     Workspace* currentWorkspace_;                      ///< Current workspace
};

} // namespace

#endif // VRN_SCREENSHOTMENUENTITY_H
