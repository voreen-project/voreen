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

#ifndef VRN_VOREENQtWORKSPACEHANDLER_H
#define VRN_VOREENQtWORKSPACEHANDLER_H

#include "voreen/core/network/processornetwork.h"
#include "voreen/core/network/workspace.h"
#include "voreen/core/network/networkevaluator.h"

#include "voreen/qt/voreenqtapi.h"

#include <QObject>

#include <set>

namespace voreen {

class VoreenQtMainWindow;
class NetworkEvaluator;
class PropertyListWidget;
class NetworkEditor;
class ProcessorNetwork;
class ProcessorListWidget;

/**
 * Singleton that handles the global workspace used in the VoreenQt application.
 *
 * Provides functions for loading and saving of workspaces, and propagates changes to registered entities.
 *
 * @see VoreenQtMainWindow::initialize
 */
class VRN_QT_API VoreenQtWorkspaceHandler : public QObject, public tgt::Singleton<VoreenQtWorkspaceHandler>,
                                            public WorkspaceObserver, public NetworkEvaluatorObserver {
    Q_OBJECT

    friend class VoreenQtMainWindow;
    friend class tgt::Singleton<VoreenQtWorkspaceHandler>;

public:
    Workspace* getWorkspace() const;

    void registerNetworkEvaluator(NetworkEvaluator* evaluator);
    void registerMainWindow(VoreenQtMainWindow* mainWindow);
    void registerNetworkEditorWidget(NetworkEditor* networkEditorWidget);
    void registerPropertyListWidget(PropertyListWidget* propertyListWidget);
    void registerWorkspaceUsingWidget(UsesWorkspace* workspaceUser);
    void unregisterWorkspaceUsingWidget(UsesWorkspace* workspaceUser);

public slots:
    void newWorkspace();
    void openWorkspace();
    void openWorkspace(const QString& filename);
    bool saveWorkspace();
    bool saveWorkspaceAs();
    bool askSave();
    /** Checks if the *.vws file is readonly. */
    bool checkForReadOnly(const std::string& wsFilepath);

    bool rebuildShaders();

signals:
    void workspaceHasBeenModifiedSignal(); //emitted, if the workspace has been modified (see workspace observer)
    void guiModeChanged(int mode);

private:
    VoreenQtWorkspaceHandler();
    ~VoreenQtWorkspaceHandler();

    void propagateWorkspace(Workspace* workspace);

    void showWorkspaceErrors();
    void showNetworkErrors();

    /**
     * Called by the NetworkEvaluator (via observation).
     * Propagates the new network to the GUI components.
     */
    virtual void networkAssigned(ProcessorNetwork* newNetwork, ProcessorNetwork* previousNetwork);
    /// workspace observer functions
    void workspaceHasBeenModified() {emit workspaceHasBeenModifiedSignal();}
    void applicationModeConfigurationChanged() {workspaceHasBeenModified();}


    Workspace* workspace_;

    // The following entities are assigned to the workspace manager, but not owned by it
    NetworkEvaluator* evaluator_;
    VoreenQtMainWindow* mainWindow_;
    NetworkEditor* networkEditorWidget_;
    PropertyListWidget* propertyListWidget_;
    std::set<UsesWorkspace*> workspaceUsers_; //< all other widgets the workspace is propagated to

    bool readOnlyWorkspace_;
    QString lastUsedWorkspacePath_;

    static const std::string loggerCat_;
};

} // namespace

#define WsHndlr tgt::Singleton<VoreenQtWorkspaceHandler>::getRef()

#endif // VRN_VOREENQTWORKSPACEHANDLER_H
