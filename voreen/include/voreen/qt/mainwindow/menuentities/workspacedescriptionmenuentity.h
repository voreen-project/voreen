/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2019 University of Muenster, Germany,                        *
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

#ifndef VRN_WORKSPACEDESCRIPTIONMENUENTITY_H
#define VRN_WORKSPACEDESCRIPTIONMENUENTITY_H

#include "voreenqtmenuentity.h"
#include "voreen/core/network/workspace.h"

#include <QTextEdit>

namespace voreen {

class VolumeViewer;

/**
 * Menu entity for the workspace description
 */
class WorkspaceDescriptionMenuEntity : public QObject, public UsesWorkspace, public VoreenQtMenuEntity {
    Q_OBJECT
    friend class VoreenQtMainWindow;

public:
    WorkspaceDescriptionMenuEntity();
    virtual ~WorkspaceDescriptionMenuEntity();
    virtual QIcon getIcon() const {return QIcon(":/qt/icons/document-icon.png");}
    virtual std::string getName() const {return "Workspace Description";}
    virtual std::string getShortCut() const {return "";}
    virtual MenuCategory getMenuCategory() const {return MC_BASIC_SIDEBAR;}
    virtual Qt::DockWidgetAreas getAllowedDockWidgetAreas() const {return static_cast<Qt::DockWidgetAreas>(Qt::LeftDockWidgetArea | Qt::RightDockWidgetArea);}
    virtual Qt::DockWidgetArea getInitialDockWidgetArea() const {return Qt::LeftDockWidgetArea;}
    virtual bool getDefaultVisibilityInApplicationMode() const {return true;}

    const std::string getDefaultWorkspaceDescription() const {
                        return  "<font color=\"#808080\"><h4>no workspace description</h4> <p>Select \"Right click -> edit\" to edit " \
                                "the workspace description.</p><p>Simple Html commands are supported.</p></font>";}

    /// @see UsesWorkspace
    virtual void setWorkspace(Workspace* workspace);
protected:
    virtual QWidget* createWidget() const;
private:
    mutable QTextEdit* wdEditor_;       //< workspace description editor
    mutable QWidget* bottomWDWidget_;   //< widget containg the save and cancel buttons
    Workspace* currentWorkspace_;       //< current workspace
private slots:
     void stateChangedWDEdit(bool checked) const;          //switch between WDEditor enabled and disabled. Stores the WD in the workspace.
     void undoWDEdit();                              //undos changes to the description
     void showWDEditorContextMenu(const QPoint &pt); //opens the context menu and adds the edit mode action
};

} // namespace

#endif // VRN_WORKSPACEDESCRIPTIONMENUENTITY_H
