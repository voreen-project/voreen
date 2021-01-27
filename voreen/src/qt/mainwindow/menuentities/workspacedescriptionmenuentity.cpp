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

#include "voreen/qt/mainwindow/menuentities/workspacedescriptionmenuentity.h"

#include "voreen/qt/mainwindow/voreenqtmainwindow.h"

#include "voreen/qt/widgets/volumeviewer.h"

#include <QPushButton>
#include <QHBoxLayout>
#include <QVBoxLayout>
#include <QAction>
#include <QMenu>

namespace voreen {

WorkspaceDescriptionMenuEntity::WorkspaceDescriptionMenuEntity()
    : VoreenQtMenuEntity()
    , currentWorkspace_(0)
{
    WsHndlr.registerWorkspaceUsingWidget(this);
}

WorkspaceDescriptionMenuEntity::~WorkspaceDescriptionMenuEntity() {
    WsHndlr.unregisterWorkspaceUsingWidget(this);
}

QWidget* WorkspaceDescriptionMenuEntity::createWidget() const {
    //editor
    QWidget* wdWidget = new QWidget();

    QVBoxLayout* wdLayout = new QVBoxLayout(wdWidget);
        wdLayout->setMargin(0);
    wdEditor_ = new QTextEdit();
    wdEditor_->setContextMenuPolicy(Qt::CustomContextMenu);
    connect(wdEditor_,SIGNAL(customContextMenuRequested(const QPoint&)),this,SLOT(showWDEditorContextMenu(const QPoint &)));
    wdLayout->addWidget(wdEditor_);
    //bottom
    bottomWDWidget_ = new QWidget();
    QHBoxLayout* bottomWDLayout = new QHBoxLayout(bottomWDWidget_);
    QPushButton* saveWD = new QPushButton("Save Changes");
    connect(saveWD,SIGNAL(clicked(bool)),this,SLOT(stateChangedWDEdit(bool)));
    QPushButton* cancelWD = new QPushButton("Cancel");
    connect(cancelWD,SIGNAL(clicked(bool)),this,SLOT(undoWDEdit()));
    bottomWDLayout->setContentsMargins(4,0,4,6);
    bottomWDLayout->addStretch(1);
    bottomWDLayout->addWidget(saveWD);
    bottomWDLayout->addStretch(1);
    bottomWDLayout->addWidget(cancelWD);
    bottomWDLayout->addStretch(1);
    wdLayout->addWidget(bottomWDWidget_);

    stateChangedWDEdit(false);

    return wdWidget;
}

void WorkspaceDescriptionMenuEntity::setWorkspace(Workspace* workspace) {
    currentWorkspace_ = workspace;
    if(workspace) {
        //workspace description
        if(!currentWorkspace_->getDescription().empty()) {
            wdEditor_->setPlainText(QString(currentWorkspace_->getDescription().c_str()));
        } else {
            wdEditor_->setPlainText(QString::fromStdString(getDefaultWorkspaceDescription()));
        }
        stateChangedWDEdit(false);
    }

}
//---------------------------------------------------------------------------
//          Slots
//---------------------------------------------------------------------------
void WorkspaceDescriptionMenuEntity::stateChangedWDEdit(bool checked) const {
    if(checked) {
        wdEditor_->setReadOnly(false);
        bottomWDWidget_->setVisible(true);
        QString tmp = QString::fromStdString(getDefaultWorkspaceDescription());
        if(currentWorkspace_ && !currentWorkspace_->getDescription().empty())
            tmp = currentWorkspace_->getDescription().c_str();
        wdEditor_->setPlainText(tmp);
    } else {// unchecked
        wdEditor_->setReadOnly(true);
        bottomWDWidget_->setVisible(false);
        QString tmp = wdEditor_->toPlainText();
        if(currentWorkspace_ && (getDefaultWorkspaceDescription().compare(tmp.toStdString()) != 0) && (currentWorkspace_->getDescription().compare(tmp.toStdString()) != 0))
            currentWorkspace_->setDescription(tmp.toStdString());
        wdEditor_->setHtml(tmp);
    }
}

void WorkspaceDescriptionMenuEntity::undoWDEdit() {
    if(currentWorkspace_) {
        if(!currentWorkspace_->getDescription().empty())
            wdEditor_->setPlainText(currentWorkspace_->getDescription().c_str());
        else
            wdEditor_->setPlainText(QString::fromStdString(getDefaultWorkspaceDescription()));
    }
    stateChangedWDEdit(false);
}

void WorkspaceDescriptionMenuEntity::showWDEditorContextMenu(const QPoint &pt) {
    QMenu *menu = wdEditor_->createStandardContextMenu();
    QAction* editAct = new QAction("edit",menu);
    editAct->setCheckable(true);
    editAct->setChecked(!wdEditor_->isReadOnly());
    connect(editAct,SIGNAL(toggled(bool)),this,SLOT(stateChangedWDEdit(bool)));
    menu->addSeparator();
    menu->addAction(editAct);
    menu->exec(wdEditor_->mapToGlobal(pt));
    delete menu;
}

} //namespace
