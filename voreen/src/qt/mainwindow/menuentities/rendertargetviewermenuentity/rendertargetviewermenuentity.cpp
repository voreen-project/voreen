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

#include "voreen/qt/mainwindow/menuentities/rendertargetviewermenuentity/rendertargetviewermenuentity.h"
#include "voreen/qt/mainwindow/menuentities/rendertargetviewermenuentity/rendertargetviewer.h"

#include "voreen/qt/widgets/voreentoolwindow.h"
#include "voreen/qt/mainwindow/voreenqtmainwindow.h"
#include "voreen/qt/networkeditor/networkeditor.h"

#include <QLayout>

namespace voreen {

RenderTargetViewerMenuEntity::RenderTargetViewerMenuEntity()
    : VoreenQtMenuEntity()
    , renderTargetViewer_(0)
{}

RenderTargetViewerMenuEntity::~RenderTargetViewerMenuEntity() {
}

QWidget* RenderTargetViewerMenuEntity::createWidget() const {
    renderTargetViewer_ = new RenderTargetViewer();
    renderTargetViewer_->setEvaluator(mainWindow_->getNetworkEvaluator());
    renderTargetViewer_->setMinimumSize(200, 200);
    QObject::connect(mainWindow_->getNetworkEditor(), SIGNAL(processorsSelected(const QList<Processor*>&)),
                renderTargetViewer_, SLOT(processorsSelected(const QList<Processor*>&)));
    WsHndlr.registerWorkspaceUsingWidget(renderTargetViewer_);
    return renderTargetViewer_;
}

void RenderTargetViewerMenuEntity::initialize() {
    VoreenQtMenuEntity::initialize();
    renderTargetViewer_->initialize();
    toolWindow_->setMouseTracking(true);
    toolWindow_->widget()->setContentsMargins(0,0,0,0);
    toolWindow_->widget()->layout()->setContentsMargins(0,0,0,0);
    toolWindow_->resize(500, 500);
}

void RenderTargetViewerMenuEntity::deinitialize() {
    if(renderTargetViewer_) {
        if(VoreenQtWorkspaceHandler::isInited())
            WsHndlr.unregisterWorkspaceUsingWidget(renderTargetViewer_);
        delete renderTargetViewer_;
        renderTargetViewer_ = 0;
    }
    VoreenQtMenuEntity::deinitialize();
}

} //namespace
