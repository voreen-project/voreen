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

#include "voreen/qt/mainwindow/menuentities/animationmenuentity.h"

#include "voreen/qt/mainwindow/voreenqtmainwindow.h"

#include "voreen/qt/animation/animationeditor.h"

namespace voreen {

AnimationMenuEntity::AnimationMenuEntity()
    : VoreenQtMenuEntity()
    , editor_(0)
{}

AnimationMenuEntity::~AnimationMenuEntity() {
}

QWidget* AnimationMenuEntity::createWidget() const {
    editor_ = new AnimationEditor(mainWindow_->getNetworkEvaluator(), WsHndlr.getWorkspace(), mainWindow_);
    WsHndlr.registerWorkspaceUsingWidget(editor_);
    editor_->resize(925, 400);
    return editor_;
}

void AnimationMenuEntity::deinitialize() {
    if(editor_)
        WsHndlr.unregisterWorkspaceUsingWidget(editor_);
    VoreenQtMenuEntity::deinitialize();
}

} //namespace
