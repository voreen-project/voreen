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

#include "voreen/qt/mainwindow/menuentities/propertylistmenuentity.h"

#include "voreen/qt/widgets/voreentoolwindow.h"
#include "voreen/qt/widgets/propertylistwidget.h"
#include "voreen/qt/mainwindow/voreenqtmainwindow.h"
#include "voreen/qt/networkeditor/networkeditor.h"

#include <QObject>

namespace voreen {

PropertyListMenuEntity::PropertyListMenuEntity()
    : VoreenQtMenuEntity()
    , propertyListWidget_(0)
{}

PropertyListMenuEntity::~PropertyListMenuEntity() {
}

QWidget* PropertyListMenuEntity::createWidget() const {
    propertyListWidget_ = new PropertyListWidget(0);
    WsHndlr.registerWorkspaceUsingWidget(propertyListWidget_);
    QObject::connect(mainWindow_->getNetworkEditor(), SIGNAL(processorsSelected(const QList<Processor*>&)),
                    propertyListWidget_, SLOT(processorsSelected(const QList<Processor*>&)));
    return propertyListWidget_;
}

void PropertyListMenuEntity::deinitialize() {
    if(propertyListWidget_)
        WsHndlr.unregisterWorkspaceUsingWidget(propertyListWidget_);
    VoreenQtMenuEntity::deinitialize();
}

} //namespace
