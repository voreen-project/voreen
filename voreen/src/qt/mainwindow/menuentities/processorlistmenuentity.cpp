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

#include "voreen/qt/mainwindow/menuentities/processorlistmenuentity.h"

#include "voreen/qt/widgets/voreentoolwindow.h"
#include "voreen/qt/widgets/processorlistwidget.h"
#include "voreen/qt/mainwindow/voreenqtmainwindow.h"
#include "voreen/qt/networkeditor/networkeditor.h"

#include <QObject>

namespace voreen {

ProcessorListMenuEntity::ProcessorListMenuEntity()
    : VoreenQtMenuEntity()
    , processorListWidget_(0)
{}

ProcessorListMenuEntity::~ProcessorListMenuEntity() {
}

QWidget* ProcessorListMenuEntity::createWidget() const {
    processorListWidget_ = new ProcessorListWidget(0);
    processorListWidget_->setMinimumSize(200, 200);
    QObject::connect(mainWindow_->getNetworkEditor(), SIGNAL(processorsSelected(const QList<Processor*>&)),
                    processorListWidget_, SLOT(processorsSelected(const QList<Processor*>&)));
    QObject::connect(processorListWidget_, SIGNAL(processorAdded(QString)), mainWindow_->getNetworkEditor(), SLOT(processorAdded(QString)));

    return processorListWidget_;
}

} //namespace
