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

#include "voreen/qt/networkeditor/styles/nwestyle_material_dark.h"

//qt
#include <QBrush>
#include <QColor>

namespace voreen{

/*********************************************************************
 *                       General Color Defines
 ********************************************************************/
//editor
QBrush NWEStyle_Material_Dark::getBackgroundBrush() const {
    return QBrush(QColor(40, 40, 40, 255));
}

//general graphicsitem
QColor NWEStyle_Material_Dark::getSelectionColor() const {
    return QColor(90, 90, 90, 255);
}
QColor NWEStyle_Material_Dark::getHoverColor() const {
    return QColor(80, 80, 80, 255);
}

//port
//processor
QColor NWEStyle_Material_Dark::getProcessorColor1() const {
    return QColor(60, 60, 60, 255);
}
//portarrow
QColor NWEStyle_Material_Dark::getPortArrowColor() const {
    return QColor(160,160,160,255);
}

QColor NWEStyle_Material_Dark::getPropertyLinkArrowColor() const {
    return QColor(245,245,245,255);
}
QColor NWEStyle_Material_Dark::getPortSizeLinkArrowColor() const {
    return QColor(255, 255, 255, 255);
}

/*********************************************************************
 *                       General Functions
 ********************************************************************/
NWEStyle_Material_Dark::NWEStyle_Material_Dark(NetworkEditor* networkeditor)
    : NWEStyle_Material(networkeditor)
{}

NWEStyle_Material_Dark::~NWEStyle_Material_Dark()
{}

void NWEStyle_Material_Dark::ProcessorGI_initializePaintSettings(ProcessorGraphicsItem* item) {
    //set text
    item->getOrCreateNameLabel()->setPlainText(item->getGuiName());
    item->getOrCreateNameLabel()->setFont(QFont("Helvetica", nwe_->getProcessorFontSizeScale()/10.f));
    item->getOrCreateNameLabel()->setDefaultTextColor(Qt::white);
}

} // namespace voreen
