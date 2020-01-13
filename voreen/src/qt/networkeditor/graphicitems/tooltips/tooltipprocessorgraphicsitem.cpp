/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2020 University of Muenster, Germany,                        *
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

#include "voreen/qt/networkeditor/graphicitems/tooltips/tooltipprocessorgraphicsitem.h"
#include "voreen/qt/networkeditor/styles/nwestyle_base.h"

//gi
#include "voreen/qt/networkeditor/graphicitems/core/processorgraphicsitem.h"
//core
#include "voreen/core/processors/processor.h"

#include <QGraphicsTextItem>

namespace voreen {

ToolTipProcessorGraphicsItem::ToolTipProcessorGraphicsItem(ProcessorGraphicsItem* item)
    : ToolTipBaseGraphicsItem(item)
{
}

//---------------------------------------------------------------------------------------------------------------
//                  nwebasegraphicsitem functions
//---------------------------------------------------------------------------------------------------------------
QRectF ToolTipProcessorGraphicsItem::boundingRect() const {
    return currentStyle()->ToolTipProcessorGI_boundingRect(this);
}

QPainterPath ToolTipProcessorGraphicsItem::shape() const {
    return currentStyle()->ToolTipProcessorGI_shape(this);
}

void ToolTipProcessorGraphicsItem::initializePaintSettings() {
    currentStyle()->ToolTipProcessorGI_initializePaintSettings(this);
}

void ToolTipProcessorGraphicsItem::mainPaint(QPainter* painter, const QStyleOptionGraphicsItem* option, QWidget* widget, NWEItemSettings& setting) {
    currentStyle()->ToolTipProcessorGI_paint(this,painter,option,widget,setting);
}

//---------------------------------------------------------------------------------------------------------------
//                  tooltipbasegraphicsitem functions
//---------------------------------------------------------------------------------------------------------------
bool ToolTipProcessorGraphicsItem::updateToolTip(qreal x, qreal y) {
    //prepare size changes
    prepareGeometryChange();
    //set text
    Processor* processor = dynamic_cast<ProcessorGraphicsItem*>(parent())->getProcessor();
    std::stringstream strstr;
    strstr  << "<center><font><b>" << processor->getGuiName() << "</b></font></center>"
            << "Type: " << processor->getClassName() << "<br>";
    switch(processor->getProcessorState()){
    case Processor::PROCESSOR_STATE_NOT_INITIALIZED:
        strstr << "Status: <font color=\"#FF0000\">Not Initialized</font>";
        break;
    case Processor::PROCESSOR_STATE_NOT_READY:
        strstr << "Status: <font color=\"#FF9900\">" << processor->getNotReadyErrorMessage() << "</font>";
        break;
    case Processor::PROCESSOR_STATE_READY:
        strstr << "Status: <font color=\"#008800\">All Fine</font>";
        return false;
        break;
    default:
        tgtAssert(false,"Shouldn't get here");
        break;
    }
    toolTipText_->setHtml(QString::fromStdString(strstr.str()));
    //set position
    return ToolTipBaseGraphicsItem::updateToolTip(x,y);
}


} // namespace voreen
