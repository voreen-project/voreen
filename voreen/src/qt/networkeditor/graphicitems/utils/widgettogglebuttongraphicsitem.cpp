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

#include "voreen/qt/networkeditor/graphicitems/utils/widgettogglebuttongraphicsitem.h"
#include "voreen/qt/networkeditor/styles/nwestyle_base.h"

#include "voreen/qt/networkeditor/networkeditor.h"
#include "voreen/qt/networkeditor/graphicitems/core/processorgraphicsitem.h"
#include "voreen/core/processors/processorwidget.h"

#include <QGraphicsSceneMouseEvent>
#include <QAction>

namespace voreen {

WidgetToggleButtonGraphicsItem::WidgetToggleButtonGraphicsItem(NWEBaseGraphicsItem* parent)
    : NWEBaseGraphicsItem(parent->getNetworkEditor())
    , widgets_(QList<ProcessorWidget*>())
{
    setParent(parent);
    setParentItem(parent);
    setAcceptHoverEvents(true);
    setToolTip(tr("show/hide processor widgets"));
}

WidgetToggleButtonGraphicsItem::~WidgetToggleButtonGraphicsItem() {}

//---------------------------------------------------------------------------------------------------------------
//                  getter and setter
//---------------------------------------------------------------------------------------------------------------
void WidgetToggleButtonGraphicsItem::setProcessorWidget(ProcessorWidget* widget) {
    widgets_.clear();

    if (widget)
        widgets_.push_back(widget);

    setContextMenuActions();
}

void WidgetToggleButtonGraphicsItem::setProcessorWidgets(const QList<ProcessorWidget*>& widgets) {
    widgets_ = widgets;
    setContextMenuActions();
}

QList<ProcessorWidget*> WidgetToggleButtonGraphicsItem::getWidgets() {
    return widgets_;
}

//---------------------------------------------------------------------------------------------------------------
//                  nwebasegraphicsitem functions
//---------------------------------------------------------------------------------------------------------------
QRectF WidgetToggleButtonGraphicsItem::boundingRect() const {
    return currentStyle()->WidgetToggleButtonGI_boundingRect(this);
}

QPainterPath WidgetToggleButtonGraphicsItem::shape() const {
    return currentStyle()->WidgetToggleButtonGI_shape(this);
}

void WidgetToggleButtonGraphicsItem::initializePaintSettings() {
    currentStyle()->WidgetToggleButtonGI_initializePaintSettings(this);
}

void WidgetToggleButtonGraphicsItem::mainPaint(QPainter* painter, const QStyleOptionGraphicsItem* option, QWidget* widget, NWEItemSettings& setting) {
    currentStyle()->WidgetToggleButtonGI_paint(this, painter, option, widget, setting);
}

void WidgetToggleButtonGraphicsItem::setContextMenuActions() {
    NWEBaseGraphicsItem::setContextMenuActions();
    foreach(ProcessorWidget* widget, widgets_) {
        QAction* action = new QAction(QString::fromStdString(widget->getProcessor()->getGuiName() + " Widget"), this);
        action->setCheckable(true);
        connect(action, SIGNAL(triggered()),dynamic_cast<PortOwnerGraphicsItem*>(parent()), SLOT(toggleProcessorWidget()));
        if (widget->isVisible())
            action->setChecked(true);
        customContextMenu_->addAction(action);
    }
}
//---------------------------------------------------------------------------------------------------------------
//                  events
//---------------------------------------------------------------------------------------------------------------
void WidgetToggleButtonGraphicsItem::mousePressEvent(QGraphicsSceneMouseEvent* event) {
    if (event->button() == Qt::LeftButton) {
        emit pressed();
        update();
    } else
        QGraphicsItem::mousePressEvent(event);
}

} // namespace
