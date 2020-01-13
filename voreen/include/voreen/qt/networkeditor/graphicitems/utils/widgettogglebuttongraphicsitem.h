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

#ifndef VRN_WIDGETTOGGLEBUTTONGRAPHICSITEM_H
#define VRN_WIDGETTOGGLEBUTTONGRAPHICSITEM_H

#include "../nwebasegraphicsitem.h"
#include "voreen/qt/networkeditor/editor_settings.h"

namespace voreen {

class ProcessorWidget;

class WidgetToggleButtonGraphicsItem : public NWEBaseGraphicsItem {
Q_OBJECT
public:
    //constructor and destructor
    WidgetToggleButtonGraphicsItem(NWEBaseGraphicsItem* parent);
    virtual ~WidgetToggleButtonGraphicsItem();

    //---------------------------------------------------------------------------------
    //      nwebasegraphicsitem functions
    //---------------------------------------------------------------------------------
    void updateNWELayerAndCursor() {};                          ///update selectable etc...
    int type() const {return UserTypesWidgetToggleButtonGraphicsItem;}   ///needed fpr select events

    //style
    QRectF boundingRect() const;
    QPainterPath shape() const;
    void mainPaint(QPainter* painter, const QStyleOptionGraphicsItem* option, QWidget* widget, NWEItemSettings& setting);
protected:
    void initializePaintSettings();

public:
    //handling child elements
    void layoutChildItems() {};
protected:
    void createChildItems() {};
    void deleteChildItems() {};

    void setContextMenuActions();

    //---------------------------------------------------------------------------------
    //      widgettogglebutton functions
    //---------------------------------------------------------------------------------
public:
    //getter and setter
    void setProcessorWidget(ProcessorWidget* widget);
    void setProcessorWidgets(const QList<ProcessorWidget*>& widgets);

    QList<ProcessorWidget*> getWidgets();

signals:
    /// This signal will be emitted if this item was been clicked upon
    void pressed();

protected:
    void mousePressEvent(QGraphicsSceneMouseEvent* event);

private:
    //members
    QList<ProcessorWidget*> widgets_;
};

} // namespace

#endif // VRN_WIDGETTOGGLEBUTTONGRAPHICSITEM_H
