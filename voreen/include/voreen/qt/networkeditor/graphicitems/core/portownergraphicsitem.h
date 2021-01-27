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

#ifndef VRN_PORTOWNERGRAPHICSITEM_H
#define VRN_PORTOWNERGRAPHICSITEM_H

#include "propertyownergraphicsitem.h"
#include "../utils/widgettogglebuttongraphicsitem.h"
#include "../utils/propertylistbuttongraphicsitem.h"
#include "../utils/renamabletextgraphicsitem.h"


namespace voreen {

    class ProgressBarGraphicsItem;
    class PortGraphicsItem;
    class PortOwnerLinkArrowGraphicsItem;

/**
 * Qt graphic item class, which collects the shared functions of processor and aggregation
 */
class PortOwnerGraphicsItem : public PropertyOwnerGraphicsItem {
Q_OBJECT
public:
    //constructor + destructor
    PortOwnerGraphicsItem(PropertyOwner* propertyOwner, NetworkEditor* networkEditor);
    virtual ~PortOwnerGraphicsItem();

    //---------------------------------------------------------------------------------
    //      nwebasegraphicsitem functions
    //---------------------------------------------------------------------------------
    void updateNWELayerAndCursor();

    void layoutChildItems();    ///layout ports
protected:
    void createChildItems();
    void deleteChildItems();

    virtual void setContextMenuActions();
    //---------------------------------------------------------------------------------
    //      NEEDS TO BE IMPLEMENTED
    //---------------------------------------------------------------------------------
public:
    virtual QList<Processor*> getProcessors() const = 0;
    virtual void saveMeta() = 0;
    virtual void loadMeta() = 0;
protected:
    virtual void registerPorts() = 0;

    //---------------------------------------------------------------------------------
    //      portownergraphicsitem functions
    //---------------------------------------------------------------------------------
public:
    //getter and setter
    QList<Port*> getPorts() const;
    QList<Port*> getInports() const;
    QList<Port*> getOutports() const;
    QList<CoProcessorPort*> getCoProcessorInports() const;
    QList<CoProcessorPort*> getCoProcessorOutports() const;
    QList<PortGraphicsItem*> getPortGraphicsItems() const;
    QList<PortOwnerLinkArrowGraphicsItem*> getPortOwnerLinkArrows() const;
    PortGraphicsItem* getPortGraphicsItem(const Port* port) const;
    void setVisible(bool visible);

    ProgressBarGraphicsItem* getProgressBar() const;
    PropertyListButtonGraphicsItem& getPropertyListButton();
    WidgetToggleButtonGraphicsItem& getWidgetToggleButton();
        //label
    const QGraphicsTextItem* getOrCreateNameLabel() const;
    QGraphicsTextItem* getOrCreateNameLabel();
public slots:
    void enterRenameMode();
protected slots:
    virtual void renameFinished() = 0;
    void nameChanged();

public:
        //link arrows
    void addGraphicalLinkArrow(PortOwnerGraphicsItem* po);
    void removeGraphicalLinkArrow(PortOwnerGraphicsItem* po);
    bool isPortOwnerLinkNeeded(PortOwnerGraphicsItem* po, NetworkEditorLayer nwel = NetworkEditorLayerGeneralLinking);

    //events
    void mousePressEvent(QGraphicsSceneMouseEvent* event);
    void mouseMoveEvent(QGraphicsSceneMouseEvent* event);
    void mouseReleaseEvent(QGraphicsSceneMouseEvent* event);
protected:
    //move events
    QVariant itemChange(GraphicsItemChange change, const QVariant& value);

public slots:
    virtual void toggleProcessorWidget() {};
    virtual void togglePropertyList() {};

signals:
    void openPropertyLinkDialog(PortOwnerGraphicsItem* src, PortOwnerGraphicsItem* dest);

public:
    //port help functions
    int getNumberOfVisiblePorts(QList<Port*> ports) const;
    int getNumberOfVisiblePorts(QList<CoProcessorPort*> ports) const;
    qreal getMinimumWidthForPorts() const;
    qreal getMinimumHeightForPorts() const;
protected:
    //members
    QList<PortGraphicsItem*> portItems_;
    QList<Port*> inports_;
    QList<Port*> outports_;
    QList<CoProcessorPort*> coInports_;
    QList<CoProcessorPort*> coOutports_;
    QList<PortOwnerLinkArrowGraphicsItem*> linkArrows_;

    RenamableTextGraphicsItem nameLabel_;
    ProgressBarGraphicsItem* progressBar_;
    WidgetToggleButtonGraphicsItem widgetToggleButton_;
    PropertyListButtonGraphicsItem propertyListButton_;

    QPointF clickPosition_;
    bool clickedPortOwner_;
    PortOwnerLinkArrowGraphicsItem* currentLinkArrow_;
};


} // namespace

#endif // VRN_PORTOWNERGraphicsItem_H
