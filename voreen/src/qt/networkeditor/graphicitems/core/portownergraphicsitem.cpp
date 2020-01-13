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

#include "voreen/qt/networkeditor/graphicitems/core/portownergraphicsitem.h"

//core
#include "voreen/core/properties/cameraproperty.h"
#include "voreen/core/ports/port.h"
#include "voreen/core/ports/coprocessorport.h"
#include "voreen/core/ports/renderport.h"

//gi
#include "voreen/qt/networkeditor/graphicitems/core/portgraphicsitem.h"
#include "voreen/qt/networkeditor/graphicitems/core/propertygraphicsitem.h"
#include "voreen/qt/networkeditor/graphicitems/connections/portownerlinkarrowgraphicsitem.h"
#include "voreen/qt/networkeditor/graphicitems/connections/portarrowgraphicsitem.h"
#include "voreen/qt/networkeditor/graphicitems/connections/propertylinkarrowgraphicsitem.h"
#include "voreen/qt/networkeditor/graphicitems/utils/progressbargraphicsitem.h"

#include <QGraphicsSceneEvent>
#include <QTextCursor>
#include <QIcon>
#include <QAction>
#include <QApplication>

namespace voreen {

PortOwnerGraphicsItem::PortOwnerGraphicsItem(PropertyOwner* propertyOwner, NetworkEditor* networkEditor)
    : PropertyOwnerGraphicsItem(propertyOwner, networkEditor)
    , nameLabel_("",this)
    , progressBar_(0)
    , widgetToggleButton_(this)
    , propertyListButton_(this)
    , currentLinkArrow_(0)
{
    tgtAssert(networkEditor, "passed null pointer");
    setFlag(QGraphicsItem::ItemSendsGeometryChanges);

    setContextMenuActions();

    setParent(networkEditor);
    networkEditor->scene()->addItem(this);

    setZValue(ZValuesPortOwnerGraphicsItem);

    QObject::connect(&nameLabel_, SIGNAL(renameFinished()), this, SLOT(renameFinished()));
    QObject::connect(&nameLabel_, SIGNAL(textChanged()), this, SLOT(nameChanged()));
    QObject::connect(&propertyListButton_, SIGNAL(pressed()), this, SLOT(togglePropertyList()));
    QObject::connect(&widgetToggleButton_, SIGNAL(pressed()), this, SLOT(toggleProcessorWidget()));
    widgetToggleButton_.hide();
    propertyListButton_.hide();
}

PortOwnerGraphicsItem::~PortOwnerGraphicsItem() {
  deleteChildItems();
}

//---------------------------------------------------------------------------------------------------------------
//                  getter and setter
//---------------------------------------------------------------------------------------------------------------
PropertyListButtonGraphicsItem& PortOwnerGraphicsItem::getPropertyListButton() {
    return propertyListButton_;
}

WidgetToggleButtonGraphicsItem& PortOwnerGraphicsItem::getWidgetToggleButton() {
    return widgetToggleButton_;
}

const QGraphicsTextItem* PortOwnerGraphicsItem::getOrCreateNameLabel() const {
    return &nameLabel_;
}

ProgressBarGraphicsItem* PortOwnerGraphicsItem::getProgressBar() const {
    return progressBar_;
}

QGraphicsTextItem* PortOwnerGraphicsItem::getOrCreateNameLabel() {
    return &nameLabel_;
}

void PortOwnerGraphicsItem::enterRenameMode() {
    nameLabel_.setTextInteractionFlags(Qt::TextEditorInteraction);
    nameLabel_.setFlag(QGraphicsItem::ItemIsFocusable, true);
    nameLabel_.setFocus();
    QTextCursor cursor = nameLabel_.textCursor();
    cursor.select(QTextCursor::Document);
    nameLabel_.setTextCursor(cursor);
}

void PortOwnerGraphicsItem::nameChanged() {
    if (!networkEditor_ || !networkEditor_->getProcessorNetwork())
        return;

    prepareGeometryChange();
    layoutChildItems();
    if (scene())
        scene()->invalidate();
}

 QList<Port*> PortOwnerGraphicsItem::getPorts() const {
    QList<Port*> result = getInports();
    result += getOutports();

    foreach (CoProcessorPort* p, getCoProcessorInports() + getCoProcessorOutports())
        result.append(p);

    return result;
}

QList<PortGraphicsItem*> PortOwnerGraphicsItem::getPortGraphicsItems() const {
    return portItems_;
}

PortGraphicsItem* PortOwnerGraphicsItem::getPortGraphicsItem(const Port* port) const {
    tgtAssert(port, "passed null pointer");

    foreach (PortGraphicsItem* portItem, portItems_) {
        if (portItem->getPort() == port)
            return portItem;
    }

    return 0;
}

QList<Port*> PortOwnerGraphicsItem::getInports() const {
    return inports_;
}

QList<Port*> PortOwnerGraphicsItem::getOutports() const {
    return outports_;
};

QList<CoProcessorPort*> PortOwnerGraphicsItem::getCoProcessorInports() const {
    return coInports_;
};

QList<CoProcessorPort*> PortOwnerGraphicsItem::getCoProcessorOutports() const {
    return coOutports_;
};

QList<PortOwnerLinkArrowGraphicsItem*> PortOwnerGraphicsItem::getPortOwnerLinkArrows() const {
    return linkArrows_;
}

//link arrows
void PortOwnerGraphicsItem::addGraphicalLinkArrow(PortOwnerGraphicsItem* po) {
    foreach(PortOwnerLinkArrowGraphicsItem* arrow, linkArrows_){
        if(arrow->getDestinationItem() == &(po->getPropertyListButton())){
            return;
        }
    }
    linkArrows_.append(new PortOwnerLinkArrowGraphicsItem(&getPropertyListButton(),&(po->getPropertyListButton())));
    po->addGraphicalLinkArrow(this);
}

void PortOwnerGraphicsItem::removeGraphicalLinkArrow(PortOwnerGraphicsItem* po) {
    foreach(PortOwnerLinkArrowGraphicsItem* arrow, linkArrows_){
        if(arrow->getDestinationItem() == &(po->getPropertyListButton())){
            PortOwnerLinkArrowGraphicsItem* save = arrow;
            linkArrows_.removeOne(arrow);
            po->removeGraphicalLinkArrow(this);
            delete save;
            return;
        }
    }
}

bool PortOwnerGraphicsItem::isPortOwnerLinkNeeded(PortOwnerGraphicsItem* po, NetworkEditorLayer nwel) {
    foreach(PropertyGraphicsItem* prop, getPropertyList()->getAllPropertyItems()) {
        switch(nwel){
        case NetworkEditorLayerDataFlow:
        case NetworkEditorLayerLinking:
        case NetworkEditorLayerGeneralLinking:
        case NetworkEditorLayerPortSizeLinking:
            foreach(PropertyLinkArrowGraphicsItem* arrow, prop->getSourceLinkList()) {
                if(dynamic_cast<PortOwnerGraphicsItem*>(arrow->getDestinationItem()->getPropertyOwnerItem()) == po)
                    return true;
                else if (dynamic_cast<PortGraphicsItem*>(arrow->getDestinationItem()->getPropertyOwnerItem()) &&
                         dynamic_cast<PortGraphicsItem*>(arrow->getDestinationItem()->getPropertyOwnerItem())->getPortOwner() == po)
                    return true;
            }
            foreach(PropertyLinkArrowGraphicsItem* arrow, prop->getDestinationLinkList()) {
                if(dynamic_cast<PortOwnerGraphicsItem*>(arrow->getSourceItem()->getPropertyOwnerItem()) == po)
                    return true;
                else if (dynamic_cast<PortGraphicsItem*>(arrow->getSourceItem()->getPropertyOwnerItem()) &&
                         dynamic_cast<PortGraphicsItem*>(arrow->getSourceItem()->getPropertyOwnerItem())->getPortOwner() == po)
                    return true;
            }
            break;
        case NetworkEditorLayerCameraLinking:
            if (prop->getProperty()->getClassName() != "CameraProperty") break;
            foreach(PropertyLinkArrowGraphicsItem* arrow, prop->getSourceLinkList()) {
                if(arrow->getDestinationItem()->getProperty()->getClassName() == "CameraProperty") {
                    if(dynamic_cast<PortOwnerGraphicsItem*>(arrow->getDestinationItem()->getPropertyOwnerItem()) == po)
                        return true;
                    else if (dynamic_cast<PortGraphicsItem*>(arrow->getDestinationItem()->getPropertyOwnerItem()) &&
                             dynamic_cast<PortGraphicsItem*>(arrow->getDestinationItem()->getPropertyOwnerItem())->getPortOwner() == po)
                        return true;
                }
            }
            foreach(PropertyLinkArrowGraphicsItem* arrow, prop->getDestinationLinkList()) {
                if(arrow->getSourceItem()->getProperty()->getClassName() == "CameraProperty") {
                    if(dynamic_cast<PortOwnerGraphicsItem*>(arrow->getSourceItem()->getPropertyOwnerItem()) == po)
                        return true;
                    else if (dynamic_cast<PortGraphicsItem*>(arrow->getSourceItem()->getPropertyOwnerItem()) &&
                             dynamic_cast<PortGraphicsItem*>(arrow->getSourceItem()->getPropertyOwnerItem())->getPortOwner() == po)
                        return true;
                }
            }
            break;
        default:
            tgtAssert(false,"should not get here");
            break;
        }
    }
    return false;
}

void PortOwnerGraphicsItem::setVisible(bool visible) {
    NWEBaseGraphicsItem::setVisible(visible);
    foreach(PortGraphicsItem* port, this->getPortGraphicsItems())
        port->setVisible(visible);
    if(!visible)
        foreach(PortOwnerLinkArrowGraphicsItem* arrow, linkArrows_)
            arrow->setVisible(false);

}

//---------------------------------------------------------------------------------------------------------------
//                  newbasegraphicsitem functions
//---------------------------------------------------------------------------------------------------------------
void PortOwnerGraphicsItem::setContextMenuActions() {
    NWEBaseGraphicsItem::setContextMenuActions();
    QAction* action = new QAction(QIcon(":/qt/icons/rename.png"), tr("Rename"), this);
    action->setShortcut(Qt::Key_F2);
    connect(action,SIGNAL(triggered()),this,SLOT(enterRenameMode()));
    customContextMenu_->addAction(action);
}

void PortOwnerGraphicsItem::createChildItems() {
    QList<Port*> ports = getPorts();
    foreach (Port* p, ports) {
        PortGraphicsItem* portItem = new PortGraphicsItem(p, this);
        propertyList_.addSubPropertyListItem(portItem->getPropertyList());
        portItems_.append(portItem);
    }
    propertyList_.hide();
    layoutChildItems();
}

void PortOwnerGraphicsItem::deleteChildItems() {
    while(!portItems_.isEmpty()) {
        PortGraphicsItem* p = portItems_.back();
        getPropertyList()->removeSubPropertyListItem(p->getPropertyList());
        portItems_.pop_back();
        delete p;
    }
}

void PortOwnerGraphicsItem::layoutChildItems() {
    if (!networkEditor_ || !networkEditor_->updatesEnabled())
        return;

    // regular inports
    int numberOfPorts = getNumberOfVisiblePorts(inports_);
    for (int i = 0, j = 0; i < inports_.size(); ++i) {
        PortGraphicsItem* inportItem = getPortGraphicsItem(inports_[i]);
        if (inportItem->isVisible()) {
            int x = (j + 1) * boundingRect().width() / (numberOfPorts + 1) - inportItem->boundingRect().width() / 2.0;
            int y = 0;
            inportItem->setPos(x,y);
            inportItem->layoutChildItems();
            ++j;
        }
    }

    // coprocessor inports
    numberOfPorts = getNumberOfVisiblePorts(getCoProcessorInports());
    for (int i = 0, j = 0; i < getCoProcessorInports().size(); ++i) {
        CoProcessorPort* inport = getCoProcessorInports()[i];
        PortGraphicsItem* inportItem = getPortGraphicsItem(inport);
        if (inportItem->isVisible()) {
            ++j;
            int x = 0;
            int y = j * boundingRect().height() / (numberOfPorts + 1) - inportItem->boundingRect().height() / 2.0;
            inportItem->setPos(x, y);
            inportItem->layoutChildItems();
        }
    }

    // regular outports
    numberOfPorts = getNumberOfVisiblePorts(getOutports());
    for (int i = 0, j = 0; i < getOutports().size(); ++i) {
        Port* outport = getOutports()[i];
        PortGraphicsItem* outportItem = getPortGraphicsItem(outport);
        if (outportItem->isVisible()) {
            ++j;
            int x = j * boundingRect().width() / (numberOfPorts + 1) - outportItem->boundingRect().width() / 2.0;
            int y = boundingRect().height() - (outportItem->boundingRect().height()); // /2.f
            outportItem->setPos(x, y);
            outportItem->layoutChildItems();
        }
    }

    // coprocessor outports
    numberOfPorts = getNumberOfVisiblePorts(getCoProcessorOutports());
    for (int i = 0, j = 0; i < getCoProcessorOutports().size(); ++i) {
        Port* outport = getCoProcessorOutports()[i];
        PortGraphicsItem* outportItem = getPortGraphicsItem(outport);
        if (outportItem->isVisible()) {
            ++j;
            int x = boundingRect().width() - (outportItem->boundingRect().width());
            int y = j * boundingRect().height() / (numberOfPorts + 1) - outportItem->boundingRect().height() / 2.0;
            outportItem->setPos(x, y);
            outportItem->layoutChildItems();
        }
    }

    // the text label should be in the middle of the processoritem
    qreal x = (boundingRect().width() / 2.0) - (nameLabel_.boundingRect().width() / 2.0);
    qreal y = (boundingRect().height() / 2.0) - (nameLabel_.boundingRect().height() / 2.0);
    nameLabel_.setPos(x, y);

    // the propertylist item should be centered on the bottom of the proecssoritem
    x = (boundingRect().width() / 2.0) - (propertyList_.boundingRect().width() / 2.0);
    y = boundingRect().height();
    propertyList_.setPos(mapToScene(QPointF(x, y)));
    propertyList_.layoutChildItems();//update();

    propertyListButton_.setPos(buttonsOffsetX, buttonsOffsetY);
    widgetToggleButton_.setPos(boundingRect().width() - widgetToggleButton_.boundingRect().width() - buttonsOffsetX, buttonsOffsetY);
}

 void PortOwnerGraphicsItem::updateNWELayerAndCursor() {
    switch(currentCursorMode()) {
    case NetworkEditorCursorSelectMode:
        switch (currentLayer()) {
        case NetworkEditorLayerDataFlow:
            setFlag(ItemIsSelectable);
            setAcceptHoverEvents(true);
            nameLabel_.setAcceptHoverEvents(true);
            propertyList_.setVisible(false);
            propertyList_.setIsVisibleInEditor(false);
            propertyListButton_.setVisible(false);
            setFlag(ItemIsMovable, true);
            break;
        case NetworkEditorLayerGeneralLinking:
            setFlag(ItemIsSelectable);
            setAcceptHoverEvents(true);
            nameLabel_.setAcceptHoverEvents(true);
            propertyListButton_.setVisible(true);
            setFlag(ItemIsMovable, false);
            propertyList_.setVisible(false);
            propertyList_.setIsVisibleInEditor(false);
            break;
        case NetworkEditorLayerCameraLinking:
            setFlag(ItemIsSelectable,false);
            setAcceptHoverEvents(false);
            nameLabel_.setAcceptHoverEvents(false);
            setFlag(ItemIsMovable, false);
            propertyList_.setVisible(false);
            propertyList_.setIsVisibleInEditor(false);
            propertyListButton_.setVisible(false);
            foreach (Processor* processor, getProcessors()){
                if(!processor->getPropertiesByType<CameraProperty>().empty()){
                    propertyListButton_.setVisible(true);
                    setFlag(ItemIsSelectable);
                    setAcceptHoverEvents(true);
                    nameLabel_.setAcceptHoverEvents(true);
                    break;
                }
            }
            break;
        case NetworkEditorLayerPortSizeLinking:
            setFlag(ItemIsSelectable, false);
            setAcceptHoverEvents(false);
            nameLabel_.setAcceptHoverEvents(false);
            setFlag(ItemIsMovable, false);
            propertyList_.setVisible(false);
            propertyList_.setIsVisibleInEditor(false);
            propertyListButton_.setVisible(false);
            foreach (Processor* processor, getProcessors()){
                foreach(Port* port, processor->getPorts()) {
                    if(!port->getPropertiesByType<RenderSizeOriginProperty>().empty() ||
                       !port->getPropertiesByType<RenderSizeReceiveProperty>().empty()){
                        setFlag(ItemIsSelectable);
                        setAcceptHoverEvents(true);
                        nameLabel_.setAcceptHoverEvents(true);
                        break;
                    }
                }
            }
            break;
        default:
            tgtAssert(false, "shouldn't get here");
        }
        break;
    case NetworkEditorCursorMoveMode:
        switch (currentLayer()) {
        case NetworkEditorLayerDataFlow:
            setFlag(ItemIsSelectable, false);
            setAcceptHoverEvents(false);
            nameLabel_.setAcceptHoverEvents(false);
            propertyList_.setVisible(false);
            propertyList_.setIsVisibleInEditor(false);
            propertyListButton_.setVisible(false);
            setFlag(ItemIsMovable, false);
            break;
        case NetworkEditorLayerGeneralLinking:
            setFlag(ItemIsSelectable, false);
            setAcceptHoverEvents(false);
            nameLabel_.setAcceptHoverEvents(false);
            propertyList_.setVisible(false);
            propertyList_.setIsVisibleInEditor(false);
            propertyListButton_.setVisible(true);
            setFlag(ItemIsMovable, false);
            break;
        case NetworkEditorLayerCameraLinking:
            setFlag(ItemIsSelectable,false);
            setAcceptHoverEvents(false);
            nameLabel_.setAcceptHoverEvents(false);
            setFlag(ItemIsMovable, false);
            propertyList_.setVisible(false);
            propertyList_.setIsVisibleInEditor(false);
            propertyListButton_.setVisible(false);
            foreach (Processor* processor, getProcessors()){
                if(!processor->getPropertiesByType<CameraProperty>().empty()){
                    propertyListButton_.setVisible(true);
                    break;
                }
            }
            break;
        case NetworkEditorLayerPortSizeLinking:
            setFlag(ItemIsSelectable, false);
            setAcceptHoverEvents(false);
            nameLabel_.setAcceptHoverEvents(false);
            propertyList_.setVisible(false);
            propertyList_.setIsVisibleInEditor(false);
            setFlag(ItemIsMovable, false);
            propertyListButton_.setVisible(false);
            break;
        default:
            tgtAssert(false, "shouldn't get here");
        }
        break;
    default:
        tgtAssert(false,"should not get here");
        break;
    }

    foreach (PortGraphicsItem* portItem, portItems_) {
        portItem->updateNWELayerAndCursor();
    }
    foreach (PortOwnerLinkArrowGraphicsItem* arrow, linkArrows_) {
        arrow->updateNWELayerAndCursor();
    }

 }

//---------------------------------------------------------------------------------------------------------------
//                  help functions
//---------------------------------------------------------------------------------------------------------------
int PortOwnerGraphicsItem::getNumberOfVisiblePorts(QList<Port*> ports) const {
    int result = 0;
    foreach (Port* p, ports) {
        PortGraphicsItem* portItem = getPortGraphicsItem(p);
        if (portItem->isVisible())
            result++;
    }
    return result;
}

int PortOwnerGraphicsItem::getNumberOfVisiblePorts(QList<CoProcessorPort*> ports) const {
    int result = 0;
    foreach (CoProcessorPort* p, ports) {
        PortGraphicsItem* portItem = getPortGraphicsItem(p);
        if (portItem->isVisible())
            result++;
    }
    return result;
}

qreal PortOwnerGraphicsItem::getMinimumWidthForPorts() const {
    int numberOfPorts = qMax(inports_.size(), outports_.size());
    if (numberOfPorts == 0)
        return 0.0;
    else {
        PortGraphicsItem* portItem = 0;
        if (inports_.size() == numberOfPorts)
            portItem = getPortGraphicsItem(inports_[0]);   //multiconnections?
        else
            portItem = getPortGraphicsItem(outports_[0]);
        return numberOfPorts * (portItem->boundingRect().width() + 5.0);
    }
}

qreal PortOwnerGraphicsItem::getMinimumHeightForPorts() const {
    int numberOfPorts = qMax(coInports_.size(), coOutports_.size());
    if (numberOfPorts == 0)
        return 0.0;
    else {
        PortGraphicsItem* portItem = 0;
        if (coInports_.size() == numberOfPorts)                    //multiconnections?
            portItem = getPortGraphicsItem(coInports_[0]);
        else
            portItem = getPortGraphicsItem(coOutports_[0]);

        return numberOfPorts * (portItem->boundingRect().height() + 5.0);
    }
}

//---------------------------------------------------------------------------------------------------------------
//                  events
//---------------------------------------------------------------------------------------------------------------
void PortOwnerGraphicsItem::mousePressEvent(QGraphicsSceneMouseEvent* event) {
    if(event->buttons() == Qt::LeftButton){
        clickedPortOwner_ = true;
        clickPosition_ = event->pos();
    } else {
        delete currentLinkArrow_;
        currentLinkArrow_ = 0;
        clickedPortOwner_ = false;
    }
    QGraphicsItem::mousePressEvent(event);
}

void PortOwnerGraphicsItem::mouseMoveEvent(QGraphicsSceneMouseEvent* event) {
    if (currentLayer() == NetworkEditorLayerGeneralLinking ||
        currentLayer() == NetworkEditorLayerCameraLinking) {
        if (event->modifiers() == Qt::ControlModifier) {
            setFlag(ItemIsMovable, true);

            delete currentLinkArrow_;
            currentLinkArrow_ = 0;
            clickedPortOwner_ = false;
//            emit endedArrow();
        }
        else {
            setFlag(ItemIsMovable, false);
            const QLineF line(event->pos(), clickPosition_);
            if (clickedPortOwner_ && (currentLinkArrow_ || (line.length() > minimumDistanceToStartDrawingArrow))) {
                if (currentLinkArrow_ == 0) {
                    currentLinkArrow_ = new PortOwnerLinkArrowGraphicsItem(&getPropertyListButton());
                    currentLinkArrow_->setZValue(ZValuesSelectedPortOwnerLinkArrowGraphicsItem);
                    currentLinkArrow_->setIsHovered(true);
                    currentLinkArrow_->updateNWELayerAndCursor();
                    //emit startedArrow();
                }
                QPointF pos = mapToScene(event->pos());
                currentLinkArrow_->setDestinationPoint(pos);
            }
        }
    }
    else {
        if(!QApplication::overrideCursor() || QApplication::overrideCursor()->shape() != Qt::SizeAllCursor) {
            QApplication::restoreOverrideCursor(); //probably not needed
            QApplication::setOverrideCursor(Qt::SizeAllCursor);
        }
    }
    QGraphicsItem::mouseMoveEvent(event);
}

void PortOwnerGraphicsItem::mouseReleaseEvent(QGraphicsSceneMouseEvent* event) {

    if (currentLayer() == NetworkEditorLayerGeneralLinking ||
        currentLayer() == NetworkEditorLayerCameraLinking ) {
        if (currentLinkArrow_) {
            QGraphicsItem* item = 0;
            QList<QGraphicsItem*> items = scene()->items(event->scenePos());
            foreach (QGraphicsItem* i, items) {
                if (i->type() != UserTypesPortOwnerLinkArrowGraphicsItem) {
                    item = i;
                    break;
                }
            }

            if (item) {
                switch (item->type()) {
                case UserTypesPortGraphicsItem:
                case QGraphicsTextItem::Type:
                case UserTypesRenamableTextGraphicsItem:
                case UserTypesWidgetToggleButtonGraphicsItem:
                case UserTypesPropertyListButtonGraphicsItem:
                    item = item->parentItem();
                case UserTypesProcessorGraphicsItem:
                    {
                    PortOwnerGraphicsItem* dest = dynamic_cast<PortOwnerGraphicsItem*>(item);
                    if(dest){
                        if(currentLayer() == NetworkEditorLayerGeneralLinking)
                            emit openPropertyLinkDialog(this, dest);
                        else if(!dest->getProcessors()[0]->getPropertiesByType<CameraProperty>().empty())
                            emit openPropertyLinkDialog(this, dest);
                    }
                    }
                    break;
                }
            }

            delete currentLinkArrow_;
            currentLinkArrow_ = 0;
            //emit endedArrow();
        }
    } else {
        QApplication::restoreOverrideCursor();
    }

    clickedPortOwner_ = false;

    QGraphicsItem::mouseReleaseEvent(event);
}

QVariant PortOwnerGraphicsItem::itemChange(GraphicsItemChange change, const QVariant &value) {
    if (change == ItemPositionChange && scene()) {
        //update port arrows
        foreach (PortGraphicsItem* pitem, portItems_){
            if(pitem->isOutport())
                foreach(PortArrowGraphicsItem* arrow, pitem->getArrowList())
                    arrow->prepareGeometryChange();
            else
                foreach (PortGraphicsItem* oitem, pitem->getConnectedOutportGraphicsItems()){
                    foreach(PortArrowGraphicsItem* arrow, oitem->getArrowList())
                        arrow->prepareGeometryChange();
                }
        }
        //update propertylist
        if(getPropertyList()->isVisible()){
            getPropertyList()->setPosition(mapToScene((boundingRect().width() / 2.0) - (propertyList_.boundingRect().width() / 2.0), boundingRect().height()));
        }

        foreach(PropertyGraphicsItem* item, getPropertyList()->getAllPropertyItems()){
            foreach(PropertyLinkArrowGraphicsItem* arrow, item->getDestinationLinkList()) {
                foreach(PropertyLinkArrowGraphicsItem* subArrow, arrow->getSourceItem()->getDestinationLinkList())
                    subArrow->prepareGeometryChange();
                foreach(PropertyLinkArrowGraphicsItem* subArrow, arrow->getSourceItem()->getSourceLinkList())
                    subArrow->prepareGeometryChange();
            }
            foreach(PropertyLinkArrowGraphicsItem* arrow, item->getSourceLinkList()) {
                foreach(PropertyLinkArrowGraphicsItem* subArrow, arrow->getDestinationItem()->getDestinationLinkList())
                    subArrow->prepareGeometryChange();
                foreach(PropertyLinkArrowGraphicsItem* subArrow, arrow->getDestinationItem()->getSourceLinkList())
                    subArrow->prepareGeometryChange();
            }
        }

        //update port owner link arrow
        foreach(PortOwnerLinkArrowGraphicsItem* arrow, linkArrows_){
            arrow->prepareGeometryChange();
        }
    }
    return QGraphicsItem::itemChange(change, value);
 }

} // namespace
