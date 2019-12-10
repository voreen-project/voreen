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

#include "voreen/qt/networkeditor/graphicitems/core/portgraphicsitem.h"
#include "voreen/qt/networkeditor/styles/nwestyle_base.h"

//core
#include "voreen/core/ports/port.h"
#include "voreen/core/ports/coprocessorport.h"
#include "voreen/core/ports/renderport.h"
#include "voreen/core/ports/geometryport.h"
#include "voreen/core/ports/volumeport.h"
#include "tgt/assert.h"
#include "voreen/core/voreenapplication.h"
#include "voreen/core/properties/link/linkevaluatorhelper.h"

//gi
#include "voreen/qt/networkeditor/graphicitems/core/portownergraphicsitem.h"
#include "voreen/qt/networkeditor/graphicitems/core/propertygraphicsitem.h"
#include "voreen/qt/networkeditor/graphicitems/connections/portarrowgraphicsitem.h"
#include "voreen/qt/networkeditor/graphicitems/connections/portsizelinkarrowgraphicsitem.h"
#include "voreen/qt/networkeditor/graphicitems/tooltips/tooltipportgraphicsitem.h"

#include <QMessageBox>
#include <QList>
#include <QPainter>
#include <QStyleOptionGraphicsItem>
#include <QGraphicsSceneMouseEvent>
#include <QFileDialog>
#include <QUrl>
#include <QDesktopServices>
#include <QAction>

namespace voreen {

PortGraphicsItem::PortGraphicsItem(Port* port, PortOwnerGraphicsItem* parent)
    : PropertyOwnerGraphicsItem(port, parent->getNetworkEditor())
    , PortWidget(port)
    , currentPortArrow_(0)
    , currentPortSizeLinkArrow_(0)
    , portOwner_(parent)
    , port_(port)
{
    tgtAssert(port, "passed null pointer");
    tgtAssert(parent, "passed null pointer");

    setParent(parent);
    setParentItem(parent);
    setGuiName(QString::fromStdString(port->getGuiName()));
    getPropertyList()->updateHeaderLabelItem();

    setContextMenuActions();

    setFlag(ItemIsSelectable);
    setFlag(ItemIsPanel);
    setAcceptHoverEvents(true);
    setZValue(ZValuesPortGraphicsItem);
    setTransformOriginPoint(boundingRect().center());

    setToolTipGraphicsItem(new ToolTipPortGraphicsItem(this));

    connect(this,SIGNAL(startedArrow()),getNetworkEditor(),SLOT(toggleToolTip()));
    connect(this,SIGNAL(endedArrow()),getNetworkEditor(),SLOT(toggleToolTip()));

    createChildItems();
}

PortGraphicsItem::~PortGraphicsItem() {
    deleteChildItems();
    foreach (PortGraphicsItem* portItem, connectedOutportGraphicsItems_){
        portItem->removeGraphicalConnection(this);
    }
}

//---------------------------------------------------------------------------------------------------------------
//                  getter and setter
//---------------------------------------------------------------------------------------------------------------
Port* PortGraphicsItem::getPort() const {
    return port_;
}

PortOwnerGraphicsItem* PortGraphicsItem::getPortOwner() const {
    return portOwner_;
}

QList<PortGraphicsItem*> PortGraphicsItem::getConnectedOutportGraphicsItems() const{
    return connectedOutportGraphicsItems_;
}

const QList<PortArrowGraphicsItem*>& PortGraphicsItem::getArrowList() const {
    return arrowList_;
}

QColor PortGraphicsItem::getColor() const {
    tgt::col3 portCol = port_->getColorHint();
    QColor result = QColor(portCol.x, portCol.y, portCol.z);
    result.setAlpha(100);
    if(currentLayer() == NetworkEditorLayerPortSizeLinking){
        if(!getPort()->getPropertiesByType<RenderSizeOriginProperty>().empty())
            result = QColor(255,0,0,255);
        if(!getPort()->getPropertiesByType<RenderSizeReceiveProperty>().empty())
            result = QColor(0,255,0,255);
    }
    return result;
}

bool PortGraphicsItem::isOutport() const {
    return port_->isOutport();
}

void PortGraphicsItem::setVisible(bool visible) {
    NWEBaseGraphicsItem::setVisible(visible);

    foreach(PortArrowGraphicsItem* arrow, arrowList_)
        arrow->setVisible(visible);

    foreach(PortSizeLinkArrowGraphicsItem* arrow, linkArrowList_)
        arrow->setVisible(visible);

}

//---------------------------------------------------------------------------------------------------------------
//                  nwebasegraphicsitem functions
//---------------------------------------------------------------------------------------------------------------
QRectF PortGraphicsItem::boundingRect() const {
    return currentStyle()->PortGI_boundingRect(this);
}

QPainterPath PortGraphicsItem::shape() const {
    return currentStyle()->PortGI_shape(this);
}

void PortGraphicsItem::initializePaintSettings() {
    currentStyle()->PortGI_initializePaintSettings(this);
}

void PortGraphicsItem::mainPaint(QPainter* painter, const QStyleOptionGraphicsItem* option, QWidget* widget, NWEItemSettings& setting) {
    currentStyle()->PortGI_paint(this, painter, option, widget, setting);
}

void PortGraphicsItem::updateFromPort() {
    if (!isOutport()) {

        // Update the port widget itself.
        update();

        // Update connected arrow since it's appearance might have changed.
        // Due to voreen architecture, we need to find the arrows owned by 
        // possibly any connected port to this one.
        foreach(PortGraphicsItem* port, connectedOutportGraphicsItems_) {
            foreach(PortArrowGraphicsItem* arrow, port->getArrowList()) {
                if (arrow->getDestinationItem() == this)
                    arrow->update();
            }
        }
    }
}

void PortGraphicsItem::setContextMenuActions() {
    NWEBaseGraphicsItem::setContextMenuActions();
    if(dynamic_cast<RenderPort*>(port_) ||
       dynamic_cast<GeometryPort*>(port_) ||
       dynamic_cast<VolumePort*>(port_)) {
        QAction* action = new QAction(QString("Save current Port Data"),this);
        connect(action,SIGNAL(triggered()),this,SLOT(savePortData()));
        customContextMenu_->addAction(action);
    }
}

void PortGraphicsItem::updateNWELayerAndCursor() {
    switch(currentCursorMode()) {
    case NetworkEditorCursorSelectMode:
        switch (currentLayer()) {
        case NetworkEditorLayerDataFlow:
            setFlag(ItemIsSelectable);
            setAcceptHoverEvents(true);
            break;
        case NetworkEditorLayerGeneralLinking:
        case NetworkEditorLayerCameraLinking:
            setFlag(ItemIsSelectable, false);
            setAcceptHoverEvents(false);
            break;
        case NetworkEditorLayerPortSizeLinking:
            if(!port_->getPropertiesByType<RenderSizeOriginProperty>().empty() || !getPort()->getPropertiesByType<RenderSizeReceiveProperty>().empty()){
                setFlag(ItemIsSelectable);
                setAcceptHoverEvents(true);
            } else {
                setFlag(ItemIsSelectable, false);
                setAcceptHoverEvents(false);
            }
            break;
        default:
            tgtAssert(false, "shouldn't get here");
            break;
        }
        break;
    case NetworkEditorCursorMoveMode:
        switch (currentLayer()) {
        case NetworkEditorLayerDataFlow:
            setFlag(ItemIsSelectable,false);
            setAcceptHoverEvents(false);
            break;
        case NetworkEditorLayerGeneralLinking:
        case NetworkEditorLayerCameraLinking:
        case NetworkEditorLayerPortSizeLinking:
            setFlag(ItemIsSelectable, false);
            setAcceptHoverEvents(false);
            break;
        default:
            tgtAssert(false, "shouldn't get here");
            break;
        }
        break;
    default:
        tgtAssert(false,"should not get here");
        break;
    }

    foreach(PortArrowGraphicsItem* item, arrowList_)
        item->updateNWELayerAndCursor();
    foreach(PortSizeLinkArrowGraphicsItem* item, linkArrowList_)
        item->updateNWELayerAndCursor();
}

void PortGraphicsItem::layoutChildItems() {
    foreach (PortArrowGraphicsItem* pa, arrowList_){
        pa->update();
    }
}

void PortGraphicsItem::createChildItems() {
    layoutChildItems();
}

void PortGraphicsItem::deleteChildItems() {
    foreach (PortArrowGraphicsItem* pa, arrowList_){
        removeGraphicalConnection(pa->getDestinationItem());
    }
}

//---------------------------------------------------------------------------------------------------------------
//                  portgraphicsitem functions
//---------------------------------------------------------------------------------------------------------------
PortArrowGraphicsItem* PortGraphicsItem::setCurrentPortArrow(PortArrowGraphicsItem* arrow) {
    delete currentPortArrow_;
    emit startedArrow();
    currentPortArrow_ = arrow;
    currentPortArrow_->setSelected(false);
    PortGraphicsItem* destPortItem = arrow->getDestinationItem();
    currentPortArrow_->setDestinationItem(0);
    portOwner_->getNetworkEditor()->getProcessorNetwork()->disconnectPorts(this->getPort(),destPortItem->getPort());
    foreach (PortArrowGraphicsItem* arrowItem, arrowList_) {
        if(arrowItem == arrow){
            arrowList_.removeOne(arrowItem);
            destPortItem->removeConnectedOutportGraphicsItem(this);
            return currentPortArrow_;
        }
    }
    //should not get here
    return currentPortArrow_;
}

PortSizeLinkArrowGraphicsItem* PortGraphicsItem::setCurrentPortSizeLinkArrow(PortSizeLinkArrowGraphicsItem* arrow) {
    delete currentPortSizeLinkArrow_;
    emit startedArrow();
    currentPortSizeLinkArrow_ = arrow;
    currentPortSizeLinkArrow_->setSelected(false);
    PropertyGraphicsItem* destPropItem = arrow->getDestinationItem();
    currentPortSizeLinkArrow_->setDestinationItem(0);
    currentPortSizeLinkArrow_->setDestinationPort(0);

    getNetworkEditor()->getProcessorNetwork()->removePropertyLink(const_cast<Property*>(currentPortSizeLinkArrow_->getSourceItem()->getProperty()),
                                                                  const_cast<Property*>(destPropItem->getProperty()));

    foreach (PortSizeLinkArrowGraphicsItem* arrowItem, linkArrowList_) {
        if(arrowItem == arrow){
            linkArrowList_.removeOne(arrowItem);
            //destPortItem->removeConnectedOutportGraphicsItem(this);
            return currentPortSizeLinkArrow_;
        }
    }
    //should not get here
    return currentPortSizeLinkArrow_;
}

void PortGraphicsItem::addGraphicalConnection(PortGraphicsItem* port) {
    tgtAssert(getPort()->isOutport(),"PortGraphicsItem::addGraphicalConnection can just be called on an outport");
    PortArrowGraphicsItem* arrow = new PortArrowGraphicsItem(this);
    arrow->setDestinationItem(port);
    arrowList_.append(arrow);
    port->addConnectedOutportGraphicsItem(this);
}

void PortGraphicsItem::addConnectedOutportGraphicsItem(PortGraphicsItem* port) {
    tgtAssert(!getPort()->isOutport(),"PortGraphicsItem::addConnectedoutportGraphicsItem can just be called on an inport");
    connectedOutportGraphicsItems_.append(port);
}

void PortGraphicsItem::removeGraphicalConnection(PortGraphicsItem* port) {
    tgtAssert(getPort()->isOutport(),"PortGraphicsItem::removeGraphicalConnection can just be called on an outport");
    foreach (PortArrowGraphicsItem* arrowItem, arrowList_) {
        if(arrowItem->getDestinationItem() == port){
            arrowList_.removeOne(arrowItem);
            port->removeConnectedOutportGraphicsItem(this);
            delete arrowItem;
            return;
        }
    }
}

void PortGraphicsItem::removeConnectedOutportGraphicsItem(PortGraphicsItem* port) {
    tgtAssert(!getPort()->isOutport(),"PortGraphicsItem::addConnectedoutportGraphicsItem can just be called on an inport");
    connectedOutportGraphicsItems_.removeOne(port);
}

void PortGraphicsItem::addPortSizeLinkArrow(PropertyGraphicsItem* src, PropertyGraphicsItem* dst) {
    PortSizeLinkArrowGraphicsItem* arrow = new PortSizeLinkArrowGraphicsItem(src, this);
    arrow->setDestinationItem(dst);
    arrow->setDestinationPort(dynamic_cast<PortGraphicsItem*>(dst->getPropertyOwnerItem()));
    linkArrowList_.append(arrow);
}

void PortGraphicsItem::removePortSizeLinkArrow(PropertyGraphicsItem* src, PropertyGraphicsItem* dst) {
    foreach(PortSizeLinkArrowGraphicsItem* arrow, linkArrowList_)
        if(arrow->getSourceItem() == src && arrow->getDestinationItem() == dst){
            linkArrowList_.removeOne(arrow);
            delete arrow;
            break;
        }
}

void PortGraphicsItem::update(const QRectF& rect) {
    if (isOutport()) {
        foreach (PortArrowGraphicsItem* arrow, arrowList_)
            arrow->update();
    }
    else {
        //foreach (PortGraphicsItem* port, connectedPorts_)
        //    port->update(rect);
    }

    QGraphicsItem::update(rect);
}

//---------------------------------------------------------------------------------------------------------------
//                  slots
//---------------------------------------------------------------------------------------------------------------
void PortGraphicsItem::savePortData() {
    QFileDialog filedialog(getNetworkEditor());
    filedialog.setWindowTitle(tr("Save current Port Data"));
    filedialog.setDirectory(VoreenApplication::app()->getUserDataPath("screenshots").c_str());

    QStringList filter;

    struct tm* Tm;
    time_t currentTime = time(NULL);
    Tm = localtime(&currentTime);
    std::stringstream timestamp;

    if(dynamic_cast<RenderPort*>(port_)) {
        filedialog.setDefaultSuffix(tr("png"));
        filter << tr("PNG image (*.png)");
        filter << tr("JPEG image (*.jpg)");
        filter << tr("Windows Bitmap (*.bmp)");
        filter << tr("TIFF image (*.tif)");
        timestamp << "screenshot " << (Tm->tm_year+1900) << "-" << (Tm->tm_mon+1) << "-" << Tm->tm_mday << "-" << Tm->tm_hour << "-" << Tm->tm_min << "-" << Tm->tm_sec;
        timestamp << ".png";
    } else if(dynamic_cast<VolumePort*>(port_)) {
        filedialog.setDefaultSuffix(tr("vvd"));
        filter << tr("Voreen Volume Data (*.vvd)");
        timestamp << "volume " << (Tm->tm_year+1900) << "-" << (Tm->tm_mon+1) << "-" << Tm->tm_mday << "-" << Tm->tm_hour << "-" << Tm->tm_min << "-" << Tm->tm_sec;
        timestamp << ".vvd";
    } else if (dynamic_cast<GeometryPort*>(port_)){
        filedialog.setDefaultSuffix(tr("xml"));
        filter << tr("XML File (*.xml)");
        timestamp << "geometry " << (Tm->tm_year+1900) << "-" << (Tm->tm_mon+1) << "-" << Tm->tm_mday << "-" << Tm->tm_hour << "-" << Tm->tm_min << "-" << Tm->tm_sec;
        timestamp << ".xml";
    } else
        tgtAssert(false,"unknown port type");

    filedialog.setNameFilters(filter);
    filedialog.setAcceptMode(QFileDialog::AcceptSave);
    filedialog.selectFile(tr(timestamp.str().c_str()));
    filedialog.setOption(QFileDialog::DontUseNativeDialog);

    QList<QUrl> urls;
    urls << QUrl::fromLocalFile(VoreenApplication::app()->getUserDataPath("screenshots").c_str());
    for (auto f : QStandardPaths::standardLocations(QStandardPaths::DesktopLocation))
        urls << QUrl::fromLocalFile(f);
    for (auto f : QStandardPaths::standardLocations(QStandardPaths::HomeLocation))
        urls << QUrl::fromLocalFile(f);
    filedialog.setSidebarUrls(urls);

    QStringList fileList;
    if (filedialog.exec())
        fileList = filedialog.selectedFiles();
    if (fileList.empty())
        return;

    QString file = fileList.at(0);

    if(dynamic_cast<RenderPort*>(port_)) {
        dynamic_cast<RenderPort*>(port_)->saveToImage(file.toStdString());
    } else if(dynamic_cast<VolumePort*>(port_)) {
        dynamic_cast<VolumePort*>(port_)->saveData(file.toStdString());
    } else if (dynamic_cast<GeometryPort*>(port_)){
        dynamic_cast<GeometryPort*>(port_)->saveData(file.toStdString());
    } else
        tgtAssert(false,"unknown port type");
}


//---------------------------------------------------------------------------------------------------------------
//                  events
//---------------------------------------------------------------------------------------------------------------
void PortGraphicsItem::mousePressEvent(QGraphicsSceneMouseEvent* event) {
    if(event->buttons() == Qt::LeftButton){
        if(currentLayer() == NetworkEditorLayerDataFlow) {
            if (isOutport()) {
                currentPortArrow_ = new PortArrowGraphicsItem(this);
                currentPortArrow_->setDestinationPoint(event->scenePos());
                currentPortArrow_->setIsHovered(true);
                emit startedArrow();
            }
        } else if (currentLayer() == NetworkEditorLayerPortSizeLinking) {
            PropertyGraphicsItem* propItem = 0;

            foreach(PropertyGraphicsItem* item, getPropertyList()->getPropertyItems())
                if(item->getProperty()->getClassName() == "RenderSizeOriginProperty" || item->getProperty()->getClassName() == "RenderSizeReceiveProperty") {
                    propItem = item;
                    break;
                }

            if(propItem) {
                currentPortSizeLinkArrow_ = new PortSizeLinkArrowGraphicsItem(propItem, this);
                currentPortSizeLinkArrow_->setDestinationPoint(event->scenePos());
                currentPortSizeLinkArrow_->setIsHovered(true);
                emit startedArrow();
            }
        }
    } else {
        delete currentPortArrow_;
        currentPortArrow_ = 0;
        delete currentPortSizeLinkArrow_;
        currentPortSizeLinkArrow_ = 0;
    }
    NWEBaseGraphicsItem::mousePressEvent(event);
}

void PortGraphicsItem::mouseMoveEvent(QGraphicsSceneMouseEvent* event) {
    if(currentPortArrow_ || currentPortSizeLinkArrow_) {
        if(currentLayer() == NetworkEditorLayerDataFlow) {
            if (isOutport()) {
                tgtAssert(currentPortArrow_ != 0, "currentPortArrow was not set");
                currentPortArrow_->setDestinationPoint(event->scenePos());
                currentPortArrow_->setColorConnectableMode(ConnectionBaseGraphicsItem::CCM_DEFAULT);
                QList<QGraphicsItem*> items = scene()->items(event->scenePos());
                bool properItemFound = false;
                foreach(QGraphicsItem* item, items) {
                    if (item && (item != this)) {
                        switch(item->type()) {
                        case UserTypesPortGraphicsItem: {
                            Port* destPort = dynamic_cast<PortGraphicsItem*>(item)->getPort();
                            if (getPort()->testConnectivity(destPort) == Port::LEGAL_CONNECTIVITY)
                                currentPortArrow_->setColorConnectableMode(ConnectionBaseGraphicsItem::CCM_YES);
                            else if (getPort()->testConnectivity(destPort) == Port::UNDEFINED_LOOP_CONNECTIVITY)
                                currentPortArrow_->setColorConnectableMode(ConnectionBaseGraphicsItem::CCM_MAYBE);
                            else
                                currentPortArrow_->setColorConnectableMode(ConnectionBaseGraphicsItem::CCM_NO);
                            properItemFound = true;
                        } break;
                        case UserTypesPropertyListButtonGraphicsItem:
                        case UserTypesWidgetToggleButtonGraphicsItem:
                        case UserTypesRenamableTextGraphicsItem:
                            item = item->parentItem();
                            if (!item || (item->type() != UserTypesProcessorGraphicsItem))
                                break;
                        case UserTypesProcessorGraphicsItem: {
                        //case UserTypesAggregationGraphicsItem:
                            PortOwnerGraphicsItem* pOwner = dynamic_cast<PortOwnerGraphicsItem*>(item);
                            if (!pOwner)
                                return;
                            currentPortArrow_->setColorConnectableMode(ConnectionBaseGraphicsItem::CCM_NO);
                            bool isNext = true;
                            if(dynamic_cast<CoProcessorPort*>(getPort())) {
                                foreach(CoProcessorPort* destPort, pOwner->getCoProcessorInports()){
                                    Port::ConnectivityState cState = getPort()->testConnectivity(destPort);
                                    if (cState != Port::ILLEGAL_CONNECTIVITY) {
                                        if(isNext) {
                                            if(cState == Port::LEGAL_CONNECTIVITY)
                                                currentPortArrow_->setColorConnectableMode(ConnectionBaseGraphicsItem::CCM_YES);
                                            else
                                                currentPortArrow_->setColorConnectableMode(ConnectionBaseGraphicsItem::CCM_MAYBE);
                                        } else {
                                            currentPortArrow_->setColorConnectableMode(ConnectionBaseGraphicsItem::CCM_MAYBE);
                                        }
                                        break;
                                    } else {
                                        if(!destPort->isConnected() || destPort->allowMultipleConnections())
                                            isNext = false;
                                    }
                                }
                            } else {
                                foreach(Port* destPort, pOwner->getInports()){
                                    Port::ConnectivityState cState = getPort()->testConnectivity(destPort);
                                    if (cState != Port::ILLEGAL_CONNECTIVITY) {
                                        if(isNext) {
                                            if(cState == Port::LEGAL_CONNECTIVITY)
                                                currentPortArrow_->setColorConnectableMode(ConnectionBaseGraphicsItem::CCM_YES);
                                            else
                                                currentPortArrow_->setColorConnectableMode(ConnectionBaseGraphicsItem::CCM_MAYBE);
                                        } else {
                                            currentPortArrow_->setColorConnectableMode(ConnectionBaseGraphicsItem::CCM_MAYBE);
                                        }
                                        break;
                                    } else {
                                        if(!destPort->isConnected() || destPort->allowMultipleConnections())
                                            isNext = false;
                                    }
                                }
                            }
                            properItemFound = true;
                        } break;
                        default:
                            break;
                        }
                    }
                    if(properItemFound) break;
                } //end for
            }
        } else if(currentLayer() == NetworkEditorLayerPortSizeLinking) {
            tgtAssert(currentPortSizeLinkArrow_ != 0, "currentPortSizeLinkArrow was not set");

            currentPortSizeLinkArrow_->setDestinationPoint(event->scenePos());
            currentPortSizeLinkArrow_->setColorConnectableMode(ConnectionBaseGraphicsItem::CCM_DEFAULT);

            QList<QGraphicsItem*> items = scene()->items(event->scenePos());
            bool properItemFound = false;
            foreach(QGraphicsItem* item, items) {
                if (item && (item != this)) {
                    switch(item->type()) {
                    case UserTypesPortGraphicsItem: {
                        PortGraphicsItem* destPort = dynamic_cast<PortGraphicsItem*>(item);
                        PropertyGraphicsItem* propItem = 0;
                        foreach(PropertyGraphicsItem* prop, destPort->getPropertyList()->getPropertyItems())
                            if(prop->getProperty()->getClassName() == "RenderSizeOriginProperty" ||
                                prop->getProperty()->getClassName() == "RenderSizeReceiveProperty") {
                            propItem = prop;
                            break;
                        }
                        if(propItem){
                            if (currentPortSizeLinkArrow_->getSourceItem()->getProperty()->getClassName() == propItem->getProperty()->getClassName())
                                currentPortSizeLinkArrow_->setColorConnectableMode(ConnectionBaseGraphicsItem::CCM_MAYBE);
                            else
                            if (currentPortSizeLinkArrow_->getSourceItem()->getProperty()->getClassName() == "RenderSizeOriginProperty"){
                                if(propItem->getDestinationLinkList().empty())
                                    currentPortSizeLinkArrow_->setColorConnectableMode(ConnectionBaseGraphicsItem::CCM_YES);
                                else
                                    currentPortSizeLinkArrow_->setColorConnectableMode(ConnectionBaseGraphicsItem::CCM_MAYBE);
                            } else
                                currentPortSizeLinkArrow_->setColorConnectableMode(ConnectionBaseGraphicsItem::CCM_NO);
                        }
                        properItemFound = true;
                    } break;
                    case UserTypesPropertyListButtonGraphicsItem:
                    case UserTypesWidgetToggleButtonGraphicsItem:
                    case UserTypesRenamableTextGraphicsItem:
                        item = item->parentItem();
                    case UserTypesProcessorGraphicsItem: {
                    //case UserTypesAggregationGraphicsItem:
                        PortOwnerGraphicsItem* pOwner = dynamic_cast<PortOwnerGraphicsItem*>(item);
                        currentPortSizeLinkArrow_->setColorConnectableMode(ConnectionBaseGraphicsItem::CCM_NO);
                        if(getPort()->isOutport()){ //outport can not be connected with processor
                            properItemFound = true;
                            break;
                        }
                        foreach(Port* dPort, pOwner->getOutports()){
                            PortGraphicsItem* destPort = pOwner->getPortGraphicsItem(dPort);
                            foreach(PropertyGraphicsItem* prop, destPort->getPropertyList()->getPropertyItems()) {
                                if(prop->getProperty()->getClassName() == "RenderSizeReceiveProperty" && prop->getDestinationLinkList().empty()) {
                                    currentPortSizeLinkArrow_->setColorConnectableMode(ConnectionBaseGraphicsItem::CCM_YES);
                                    properItemFound = true;
                                    break;
                                }
                            }
                            if(properItemFound) break;
                        }
                        properItemFound = true;
                    } break;
                    default:
                        break;
                    }
                    if(properItemFound) break;
                } //if item != this
            } //end for
        } //if layer
    }
    NWEBaseGraphicsItem::mouseMoveEvent(event);
}

void PortGraphicsItem::mouseReleaseEvent(QGraphicsSceneMouseEvent* event) {
    if (currentPortArrow_) {
        QGraphicsItem* item = scene()->itemAt(event->scenePos(), QTransform());
        QList<QGraphicsItem*> items = scene()->items(event->scenePos());
        bool properItemFound = false;
        foreach(QGraphicsItem* item, items) {
            if (item && (item != this)){
                switch(item->type()) {
                case UserTypesPortGraphicsItem: {
                    Port* destPort = dynamic_cast<PortGraphicsItem*>(item)->getPort();
                    Port::ConnectivityState cState = getPort()->testConnectivity(destPort);
                    bool allowConnection = cState != Port::ILLEGAL_CONNECTIVITY;
                    if(cState == Port::UNDEFINED_LOOP_CONNECTIVITY) {
                        allowConnection = QMessageBox::warning(0, "Undefined loop", "The proposed port connection creates a loop within the data flow. "
                                                                                    "The resulting behaviour is potentially undefined. Would you like to add the connection anyway?",
                                                               QMessageBox::Yes, QMessageBox::No) == QMessageBox::Yes;
                    }
                    if (allowConnection)
                        portOwner_->getNetworkEditor()->getProcessorNetwork()->connectPorts(getPort(), destPort);
                    properItemFound = true;
                } break;
                case UserTypesPropertyListButtonGraphicsItem:
                case UserTypesWidgetToggleButtonGraphicsItem:
                case QGraphicsTextItem::Type:
                case UserTypesRenamableTextGraphicsItem:
                    item = item->parentItem();
                case UserTypesProcessorGraphicsItem: {
                //case UserTypesAggregationGraphicsItem:
                    PortOwnerGraphicsItem* pOwner = dynamic_cast<PortOwnerGraphicsItem*>(item);
                    if(dynamic_cast<CoProcessorPort*>(getPort())) {
                        foreach(CoProcessorPort* destPort, pOwner->getCoProcessorInports())
                            if (getPort()->testConnectivity(destPort) == Port::LEGAL_CONNECTIVITY) {
                                portOwner_->getNetworkEditor()->getProcessorNetwork()->connectPorts(getPort(), destPort);
                                break;
                            }
                    } else {
                        foreach(Port* destPort, pOwner->getInports())
                            if (getPort()->testConnectivity(destPort) == Port::LEGAL_CONNECTIVITY) {
                                portOwner_->getNetworkEditor()->getProcessorNetwork()->connectPorts(getPort(), destPort);
                                break;
                            }
                    }
                    properItemFound = true;
                } break;
                default:
                    break;
                }
            }
            if(properItemFound) break;
        } //end for
        currentPortArrow_->scene()->removeItem(currentPortArrow_);
        delete currentPortArrow_;
        currentPortArrow_ = 0;
        emit endedArrow();
    } else if (currentPortSizeLinkArrow_) {
        QList<QGraphicsItem*> items = scene()->items(event->scenePos());
        bool properItemFound = false;
        foreach(QGraphicsItem* item, items) {
            if (item && (item != this)) {
                switch(item->type()) {
                case UserTypesPortGraphicsItem: {
                    PortGraphicsItem* destPort = dynamic_cast<PortGraphicsItem*>(item);
                    PropertyGraphicsItem* propItem = 0;
                    foreach(PropertyGraphicsItem* prop, destPort->getPropertyList()->getPropertyItems())
                        if(prop->getProperty()->getClassName() == "RenderSizeOriginProperty" ||
                            prop->getProperty()->getClassName() == "RenderSizeReceiveProperty") {
                            propItem = prop;
                            break;
                        }

                    if(propItem && (currentPortSizeLinkArrow_->getSourceItem()->getProperty()->getClassName() == propItem->getProperty()->getClassName() ||
                                    currentPortSizeLinkArrow_->getSourceItem()->getProperty()->getClassName() == "RenderSizeOriginProperty"))
                        getNetworkEditor()->getProcessorNetwork()->createPropertyLink(const_cast<Property*>(currentPortSizeLinkArrow_->getSourceItem()->getProperty()),
                                                                                      const_cast<Property*>(propItem->getProperty()),
                                                                                      LinkEvaluatorHelper::createEvaluator("LinkEvaluatorRenderSize"));
                    properItemFound = true;
                } break;
                case UserTypesPropertyListButtonGraphicsItem:
                case UserTypesWidgetToggleButtonGraphicsItem:
                case QGraphicsTextItem::Type:
                case UserTypesRenamableTextGraphicsItem:
                    item = item->parentItem();
                case UserTypesProcessorGraphicsItem: {
                //case UserTypesAggregationGraphicsItem:
                    PortOwnerGraphicsItem* pOwner = dynamic_cast<PortOwnerGraphicsItem*>(item);
                    if(getPort()->isOutport()){ //outport can not be connected with processor
                        properItemFound = true;
                        break;
                    }
                    foreach(Port* dPort, pOwner->getOutports()){
                        PortGraphicsItem* destPort = pOwner->getPortGraphicsItem(dPort);
                        foreach(PropertyGraphicsItem* prop, destPort->getPropertyList()->getPropertyItems()) {
                            if(prop->getProperty()->getClassName() == "RenderSizeReceiveProperty" && prop->getDestinationLinkList().empty()) {
                                getNetworkEditor()->getProcessorNetwork()->createPropertyLink(const_cast<Property*>(currentPortSizeLinkArrow_->getSourceItem()->getProperty()),
                                                                                      const_cast<Property*>(prop->getProperty()),
                                                                                      LinkEvaluatorHelper::createEvaluator("LinkEvaluatorRenderSize"));
                                properItemFound = true;
                                break;
                            }
                        }
                        if(properItemFound) break;
                    }
                    properItemFound = true;
                } break;
                default:
                    break;
                }
                if(properItemFound) break;
            } //if item != this
        } //end for
        delete currentPortSizeLinkArrow_;
        currentPortSizeLinkArrow_ = 0;
        emit endedArrow();
    }
    NWEBaseGraphicsItem::mouseReleaseEvent(event);
}

} // namespace
