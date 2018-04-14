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

#include "voreen/qt/networkeditor/graphicitems/core/processorgraphicsitem.h"

#include "voreen/qt/networkeditor/styles/nwestyle_base.h"
#include "voreen/qt/networkeditor/editor_settings.h"
//core
#include "voreen/core/processors/processor.h"
#include "voreen/core/processors/processorwidget.h"
#include "voreen/core/ports/port.h"
//gi
#include "voreen/qt/networkeditor/graphicitems/utils/progressbargraphicsitem.h"
#include "voreen/qt/networkeditor/graphicitems/core/portgraphicsitem.h"
#include "voreen/qt/networkeditor/graphicitems/connections/portarrowgraphicsitem.h"
#include "voreen/qt/networkeditor/graphicitems/tooltips/tooltipprocessorgraphicsitem.h"
//meta
#include "voreen/core/datastructures/meta/positionmetadata.h"

//#include <QStyleOptionGraphicsItem>
#include <QGraphicsSceneMouseEvent>

namespace voreen {


ProcessorGraphicsItem::ProcessorGraphicsItem(Processor* processor, NetworkEditor* networkEditor)
    : PortOwnerGraphicsItem(processor, networkEditor)
    , processor_(processor)
{
    tgtAssert(processor_ != 0, "passed null pointer");

    setGuiName(QString::fromStdString(processor_->getID()));

    //set text
    initializePaintSettings();
    paintHasBeenInitialized_ = true;
    getPropertyList()->updateHeaderLabelItem();

    registerPorts();

    if (processor_->usesExpensiveComputation()) {
        QRectF boundRect = boundingRect();
        QPointF center(boundRect.x() + boundRect.width() / 2.0, boundRect.y() + boundRect.height() * 0.7f); //0.775
        qreal width = boundRect.width() * 0.8;
        qreal height = 8;
        progressBar_ = new ProgressBarGraphicsItem(this, center, width, height);
        processor_->addProgressBar(progressBar_);
    }
    //connect processor widget button with processor widget
    dynamic_cast<Observable<ProcessorObserver>* >(processor_)->addObserver(this);
    //HACK: in new workspaces the processor widget has already been created
    if(processor_->getProcessorWidget())
        processorWidgetCreated(processor_);

    setToolTipGraphicsItem(new ToolTipProcessorGraphicsItem(this));

}

ProcessorGraphicsItem::~ProcessorGraphicsItem() {
    if(progressBar_) {
        processor_->removeProgressBar(progressBar_);
    }
}

//---------------------------------------------------------------------------------------------------------------
//                  nwebasegraphicsitem functions
//---------------------------------------------------------------------------------------------------------------
QRectF ProcessorGraphicsItem::boundingRect() const {
    return currentStyle()->ProcessorGI_boundingRect(this);
}

QPainterPath ProcessorGraphicsItem::shape() const {
    return currentStyle()->ProcessorGI_shape(this);
}

void ProcessorGraphicsItem::initializePaintSettings(){
    currentStyle()->ProcessorGI_initializePaintSettings(this);
}

void ProcessorGraphicsItem::mainPaint(QPainter* painter, const QStyleOptionGraphicsItem* option, QWidget* widget, NWEItemSettings& setting) {
    currentStyle()->ProcessorGI_paint(this, painter, option, widget, setting);
}

void ProcessorGraphicsItem::layoutChildItems() {
    if (progressBar_) {
        QRectF br = boundingRect();
        QPointF center(br.x() + br.width() / 2.0, br.y() + br.height() * 0.7f); //0.775
        qreal width = br.width() * 0.8;
        qreal height = 8;

        progressBar_->resize(center, width, height);
    }
    PortOwnerGraphicsItem::layoutChildItems();
}

//---------------------------------------------------------------------------------------------------------------
//                  portownergraphicsitem functions
//---------------------------------------------------------------------------------------------------------------
void ProcessorGraphicsItem::registerPorts() {
    inports_ = stdVectorToQList<Port*>(processor_->getInports());
    outports_ = stdVectorToQList<Port*>(processor_->getOutports());
    coInports_ = stdVectorToQList<CoProcessorPort*>(processor_->getCoProcessorInports());
    coOutports_ = stdVectorToQList<CoProcessorPort*>(processor_->getCoProcessorOutports());

    createChildItems();
}

void ProcessorGraphicsItem::renameFinished() {
    nameLabel_.setTextInteractionFlags(Qt::NoTextInteraction);
    nameLabel_.setFlag(QGraphicsItem::ItemIsFocusable, false);
    nameLabel_.setFlag(QGraphicsItem::ItemIsSelectable, false);

    // name unchanged
    if (nameLabel_.toPlainText().toStdString() == processor_->getID())
        return;

    try {
        ProcessorNetwork* network = networkEditor_->getProcessorNetwork();
        std::string name = network->generateUniqueProcessorName(nameLabel_.toPlainText().toStdString());
        network->setProcessorName(processor_, name);
        setGuiName(QString::fromStdString(name));
        nameLabel_.setPlainText(QString::fromStdString(name));
        getPropertyList()->updateHeaderLabelItem();
    }
    catch (VoreenException& e) {
        LWARNINGC("voreen.qt.ProcessorGraphicsItem", e.what());
        nameLabel_.setPlainText(QString::fromStdString(processor_->getID()));
    }

    nameChanged();
}

QList<Processor*> ProcessorGraphicsItem::getProcessors() const {
    QList<Processor*> result;
    result.append(processor_);
    return result;
}

void ProcessorGraphicsItem::toggleProcessorWidget() {
    processor_->getProcessorWidget()->setVisible(!processor_->getProcessorWidget()->isVisible());
    widgetToggleButton_.update();
}

void ProcessorGraphicsItem::togglePropertyList() {
    if (!networkEditor_ || !networkEditor_->updatesEnabled())
        return;

    prepareGeometryChange();
    if (propertyList_.getIsVisibleInEditor()){
        propertyList_.setVisible(false);
        propertyList_.setIsVisibleInEditor(false);
    } else {
        //layoutChildItems();
        propertyList_.setIsVisibleInEditor(true);
        if(currentLayer() == NetworkEditorLayerGeneralLinking){
            propertyList_.setPropertyVisibleMode(PropertyListGraphicsItem::SHOW_LINKED_PROPERTIES);
            if(getNetworkEditor()->cameraLinksHidden() && getNetworkEditor()->portSizeLinksHidden())
                propertyList_.setPropertyVisibleModifier(PropertyListGraphicsItem::HIDE_CAMERA_AND_SIZE_PROPERTIES);
            else if(getNetworkEditor()->cameraLinksHidden())
                propertyList_.setPropertyVisibleModifier(PropertyListGraphicsItem::HIDE_CAMERA_PROPERTIES);
            else if(getNetworkEditor()->portSizeLinksHidden())
                propertyList_.setPropertyVisibleModifier(PropertyListGraphicsItem::HIDE_SIZE_PROPERTIES);
            else
                propertyList_.setPropertyVisibleModifier(PropertyListGraphicsItem::HIDE_NO_PROPERTIES);
        }
        else if (currentLayer() == NetworkEditorLayerCameraLinking) {
            propertyList_.setPropertyVisibleMode(PropertyListGraphicsItem::SHOW_CAMERA_PROPERTIES);
            propertyList_.setPropertyVisibleModifier(PropertyListGraphicsItem::HIDE_NO_PROPERTIES);
        }
        propertyList_.setHeaderMode(PropertyListGraphicsItem::NO_HEADER);
        propertyList_.setLinkArrowMode(PropertyListGraphicsItem::SHOW_ALL_ARROWS);
        propertyList_.setVisible(true);
        if(propertyList_.getAllVisiblePropertyItems().empty()) {
            propertyList_.setVisible(false);
            propertyList_.setIsVisibleInEditor(false);
        }
    }
    layoutChildItems();

    networkEditor_->resetCachedContent();
}

void ProcessorGraphicsItem::saveMeta() {
    PositionMetaData* meta = new PositionMetaData(static_cast<int>(x()), static_cast<int>(y()));

    processor_->getMetaDataContainer().addMetaData("ProcessorGraphicsItem", meta);
}

void ProcessorGraphicsItem::loadMeta() {
    if (processor_->getMetaDataContainer().hasMetaData("ProcessorGraphicsItem")) {
        PositionMetaData* meta = dynamic_cast<PositionMetaData*>(processor_->getMetaDataContainer().getMetaData("ProcessorGraphicsItem"));
        setPos(meta->getX(), meta->getY());
    }
    layoutChildItems();
}

//---------------------------------------------------------------------------------------------------------------
//                  processorobserver functions
//---------------------------------------------------------------------------------------------------------------
void ProcessorGraphicsItem::processorWidgetCreated(const Processor* processor) {
    widgetToggleButton_.show();
    widgetToggleButton_.setProcessorWidget(processor->getProcessorWidget());
}

void ProcessorGraphicsItem::processorWidgetDeleted(const Processor*) {
    widgetToggleButton_.hide();
    widgetToggleButton_.setProcessorWidget(0);
}

void ProcessorGraphicsItem::propertiesChanged(const PropertyOwner*) {
    //RootGraphicsItem::propertiesChanged();
}

void ProcessorGraphicsItem::portsChanged(const Processor*) {
    std::vector<std::pair<std::string,PortGraphicsItem*> > connections;
    foreach (PortGraphicsItem* ownPort, portItems_) {
        foreach (PortArrowGraphicsItem* connection, ownPort->getArrowList()) {
            connections.push_back(std::pair<std::string,PortGraphicsItem*>(ownPort->getPort()->getID(), connection->getDestinationItem()));
        }
        const QList<PortArrowGraphicsItem*> connectedPorts = ownPort->getArrowList();
        foreach (PortGraphicsItem* externalPort, ownPort->getConnectedOutportGraphicsItems()) {
            connections.push_back(std::pair<std::string,PortGraphicsItem*>(ownPort->getPort()->getID(), externalPort));
        }
    }
    deleteChildItems();
    registerPorts();
    for (size_t i=0; i<connections.size(); i++) {
        foreach (PortGraphicsItem* p, portItems_) {
            if (p->getPort()->getID() == connections[i].first) {
                if(p->isOutport())
                    p->addGraphicalConnection(connections[i].second);
                else
                    connections[i].second->addGraphicalConnection(p);
            }
        }
    }
    if (scene())
        scene()->invalidate();
}

void ProcessorGraphicsItem::stateChanged(const Processor* processor) {
    update();
}

//---------------------------------------------------------------------------------------------------------------
//                  processor functions
//---------------------------------------------------------------------------------------------------------------
Processor* ProcessorGraphicsItem::getProcessor() const {
    return processor_;
}

//---------------------------------------------------------------------------------------------------------------
//                  events
//---------------------------------------------------------------------------------------------------------------
void ProcessorGraphicsItem::mouseDoubleClickEvent(QGraphicsSceneMouseEvent* event) {
    if(event->buttons() != 1)
        return;
    switch(currentLayer()){
    case NetworkEditorLayerDataFlow:
        // toggle visibility of processor widget
        if (processor_->getProcessorWidget())
            toggleProcessorWidget();
        break;
    case NetworkEditorLayerGeneralLinking:
    case NetworkEditorLayerCameraLinking:
        togglePropertyList();
        break;
    default:
        break;
    }
    event->accept();
}

} // namespace voreen
