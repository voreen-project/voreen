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

//super class
#include "voreen/qt/networkeditor/networkeditor.h"
#include "voreen/qt/networkeditor/editor_settings.h"
//general includes
#include "voreen/core/network/networkevaluator.h"
#include "voreen/core/network/processornetwork.h"
#include "voreen/core/voreenapplication.h"
#include "voreen/qt/voreenapplicationqt.h"

//meta informations
#include "voreen/core/datastructures/meta/serializablevectormetadata.h"
#include "voreen/qt/networkeditor/meta/textboxmetadata.h"
#include "voreen/core/datastructures/meta/positionmetadata.h"

//processors
#include "voreen/core/processors/processor.h"
#include "voreen/core/processors/renderprocessor.h"

//properties
#include "voreen/core/properties/cameraproperty.h"
#include "voreen/core/properties/link/linkevaluatorid.h"
#include "voreen/core/datastructures/callback/memberfunctioncallback.h"

//dialogs
#include "voreen/qt/networkeditor/dialogs/propertylinkdialog.h"

//styles
#include "voreen/qt/networkeditor/styles/nwestyle_base.h"
#include "voreen/qt/networkeditor/styles/nwestyle_classic.h"
#include "voreen/qt/networkeditor/styles/nwestyle_classic_print.h"
#include "voreen/qt/networkeditor/styles/nwestyle_material.h"

//graph layouts
#include "voreen/qt/networkeditor/graphlayouts/nwegl_base.h"
#include "voreen/qt/networkeditor/graphlayouts/nwegl_sugiyama.h"

//graphic items
    //core
#include "voreen/qt/networkeditor/graphicitems/core/processorgraphicsitem.h"
#include "voreen/qt/networkeditor/graphicitems/core/portgraphicsitem.h"
#include "voreen/qt/networkeditor/graphicitems/core/portownergraphicsitem.h"
#include "voreen/qt/networkeditor/graphicitems/core/propertygraphicsitem.h"
    //connections
#include "voreen/qt/networkeditor/graphicitems/connections/portarrowgraphicsitem.h"
#include "voreen/qt/networkeditor/graphicitems/connections/propertylinkarrowgraphicsitem.h"
#include "voreen/qt/networkeditor/graphicitems/connections/portownerlinkarrowgraphicsitem.h"
#include "voreen/qt/networkeditor/graphicitems/connections/portsizelinkarrowgraphicsitem.h"
    // utils
#include "voreen/qt/networkeditor/graphicitems/utils/widgettogglebuttongraphicsitem.h"
    //textboxes
#include "voreen/qt/networkeditor/graphicitems/textboxes/frameboxgraphicsitem.h"
#include "voreen/qt/networkeditor/graphicitems/textboxes/textboxbasegraphicsitem.h"
#include "voreen/qt/networkeditor/graphicitems/textboxes/textboxgraphicsitem.h"

#include "voreen/qt/networkeditor/graphicitems/utils/progressbargraphicsitem.h"

// editor buttons
#include <QToolButton>
#include <QBoxLayout>
#include <QButtonGroup>
#include <QMessageBox>
#include <QTextStream>
#include <QComboBox>
#include <QSpinBox>
#include <QSizePolicy>
#include <QLabel>
#include <QCheckBox>
#include <QPushButton>
#include <QTextEdit>
#include <QGraphicsTextItem>
#include <QApplication>
#include <QClipboard>
#include <QMenu>
#include <QList>
#include <QMimeData>
//events
#include <QWheelEvent>
#include <QScrollBar>


namespace voreen {

const std::string NetworkEditor::loggerCat_("voreenve.NetworkEditor");

NetworkEditor::NetworkEditor(QWidget* parent, NetworkEvaluator* evaluator)
    : QGraphicsView(parent), ProcessorNetworkObserver()
    //general members
    , workspace_(0)
    , evaluator_(evaluator)
    //properties
    , scaleProcessorFontSizeProperty_("scaleProcessorFontSize","Scale Processor Fonts:",100,50,150)
    , networkEditorStyleProperty_("networkEditorStyleProperty","Style:")
    , networkEditorGraphLayoutsProperty_("networkEditorGraphLayoutsProperty","Graph Layout:")
    , sugiShiftXProperty_("sugiShiftXProperty", "Processor Gap:",300.f,100.f,1000.f)
    , sugiOverlapProperty_("sugiOverlapProperty", "Overlapping", false)
    , sugiMedianProperty_("sugiMedianProperty", "Median", true)
    , sugiPortFlushProperty_("sugiPortFlushProperty", "Port Alignment", true)
    , showDocumentationProperty_("showDocumentationProperty","should not be visible",true)
    //style, layer and cursor
    , currentLayer_(NetworkEditorLayerUndefined)
    , currentCursorMode_(NetworkEditorCursorSelectMode)
    , currentStyle_(0)
    , currentGraphLayout_(0)
    , currentToolTipMode_(true)
    //editor buttons
    , layerButtonContainer_(0), dataFlowLayerButton_(0), linkingLayerButton_(0)
    , linkingLayerButtonContainer_(0), generalLinkingButton_(0), cameraLinkingButton_(0), portSizeLinkingButton_(0)
    , generalLinkingLayerButtonContainer_(0), hideCameraLinksButton_(0), hidePortSizeLinksButton_(0), removeAllPropertyLinksButton_(0)
    , cameraLinkingLayerButtonContainer_(0), linkCamerasAutoButton_(0), linkCamerasButton_(0), removeAllCameraLinksButton_(0)
    , portSizeLinkingLayerButtonContainer_(0), linkPortSizeAutoButton_(0), linkPortSizeButton_(0), removeAllPortSizeLinksButton_(0)
    , stopButtonContainer_(0), stopNetworkEvaluatorButton_(0), networkEvaluatorIsLockedByButton_(false)
    , navigationButtonContainer_(0), selectCursorButton_(0), moveCursorButton_(0)
    , layoutButtonContainer_(0), centerViewButton_(0), graphLayoutButton_(0)
    //scale
    , needsScale_(false)
    //navigation
    , translateScene_(false)
    , translateSceneVector_(0.0,0.0)
    , lastTranslateCenter_(0.0,0.0)
    //drag&drop
    , selectedItem_(nullptr)
{
    tgtAssert(evaluator_ != 0, "passed null pointer");

    //create scene
    setScene(new QGraphicsScene(this));
    scene()->setSceneRect(sceneRectSpacing.x(), sceneRectSpacing.y(), sceneRectSpacing.width(), sceneRectSpacing.height());
    scene()->setItemIndexMethod(QGraphicsScene::NoIndex);
    connect(scene(), SIGNAL(changed(QList<QRectF>)), this, SIGNAL(networkEditor_visibleSceneSizeChanged_Signal()));

    //network editor properties
    addProperty(scaleProcessorFontSizeProperty_);
        scaleProcessorFontSizeProperty_.setGroupID("style");
    addProperty(networkEditorStyleProperty_);
        networkEditorStyleProperty_.addOption("first","Classic",NWESTYLE_CLASSIC);
        networkEditorStyleProperty_.addOption("second","Classic Print",NWESTYLE_CLASSIC_PRINT);
        networkEditorStyleProperty_.addOption("third","Material",NWESTYLE_MATERIAL);
        networkEditorStyleProperty_.setGroupID("style");
    addProperty(networkEditorGraphLayoutsProperty_);
        networkEditorGraphLayoutsProperty_.addOption("first","Sugiyama",NWEGL_SUGIYAMA);
        networkEditorGraphLayoutsProperty_.setGroupID("layout");
    addProperty(sugiShiftXProperty_);
        sugiShiftXProperty_.setGroupID("layout");
    addProperty(sugiOverlapProperty_);
        sugiOverlapProperty_.setGroupID("layout");
    addProperty(sugiMedianProperty_);
        sugiMedianProperty_.setGroupID("layout");
    addProperty(sugiPortFlushProperty_);
        sugiPortFlushProperty_.setGroupID("layout");
    //setPropertyGroupGuiName("lod","lod");
    setPropertyGroupGuiName("style","Graphics Item Style");
    setPropertyGroupGuiName("layout","Graph Layouting");
    addProperty(showDocumentationProperty_);
    showDocumentationProperty_.setVisibleFlag(false);
    //deserialize editor settings
    std::string filename = VoreenApplication::app()->getUserDataPath("networkeditor_settings.xml");
    if (!deserializeSettings(this, filename))
        LWARNING("Failed to deserialize networkeditor settings! Ignore on first start of voreen.");
    //set current Style
    updateStyle();
    //set current graph layout
    switch(networkEditorGraphLayoutsProperty_.getValue()){
    case NWEGL_SUGIYAMA: {
        currentGraphLayout_ = new NWEGL_Sugiyama();
        qreal shift = 300.f; bool overlap = false, median = true, portflush = true;
        shift = sugiShiftXProperty_.get();
        overlap = sugiOverlapProperty_.get();
        median = sugiMedianProperty_.get();
        portflush = sugiPortFlushProperty_.get();
        static_cast<NWEGL_Sugiyama*>(currentGraphLayout_)->setSortParameter(shift,overlap,median,portflush);
        } break;
    default:
        tgtAssert(false,"Unknown NetworkEditorGraphLayout!!!");
        LERROR("Unknown NetworkEditorGraphLayout!!!");
        break;
    }
    //link properties
    ON_PROPERTY_CHANGE(networkEditorStyleProperty_,NetworkEditor,updateStyle);
    ON_PROPERTY_CHANGE(scaleProcessorFontSizeProperty_,NetworkEditor,processorFontOnChange);
    ON_PROPERTY_CHANGE(sugiShiftXProperty_,NetworkEditor,updateGraphLayout);
    ON_PROPERTY_CHANGE(sugiOverlapProperty_,NetworkEditor,updateGraphLayout);
    ON_PROPERTY_CHANGE(sugiMedianProperty_,NetworkEditor,updateGraphLayout);
    ON_PROPERTY_CHANGE(sugiPortFlushProperty_,NetworkEditor,updateGraphLayout);

    //paint settings
    setBackgroundBrush(currentStyle_->getBackgroundBrush());
    setCacheMode(CacheBackground);
    setRenderHint(QPainter::Antialiasing);
    setTransformationAnchor(QGraphicsView::AnchorUnderMouse);
    setResizeAnchor(QGraphicsView::AnchorViewCenter);
    setMouseTracking(true);
    setDragMode(QGraphicsView::RubberBandDrag);
    setMinimumSize(QSize(400, 400));
    setVerticalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
    setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);

    //initilize buttons
    initilizeEditorButtons();

    //context menu
    createContextMenuActions();

    // Signal that this widget accepts drop events (e.g., dropping a processor from the list onto the canvas).
    setAcceptDrops(true);
}

NetworkEditor::~NetworkEditor() {
    //store property settings
    std::string filename = VoreenApplication::app()->getUserDataPath("networkeditor_settings.xml");
    if(!serializeSettings(this, filename))
        LWARNING("Failed to save networkeditor settings");
    //delete all created items/scenes
    delete scene();
    setScene(0);
    delete currentStyle_;
    delete currentGraphLayout_;
}
//---------------------------------------------------------------------------------------------------------------
//                  general members and functions
//---------------------------------------------------------------------------------------------------------------
void NetworkEditor::setWorkspace(Workspace* workspace)  {
    //hideTooltip();

    if (workspace_ && workspace_->getProcessorNetwork() && workspace_ != workspace)
        workspace_->getProcessorNetwork()->removeObserver(this);

    resetScene();
    //linkMap_.clear();

    //processorNetwork_ = network;
    workspace_ = workspace;
    ProcessorNetwork* processorNetwork = getProcessorNetwork();
    if (processorNetwork)
        processorNetwork->addObserver(this);

    generateGraphicsItems();

    resetMatrix();
    resetTransform();

    // do not scale immediately when we might be in visualization mode and this window might
    // not have its final size yet
    if (isVisible()) {
        scaleView();
        needsScale_ = false;
    } else {
        needsScale_ = true;
    }

    dataFlowLayerButton_->setChecked(true);
    setLayer(NetworkEditorLayerDataFlow);

    // set state of cam auto-linking button according to stored meta data, if present
    BoolMetaData* autoMeta = 0;
    if (processorNetwork)
        autoMeta = dynamic_cast<BoolMetaData*>(processorNetwork->getMetaDataContainer().getMetaData("autoLinkCameras"));
    if (autoMeta)
        linkCamerasAutoButton_->setChecked(autoMeta->getValue());
    else
        linkCamerasAutoButton_->setChecked(true);
    autoMeta = 0;
    if (processorNetwork)
        autoMeta = dynamic_cast<BoolMetaData*>(processorNetwork->getMetaDataContainer().getMetaData("autoLinkPortSize"));
    if (autoMeta)
        linkPortSizeAutoButton_->setChecked(autoMeta->getValue());
    else
        linkPortSizeAutoButton_->setChecked(true);

    // camera auto linking
    if (processorNetwork && linkCamerasAutoButton_->isChecked()) {
        int numCreated = processorNetwork->createPropertyLinksWithinSubNetwork<CameraProperty>(
            processorNetwork->getProcessors(), std::vector<std::string>(), new LinkEvaluatorCameraId());
        if (numCreated)
            LINFO("Created " << numCreated << " camera property links.");
    }

    // render size auto linking
    if (processorNetwork && linkPortSizeAutoButton_->isChecked()) {
        int numCreated = processorNetwork->createRenderSizeLinksWithinSubNetwork(processorNetwork->getProcessors(), false);
        if (numCreated)
            LINFO("Created " << numCreated << " render size links.");
    }

    //load documentation items
    if(workspace) {
        deserializeTextItems();
    }
}

Workspace* NetworkEditor::getWorkspace() const {
    return workspace_;
}

ProcessorNetwork* NetworkEditor::getProcessorNetwork() const {
    if (workspace_)
        return workspace_->getProcessorNetwork();
    else
        return 0;
}

float NetworkEditor::getProcessorFontSizeScale() const {
    return static_cast<float>(scaleProcessorFontSizeProperty_.get());
}

void NetworkEditor::createContextMenuActions() {
    copyAction_ = new QAction(QIcon(":/qt/icons/edit-copy.png"), tr("Copy Items"), this);
    pasteAction_ = new QAction(QIcon(":/qt/icons/edit-paste.png"), tr("Paste Items"), this);
    deleteAction_ = new QAction(QIcon(":/qt/icons/eraser.png"), tr("Delete Items"), this);
    editLinkAction_ = new QAction(tr("Edit Link"), this);
    deleteAllLinksAction_ = new QAction(QIcon(":/qt/icons/eraser.png"), tr("Delete All Property Links"), this);
    deleteAllCameraLinksAction_ = new QAction(QIcon(":/qt/icons/eraser.png"), tr("Delete All Camera Links"), this);
    deleteAllPortSizeLinksAction_ = new QAction(QIcon(":/qt/icons/eraser.png"), tr("Delete All Render Size Links"), this);
    createAllCameraLinksAction_ = new QAction(tr("Create All Camera Links"), this);
    createAllPortSizeLinksAction_ = new QAction(tr("Create All Render Size Links"), this);
    deleteInnerLinksAction_ = new QAction(QIcon(":/qt/icons/eraser.png"), tr("Delete Inner Property Links"), this);
    deleteInnerCameraLinksAction_ = new QAction(QIcon(":/qt/icons/eraser.png"), tr("Delete Inner Camera Links"), this);
    deleteInnerPortSizeLinksAction_ = new QAction(QIcon(":/qt/icons/eraser.png"), tr("Delete Inner Render Size Links"), this);
    createInnerCameraLinksAction_ = new QAction(tr("Create Inner Camera Links"), this);
    createInnerPortSizeLinksAction_ = new QAction(tr("Create Inner Render Size Links"), this);
    deletePortOwnerLinksAction_ = new QAction(QIcon(":/qt/icons/eraser.png"), tr("Delete Processor Property Links"), this);
    deletePortOwnerCameraLinksAction_ = new QAction(QIcon(":/qt/icons/eraser.png"), tr("Delete Processor Camera Links"), this);
    deletePortOwnerPortSizeLinksAction_ = new QAction(QIcon(":/qt/icons/eraser.png"), tr("Delete Processor Render Size Links"), this);
    createPortOwnerCameraLinksAction_ = new QAction(tr("Create Processor Camera Links"), this);
    createPortOwnerPortSizeLinksAction_ = new QAction(tr("Create Processor Render Size Links"), this);
    sortSubNetworkAction_ = new QAction(QIcon(":/qt/icons/sortGraph.png"), tr("Sort Selected Processors"), this);
    createNewTextNoteAction_ = new QAction(QIcon(":/qt/icons/textnote-icon.png"), tr("Add Text Note"), this);
    createNewTextFrameAction_ = new QAction(QIcon(":/qt/icons/frame-gray.png"), tr("Add Colored Frame"), this);

    connect(copyAction_, SIGNAL(triggered()), this, SLOT(copyActionSlot()));
    connect(pasteAction_, SIGNAL(triggered()), this, SLOT(pasteActionSlot()));
    connect(deleteAction_, SIGNAL(triggered()), this, SLOT(deleteActionSlot()));
    connect(editLinkAction_, SIGNAL(triggered()), this, SLOT(editLinkActionSlot()));
    connect(deleteAllLinksAction_, SIGNAL(triggered()), this, SLOT(deleteLinksActionSlot()));
    connect(deleteAllCameraLinksAction_, SIGNAL(triggered()), this, SLOT(deleteCameraLinksActionSlot()));
    connect(deleteAllPortSizeLinksAction_, SIGNAL(triggered()), this, SLOT(deletePortSizeLinksActionSlot()));
    connect(createAllCameraLinksAction_, SIGNAL(triggered()), this, SLOT(createCameraLinksActionSlot()));
    connect(createAllPortSizeLinksAction_, SIGNAL(triggered()), this, SLOT(createPortSizeLinksActionSlot()));
    connect(deleteInnerLinksAction_, SIGNAL(triggered()), this, SLOT(deleteLinksActionSlot()));
    connect(deleteInnerCameraLinksAction_, SIGNAL(triggered()), this, SLOT(deleteCameraLinksActionSlot()));
    connect(deleteInnerPortSizeLinksAction_, SIGNAL(triggered()), this, SLOT(deletePortSizeLinksActionSlot()));
    connect(createInnerCameraLinksAction_, SIGNAL(triggered()), this, SLOT(createCameraLinksActionSlot()));
    connect(createInnerPortSizeLinksAction_, SIGNAL(triggered()), this, SLOT(createPortSizeLinksActionSlot()));
    connect(deletePortOwnerLinksAction_, SIGNAL(triggered()), this, SLOT(deleteLinksActionSlot()));
    connect(deletePortOwnerCameraLinksAction_, SIGNAL(triggered()), this, SLOT(deleteCameraLinksActionSlot()));
    connect(deletePortOwnerPortSizeLinksAction_, SIGNAL(triggered()), this, SLOT(deletePortSizeLinksActionSlot()));
    connect(createPortOwnerCameraLinksAction_, SIGNAL(triggered()), this, SLOT(createCameraLinksActionSlot()));
    connect(createPortOwnerPortSizeLinksAction_, SIGNAL(triggered()), this, SLOT(createPortSizeLinksActionSlot()));
    connect(sortSubNetworkAction_, SIGNAL(triggered()), this, SLOT(sortSubNetwork()));
    connect(createNewTextNoteAction_, SIGNAL(triggered()), this, SLOT(createNewTextNoteSlot()));
    connect(createNewTextFrameAction_, SIGNAL(triggered()), this, SLOT(createNewTextFrameSlot()));

    copyAction_->setShortcut(Qt::CTRL + Qt::Key_C);
    pasteAction_->setShortcut(Qt::CTRL + Qt::Key_V);
    //replaceAction_
    deleteAction_->setShortcut(Qt::Key_Delete);

}

//---------------------------------------------------------------------------------------------------------------
//                  ProcessorNetworkObserver functions
//---------------------------------------------------------------------------------------------------------------
void NetworkEditor::networkChanged() {}

void NetworkEditor::processorAdded(const Processor* processor) {
    if (processorItemMap_.find(const_cast<Processor*>(processor)) == processorItemMap_.end()) {
        ProcessorGraphicsItem* item = createProcessorGraphicsItem(const_cast<Processor*>(processor));
        if (linkCamerasAutoButton_->isChecked())
            linkCamerasOfProcessor(processor);
    }
}

void NetworkEditor::processorAdded(QString id) {
    if (!getProcessorNetwork())
        return;

    if (evaluator_->isLocked() && !networkEvaluatorIsLockedByButton_) {
        QMessageBox::information(this, tr("Network Locked"), tr("The network is being evaluated, so dropping processors is not allowed"));
        return;
    }

    Processor* proc = dynamic_cast<Processor*>(VoreenApplication::app()->createSerializableType(id.toStdString()));
    if (!proc) {
        LDEBUG("Drop event contains no valid processor");
        return;
    }

    std::string processorName = proc->getClassName();
    tgtAssert(!processorName.empty(), "Processor class name is empty");

    QPoint position = viewport()->rect().center();
    while (itemAt(position) != nullptr) {
        position = QPoint(position.x() + 10, position.y() + 10);
    }

    getProcessorNetwork()->addProcessor(proc, processorName);

    tgtAssert(processorItemMap_.contains(proc), "processorItemMap didn't contain the processor");
    ProcessorGraphicsItem* item = processorItemMap_[proc];
    tgtAssert(item, "no ProcessorGraphicsItem for proc");
    item->setPos(mapToScene(position));
    item->saveMeta();
    // make sure that the added processor is initialized
    evaluator_->initializeNetwork();
    scene()->clearSelection();
}

void NetworkEditor::processorRemoved(const Processor* processor) {
    if (processorItemMap_.find(const_cast<Processor*>(processor)) != processorItemMap_.end()) {
        Processor* nonConst = const_cast<Processor*>(processor);
        ProcessorGraphicsItem* item = processorItemMap_.value(nonConst);
        setUpdatesEnabled(false);
        scene()->removeItem(item);
        delete item;
        setUpdatesEnabled(true);
        processorItemMap_.remove(nonConst);
    }
}

void NetworkEditor::propertyLinkAdded(const PropertyLink* link) {
    createLinkArrowForPropertyLink(link);
}

void NetworkEditor::propertyLinkRemoved(const PropertyLink* link) {
    PropertyGraphicsItem* srcProp = 0, *dstProp = 0;
    PortOwnerGraphicsItem* srcPOItem = 0, *dstPOItem = 0;

    //find src items
    foreach(ProcessorGraphicsItem* item, processorItemMap_) {
        if(item->getPropertyList()->getPropertyItem(link->getSourceProperty())){
            srcProp = item->getPropertyList()->getPropertyItem(link->getSourceProperty());
            srcPOItem = item;
            break;
        }
        foreach(PortGraphicsItem* portItem, item->getPortGraphicsItems()){
            if(portItem->getPropertyList()->getPropertyItem(link->getSourceProperty())){
                srcProp = portItem->getPropertyList()->getPropertyItem(link->getSourceProperty());
                srcPOItem = portItem->getPortOwner();
                break;
            }
        }
        if(srcProp)
            break;
    }

    //find dst items
    foreach(ProcessorGraphicsItem* item, processorItemMap_) {
        if(item->getPropertyList()->getPropertyItem(link->getDestinationProperty())){
            dstProp = item->getPropertyList()->getPropertyItem(link->getDestinationProperty());
            dstPOItem = item;
            break;
        }
        foreach(PortGraphicsItem* portItem, item->getPortGraphicsItems()){
            if(portItem->getPropertyList()->getPropertyItem(link->getDestinationProperty())){
                dstProp = portItem->getPropertyList()->getPropertyItem(link->getDestinationProperty());
                dstPOItem = portItem->getPortOwner();
                break;
            }
        }
        if(dstProp)
            break;
    }

    srcProp->removeGraphicalLink(dstProp);

    //remove portsizelinkarrow if exists
    if(((srcProp->getProperty()->getClassName() == "RenderSizeOriginProperty") ||
        (srcProp->getProperty()->getClassName() == "RenderSizeReceiveProperty")) &&
       ((dstProp->getProperty()->getClassName() == "RenderSizeOriginProperty") ||
        (dstProp->getProperty()->getClassName() == "RenderSizeReceiveProperty")))
        dynamic_cast<PortGraphicsItem*>(srcProp->getPropertyOwnerItem())->removePortSizeLinkArrow(srcProp, dstProp);


    //remove portownerlinkarrow if no longer needed
    if(!srcPOItem->isPortOwnerLinkNeeded(dstPOItem))
        srcPOItem->removeGraphicalLinkArrow(dstPOItem);
    else if(currentLayer_ == NetworkEditorLayerCameraLinking) {//test if it is still visible in camera mode
        if(!srcPOItem->isPortOwnerLinkNeeded(dstPOItem, currentLayer_)) {
            foreach(PortOwnerLinkArrowGraphicsItem* arrow, srcPOItem->getPortOwnerLinkArrows())
                if(arrow->getDestinationItem()->parent() == dstPOItem)
                    arrow->setVisible(false);
            foreach(PortOwnerLinkArrowGraphicsItem* arrow, dstPOItem->getPortOwnerLinkArrows())
                if(arrow->getDestinationItem()->parent() == srcPOItem)
                    arrow->setVisible(false);
        }
    }

    srcPOItem->getPropertyList()->setVisible(srcPOItem->getPropertyList()->isVisible());
    dstPOItem->getPropertyList()->setVisible(dstPOItem->getPropertyList()->isVisible());
}

void NetworkEditor::processorRenamed(const Processor* processor, const std::string& /*prevName*/) {
    if (processorItemMap_.find(const_cast<Processor*>(processor)) != processorItemMap_.end()) {
        Processor* nonConst = const_cast<Processor*>(processor);
        ProcessorGraphicsItem* item = processorItemMap_.value(nonConst);
    }
}

void NetworkEditor::portConnectionAdded(const Port* outport, const Port* inport) {
    if (!getProcessorNetwork())
        return;

    getPortGraphicsItem(outport)->addGraphicalConnection(getPortGraphicsItem(inport));

    // render size linking
    if (linkPortSizeAutoButton_->isChecked()) {
        const RenderPort* renderOutport = dynamic_cast<const RenderPort*>(outport);
        const RenderPort* renderInport = dynamic_cast<const RenderPort*>(inport);
        if (renderOutport && renderInport) {
            int numCreated = getProcessorNetwork()->createRenderSizeLinksOverConnection(
                const_cast<RenderPort*>(renderOutport), const_cast<RenderPort*>(renderInport), false);
            if (numCreated)
                LINFO("Created " << numCreated << " render size links.");
        }
    }
}

void NetworkEditor::portConnectionRemoved(const Port* outport, const Port* inport) {
    if (!getProcessorNetwork())
        return;

    getPortGraphicsItem(outport)->removeGraphicalConnection(getPortGraphicsItem(inport));

    // render size linking
    if (linkPortSizeAutoButton_->isChecked()) {
        const RenderPort* renderOutport = dynamic_cast<const RenderPort*>(outport);
        const RenderPort* renderInport = dynamic_cast<const RenderPort*>(inport);
        if (renderOutport && renderInport) {
            int numRemoved = getProcessorNetwork()->removeRenderSizeLinksOverConnection(
                const_cast<RenderPort*>(renderOutport), const_cast<RenderPort*>(renderInport));
            if (numRemoved)
                LINFO("Removed " << numRemoved << " render size links.");
        }
    }

}

//---------------------------------------------------------------------------------------------------------------
//                  scene transformations
//---------------------------------------------------------------------------------------------------------------
void NetworkEditor::scaleView() {
    if (!getProcessorNetwork())
        return;

    QRectF visibleRect;
    // scene()->itemsBoundingRect() will consider all invisible (i.e. properties) as well, so we have to to it by ourselves
    foreach (QGraphicsItem* item, scene()->items()) {
        if (item->isVisible()) {
            QRectF iRect = item->mapRectToScene(item->boundingRect());
            visibleRect = visibleRect.united(iRect);
        }
    }

    //set scale
    FloatMetaData* metaFactor = dynamic_cast<FloatMetaData*>(getProcessorNetwork()->getMetaDataContainer().getMetaData("ZoomFactor"));
    if (metaFactor) {
        QGraphicsView::scale(metaFactor->getValue(),metaFactor->getValue());
    }
    else {
        qreal scaleFactor;
        if (visibleRect.isEmpty()) {
            scaleFactor = 1.0;
            QGraphicsView::scale(scaleFactor,scaleFactor);
        } else {
            QSizeF vps = mapToScene(viewport()->rect()).boundingRect().size();
            QSizeF vs = visibleRect.size();
            scaleFactor = std::min(vps.width()/vs.width(), vps.height()/vs.height());
            QGraphicsView::scale(scaleFactor-0.1,scaleFactor-0.1);
        }
        getProcessorNetwork()->getMetaDataContainer().addMetaData("ZoomFactor", new FloatMetaData(transform().m11()));
    }

    //set center
    Vec2MetaData* metaCenter = dynamic_cast<Vec2MetaData*>(getProcessorNetwork()->getMetaDataContainer().getMetaData("ZoomCenter"));
    if (metaCenter) {
        centerOn(metaCenter->getValue().x, metaCenter->getValue().y);
    }
    else {
        QPointF center;
        // set center
        if (visibleRect.isEmpty()) {
            center = QPointF(0.0,0.0);
        } else {
            center = QPointF(visibleRect.x() + visibleRect.width() / 2.f, visibleRect.y() + visibleRect.height() / 2.f);
        }
        centerOn(center);
        getProcessorNetwork()->getMetaDataContainer().addMetaData("ZoomCenter", new Vec2MetaData(tgt::vec2(center.x(),center.y())));
    }

    invalidateScene(QRectF(), QGraphicsScene::ForegroundLayer);
}

void NetworkEditor::scale(qreal sx, qreal sy) {
    tgtAssert(sx == sy, "no rectangular zoom performed");
    QPolygonF sv = mapToScene(viewport()->rect());
    QRectF vr;
    foreach (QGraphicsItem* item, scene()->items()) {
        if (item->isVisible()) {
            QRectF iRect = item->mapRectToScene(item->boundingRect());
            vr = vr.united(iRect);
        }
    }
    //QRectF sr = scene()->itemsBoundingRect();
    vr.setCoords(vr.left()-vr.width()/2.f,vr.top()-vr.height()/2.f,vr.right()-vr.width()/2.f,vr.bottom()-vr.height()/2.f);

    // max. zoom
    if((sx > 1.f && sv.boundingRect().size().width() < 400) || (sy > 1.f && sv.boundingRect().size().height() < 600))
        return;

    // min. zoom
    if(sx < 1.f && sv.boundingRect().size().width() > 3*vr.size().width() && sv.boundingRect().size().height() > 3*vr.size().height())
        return;

    QGraphicsView::scale(sx, sy);

    if (getProcessorNetwork()) {
        if (getProcessorNetwork()->getMetaDataContainer().hasMetaData("ZoomFactor")) {
            MetaDataBase* base = getProcessorNetwork()->getMetaDataContainer().getMetaData("ZoomFactor");
            FloatMetaData* meta = dynamic_cast<FloatMetaData*>(base);
            meta->setValue(transform().m11());
        } else {
            getProcessorNetwork()->getMetaDataContainer().addMetaData("ZoomFactor", new FloatMetaData(transform().m11()));
        }
        QPointF center = mapToScene(viewport()->rect().center());
        if (getProcessorNetwork()->getMetaDataContainer().hasMetaData("ZoomCenter")) {
            MetaDataBase* base = getProcessorNetwork()->getMetaDataContainer().getMetaData("ZoomCenter");
            Vec2MetaData* meta = dynamic_cast<Vec2MetaData*>(base);
            meta->setValue(tgt::vec2(center.x(),center.y()));
        } else {
            getProcessorNetwork()->getMetaDataContainer().addMetaData("ZoomCenter", new Vec2MetaData(tgt::vec2(center.x(),center.y())));
        }
    }
}

QSize NetworkEditor::sizeHint() const {
    return QSize(400, 600);
}

//---------------------------------------------------------------------------------------------------------------
//                  create and handle graphicitems
//---------------------------------------------------------------------------------------------------------------
void NetworkEditor::resetScene() {
    emit update();
    foreach(ProcessorGraphicsItem* item, processorItemMap_) {
        delete item;
    }
    foreach(TextBoxBaseGraphicsItem* item, documentationItems_){
        delete item;
    }
    processorItemMap_.clear();
    documentationItems_.clear();

    // deletion is necessary because the QGraphicsScene's sceneRect will consider all added items
    setUpdatesEnabled(false);
    scene()->blockSignals(true);
    scene()->clear();
    delete scene();
    setScene(new QGraphicsScene(this));
    scene()->setSceneRect(sceneRectSpacing.x(), sceneRectSpacing.y(), sceneRectSpacing.width(), sceneRectSpacing.height());
    scene()->setItemIndexMethod(QGraphicsScene::NoIndex);
    connect(scene(), SIGNAL(changed(QList<QRectF>)), this, SIGNAL(networkEditor_visibleSceneSizeChanged_Signal()));
    setUpdatesEnabled(true);
}

void NetworkEditor::generateGraphicsItems() {
    if (!getProcessorNetwork())
        return;

    foreach (Processor* proc, getProcessorNetwork()->getProcessors())
        createProcessorGraphicsItem(proc);

    foreach (Processor* proc, getProcessorNetwork()->getProcessors()) {
        std::vector<Port*> outports = proc->getOutports();
        std::vector<CoProcessorPort*> coprocessoroutports = proc->getCoProcessorOutports();
        // append coprocessoroutports to outports because we can handle them identically
        outports.insert(outports.end(), coprocessoroutports.begin(), coprocessoroutports.end());

        foreach (Port* port, outports) {
            std::vector<const Port*> connectedPorts = port->getConnected();

            foreach (const Port* connectedPort, connectedPorts)
                getPortGraphicsItem(port)->addGraphicalConnection(getPortGraphicsItem(connectedPort));
        }
    }

    foreach (PropertyLink* link, getProcessorNetwork()->getPropertyLinks()) {
        createLinkArrowForPropertyLink(link);
    }
}


ProcessorGraphicsItem* NetworkEditor::createProcessorGraphicsItem(Processor* processor) {
    ProcessorGraphicsItem* result = new ProcessorGraphicsItem(processor, this);

    processorItemMap_[processor] = result;
    result->updateNWELayerAndCursor();
    result->loadMeta();

    connect(result, SIGNAL(openPropertyLinkDialog(PortOwnerGraphicsItem*, PortOwnerGraphicsItem*)), this, SLOT(openPropertyLinkDialog(PortOwnerGraphicsItem*, PortOwnerGraphicsItem*)));
    //connect(result, SIGNAL(startedArrow()), this, SLOT(disableTooltips()));
    //connect(result, SIGNAL(endedArrow()), this, SLOT(enableTooltips()));

    return result;
}

void NetworkEditor::createLinkArrowForPropertyLink(const PropertyLink* link) {
    if (link == 0)
        return;

    PortOwnerGraphicsItem* srcPOItem = 0;
    PropertyGraphicsItem* srcPropItem = 0;
    //source prop belongs to processor
    if (dynamic_cast<Processor*>(link->getSourceProperty()->getOwner())){
        Processor* srcProc = static_cast<Processor*>(link->getSourceProperty()->getOwner());

        if (processorItemMap_.contains(srcProc))
            srcPOItem = processorItemMap_[srcProc];

        tgtAssert(srcPOItem, "source item was not found");

        srcPropItem = srcPOItem->getPropertyList()->getPropertyItem(link->getSourceProperty());
    }
    //source prop belongs to port
    else if (dynamic_cast<Port*>(link->getSourceProperty()->getOwner())){
        Port* srcPort = static_cast<Port*>(link->getSourceProperty()->getOwner());

        PropertyOwnerGraphicsItem* srcPortItem = 0;
        srcPortItem = getPortGraphicsItem(srcPort);
        tgtAssert(srcPortItem, "source item was not found");

        srcPOItem = dynamic_cast<PortOwnerGraphicsItem*>(srcPortItem->parent());
        srcPropItem = srcPortItem->getPropertyList()->getPropertyItem(link->getSourceProperty());

    } else {
        LWARNING("createLinkArrowForPropertyLink(): unknown link source property owner");
        return;
    }
    tgtAssert(srcPropItem, "no src prop gi");

    PortOwnerGraphicsItem* dstPOItem = 0;
    PropertyGraphicsItem* dstPropItem = 0;
    //destination prop belongs to processor
    if (dynamic_cast<Processor*>(link->getDestinationProperty()->getOwner())){
        Processor* dstProc = static_cast<Processor*>(link->getDestinationProperty()->getOwner());

        if (processorItemMap_.contains(dstProc))
            dstPOItem = processorItemMap_[dstProc];

        tgtAssert(dstPOItem, "destination item was not found");

        dstPropItem = dstPOItem->getPropertyList()->getPropertyItem(link->getDestinationProperty());
    }
    //source prop belongs to port
    else if (dynamic_cast<Port*>(link->getDestinationProperty()->getOwner())){
        Port* dstPort = static_cast<Port*>(link->getDestinationProperty()->getOwner());

        PropertyOwnerGraphicsItem* dstPortItem = 0;
        dstPortItem = getPortGraphicsItem(dstPort);
        tgtAssert(dstPortItem, "destination item was not found");

        dstPOItem = dynamic_cast<PortOwnerGraphicsItem*>(dstPortItem->parent());
        dstPropItem = dstPortItem->getPropertyList()->getPropertyItem(link->getDestinationProperty());

    } else {
        LWARNING("createLinkArrowForPropertyLink(): unknown link destination property owner");
        return;
    }
    tgtAssert(dstPropItem, "no dst prop gi");

    srcPropItem->addGraphicalLink(dstPropItem);

    //create portsizelinkarrow if needed
    if(link->getLinkEvaluator()->getClassName() == "LinkEvaluatorRenderSize")
        dynamic_cast<PortGraphicsItem*>(srcPropItem->getPropertyOwnerItem())->addPortSizeLinkArrow(srcPropItem, dstPropItem);

    //create portownerlinkarrow if needed
    if(srcPOItem->isPortOwnerLinkNeeded(dstPOItem))
        srcPOItem->addGraphicalLinkArrow(dstPOItem);

    if(currentLayer_ == NetworkEditorLayerCameraLinking) {//test if it is still visible in camera mode
        if(srcPOItem->isPortOwnerLinkNeeded(dstPOItem, currentLayer_)) {
            foreach(PortOwnerLinkArrowGraphicsItem* arrow, srcPOItem->getPortOwnerLinkArrows())
                if(arrow->getDestinationItem()->parent() == dstPOItem)
                    arrow->setVisible(true);
            foreach(PortOwnerLinkArrowGraphicsItem* arrow, dstPOItem->getPortOwnerLinkArrows())
                if(arrow->getDestinationItem()->parent() == srcPOItem)
                    arrow->setVisible(true);
        }
    }
}

void NetworkEditor::removeItems(QList<QGraphicsItem*> items) {
    //hideTooltip();

    if (evaluator_->isLocked() && !networkEvaluatorIsLockedByButton_) {
        QMessageBox::information(this, tr("Network Locked"), tr("Network is being evaluated, so no delete operation is allowed"));
        return;
    }

    // make sure the evaluator does not operate on a temporarily inconsistent network
    evaluator_->lock();

    // sort items by their type to not delete a port/arrow-item
    // that has already been deleted indirectly with the guiitem
    // so at first kick out the ports:
    foreach (QGraphicsItem* item, items) {
        switch (item->type()) {
        case UserTypesPortGraphicsItem:
        case QGraphicsTextItem::Type:
        case UserTypesRenamableTextGraphicsItem:
        case UserTypesPropertyGraphicsItem:
            items.removeOne(item);
            break;
        case UserTypesTextBoxBaseGraphicsItem:
        case UserTypesTextBoxGraphicsItem:
        case UserTypesFrameBoxGraphicsItem:
            TextBoxBaseGraphicsItem* tItem = dynamic_cast<TextBoxBaseGraphicsItem*>(item);
            items.removeOne(item);
            removeTextBoxGraphicsItem(tItem);
            break;
        }
    }

    switch(currentLayer_){
    case NetworkEditorLayerDataFlow:
        removeItemsInDataFlow(items);
        break;
    case NetworkEditorLayerGeneralLinking:
        removeItemsInGeneralLinking(items);
        break;
    case NetworkEditorLayerCameraLinking:
        removeItemsInCameraLinking(items);
        break;
    case NetworkEditorLayerPortSizeLinking:
        removeItemsInPortSizeLinking(items);
        break;
    default:
        tgtAssert(false,"shouldn't get here");
        break;
    }

    resetCachedContent();
    scene()->clearSelection();

    // unlock evaluator (update is triggered automatically by network invalidations)
    // only unlock, if the appropriate button is not checked
    if (!networkEvaluatorIsLockedByButton_)
        evaluator_->unlock();

    if(workspace_) workspace_->setModified(true);
}


void NetworkEditor::removeItemsInDataFlow(QList<QGraphicsItem*> items){
    // next delete port arrows
    foreach (QGraphicsItem* item, items) {
        if (item->type() == UserTypesPortArrowGraphicsItem) {
            PortArrowGraphicsItem* arrow = qgraphicsitem_cast<PortArrowGraphicsItem*>(item);
            items.removeOne(item);
            removePortArrowGraphicsItem(arrow);
        }
    }

    // finally delete processor items
    foreach (QGraphicsItem* item, items) {
        if ((item->type() == UserTypesProcessorGraphicsItem) ) {
            PortOwnerGraphicsItem* poItem = static_cast<PortOwnerGraphicsItem*>(item);
            items.removeOne(item);
            removePortOwnerGraphicsItem(poItem);
        }
    }
}

void NetworkEditor::removeItemsInGeneralLinking(QList<QGraphicsItem*> items){
    foreach (QGraphicsItem* item, items) {
        switch (item->type()) {
        case UserTypesPortArrowGraphicsItem:
        case UserTypesPortSizeLinkArrowGraphicsItem:
            items.removeOne(item);
        default:
            break;
        }
    }

    //delete property link arrows
    QList<PropertyLinkArrowGraphicsItem*> linkList;
    foreach (QGraphicsItem* item, items) {
        if (item->type() == UserTypesPropertyLinkArrowGraphicsItem) {
            linkList.append(qgraphicsitem_cast<PropertyLinkArrowGraphicsItem*>(item));
            items.removeOne(item);
        }
    }
    while(!linkList.empty()){
        PropertyLinkArrowGraphicsItem* item = linkList.first();
        linkList.removeOne(item);
        //remove dual link if exists
        foreach (PropertyLinkArrowGraphicsItem* dual, linkList) {
            if(dual->getDestinationItem() == item->getSourceItem() && dual->getSourceItem() == item->getDestinationItem()) {
                linkList.removeOne(dual);
            } else
            if(dual->getSourceItem() == item->getSourceItem() && dual->getDestinationItem() == item->getDestinationItem()) {
                linkList.removeOne(dual);
            }
        }
        removePropertyLinkArrowGraphicsItem(item);
    }

    //delete links between processors
    QList<ProcessorGraphicsItem*> processorList;
    foreach (QGraphicsItem* item, items) {
        if (item->type() == UserTypesProcessorGraphicsItem) {
            processorList.append(qgraphicsitem_cast<ProcessorGraphicsItem*>(item));
            items.removeOne(item);
        }
    }
    if(processorList.size() == 1){
        foreach(PropertyGraphicsItem* prop, processorList.first()->getPropertyList()->getAllPropertyItems()){
            foreach(PropertyLinkArrowGraphicsItem* arrow, prop->getSourceLinkList())
                removePropertyLinkArrowGraphicsItem(arrow);
            foreach(PropertyLinkArrowGraphicsItem* arrow, prop->getDestinationLinkList())
                removePropertyLinkArrowGraphicsItem(arrow);
        }
    } else {
        foreach(ProcessorGraphicsItem* proc1, processorList){
            foreach(ProcessorGraphicsItem* proc2, processorList){
                foreach(PropertyGraphicsItem* prop1, proc1->getPropertyList()->getAllPropertyItems()){
                    foreach(PropertyLinkArrowGraphicsItem* arrow, prop1->getSourceLinkList()){
                        foreach(PropertyGraphicsItem* prop2, proc2->getPropertyList()->getAllPropertyItems()){
                            if(arrow->getDestinationItem() == prop2){
                                removePropertyLinkArrowGraphicsItem(arrow);
                                break;
                            }
                        }
                    }
                }
            }
        }
    }
}

void NetworkEditor::removeItemsInCameraLinking(QList<QGraphicsItem*> items) {
    foreach (QGraphicsItem* item, items) {
        switch (item->type()) {
        case UserTypesPortArrowGraphicsItem:
        case UserTypesPortSizeLinkArrowGraphicsItem:
            items.removeOne(item);
        default:
            break;
        }
    }

    //delete property link arrows
    QList<PropertyLinkArrowGraphicsItem*> linkList;
    foreach (QGraphicsItem* item, items) {
        if (item->type() == UserTypesPropertyLinkArrowGraphicsItem) {
            linkList.append(qgraphicsitem_cast<PropertyLinkArrowGraphicsItem*>(item));
            items.removeOne(item);
        }
    }
    while(!linkList.empty()){
        PropertyLinkArrowGraphicsItem* item = linkList.first();
        linkList.removeOne(item);
        //remove dual link if exists
        foreach (PropertyLinkArrowGraphicsItem* dual, linkList) {
            if(dual->getSourceItem() == item->getSourceItem() && dual->getDestinationItem() == item->getDestinationItem()) {
                linkList.removeOne(dual);
            } else
            if(dual->getDestinationItem() == item->getSourceItem() && dual->getSourceItem() == item->getDestinationItem()) {
                linkList.removeOne(dual);
            }
        }
        removePropertyLinkArrowGraphicsItem(item);
    }

    //delete links between processors
    QList<ProcessorGraphicsItem*> processorList;
    foreach (QGraphicsItem* item, items) {
        if (item->type() == UserTypesProcessorGraphicsItem) {
            processorList.append(qgraphicsitem_cast<ProcessorGraphicsItem*>(item));
            items.removeOne(item);
        }
    }
    if(processorList.size() == 1){
        foreach(PropertyGraphicsItem* prop, processorList.first()->getPropertyList()->getAllPropertyItems()){
            foreach(PropertyLinkArrowGraphicsItem* arrow, prop->getSourceLinkList())
                if(arrow->getSourceItem()->getProperty()->getClassName() == "CameraProperty" &&
                   arrow->getDestinationItem()->getProperty()->getClassName() == "CameraProperty")
                    removePropertyLinkArrowGraphicsItem(arrow);
            foreach(PropertyLinkArrowGraphicsItem* arrow, prop->getDestinationLinkList())
                if(arrow->getSourceItem()->getProperty()->getClassName() == "CameraProperty" &&
                   arrow->getDestinationItem()->getProperty()->getClassName() == "CameraProperty")
                    removePropertyLinkArrowGraphicsItem(arrow);
        }
    } else {
        foreach(ProcessorGraphicsItem* proc1, processorList){
            foreach(ProcessorGraphicsItem* proc2, processorList){
                foreach(PropertyGraphicsItem* prop1, proc1->getPropertyList()->getAllPropertyItems()){
                    if(prop1->getProperty()->getClassName() == "CameraProperty"){
                        foreach(PropertyLinkArrowGraphicsItem* arrow, prop1->getSourceLinkList()){
                            foreach(PropertyGraphicsItem* prop2, proc2->getPropertyList()->getAllPropertyItems()){
                                if(prop2->getProperty()->getClassName() == "CameraProperty" &&
                                    arrow->getDestinationItem() == prop2){
                                    removePropertyLinkArrowGraphicsItem(arrow);
                                    break;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

void NetworkEditor::removeItemsInPortSizeLinking(QList<QGraphicsItem*> items) {
    foreach (QGraphicsItem* item, items) {
        switch (item->type()) {
        case UserTypesPortArrowGraphicsItem:
            items.removeOne(item);
        default:
            break;
        }
    }

    //delete port size link arrows
    foreach (QGraphicsItem* item, items) {
        if (item->type() == UserTypesPortSizeLinkArrowGraphicsItem) {
            PortSizeLinkArrowGraphicsItem* arrow = qgraphicsitem_cast<PortSizeLinkArrowGraphicsItem*>(item);
            items.removeOne(item);
            removePortSizeLinkArrowGraphicsItem(arrow);
        }
    }

    //delete links between processors
    QList<ProcessorGraphicsItem*> processorList;
    foreach (QGraphicsItem* item, items) {
        if (item->type() == UserTypesProcessorGraphicsItem) {
            processorList.append(qgraphicsitem_cast<ProcessorGraphicsItem*>(item));
            items.removeOne(item);
        }
    }
    if(processorList.size() == 1){
        foreach(PropertyGraphicsItem* prop, processorList.first()->getPropertyList()->getAllPropertyItems()){
            foreach(PropertyLinkArrowGraphicsItem* arrow, prop->getSourceLinkList())
                if(arrow->getSourceItem()->getProperty()->getLink(arrow->getDestinationItem()->getProperty())
                   ->getLinkEvaluator()->getClassName() == "LinkEvaluatorRenderSize")
                   removePropertyLinkArrowGraphicsItem(arrow);
            foreach(PropertyLinkArrowGraphicsItem* arrow, prop->getDestinationLinkList())
                if(arrow->getSourceItem()->getProperty()->getLink(arrow->getDestinationItem()->getProperty())
                   ->getLinkEvaluator()->getClassName() == "LinkEvaluatorRenderSize")
                    removePropertyLinkArrowGraphicsItem(arrow);
        }
    } else {
        foreach(ProcessorGraphicsItem* proc1, processorList){
            foreach(ProcessorGraphicsItem* proc2, processorList){
                foreach(PropertyGraphicsItem* prop1, proc1->getPropertyList()->getAllPropertyItems()){
                    if(prop1->getProperty()->getClassName() == "RenderSizeOriginProperty" ||
                       prop1->getProperty()->getClassName() == "RenderSizeReceiveProperty"){
                        foreach(PropertyLinkArrowGraphicsItem* arrow, prop1->getSourceLinkList()){
                            foreach(PropertyGraphicsItem* prop2, proc2->getPropertyList()->getAllPropertyItems()){
                                if(arrow->getDestinationItem() == prop2 &&
                                  (prop2->getProperty()->getClassName() == "RenderSizeOriginProperty" ||
                                   prop2->getProperty()->getClassName() == "RenderSizeReceiveProperty")){
                                    removePropertyLinkArrowGraphicsItem(arrow);
                                    break;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

void NetworkEditor::removePortOwnerGraphicsItem(PortOwnerGraphicsItem* poItem) {
    if (!getProcessorNetwork())
        return;

    if (poItem->type() == UserTypesProcessorGraphicsItem)
        foreach (Processor* processor, poItem->getProcessors())
            getProcessorNetwork()->removeProcessor(processor);
}

void NetworkEditor::removePortArrowGraphicsItem(PortArrowGraphicsItem* arrow) {
    if (!getProcessorNetwork())
        return;

    if (arrow->getDestinationItem() != 0) {
        getProcessorNetwork()->disconnectPorts(arrow->getSourceItem()->getPort(), arrow->getDestinationItem()->getPort());
    }
}

void NetworkEditor::removePropertyLinkArrowGraphicsItem(PropertyLinkArrowGraphicsItem* arrow) {
    if (!getProcessorNetwork())
        return;

    if (arrow->getDestinationItem() != 0) {
        const Property* src = arrow->getSourceItem()->getProperty();
        const Property* dst = arrow->getDestinationItem()->getProperty();
        getProcessorNetwork()->removePropertyLink(src->getLink(dst));
        //remove other link if exists
        if(dst->getLink(src))
            getProcessorNetwork()->removePropertyLink(dst->getLink(src));
    }
}

void NetworkEditor::removePortSizeLinkArrowGraphicsItem(PortSizeLinkArrowGraphicsItem* arrow) {
    if (!getProcessorNetwork())
        return;

    if (arrow->getDestinationItem()) {
        const Property* src = arrow->getSourceItem()->getProperty();
        const Property* dst = arrow->getDestinationItem()->getProperty();
        getProcessorNetwork()->removePropertyLink(src->getLink(dst));
        //remove other link if exists
        if(dst->getLink(src))
            getProcessorNetwork()->removePropertyLink(dst->getLink(src));
    }
}

void NetworkEditor::removePortOwnerLinkArrowGraphicsItem(PortOwnerLinkArrowGraphicsItem* arrow) {
    if (!getProcessorNetwork())
        return;

    if(arrow->getDestinationItem()) {
        LERROR("Not Implemented Yet!!!!!");
    }
}

void NetworkEditor::removeTextBoxGraphicsItem(TextBoxBaseGraphicsItem* item) {
    documentationItems_.removeOne(item);
    scene()->removeItem(item);
    delete item;
}

ProcessorGraphicsItem* NetworkEditor::getProcessorGraphicsItem(const Processor* processor) const {
    tgtAssert(processor, "null pointer passed");
    return processorItemMap_.find(const_cast<Processor*>(processor)).value();
}

PortGraphicsItem* NetworkEditor::getPortGraphicsItem(const Port* port) const {
    tgtAssert(port, "null pointer passed");

    PortGraphicsItem* portItem = 0;
    if (processorItemMap_.contains(port->getProcessor())) {
        ProcessorGraphicsItem* outportProcItem = processorItemMap_[port->getProcessor()];
        portItem = outportProcItem->getPortGraphicsItem(port);
    }

    tgtAssert(portItem, "no portgraphicsitem found");
    return portItem;
}

PropertyGraphicsItem* NetworkEditor::getPropertyGraphicsItem(const Property* prop) const {
    tgtAssert(prop, "null pointer passed");

    PropertyGraphicsItem* propItem = 0;
    foreach(ProcessorGraphicsItem* item, processorItemMap_) {
        if(item->getPropertyList()->getPropertyItem(prop))
            return item->getPropertyList()->getPropertyItem(prop);
        foreach(PortGraphicsItem* portItem, item->getPortGraphicsItems()){
            if(portItem->getPropertyList()->getPropertyItem(prop))
                return portItem->getPropertyList()->getPropertyItem(prop);
        }
    }

    tgtAssert(propItem, "no portgraphicsitem found");
    return propItem;
}

void NetworkEditor::updateSelectedItems() {
    if (!getProcessorNetwork())
        return;

    QList<QGraphicsItem*> selectedItems = scene()->selectedItems();

    QList<Processor*> selectedProcessors;
    foreach (QGraphicsItem* selectedItem, selectedItems) {
        if (selectedItem->type() == UserTypesProcessorGraphicsItem) {
            PortOwnerGraphicsItem* portOwnerItem = static_cast<PortOwnerGraphicsItem*>(selectedItem);
            selectedProcessors += portOwnerItem->getProcessors();
            portOwnerItem->saveMeta();
        }
    }

    if (selectedProcessors.isEmpty())
        getProcessorNetwork()->getMetaDataContainer().removeMetaData("ProcessorSelection");
    else {
        std::vector<Processor*> selectedProc = qListToStdVector(selectedProcessors);
        getProcessorNetwork()->getMetaDataContainer().addMetaData("ProcessorSelection", new SerializableVectorMetaData<Processor*>(selectedProc));
    }

    emit processorsSelected(selectedProcessors);
}

void NetworkEditor::selectPreviouslySelectedProcessors() {
    if (!getProcessorNetwork())
        return;

    SerializableVectorMetaData<Processor*>* selectionMetaData = dynamic_cast<SerializableVectorMetaData<Processor*>*>(getProcessorNetwork()->getMetaDataContainer().getMetaData("ProcessorSelection"));
    if (selectionMetaData) {
        std::vector<Processor*> selectionVector = selectionMetaData->getValues();
        QList<Processor*> selectedProcessors;
        for (size_t i=0; i<selectionVector.size(); i++) {
            Processor* proc = selectionVector.at(i);
            if (!proc)
                continue;

            selectedProcessors.push_back(proc);

            ProcessorGraphicsItem* item = getProcessorGraphicsItem(proc);

            // item might not exist (e.g. network was build with a voreen application using different defines)
            if (item)
                item->setSelected(true);
        }

        emit processorsSelected(selectedProcessors);
    }
}

QList<PortOwnerGraphicsItem*> NetworkEditor::getSelectedPortOwnerGraphicsItems() {
    QList<PortOwnerGraphicsItem*> portOwners;
    foreach(QGraphicsItem* item, scene()->selectedItems()) {
        if(item->type() == UserTypesProcessorGraphicsItem)
            portOwners.append(dynamic_cast<PortOwnerGraphicsItem*>(item));
    }
    return portOwners;
}

QList<ProcessorGraphicsItem*> NetworkEditor::getSelectedProcessorGraphicsItems() {
    QList<ProcessorGraphicsItem*> processors;
    foreach(QGraphicsItem* item, scene()->selectedItems()) {
        if(item->type() == UserTypesProcessorGraphicsItem )
            processors.append(dynamic_cast<ProcessorGraphicsItem*>(item));
    }
    return processors;
}

void NetworkEditor::createNewTextNoteSlot() {
    //activate documentation if deactive
    if(!documentationButton_->isChecked())
        documentationButton_->click();
    TextBoxGraphicsItem* item = new TextBoxGraphicsItem(this);
    item->setPos(mapToScene(nextTextBoxPos_));
    documentationItems_.push_back(item);
    scene()->addItem(item);
}

void NetworkEditor::createNewTextFrameSlot() {
    //activate documentation if deactive
    if(!documentationButton_->isChecked())
        documentationButton_->click();
    FrameBoxGraphicsItem* item = new FrameBoxGraphicsItem(this);
    item->setPos(mapToScene(nextTextBoxPos_));
    documentationItems_.push_back(item);
    scene()->addItem(item);
}

//---------------------------------------------------------------------------------------------------------------
//                  style, layer and cursor management
//---------------------------------------------------------------------------------------------------------------
NetworkEditorLayer NetworkEditor::getCurrentLayer() const {
    return currentLayer_;
}

NetworkEditorCursorMode NetworkEditor::getCurrentCursorMode() const {
    return currentCursorMode_;
}

NWEStyle_Base* NetworkEditor::getCurrentStyle() const {
    return currentStyle_;
}

bool NetworkEditor::getCurrentToolTipMode() const {
    return currentToolTipMode_;
}

void NetworkEditor::setLayer(NetworkEditorLayer layer) {
    if (layer == currentLayer_)
        return;

    switch(layer){
    case NetworkEditorLayerDataFlow:
        linkingLayerButtonContainer_->setVisible(false);
        generalLinkingLayerButtonContainer_->setVisible(false);
        cameraLinkingLayerButtonContainer_->setVisible(false);
        portSizeLinkingLayerButtonContainer_->setVisible(false);
        break;
    case NetworkEditorLayerLinking:
        if(generalLinkingButton_->isChecked())
            setLayer(NetworkEditorLayerGeneralLinking);
        else if(cameraLinkingButton_->isChecked())
            setLayer(NetworkEditorLayerCameraLinking);
        else if(portSizeLinkingButton_->isChecked())
            setLayer(NetworkEditorLayerPortSizeLinking);
        else
            tgtAssert(false,"no button checked");
        return;
        break;
    case NetworkEditorLayerGeneralLinking:
        linkingLayerButtonContainer_->setVisible(true);
        generalLinkingLayerButtonContainer_->setVisible(true);
        cameraLinkingLayerButtonContainer_->setVisible(false);
        portSizeLinkingLayerButtonContainer_->setVisible(false);
        break;
    case NetworkEditorLayerCameraLinking:
        linkingLayerButtonContainer_->setVisible(true);
        generalLinkingLayerButtonContainer_->setVisible(false);
        cameraLinkingLayerButtonContainer_->setVisible(true);
        portSizeLinkingLayerButtonContainer_->setVisible(false);
        break;
    case NetworkEditorLayerPortSizeLinking:
        linkingLayerButtonContainer_->setVisible(true);
        generalLinkingLayerButtonContainer_->setVisible(false);
        cameraLinkingLayerButtonContainer_->setVisible(false);
        portSizeLinkingLayerButtonContainer_->setVisible(true);
        break;
    default:
        tgtAssert(false,"Unknown Layer");
        break;
    }

    currentLayer_ = layer;

    //set layer for each processor
    foreach (ProcessorGraphicsItem* procItem, processorItemMap_.values()) {
        if (procItem)
            procItem->updateNWELayerAndCursor();
    }

    resetCachedContent(); //needed?
}

void NetworkEditor::setCursorMode(NetworkEditorCursorMode mode) {
    if (currentCursorMode_ == mode)
        return;

    currentCursorMode_ = mode;

    //set cursor mode for each processor
    foreach (ProcessorGraphicsItem* procItem, processorItemMap_.values()) {
        if (procItem)
            procItem->updateNWELayerAndCursor();
    }

    switch(currentCursorMode_){
    case NetworkEditorCursorSelectMode:
        setCursor(Qt::ArrowCursor);
        break;
    case NetworkEditorCursorMoveMode:
        setCursor(Qt::OpenHandCursor);
        scene()->clearSelection();
        emit processorsSelected(QList<Processor*>());
        break;
    default:
        tgtAssert(false, "should not get here");
        break;
    }
}

void NetworkEditor::processorFontOnChange() {
    NWEBaseGraphicsItem* base = 0;
    foreach(QGraphicsItem* item,scene()->items()){
        if((base = dynamic_cast<NWEBaseGraphicsItem*>(item))) {
            base->resetPaintInitialization();
            base->update();
        }
    }
    scene()->update();
}

void NetworkEditor::updateStyle() {
    switch(networkEditorStyleProperty_.getValue()){
    case NWESTYLE_CLASSIC:
        setStyle(new NWEStyle_Classic(this));
        break;
    case NWESTYLE_CLASSIC_PRINT:
        setStyle(new NWEStyle_Classic_Print(this));
        break;
    case NWESTYLE_MATERIAL:
        setStyle(new NWEStyle_Material(this));
        break;
    default:
        tgtAssert(false,"Unknown NetworkEditorStyle. Style not changed.");
        LERROR("Unknown NetworkEditorStyle. Style not changed.");
        break;
    }
}

void NetworkEditor::updateGraphLayout() {
  switch(networkEditorGraphLayoutsProperty_.getValue()){
    case NWEGL_SUGIYAMA: {
        qreal shift = 300.f; bool overlap = false, median = true, portflush = true;
        shift = sugiShiftXProperty_.get();
        overlap = sugiOverlapProperty_.get();
        median = sugiMedianProperty_.get();
        portflush = sugiPortFlushProperty_.get();
        static_cast<NWEGL_Sugiyama*>(currentGraphLayout_)->setSortParameter(shift,overlap,median,portflush);
        } break;
    default:
        tgtAssert(false,"Unknown NetworkEditorGraphLayout!!!");
        LERROR("Unknown NetworkEditorGraphLayout!!!");
        break;
    }
}

void NetworkEditor::setStyle(NWEStyle_Base* style) {
    if(!currentStyle_) {
        currentStyle_ = style;
        return;
    }

    delete currentStyle_;
    currentStyle_ = style;
    NWEBaseGraphicsItem* base = 0;
    foreach(QGraphicsItem* item,scene()->items()){
        if((base = dynamic_cast<NWEBaseGraphicsItem*>(item))) {
            base->resetPaintInitialization();
            base->update();
        }
    }
    //set background color
    setBackgroundBrush(currentStyle_->getBackgroundBrush());
    //set button background
    QPalette pal = generalLinkingLayerButtonContainer_->palette();
    pal.setColor(generalLinkingLayerButtonContainer_->backgroundRole(), currentStyle_->getButtonBackgroundColor());
    generalLinkingLayerButtonContainer_->setPalette(pal);
    pal = cameraLinkingLayerButtonContainer_->palette();
    pal.setColor(cameraLinkingLayerButtonContainer_->backgroundRole(), currentStyle_->getButtonBackgroundColor());
    cameraLinkingLayerButtonContainer_->setPalette(pal);
    pal = portSizeLinkingLayerButtonContainer_->palette();
    pal.setColor(portSizeLinkingLayerButtonContainer_->backgroundRole(), currentStyle_->getButtonBackgroundColor());
    portSizeLinkingLayerButtonContainer_->setPalette(pal);
    pal = linkingLayerButtonContainer_->palette();
    pal.setColor(linkingLayerButtonContainer_->backgroundRole(), currentStyle_->getButtonBackgroundColor());
    linkingLayerButtonContainer_->setPalette(pal);

    update();
}

void NetworkEditor::setToolTipMode(bool mode) {
    if(currentToolTipMode_ == mode) return;

    currentToolTipMode_ = mode;
}

void NetworkEditor::toggleToolTip() {
    currentToolTipMode_ = !currentToolTipMode_;
}

//---------------------------------------------------------------------------------------------------------------
//                  button management
//---------------------------------------------------------------------------------------------------------------
void NetworkEditor::initilizeEditorButtons() {
//general linking
    generalLinkingLayerButtonContainer_ = new QWidget(this);
    generalLinkingLayerButtonContainer_->setAutoFillBackground(true);
    QPalette pal = generalLinkingLayerButtonContainer_->palette();
    pal.setColor(generalLinkingLayerButtonContainer_->backgroundRole(), currentStyle_->getButtonBackgroundColor());
    generalLinkingLayerButtonContainer_->setPalette(pal);
    QBoxLayout* generalLinkingLayerButtonLayout = new QHBoxLayout(generalLinkingLayerButtonContainer_);
    generalLinkingLayerButtonLayout->setMargin(NWEButtonBackgroundMargin);
    //create buttons
    hideCameraLinksButton_ = new QToolButton;
        hideCameraLinksButton_->setIcon(QIcon(":/qt/icons/hide-linking-camera.png"));
        hideCameraLinksButton_->setIconSize(NWESub2ButtonSize);
        hideCameraLinksButton_->setToolTip(tr("hides all camera links"));
        hideCameraLinksButton_->setCheckable(true);
        hideCameraLinksButton_->setChecked(false);
        connect(hideCameraLinksButton_, SIGNAL(clicked()), this, SLOT(hideCameraLinks()));
        generalLinkingLayerButtonLayout->addWidget(hideCameraLinksButton_);
    hidePortSizeLinksButton_ = new QToolButton;
        hidePortSizeLinksButton_->setIcon(QIcon(":/qt/icons/hide-linking-port.png"));
        hidePortSizeLinksButton_->setIconSize(NWESub2ButtonSize);
        hidePortSizeLinksButton_->setToolTip(tr("hides all port size links"));
        hidePortSizeLinksButton_->setCheckable(true);
        hidePortSizeLinksButton_->setChecked(true);
        connect(hidePortSizeLinksButton_, SIGNAL(clicked()), this, SLOT(hidePortSizeLinks()));
        generalLinkingLayerButtonLayout->addWidget(hidePortSizeLinksButton_);
    removeAllPropertyLinksButton_ = new QToolButton;
        removeAllPropertyLinksButton_->setIcon(QIcon(":/qt/icons/linking-remove.png"));
        removeAllPropertyLinksButton_->setIconSize(NWESub2ButtonSize);
        removeAllPropertyLinksButton_->setToolTip(tr("remove all property links from the network"));
        connect(removeAllPropertyLinksButton_, SIGNAL(clicked()), this, SLOT(removeAllPropertyLinks()));
        generalLinkingLayerButtonLayout->addWidget(removeAllPropertyLinksButton_);

//camera linking
    cameraLinkingLayerButtonContainer_ = new QWidget(this);
    cameraLinkingLayerButtonContainer_->setAutoFillBackground(true);
    pal = cameraLinkingLayerButtonContainer_->palette();
    pal.setColor(cameraLinkingLayerButtonContainer_->backgroundRole(), currentStyle_->getButtonBackgroundColor());
    cameraLinkingLayerButtonContainer_->setPalette(pal);
    QBoxLayout* cameraLinkingLayerButtonLayout = new QHBoxLayout(cameraLinkingLayerButtonContainer_);
    cameraLinkingLayerButtonLayout->setMargin(NWEButtonBackgroundMargin);
    //create buttons
    linkCamerasAutoButton_ = new QToolButton;
        linkCamerasAutoButton_->setIcon(QIcon(":/qt/icons/linking-camera-auto.png"));
        linkCamerasAutoButton_->setIconSize(NWESub2ButtonSize);
        linkCamerasAutoButton_->setToolTip(tr("link cameras of processors when they are added to the network"));
        linkCamerasAutoButton_->setCheckable(true);
        linkCamerasAutoButton_->setChecked(true);
        connect(linkCamerasAutoButton_, SIGNAL(clicked()), this, SLOT(linkCamerasAutoChanged()));
        cameraLinkingLayerButtonLayout->addWidget(linkCamerasAutoButton_);
    linkCamerasButton_ = new QToolButton;
        linkCamerasButton_->setIcon(QIcon(":/qt/icons/linking-camera.png"));
        linkCamerasButton_->setIconSize(NWESub2ButtonSize);
        linkCamerasButton_->setToolTip(tr("link all cameras in the network"));
        connect(linkCamerasButton_, SIGNAL(clicked()), this, SLOT(linkCameras()));
        cameraLinkingLayerButtonLayout->addWidget(linkCamerasButton_);
    removeAllCameraLinksButton_ = new QToolButton;
        removeAllCameraLinksButton_->setIcon(QIcon(":/qt/icons/linking-remove.png"));
        removeAllCameraLinksButton_->setIconSize(NWESub2ButtonSize);
        removeAllCameraLinksButton_->setToolTip(tr("remove all camera property links from the network"));
        connect(removeAllCameraLinksButton_, SIGNAL(clicked()), this, SLOT(removeAllCameraLinks()));
        cameraLinkingLayerButtonLayout->addWidget(removeAllCameraLinksButton_);

//port size linking
    portSizeLinkingLayerButtonContainer_ = new QWidget(this);
    portSizeLinkingLayerButtonContainer_->setAutoFillBackground(true);
    pal = portSizeLinkingLayerButtonContainer_->palette();
    pal.setColor(portSizeLinkingLayerButtonContainer_->backgroundRole(), currentStyle_->getButtonBackgroundColor());
    portSizeLinkingLayerButtonContainer_->setPalette(pal);
    QBoxLayout* portSizeLinkingLayerButtonLayout = new QHBoxLayout(portSizeLinkingLayerButtonContainer_);
    portSizeLinkingLayerButtonLayout->setMargin(NWEButtonBackgroundMargin);
    //create buttons
    linkPortSizeAutoButton_ = new QToolButton;
        linkPortSizeAutoButton_->setIcon(QIcon(":/qt/icons/linking-port-auto.png"));
        linkPortSizeAutoButton_->setIconSize(NWESub2ButtonSize);
        linkPortSizeAutoButton_->setToolTip(tr("link render port sizes of processors when they are added to the network"));
        linkPortSizeAutoButton_->setCheckable(true);
        linkPortSizeAutoButton_->setChecked(true);
        connect(linkPortSizeAutoButton_, SIGNAL(clicked()), this, SLOT(linkPortSizeAutoChanged()));
        portSizeLinkingLayerButtonLayout->addWidget(linkPortSizeAutoButton_);
    linkPortSizeButton_ = new QToolButton;
        linkPortSizeButton_->setIcon(QIcon(":/qt/icons/linking-port.png"));
        linkPortSizeButton_->setIconSize(NWESub2ButtonSize);
        linkPortSizeButton_->setToolTip(tr("link all cameras in the network"));
        connect(linkPortSizeButton_, SIGNAL(clicked()), this, SLOT(linkPortSize()));
        portSizeLinkingLayerButtonLayout->addWidget(linkPortSizeButton_);
    removeAllPortSizeLinksButton_ = new QToolButton;
        removeAllPortSizeLinksButton_->setIcon(QIcon(":/qt/icons/linking-remove.png"));
        removeAllPortSizeLinksButton_->setIconSize(NWESub2ButtonSize);
        removeAllPortSizeLinksButton_->setToolTip(tr("remove all port size property links from the network"));
        connect(removeAllPortSizeLinksButton_, SIGNAL(clicked()), this, SLOT(removeAllPortSizeLinks()));
        portSizeLinkingLayerButtonLayout->addWidget(removeAllPortSizeLinksButton_);

//linking layers
    linkingLayerButtonContainer_ = new QWidget(this);
    linkingLayerButtonContainer_->setAutoFillBackground(true);
    pal = linkingLayerButtonContainer_->palette();
    pal.setColor(linkingLayerButtonContainer_->backgroundRole(), currentStyle_->getButtonBackgroundColor());
    linkingLayerButtonContainer_->setPalette(pal);
    linkingLayerButtonContainer_->setMinimumSize(NWEButtonBackgroundMargin*2+NWEMainButtonSize.width(),1);

    QBoxLayout* linkingLayerButtonLayout = new QVBoxLayout(linkingLayerButtonContainer_);
    linkingLayerButtonLayout->setMargin(NWEButtonBackgroundMargin);
    linkingLayerButtonLayout->addSpacing(NWEMainButtonSize.height()+NWEButtonBackgroundMargin+NWEMarginLayerToLinking);
    //create buttons
    generalLinkingButton_ = new QToolButton;
        generalLinkingButton_->setIcon(QIcon(":/qt/icons/linking-general.png"));
        generalLinkingButton_->setIconSize(NWESub1ButtonSize);
        generalLinkingButton_->setToolTip(tr("shows all property links"));
        generalLinkingButton_->setCheckable(true);
        generalLinkingButton_->setChecked(true);
        connect(generalLinkingButton_, SIGNAL(clicked()), this, SLOT(setLayerToGeneralLinking()));
        linkingLayerButtonLayout->addWidget(generalLinkingButton_);
    cameraLinkingButton_ = new QToolButton;
        cameraLinkingButton_->setIcon(QIcon(":/qt/icons/linking-camera.png"));
        cameraLinkingButton_->setIconSize(NWESub1ButtonSize);
        cameraLinkingButton_->setToolTip(tr("shows all camera property links"));
        cameraLinkingButton_->setCheckable(true);
        connect(cameraLinkingButton_, SIGNAL(clicked()), this, SLOT(setLayerToCameraLinking()));
        linkingLayerButtonLayout->addWidget(cameraLinkingButton_);
    portSizeLinkingButton_ = new QToolButton;
        portSizeLinkingButton_->setIcon(QIcon(":/qt/icons/linking-port.png"));
        portSizeLinkingButton_->setIconSize(NWESub1ButtonSize);
        portSizeLinkingButton_->setToolTip(tr("shows all port size property links"));
        portSizeLinkingButton_->setCheckable(true);
        connect(portSizeLinkingButton_, SIGNAL(clicked()), this, SLOT(setLayerToPortSizeLinking()));
        linkingLayerButtonLayout->addWidget(portSizeLinkingButton_);
    // add to button group, so only one can be checked at the same time
    QButtonGroup* linkingLayerButtonGroup = new QButtonGroup(this);
        linkingLayerButtonGroup->addButton(generalLinkingButton_);
        linkingLayerButtonGroup->addButton(cameraLinkingButton_);
        linkingLayerButtonGroup->addButton(portSizeLinkingButton_);

//layer buttons
    layerButtonContainer_ = new QWidget(this);
    QBoxLayout* layerButtonLayout = new QHBoxLayout(layerButtonContainer_);
    layerButtonLayout->setMargin(NWEButtonBackgroundMargin);
    //create buttons
    dataFlowLayerButton_ = new QToolButton;
        dataFlowLayerButton_->setIcon(QIcon(":/qt/icons/dataflow-mode.png"));
        dataFlowLayerButton_->setIconSize(NWEMainButtonSize);
        dataFlowLayerButton_->setToolTip(tr("Switch to data flow mode (ctrl+1)"));
        dataFlowLayerButton_->setCheckable(true);
        dataFlowLayerButton_->setShortcut(Qt::CTRL + Qt::Key_1);
        connect(dataFlowLayerButton_, SIGNAL(clicked()), this, SLOT(setLayerToDataFlow()));
        layerButtonLayout->addWidget(dataFlowLayerButton_);
    linkingLayerButton_ = new QToolButton;
        linkingLayerButton_->setIcon(QIcon(":/qt/icons/linking-mode.png"));
        linkingLayerButton_->setIconSize(NWEMainButtonSize);
        linkingLayerButton_->setToolTip(tr("Switch to linking mode (ctrl+2)"));
        linkingLayerButton_->setCheckable(true);
        linkingLayerButton_->setShortcut(Qt::CTRL + Qt::Key_2);
        connect(linkingLayerButton_, SIGNAL(clicked()), this, SLOT(setLayerToLinking()));
        layerButtonLayout->addWidget(linkingLayerButton_);
    // add to button group, so only one can be checked at the same time
    QButtonGroup* layerButtonGroup = new QButtonGroup(this);
        layerButtonGroup->addButton(dataFlowLayerButton_);
        layerButtonGroup->addButton(linkingLayerButton_);

//stop network button
    stopButtonContainer_ = new QWidget(this);
    QBoxLayout* stopButtonLayout = new QHBoxLayout(stopButtonContainer_);
    stopButtonLayout->setMargin(NWEButtonBackgroundMargin);
    stopNetworkEvaluatorButton_ = new QToolButton;
    stopNetworkEvaluatorButton_->setIcon(QIcon(":/qt/icons/player-pause.png"));
    stopNetworkEvaluatorButton_->setIconSize(NWEMainButtonSize);
    stopNetworkEvaluatorButton_->setToolTip(tr("Stop the automatic evaluation of the network"));
    connect(stopNetworkEvaluatorButton_, SIGNAL(clicked()), this, SLOT(toggleNetworkEvaluator()));
    stopNetworkEvaluatorButton_->setCheckable(true);
    stopButtonLayout->addWidget(stopNetworkEvaluatorButton_);

//navigation
    navigationButtonContainer_ = new QWidget(this);
    QHBoxLayout* navigationButtonContainerLayout = new QHBoxLayout(navigationButtonContainer_);
    navigationButtonContainerLayout->setMargin(NWEButtonBackgroundMargin);
    //create buttons
    selectCursorButton_ = new QToolButton;
        selectCursorButton_->setIcon(QIcon(":/qt/icons/cursor_arrow.svg"));
        selectCursorButton_->setIconSize(NWEMainButtonSize);
        selectCursorButton_->setToolTip(tr("Select network items"));
        selectCursorButton_->setCheckable(true);
        selectCursorButton_->setChecked(true);
        connect(selectCursorButton_, SIGNAL(clicked()), this, SLOT(setCursorSelect()));
        navigationButtonContainerLayout->addWidget(selectCursorButton_);
    moveCursorButton_ = new QToolButton;
        moveCursorButton_->setIcon(QIcon(":/qt/icons/cursor_hand.svg"));
        moveCursorButton_->setIconSize(NWEMainButtonSize);
        moveCursorButton_->setToolTip(tr("Navigate in network (shift)"));
        moveCursorButton_->setCheckable(true);
        connect(moveCursorButton_, SIGNAL(clicked()), this, SLOT(setCursorMove()));
        navigationButtonContainerLayout->addWidget(moveCursorButton_);
    QButtonGroup* navigationButtonGroup = new QButtonGroup(this);
        navigationButtonGroup->addButton(selectCursorButton_);
        navigationButtonGroup->addButton(moveCursorButton_);

//layout
    layoutButtonContainer_ = new QWidget(this);
    QHBoxLayout* layoutButtonContainerLayout = new QHBoxLayout(layoutButtonContainer_);
    layoutButtonContainerLayout->setMargin(NWEButtonBackgroundMargin);

    centerViewButton_ = new QToolButton;
        centerViewButton_->setIcon(QIcon(":/qt/icons/center.png"));
        centerViewButton_->setIconSize(NWEMainButtonSize);
        centerViewButton_->setToolTip(tr("Overview entire network"));
        connect(centerViewButton_, SIGNAL(clicked()), this, SLOT(setViewCenter()));
        layoutButtonContainerLayout->addWidget(centerViewButton_);

    graphLayoutButton_ = new QToolButton;
        graphLayoutButton_->setIcon(QIcon(":/qt/icons/sortGraph.png"));
        graphLayoutButton_->setIconSize(NWEMainButtonSize);
        graphLayoutButton_->setToolTip(tr("Lay out network graph"));
        connect(graphLayoutButton_, SIGNAL(clicked()), this, SLOT(sortEntireNetwork()));
        layoutButtonContainerLayout->addWidget(graphLayoutButton_);

//documentation
    documentationButtonContainer_ = new QWidget(this);
    QHBoxLayout* documenationButtonContainerLayout = new QHBoxLayout(documentationButtonContainer_);
    documenationButtonContainerLayout->setMargin(NWEButtonBackgroundMargin);
    documentationButton_ = new QToolButton(this);
        documentationButton_->setIcon(QIcon(":/qt/icons/docu-icons.png"));
        documentationButton_->setIconSize(NWEMainButtonSize);
        documentationButton_->setToolTip(tr("Show/hide documentation items"));
        documentationButton_->setCheckable(true);
        connect(documentationButton_, SIGNAL(clicked()), this, SLOT(toggleDocumentation()));
        documentationButton_->setChecked(showDocumentationProperty_.get());
        toggleDocumentation();
        documenationButtonContainerLayout->addWidget(documentationButton_);

    layoutEditorButtons();
}

void NetworkEditor::layoutEditorButtons() {
//layer buttons
    int x = size().width() - layerButtonContainer_->size().width() - NWEButtonBorderSpacingX;
    int y = NWEButtonBorderSpacingY;
    layerButtonContainer_->move(x,y);
//linking buttons
    x = size().width() - NWEMainButtonSize.width() - NWEButtonBorderSpacingX - NWEButtonBackgroundMargin*3;
    linkingLayerButtonContainer_->move(x,y);
//general linking buttons
    x -= generalLinkingLayerButtonContainer_->width();
    y += layerButtonContainer_->size().height() + NWEMarginLayerToLinking - NWEButtonBackgroundMargin*2;
    generalLinkingLayerButtonContainer_->move(x,y);
//camera linking buttons
    y += generalLinkingLayerButtonContainer_->height() + linkingLayerButtonContainer_->layout()->contentsMargins().top();
    cameraLinkingLayerButtonContainer_->move(x,y);
//port size linking buttons
    y += portSizeLinkingLayerButtonContainer_->height() + linkingLayerButtonContainer_->layout()->contentsMargins().top();
    portSizeLinkingLayerButtonContainer_->move(x,y);
//stop network button
    x = size().width() - stopButtonContainer_->size().width() - NWEButtonBorderSpacingX;
    y = size().height() - stopButtonContainer_->size().height() - NWEButtonBorderSpacingY;
    stopButtonContainer_->move(x,y);
//navigation buttons
    x = NWEButtonBorderSpacingX;
    y = NWEButtonBorderSpacingY;
    navigationButtonContainer_->move(x,y);
//documentation button
    x = NWEButtonBorderSpacingX;
    y = (size().height() - NWEMainButtonSize.height())/2;
    documentationButtonContainer_->move(x,y);
//layout button
    x = NWEButtonBorderSpacingX;
    y = size().height() - layoutButtonContainer_->size().height() - NWEButtonBorderSpacingY;
    layoutButtonContainer_->move(x,y);
}

bool NetworkEditor::cameraLinksHidden() {
    return hideCameraLinksButton_->isChecked();
}

bool NetworkEditor::portSizeLinksHidden() {
    return hidePortSizeLinksButton_->isChecked();
}

void NetworkEditor::setLayerToDataFlow() {
    setLayer(NetworkEditorLayerDataFlow);
}

void NetworkEditor::setLayerToLinking() {
    setLayer(NetworkEditorLayerLinking);
}

void NetworkEditor::setLayerToGeneralLinking() {
    setLayer(NetworkEditorLayerGeneralLinking);
}

void NetworkEditor::setLayerToCameraLinking() {
    setLayer(NetworkEditorLayerCameraLinking);
}

void NetworkEditor::setLayerToPortSizeLinking() {
    setLayer(NetworkEditorLayerPortSizeLinking);
}

void NetworkEditor::setCursorSelect() {
    setCursorMode(NetworkEditorCursorSelectMode);
}

void NetworkEditor::setCursorMove() {
    setCursorMode(NetworkEditorCursorMoveMode);
}

void NetworkEditor::hideCameraLinks() {
    foreach(ProcessorGraphicsItem* processorItem, processorItemMap_) {
        if(processorItem->getPropertyList()->getIsVisibleInEditor()){
            processorItem->getPropertyList()->prepareGeometryChange();
            if(cameraLinksHidden() && portSizeLinksHidden())
                processorItem->getPropertyList()->setPropertyVisibleModifier(PropertyListGraphicsItem::HIDE_CAMERA_AND_SIZE_PROPERTIES);
            else if(cameraLinksHidden())
                processorItem->getPropertyList()->setPropertyVisibleModifier(PropertyListGraphicsItem::HIDE_CAMERA_PROPERTIES);
            else if(portSizeLinksHidden())
                processorItem->getPropertyList()->setPropertyVisibleModifier(PropertyListGraphicsItem::HIDE_SIZE_PROPERTIES);
            else
                processorItem->getPropertyList()->setPropertyVisibleModifier(PropertyListGraphicsItem::HIDE_NO_PROPERTIES);
            processorItem->getPropertyList()->setVisible(true);
        }
        foreach(PortOwnerLinkArrowGraphicsItem* arrow, processorItem->getPortOwnerLinkArrows())
            arrow->setVisible(true);
    }
}

void NetworkEditor::hidePortSizeLinks() {
    foreach(ProcessorGraphicsItem* processorItem, processorItemMap_) {
        if(processorItem->getPropertyList()->getIsVisibleInEditor()){
            processorItem->getPropertyList()->prepareGeometryChange();
            if(cameraLinksHidden() && portSizeLinksHidden())
                processorItem->getPropertyList()->setPropertyVisibleModifier(PropertyListGraphicsItem::HIDE_CAMERA_AND_SIZE_PROPERTIES);
            else if(cameraLinksHidden())
                processorItem->getPropertyList()->setPropertyVisibleModifier(PropertyListGraphicsItem::HIDE_CAMERA_PROPERTIES);
            else if(portSizeLinksHidden())
                processorItem->getPropertyList()->setPropertyVisibleModifier(PropertyListGraphicsItem::HIDE_SIZE_PROPERTIES);
            else
                processorItem->getPropertyList()->setPropertyVisibleModifier(PropertyListGraphicsItem::HIDE_NO_PROPERTIES);
            processorItem->getPropertyList()->setVisible(true);
        }
        foreach(PortOwnerLinkArrowGraphicsItem* arrow, processorItem->getPortOwnerLinkArrows())
            arrow->setVisible(true);
    }
}

void NetworkEditor::removeAllPropertyLinks() {
    if (!getProcessorNetwork())
        return;

    if (QMessageBox::question(this, tr("VoreenVE"), tr("Remove all property links from the current network?"), QMessageBox::Ok | QMessageBox::Cancel) == QMessageBox::Ok) {
        int numLinks = getProcessorNetwork()->removeAllPropertyLinks();
        if (numLinks)
            LINFO("Removed " << numLinks << " property links");
    }
}

void NetworkEditor::linkCamerasAutoChanged() {
    if (!getProcessorNetwork())
        return;

    tgtAssert(getProcessorNetwork(), "No processor network");
    if (getProcessorNetwork()->getMetaDataContainer().hasMetaData("autoLinkCameras")) {
        MetaDataBase* base = getProcessorNetwork()->getMetaDataContainer().getMetaData("autoLinkCameras");
        BoolMetaData* meta = dynamic_cast<BoolMetaData*>(base);
        if (!meta) {
            LWARNING("Meta data object not of expected type 'BoolMetaData'");
            return;
        }
        meta->setValue(linkCamerasAutoButton_->isChecked());
    }
    else {
        BoolMetaData* meta = new BoolMetaData;
        meta->setValue(linkCamerasAutoButton_->isChecked());
        getProcessorNetwork()->getMetaDataContainer().addMetaData("autoLinkCameras", meta);
    }
}

void NetworkEditor::linkCameras() {
    if (!getProcessorNetwork())
        return;

    if (QMessageBox::question(this, tr("VoreenVE"), tr("Link all cameras in the current network?"), QMessageBox::Ok | QMessageBox::Cancel) == QMessageBox::Ok) {
        int numLinks = getProcessorNetwork()->createPropertyLinksWithinSubNetwork<CameraProperty>(
            getProcessorNetwork()->getProcessors(), std::vector<std::string>(), new LinkEvaluatorCameraId());
        if (numLinks)
            LINFO("Created " << numLinks << " camera property links");
    }
}

void NetworkEditor::removeAllCameraLinks() {
    if (!getProcessorNetwork())
        return;

    if (QMessageBox::question(this, tr("VoreenVE"), tr("Remove all camera links from the current network?"), QMessageBox::Ok | QMessageBox::Cancel) == QMessageBox::Ok) {
        int numRemoved = getProcessorNetwork()->removePropertyLinksFromSubNetwork<CameraProperty>(getProcessorNetwork()->getProcessors());
        if (numRemoved)
            LINFO("Removed " << numRemoved << " camera property links");
    }
}

void NetworkEditor::linkPortSizeAutoChanged() {
    if (!getProcessorNetwork())
        return;

    if (getProcessorNetwork()->getMetaDataContainer().hasMetaData("autoLinkPortSize")) {
        MetaDataBase* base = getProcessorNetwork()->getMetaDataContainer().getMetaData("autoLinkPortSize");
        BoolMetaData* meta = dynamic_cast<BoolMetaData*>(base);
        if (!meta) {
            LWARNING("Meta data object not of expected type 'BoolMetaData'");
            return;
        }
        meta->setValue(linkCamerasAutoButton_->isChecked());
    }
    else {
        BoolMetaData* meta = new BoolMetaData;
        meta->setValue(linkPortSizeAutoButton_->isChecked());
        getProcessorNetwork()->getMetaDataContainer().addMetaData("autoLinkPortSize", meta);
    }
}

void NetworkEditor::linkPortSize() {
    if (!getProcessorNetwork())
        return;

    if (QMessageBox::question(this, tr("VoreenVE"), tr("Create render size links for all RenderPorts in the current network?"), QMessageBox::Ok | QMessageBox::Cancel) == QMessageBox::Ok) {
        int numCreated = getProcessorNetwork()->createRenderSizeLinksWithinSubNetwork(getProcessorNetwork()->getProcessors(), false);
        if (numCreated)
            LINFO("Created " << numCreated << " render size links.");
    }
}

void NetworkEditor::removeAllPortSizeLinks() {
    if (!getProcessorNetwork())
        return;

    if (QMessageBox::question(this, tr("VoreenVE"), tr("Remove all render size links from the current network?"), QMessageBox::Ok | QMessageBox::Cancel) == QMessageBox::Ok) {
        int numRemoved = getProcessorNetwork()->removeRenderSizeLinksFromSubNetwork(getProcessorNetwork()->getProcessors());
        if (numRemoved)
            LINFO("Removed " << numRemoved << " render size links.");
    }
}

void NetworkEditor::toggleNetworkEvaluator() {
    if (networkEvaluatorIsLockedByButton_) {
        evaluator_->unlock();
        networkEvaluatorIsLockedByButton_ = false;
        stopNetworkEvaluatorButton_->setIcon(QIcon(":/qt/icons/player-pause.png"));
        stopNetworkEvaluatorButton_->setToolTip(tr("Stop automatic network evaluation"));
        evaluator_->process();
    }
    else {
        evaluator_->lock();
        stopNetworkEvaluatorButton_->setIcon(QIcon(":/qt/icons/player-start.png"));
        stopNetworkEvaluatorButton_->setToolTip(tr("Start automatic network evaluation"));
        networkEvaluatorIsLockedByButton_ = true;
    }
}

void NetworkEditor::sortEntireNetwork() {
    //return, if nothing to do
    if(!getProcessorNetwork() || getProcessorNetwork()->getProcessors().empty())
        return;

    tgtAssert(currentGraphLayout_, "No GraphLayout");

    //list that saves all positions of processors before repositioning
    std::vector<std::pair<QPointF,Processor*> > savePosList;
    //the current network processors
    std::vector<Processor*> processors = getProcessorNetwork()->getProcessors();
    //save position of all processors before sorting anything
    for (size_t j = 0; j < processors.size(); j++) {
        QMap<Processor*,ProcessorGraphicsItem*>::const_iterator i = processorItemMap_.constBegin();
        while (i != processorItemMap_.constEnd()) {
            if(processors[j] == i.key()) {
                QPointF pos = i.value()->pos();
                std::pair<QPointF,Processor*> savePos = std::make_pair(pos,processors[j]);
                savePosList.push_back(savePos);
                ++i;
            }
            else ++i;
        }
    }

    // empty vector to indicate a sort of the entire network
    std::vector<Processor*> selectedProc;
    currentGraphLayout_->sort(getProcessorNetwork(),&selectedProc,&processorItemMap_);

    //center view
    setViewCenter();

    //popup that asks if you want to keep changes. if not, reset position of all processors
    if (QMessageBox::information(NULL, "Network Layout", "An automatic network layout has been created. Do you want to keep the result?", QMessageBox::Yes, QMessageBox::No) == QMessageBox::No) {
        //undo changes
        for(size_t it = 0; it < savePosList.size(); it++) {
            for (size_t it2 = 0; it2 < processors.size(); it2++) {
                if (savePosList[it].second == processors[it2]) {
                    QMap<Processor*,ProcessorGraphicsItem*>::const_iterator i = processorItemMap_.constBegin();
                    while (i != processorItemMap_.constEnd()) {
                        if (i.key() == processors[it2]) {
                            QPointF pos = savePosList[it].first;
                            i.value()->setPos(pos);
                            ++i;
                        }
                        else ++i;
                    }
                }
            }
        }
        setViewCenter();
    } else {
        //save changes
        foreach (ProcessorGraphicsItem* processorItem, processorItemMap_.values()) {
            processorItem->saveMeta();
        }
    }
}

void NetworkEditor::sortSubNetwork() {
    tgtAssert(currentGraphLayout_, "No GraphLayout");

    //get selected Processors
    std::vector<Processor*> selectedProc;
    //list that saves positions of delected processors before repositioning
    std::vector<std::pair<QPointF,Processor*> > savePosList;
    //the current network processors
    std::vector<Processor*> processors = getProcessorNetwork()->getProcessors();
    foreach (QGraphicsItem* item, scene()->selectedItems()) {
        if(item->type() == UserTypesProcessorGraphicsItem) {
            ProcessorGraphicsItem* procItem = qgraphicsitem_cast<ProcessorGraphicsItem*>(item);
            Processor* processor = procItem->getProcessor();
            selectedProc.push_back(processor);
            std::pair<QPointF,Processor*> savePos = std::make_pair(procItem->pos(),processor);
            savePosList.push_back(savePos);
        }
    }
    tgtAssert(selectedProc.size() > 1, "sortSubNetwork called on less then 2 items!");

    //sort
    currentGraphLayout_->sort(getProcessorNetwork(),&selectedProc,&processorItemMap_);

    //popup that asks if you want to keep changes. if not, reset position of all processors
    if (QMessageBox::information(NULL, "Sorting Network", "The processor items have been sorted. Do you want to keep the result?", QMessageBox::Yes, QMessageBox::No) == QMessageBox::No) {
        //undo changes
        for(size_t it = 0; it < savePosList.size(); it++) {
            for (size_t it2 = 0; it2 < processors.size(); it2++) {
                if (savePosList[it].second == processors[it2]) {
                    QMap<Processor*,ProcessorGraphicsItem*>::const_iterator i = processorItemMap_.constBegin();
                    while (i != processorItemMap_.constEnd()) {
                        if (i.key() == processors[it2]) {
                            QPointF pos = savePosList[it].first;
                            i.value()->setPos(pos);
                            ++i;
                        }
                        else ++i;
                    }
                }
            }
        }
    } else {
        //save changes
        foreach (QGraphicsItem* item, scene()->selectedItems()) {
            if(item->type() == UserTypesProcessorGraphicsItem)
            {
                ProcessorGraphicsItem* procItem = qgraphicsitem_cast<ProcessorGraphicsItem*>(item);
                procItem->saveMeta();
            }
        }
    }
}


void NetworkEditor::setViewCenter() {
    resetMatrix();
    resetTransform();
    if (!getProcessorNetwork()){
        centerOn(0.0,0.0);
    } else {
        getProcessorNetwork()->getMetaDataContainer().removeMetaData("ZoomFactor");
        getProcessorNetwork()->getMetaDataContainer().removeMetaData("ZoomCenter");
        scaleView();
    }
}

void NetworkEditor::linkCamerasOfProcessor(const Processor* processor) {

    if (!getProcessorNetwork())
        return;

    // get camera properties of processor to link
    std::vector<CameraProperty*> camPropsProcessor = processor->getPropertiesByType<CameraProperty>();
    if (camPropsProcessor.empty())
        return;

    // get camera properties of network processors
    std::vector<CameraProperty*> camPropsNetwork = getProcessorNetwork()->getPropertiesByType<CameraProperty>();
    // remove processor's properties from network prop vector
    for (size_t i=0; i<camPropsProcessor.size(); ++i) {
        std::vector<CameraProperty*>::iterator iter = std::find(camPropsNetwork.begin(), camPropsNetwork.end(), camPropsProcessor[i]);
        if (iter != camPropsNetwork.end())
            camPropsNetwork.erase(iter);
    }

    if (camPropsNetwork.empty())
        return;

    int numLinks = 0;
    for (size_t i=0; i<camPropsProcessor.size(); ++i) {

        // create link with first network camera property
        if (!camPropsNetwork.front()->isLinkedWith(camPropsProcessor[i])) {
            if (getProcessorNetwork()->createPropertyLink(camPropsNetwork.front(), camPropsProcessor[i]))
                numLinks++;
        }
        if (!camPropsProcessor[i]->isLinkedWith(camPropsNetwork.front()))
            if (getProcessorNetwork()->createPropertyLink(camPropsProcessor[i], camPropsNetwork.front()))
                numLinks++;

        if (camPropsNetwork.size() == 1)
            continue;

        // create link with last network camera property
        if (!camPropsNetwork.back()->isLinkedWith(camPropsProcessor[i])) {
            if (getProcessorNetwork()->createPropertyLink(camPropsNetwork.back(), camPropsProcessor[i]))
                numLinks++;
        }
        if (!camPropsProcessor[i]->isLinkedWith(camPropsNetwork.back()))
            if (getProcessorNetwork()->createPropertyLink(camPropsProcessor[i], camPropsNetwork.back()))
                numLinks++;

        // create link with all network camera properties,
        // unless they are already (directly or indirectly) linked
        for (size_t j=0; j<camPropsNetwork.size(); ++j) {
            if (!camPropsNetwork[j]->isLinkedWith(camPropsProcessor[i], true)) {
                if (getProcessorNetwork()->createPropertyLink(camPropsNetwork[j], camPropsProcessor[i]))
                    numLinks++;
            }
            if (!camPropsProcessor[i]->isLinkedWith(camPropsNetwork[j], true))
                if (getProcessorNetwork()->createPropertyLink(camPropsProcessor[i], camPropsNetwork[j]))
                    numLinks++;
        }
    }

    LINFO("Created " << numLinks << " camera property links");
}

//---------------------------------------------------------------------------------------------------------------
//                  documentations
//---------------------------------------------------------------------------------------------------------------
void NetworkEditor::toggleDocumentation() {
    bool checked = documentationButton_->isChecked();
    showDocumentationProperty_.set(checked);
    foreach(TextBoxBaseGraphicsItem* item, documentationItems_) {
        item->setVisible(checked);
    }
}

void NetworkEditor::serializeTextItems() {
    if(workspace_ && workspace_->getProcessorNetwork()) {
        //replace all old items
        SerializableVectorMetaData<TextBoxMetaData*>*meta = new SerializableVectorMetaData<TextBoxMetaData*>(std::vector<TextBoxMetaData*>(),true);
        workspace_->getProcessorNetwork()->getMetaDataContainer().addMetaData("DocumentationGraphicsItems",meta);

        foreach(TextBoxBaseGraphicsItem* item, documentationItems_) {
            meta->addValue(new TextBoxMetaData(item));
        }
        //notify workspace change
        workspace_->setModified(true);
    }
}

void NetworkEditor::deserializeTextItems() {
    if(workspace_ && workspace_->getProcessorNetwork()) {
        SerializableVectorMetaData<TextBoxMetaData*>* meta = dynamic_cast<SerializableVectorMetaData<TextBoxMetaData*>*>
                (workspace_->getProcessorNetwork()->getMetaDataContainer().getMetaData("DocumentationGraphicsItems"));
        //check if meta exists
        if(!meta) return;
        for(size_t i = 0; i < meta->getValues().size(); i++) {
            TextBoxBaseGraphicsItem* item = NULL;
            switch(meta->getValues()[i]->getUserType()) {
            case UserTypesTextBoxGraphicsItem:
                item = new TextBoxGraphicsItem(this);
                break;
            case UserTypesFrameBoxGraphicsItem:
                item = new FrameBoxGraphicsItem(this);
                break;
            default:
                tgtAssert(false,"unknown type");
            }
            item->copyFromMeta(meta->getValues()[i]);
            item->setVisible(documentationButton_->isChecked());
            documentationItems_.push_back(item);
            scene()->addItem(item);
        }
    }
}

//---------------------------------------------------------------------------------------------------------------
//                  events
//---------------------------------------------------------------------------------------------------------------
void NetworkEditor::resizeEvent(QResizeEvent* event) {
    QGraphicsView::resizeEvent(event);

    layoutEditorButtons();
    if (isVisible() && needsScale_) {
        scaleView();
        needsScale_ = false;
    }
}

void NetworkEditor::wheelEvent(QWheelEvent *event) {
    float factor = pow(2.0, event->delta() / 360.0);
    scale(factor, factor);
}

void NetworkEditor::mousePressEvent(QMouseEvent* event) {
    // shift and left button activate translation of scene.
    if (event->button() == Qt::LeftButton && (event->modifiers() == Qt::ShiftModifier || currentCursorMode_ == NetworkEditorCursorMoveMode)) {
        translateScene_ = true;
        translateSceneVector_ = mapToScene(event->pos());
        lastTranslateCenter_ =  mapToScene(viewport()->rect().center());
        setCursor(Qt::ClosedHandCursor);
    }
    else
        if(!(event->button() & Qt::RightButton))
            QGraphicsView::mousePressEvent(event);
        else
            nextTextBoxPos_ = event->pos();
}

void NetworkEditor::mouseMoveEvent(QMouseEvent* event) {
    if (translateScene_) {
        //set new center
        translateSceneVector_ -= mapToScene(event->pos());
        lastTranslateCenter_ += translateSceneVector_;
        centerOn(lastTranslateCenter_);
        translateSceneVector_ = mapToScene(event->pos());
        //save new center
        if (getProcessorNetwork()) {
            if (getProcessorNetwork()->getMetaDataContainer().hasMetaData("ZoomCenter")) {
                MetaDataBase* base = getProcessorNetwork()->getMetaDataContainer().getMetaData("ZoomCenter");
                Vec2MetaData* meta = dynamic_cast<Vec2MetaData*>(base);
                meta->setValue(tgt::vec2(lastTranslateCenter_.x(),lastTranslateCenter_.y()));
            } else {
                getProcessorNetwork()->getMetaDataContainer().addMetaData("ZoomCenter", new Vec2MetaData(tgt::vec2(lastTranslateCenter_.x(),lastTranslateCenter_.y())));
            }
        }
    }
    QGraphicsView::mouseMoveEvent(event);
}

void NetworkEditor::mouseReleaseEvent(QMouseEvent* event) {
    if (event->button() == Qt::LeftButton){
        if(translateScene_){
            if(currentCursorMode_ == NetworkEditorCursorMoveMode)
                setCursor(Qt::OpenHandCursor);
            else
                setCursor(Qt::ArrowCursor);
            translateScene_ = false;
        } else {
            QGraphicsView::mouseReleaseEvent(event);
            if(currentCursorMode_ != NetworkEditorCursorMoveMode){
                updateSelectedItems();
            }
        }
    } else {
       nextTextBoxPos_ = event->pos();
       QGraphicsView::mouseReleaseEvent(event);
    }
}

void NetworkEditor::keyPressEvent(QKeyEvent* event) {
    QGraphicsView::keyPressEvent(event);

    if (!event->isAccepted() && (event->key() == Qt::Key_Delete || event->key() == Qt::Key_Backspace)) {
        if(QApplication::overrideCursor()) QApplication::restoreOverrideCursor(); //restore cursor to prevent getting stuck
        deleteActionSlot();
    }
    else if (!event->isAccepted() && (event->key() == Qt::Key_Plus)) {
        float factor = keyPressScaleFactor;
        scale(factor, factor);
    }
    else if (!event->isAccepted() && (event->key() == Qt::Key_Minus)) {
        float factor = 1.f / keyPressScaleFactor;
        scale(factor, factor);
    }
    else if (!event->isAccepted() && (event->key() == Qt::Key_C) && (event->modifiers() & Qt::CTRL))
        copyActionSlot();
    else if (!event->isAccepted() && (event->key() == Qt::Key_V) && (event->modifiers() & Qt::CTRL)) {
        rightClickPosition_ = mapToScene(mapFromGlobal(QCursor::pos()));
        pasteActionSlot();
    }
    else if (!event->isAccepted() && (event->key() == Qt::Key_F2) && currentLayer_ == NetworkEditorLayerDataFlow) {
        QList<QGraphicsItem*> items = scene()->selectedItems();
        if (items.count() == 1) {
            QGraphicsItem* item = items[0];
            PortOwnerGraphicsItem* poItem = dynamic_cast<PortOwnerGraphicsItem*>(item);
            if (poItem)
                poItem->enterRenameMode();
        }
    }
}

void NetworkEditor::dragEnterEvent(QDragEnterEvent* event) {
    if (event->mimeData()->hasText()) {
        event->setDropAction(Qt::CopyAction);
        event->accept();
    }
    else
        event->ignore();
}

void NetworkEditor::dragMoveEvent(QDragMoveEvent* event) {
    QGraphicsItem* item = itemAt(event->pos());

    // Remove old selection.
    if (selectedItem_) {
        if(selectedItem_->type() == UserTypesPortGraphicsItem)
            qgraphicsitem_cast<PortGraphicsItem*>(selectedItem_)->setScale(1.0);
        else
            selectedItem_->setSelected(false);
    }

    if (item) {
        // separate handling for processors and documentation boxes
        if (item->type() != UserTypesPortGraphicsItem && item->parentItem() && item->parentItem()->type() == UserTypesProcessorGraphicsItem)
            item = item->parentItem();

        if (item->type() == UserTypesPortGraphicsItem) {
            PortGraphicsItem* port = qgraphicsitem_cast<PortGraphicsItem*>(item);
            if (port && port->isOutport()) {
                // Create a temporary processor to check connectivity.
                std::unique_ptr<Processor> processor(dynamic_cast<Processor*>(VoreenApplication::app()->createSerializableType(event->mimeData()->text().toStdString())));
                if (processor && getProcessorNetwork()->canTakeOverConnections(port->getPort(), processor.get())) {
                    selectedItem_ = item;
                    selectedItem_->setScale(2.0);
                    return;
                }
            }
        } else if (item->type() == UserTypesPortArrowGraphicsItem || item->type() == UserTypesProcessorGraphicsItem) {
            selectedItem_ = item;
            selectedItem_->setSelected(true);
            return;
        }
    }

    selectedItem_ = nullptr;

    // Case of dragging a processor (or possibly text files) onto the canvas:
    // Qt >= 5.12 apparently requires accepting every single DragMoveEvent.
    if (event->mimeData()->hasText()) {
        event->accept();
    }
}

void NetworkEditor::dragLeaveEvent(QDragLeaveEvent*) {}

void NetworkEditor::dropEvent(QDropEvent* event) {

    if (!getProcessorNetwork())
        return;

    if (evaluator_->isLocked() && !networkEvaluatorIsLockedByButton_) {
        QMessageBox::information(this, tr("Network Locked"), tr("The network is being evaluated, so dropping processors is not allowed"));
        return;
    }

    QGraphicsItem* lowerItem = itemAt(event->pos());
    if (event->mimeData()->hasText()) {
        event->setDropAction(Qt::CopyAction);
        event->accept();

        Processor* proc = dynamic_cast<Processor*>(VoreenApplication::app()->createSerializableType(event->mimeData()->text().toStdString()));

        if(!proc) {
            LDEBUG("Drop event contains no valid processor");
            return;
        }

        std::string processorName = proc->getClassName();
        tgtAssert(!processorName.empty(), "Processor class name is empty");

        ProcessorGraphicsItem* item;
        QPointF p;
        bool hasBeenReplaced = false;   // check if dropped on a processor to replace it
        if (lowerItem) {
            switch (lowerItem->type()) {
                case QGraphicsTextItem::Type:
                case UserTypesRenamableTextGraphicsItem:
                case UserTypesProgressBarGraphicsItem:
                case UserTypesPropertyListButtonGraphicsItem:
                    {
                        // separate handling for processors and documentation boxes
                        if (lowerItem->parentItem() && lowerItem->parentItem()->type() == UserTypesProcessorGraphicsItem)
                            lowerItem = lowerItem->parentItem();
                        else
                            break;
                    }
                case UserTypesProcessorGraphicsItem:
                    {
                        p = lowerItem->pos();
                        if(ProcessorGraphicsItem* i = qgraphicsitem_cast<ProcessorGraphicsItem*>(lowerItem))
                            getProcessorNetwork()->replaceProcessor(i->getProcessor(), proc);
                        else
                            return;
                        updateSelectedItems();
                        hasBeenReplaced = true;
                        break;
                    }
                case UserTypesPortArrowGraphicsItem:
                    {
                        // The following cast is guaranteed to work, since the arrow is selected.
                        PortArrowGraphicsItem* selectedPortArrow = qgraphicsitem_cast<PortArrowGraphicsItem*>(selectedItem_);
                        PortGraphicsItem* srcItem = selectedPortArrow->getSourceItem();
                        PortGraphicsItem* dstItem = selectedPortArrow->getDestinationItem();
                        getProcessorNetwork()->addProcessorInConnection(srcItem->getPort(), dstItem->getPort(), proc);
                        p = mapToScene(event->pos());
                        hasBeenReplaced = true;     // not really replaced, but has been handled already
                        break;
                    }
                case UserTypesPortGraphicsItem:
                    {
                        // The following cast is guaranteed to work, since the arrow is selected.
                        PortGraphicsItem* port = qgraphicsitem_cast<PortGraphicsItem*>(lowerItem);
                        if (port->isOutport() && getProcessorNetwork()->canTakeOverConnections(port->getPort(), proc)) {
                            getProcessorNetwork()->takeOverConnections(port->getPort(), proc);
                            port->setScale(1.0); // restore scale from prior drag move event
                            p = mapToScene(event->pos());
                            hasBeenReplaced = true;
                        }
                        break;
                    }
                case UserTypesFrameBoxGraphicsItem:
                    break;
                default:
                    tgtAssert(false, "should not get here");
            }
        }

        if (!hasBeenReplaced) {
            // not replaced / not handled before
            getProcessorNetwork()->addProcessor(proc, processorName);
            p = mapToScene(event->pos());
        }
        tgtAssert(processorItemMap_.contains(proc), "processorItemMap didn't contain the processor");
        item = processorItemMap_[proc];
        tgtAssert(item, "no ProcessorGraphicsItem for proc");
        item->setPos(p);
        item->saveMeta();
        //adjustLinkArrowGraphicsItems();
        // make sure that the added processor is initialized
        evaluator_->initializeNetwork();
        //select new processor
        scene()->clearSelection();
        item->setSelected(true);
        selectedItem_ = nullptr;
        updateSelectedItems();
        //to allow direct delete. otherwise the processorlistwidget keeps the focus
        setFocus();
    }
    else {
        event->ignore();
    }
}

//---------------------------------------------------------------------------------------------------------------
//                  action slots
//---------------------------------------------------------------------------------------------------------------
bool NetworkEditor::clipboardHasValidContent() {
    QClipboard* clipboard = QApplication::clipboard();
    QString clip = clipboard->text();
    std::stringstream stream(clip.toStdString());

    XmlDeserializer d;
    d.setUseAttributes(true);
    try {
        d.read(stream);
    }
    catch (SerializationException&) {
        return false;
    }

    ProcessorNetwork* tmp = 0;
    try {
        d.deserialize("ProcessorNetwork", tmp);
    }
    catch (SerializationException& e){
        LWARNINGC("voreen.NetworkEditor", "Error during deserialization of temporary network:" << e.what());
        delete tmp;
        return false;
    }
    delete tmp;
    return true;
}

void NetworkEditor::deleteActionSlot() {
    removeItems(scene()->selectedItems());
    updateSelectedItems();

    emit processorsSelected(QList<Processor*>());
}

void NetworkEditor::copyActionSlot() {
    if (evaluator_->isLocked() && !networkEvaluatorIsLockedByButton_) {
        QMessageBox::information(this, tr("Network Locked"), tr("The network is being evaluated, so copying is not allowed"));
        return;
    }

    //stream used for serialization
    std::stringstream stream;
    XmlSerializer s;
    s.setUseAttributes(true);

    //processor and documentation items can be stored
    std::vector<Processor*> processors;
    std::vector<TextBoxMetaData> docItems;


    foreach (QGraphicsItem* item, scene()->selectedItems()) {
        switch (item->type()) {
        case UserTypesProcessorGraphicsItem:
            {
                ProcessorGraphicsItem* procItem = qgraphicsitem_cast<ProcessorGraphicsItem*>(item);
                Processor* processor = procItem->getProcessor();
                processors.push_back(processor);
                break;
            }
        case UserTypesTextBoxBaseGraphicsItem:
        case UserTypesTextBoxGraphicsItem:
        case UserTypesFrameBoxGraphicsItem:
            docItems.push_back(TextBoxMetaData(qgraphicsitem_cast<TextBoxBaseGraphicsItem*>(item)));
        }
    }
    ProcessorNetwork* copyNetwork = getProcessorNetwork()->cloneSubNetwork(processors);
    Serializer serializer(s);
    serializer.serialize("ProcessorNetwork", copyNetwork);

    serializer.serialize("DocumentationItems", docItems);

    s.write(stream);

    QString value = QString::fromStdString(stream.str());

    QClipboard* clipboard = QApplication::clipboard();
    clipboard->setText(value);

    delete copyNetwork;
}

void NetworkEditor::pasteActionSlot() {
    if (!getProcessorNetwork())
        return;

    if (evaluator_->isLocked() && !networkEvaluatorIsLockedByButton_) {
        QMessageBox::information(this, tr("Network Locked"), tr("The network is being evaluated, so pasting is not allowed"));
        return;
    }

    QClipboard* clipboard = QApplication::clipboard();
    QString clip = clipboard->text();
    std::stringstream stream(clip.toStdString());

    XmlDeserializer d;
    Deserializer deserializer(d);
    d.setUseAttributes(true);
    try{
        d.read(stream);
    } catch (SerializationException&){
        deserializer.removeLastError();
        return; //no valid string
    }

    //try to deserialize network
    ProcessorNetwork* pasteNetwork = 0;
    try{
        deserializer.deserialize("ProcessorNetwork", pasteNetwork);
    }
    catch (SerializationException&){
        pasteNetwork = 0;
        deserializer.removeLastError();
    }

    //try to deserialize docItems
    std::vector<TextBoxMetaData> docMeta;
    std::vector<TextBoxBaseGraphicsItem*> docItems;
    try{
        deserializer.deserialize("DocumentationItems", docMeta);
    }
    catch (SerializationException&){
        deserializer.removeLastError();
    }
    for(size_t i = 0; i < docMeta.size(); i++) {
        TextBoxBaseGraphicsItem* item = NULL;
        switch(docMeta[i].getUserType()) {
        case UserTypesTextBoxGraphicsItem:
            item = new TextBoxGraphicsItem(this);
            break;
        case UserTypesFrameBoxGraphicsItem:
            item = new FrameBoxGraphicsItem(this);
            break;
        default:
            tgtAssert(false,"unknown type");
        }
        item->copyFromMeta(&(docMeta[i]));
        docItems.push_back(item);
    }

    //calculate new position
    int* minX = 0; int* maxX = 0;
    int* minY = 0; int* maxY = 0;
    for (size_t i = 0; i < (pasteNetwork ? pasteNetwork->getProcessors().size() : 0); ++i) {
        if (PositionMetaData* posData = dynamic_cast<PositionMetaData*>(pasteNetwork->getProcessors()[i]->getMetaDataContainer().getMetaData("ProcessorGraphicsItem"))){
            tgt::ivec2 pos(posData->getX(),posData->getY());
            if(!minX) minX = new int(pos.x);
            else if (pos.x < *minX) *minX = pos.x;
            if(!maxX) maxX = new int(pos.x);
            else if (pos.x > *maxX) *maxX = pos.x;
            if(!minY) minY = new int(pos.y);
            else if (pos.y < *minY) *minY = pos.y;
            if(!maxY) maxY = new int(pos.y);
            else if (pos.y > *maxY) *maxY = pos.y;
        }
    }
    for (size_t i = 0; i < docItems.size(); ++i) {
        tgt::ivec2 pos(docItems[i]->pos().rx(),docItems[i]->pos().ry());
        if(!minX) minX = new int(pos.x);
        else if (pos.x < *minX) *minX = pos.x;
        if(!maxX) maxX = new int(pos.x);
        else if (pos.x > *maxX) *maxX = pos.x;
        if(!minY) minY = new int(pos.y);
        else if (pos.y < *minY) *minY = pos.y;
        if(!maxY) maxY = new int(pos.y);
        else if (pos.y > *maxY) *maxY = pos.y;
    }

    // add least one item was valid
    if(minX && minY && maxX && maxY) {
        QPointF center((*minX + *maxX)/2,(*minY + *maxY)/2);
        center = rightClickPosition_-center;
        delete minX; delete minY; delete maxX; delete maxY;

        //set new positions
        for (size_t i = 0; i < (pasteNetwork ? pasteNetwork->getProcessors().size() : 0); ++i) {
            if (PositionMetaData* posData = dynamic_cast<PositionMetaData*>(pasteNetwork->getProcessors()[i]->getMetaDataContainer().getMetaData("ProcessorGraphicsItem"))){
                tgt::ivec2 pos(posData->getX(),posData->getY());
                posData->setX(pos.x+static_cast<int>(center.rx()));
                posData->setY(pos.y+static_cast<int>(center.ry()));
            }
        }
        for (size_t i = 0; i < docItems.size(); ++i) {
            docItems[i]->setPos(docItems[i]->pos()+center);
        }

        // add pasted elements to workspace
        for(size_t i = 0; i < docItems.size(); i++) {
            docItems[i]->setVisible(documentationButton_->isChecked());
            documentationItems_.push_back(docItems[i]);
            scene()->addItem(docItems[i]);
        }
        if(pasteNetwork) {
            getProcessorNetwork()->mergeSubNetwork(pasteNetwork);
            evaluator_->initializeNetwork();
        }
        //notify workspace change
        if(pasteNetwork || !docItems.empty())
            workspace_->setModified(true);
    }
}

void NetworkEditor::openPropertyLinkDialog(PortOwnerGraphicsItem* src, PortOwnerGraphicsItem* dest) {
    if (!getProcessorNetwork())
        return;

    if (evaluator_->isLocked() && !networkEvaluatorIsLockedByButton_) {
        QMessageBox::information(this, tr("Network Locked"), tr("The network is being evaluated, so creating links is not allowed"));
        return;
    }

    //switch all propertylists off
    QList<ProcessorGraphicsItem*> list;
    foreach(ProcessorGraphicsItem* item, processorItemMap_){
        if(item != src && item != dest && item->getPropertyList()->getIsVisibleInEditor()){
            list.append(item);
            item->togglePropertyList();
        }
    }

    PropertyLinkDialog* dialog = new PropertyLinkDialog(this, src, dest);
    dialog->exec(); //deletes dialog on return

    //switch lists back on
    foreach(ProcessorGraphicsItem* item, list)
        item->togglePropertyList();

    src->getPropertyList()->setIsVisibleInEditor(!src->getPropertyList()->getIsVisibleInEditor());
    src->togglePropertyList();
    if(src != dest){
        dest->getPropertyList()->setIsVisibleInEditor(!dest->getPropertyList()->getIsVisibleInEditor());
        dest->togglePropertyList();
    }
}

void NetworkEditor::editLinkActionSlot() {
    if (evaluator_->isLocked() && !networkEvaluatorIsLockedByButton_) {
        QMessageBox::information(this, tr("Network Locked"), tr("The network is being evaluated, so editing links is not allowed"));
        return;
    }

    tgtAssert(scene()->selectedItems().size() == 1,"to many items for this slot");
    tgtAssert(dynamic_cast<PropertyLinkArrowGraphicsItem*>(scene()->selectedItems().first()) ||
              dynamic_cast<PortSizeLinkArrowGraphicsItem*>(scene()->selectedItems().first()),"wrong type of selected item");

    PropertyLinkArrowGraphicsItem* plItem = 0;
    PortSizeLinkArrowGraphicsItem* pslItem = 0;
    if((plItem = dynamic_cast<PropertyLinkArrowGraphicsItem*>(scene()->selectedItems().first()))){
        PortOwnerGraphicsItem* src = 0, *dst = 0;
        if((src = dynamic_cast<PortOwnerGraphicsItem*>(plItem->getSourceItem()->getPropertyOwnerItem()))) {}
        else
            src = dynamic_cast<PortOwnerGraphicsItem*>(plItem->getSourceItem()->getPropertyOwnerItem()->parent());
        if((dst = dynamic_cast<PortOwnerGraphicsItem*>(plItem->getDestinationItem()->getPropertyOwnerItem()))) {}
        else
            dst = dynamic_cast<PortOwnerGraphicsItem*>(plItem->getDestinationItem()->getPropertyOwnerItem()->parent());
        openPropertyLinkDialog(src,dst);
    }
    else if((pslItem = dynamic_cast<PortSizeLinkArrowGraphicsItem*>(scene()->selectedItems().first()))) {
        openPropertyLinkDialog(pslItem->getSourcePort()->getPortOwner(),pslItem->getDestinationPort()->getPortOwner());
    }
}

void NetworkEditor::deleteLinksActionSlot() {
    if (!getProcessorNetwork())
        return;
    QList<ProcessorGraphicsItem*> processors = this->getSelectedProcessorGraphicsItems();

    int number = 0;
    if(processors.size() == 0) {//do it for all
        removeAllPropertyLinks();
    }
    else if(processors.size() == 1) {
        number = getProcessorNetwork()->removePropertyLinksFromProcessor(processors.first()->getProcessor());
    }
    else {
        std::vector<Processor*> vec;
        foreach(ProcessorGraphicsItem* item, processors)
            vec.push_back(item->getProcessor());
        number = getProcessorNetwork()->removePropertyLinksFromSubNetwork(vec);
    }

    if (number)
        LINFO("Removed " << number << " property links.");
}

void NetworkEditor::deleteCameraLinksActionSlot() {
    if (!getProcessorNetwork())
        return;
    QList<ProcessorGraphicsItem*> processors = this->getSelectedProcessorGraphicsItems();

    int number = 0;
    if (processors.size() == 0) //do it for all
        removeAllCameraLinks();
    else if(processors.size() == 1) {
        number = getProcessorNetwork()->removePropertyLinksFromProcessor<CameraProperty>(processors.first()->getProcessor());
    }
    else {
        std::vector<Processor*> vec;
        foreach(ProcessorGraphicsItem* item, processors)
            vec.push_back(item->getProcessor());
       number = getProcessorNetwork()->removePropertyLinksFromSubNetwork<CameraProperty>(vec);
    }

    if (number)
        LINFO("Removed " << number << " camera property links.");
}

void NetworkEditor::deletePortSizeLinksActionSlot() {
    if (!getProcessorNetwork())
        return;
    QList<ProcessorGraphicsItem*> processors = this->getSelectedProcessorGraphicsItems();

    int number = 0;
    if (processors.size() == 0) //do it for all
        removeAllPortSizeLinks();
    else if (processors.size() == 1){
        number = getProcessorNetwork()->removeRenderSizeLinksFromProcessor(processors.first()->getProcessor());
    }
    else {
        std::vector<Processor*> vec;
        foreach(ProcessorGraphicsItem* item, processors)
            vec.push_back(item->getProcessor());
        number = getProcessorNetwork()->removeRenderSizeLinksFromSubNetwork(vec);
    }

    if (number)
        LINFO("Removed " << number << " render size links.");
}

void NetworkEditor::createCameraLinksActionSlot() {
    if (!getProcessorNetwork())
        return;
    QList<ProcessorGraphicsItem*> processors = this->getSelectedProcessorGraphicsItems();

    int number = 0;
    if (processors.size() == 0) //do it for all
        number = getProcessorNetwork()->createPropertyLinksWithinSubNetwork<CameraProperty>(getProcessorNetwork()->getProcessors(),std::vector<std::string>());
    else if (processors.size() == 1) {
        number = getProcessorNetwork()->createPropertyLinksForProcessor<CameraProperty>(processors.first()->getProcessor());
    }
    else {
        std::vector<Processor*> vec;
        foreach(ProcessorGraphicsItem* item, processors)
            vec.push_back(item->getProcessor());
        number = getProcessorNetwork()->createPropertyLinksWithinSubNetwork<CameraProperty>(vec, std::vector<std::string>());
    }

    if (number)
        LINFO("Created " << number << " camera property links.");
}

void NetworkEditor::createPortSizeLinksActionSlot() {
    if (!getProcessorNetwork())
        return;
    QList<ProcessorGraphicsItem*> processors = this->getSelectedProcessorGraphicsItems();

    int number = 0;
    if(processors.size() == 0)//do it for all
        number = getProcessorNetwork()->createRenderSizeLinksWithinSubNetwork(getProcessorNetwork()->getProcessors());
    else if (processors.size() == 1) {
        number = getProcessorNetwork()->createRenderSizeLinksForProcessor(processors.first()->getProcessor());
    }
    else {
        std::vector<Processor*> vec;
        foreach(ProcessorGraphicsItem* item, processors)
            vec.push_back(item->getProcessor());
        number = getProcessorNetwork()->createRenderSizeLinksWithinSubNetwork(vec);
    }

    if (number)
        LINFO("Created " << number << " render size links.");
}

void NetworkEditor::contextMenuEvent(QContextMenuEvent* event) {
    //create new menu
    QMenu currentMenu;
    event->accept();

    //empty action
    QAction* emptyAction = new QAction(QString("no action"),&currentMenu);
    emptyAction->setEnabled(false);

    QList<QGraphicsItem*> selectedItems = scene()->selectedItems();

    //HACK to handle right click, while an arrow moved
    if(selectedItems.size() == 1 && dynamic_cast<PortGraphicsItem*>(selectedItems.first())){
        scene()->clearSelection(); //still needs improvement :(
        selectedItems.clear();
    }

    QGraphicsItem* item = itemAt(event->pos());

    //move to parent item
    if (item) {
        switch (item->type()) {
        case UserTypesRenamableTextGraphicsItem:
        case UserTypesPropertyListButtonGraphicsItem:
        case QGraphicsProxyWidget::Type:
            item = item->parentItem();
        }
    }

    //action belonging to just one item
    if((item && item->isSelected() && scene()->selectedItems().size() == 1) ||
       (item && !item->isSelected() && item->acceptHoverEvents())) {
        scene()->clearSelection();
        selectedItems.clear();
        item->setSelected(true);
        selectedItems.append(item);

        switch (item->type()) {
        case UserTypesProcessorGraphicsItem: {
            dynamic_cast<NWEBaseGraphicsItem*>(item)->addActionsToContextMenu(&currentMenu);

            currentMenu.addSeparator();
            currentMenu.addAction(createPortOwnerCameraLinksAction_);
            currentMenu.addAction(createPortOwnerPortSizeLinksAction_);
            currentMenu.addAction(deletePortOwnerLinksAction_);
            currentMenu.addAction(deletePortOwnerCameraLinksAction_);
            currentMenu.addAction(deletePortOwnerPortSizeLinksAction_);
            currentMenu.addSeparator();
            currentMenu.addAction(deleteAction_);
            currentMenu.addAction(copyAction_);

            PortOwnerGraphicsItem* poItem = dynamic_cast<PortOwnerGraphicsItem*>(item);
            if(!poItem->getWidgetToggleButton().getWidgets().empty()) {
                currentMenu.addSeparator();
                poItem->getWidgetToggleButton().addActionsToContextMenu(&currentMenu);
            }

            } break;
        case UserTypesPortArrowGraphicsItem:
            {
            /*bool isPartOfBundle = false;
            QList<PortArrowGraphicsItem*> items = convertQList<QGraphicsItem*, PortArrowGraphicsItem*>(selectedItems);
            foreach(PortArrowGraphicsItem* arrow, items) {
                if(arrow->isBundled()) {
                    isPartOfBundle = true;
                    break;
                }
            }
            if(isPartOfBundle) {
                currentMenu.addAction(unbundleAction_);
                currentMenu.addAction(addHandleAction_);
            }*/
            currentMenu.addAction(deleteAction_);
            break;
            }
        case UserTypesPropertyLinkArrowGraphicsItem:
            {
            currentMenu.addAction(deleteAction_);
            currentMenu.addAction(editLinkAction_);
            break;
            }
        case UserTypesPortSizeLinkArrowGraphicsItem:
            {
            currentMenu.addAction(deleteAction_);
            currentMenu.addAction(editLinkAction_);
            break;
            }
        case UserTypesWidgetToggleButtonGraphicsItem:
            {
            WidgetToggleButtonGraphicsItem* wtb = dynamic_cast<WidgetToggleButtonGraphicsItem*>(item);
            if(!wtb->getWidgets().empty())
                wtb->addActionsToContextMenu(&currentMenu);
            break;
            }
        case UserTypesPortGraphicsItem:
            {
            qgraphicsitem_cast<PortGraphicsItem*>(item)->addActionsToContextMenu(&currentMenu);
            if(currentMenu.actions().empty())
                currentMenu.addAction(emptyAction);
            break;
            }
        case UserTypesTextBoxBaseGraphicsItem:
        case UserTypesTextBoxGraphicsItem:
        case UserTypesFrameBoxGraphicsItem:
            {
            qgraphicsitem_cast<TextBoxBaseGraphicsItem*>(item)->addActionsToContextMenu(&currentMenu);
            if(currentMenu.actions().empty())
                currentMenu.addAction(emptyAction);
            break;
            }
        default:
            return;
        }
    } else
    //action with multiple items
    if(item && item->isSelected()) {
        QList<PortOwnerGraphicsItem*> list = getSelectedPortOwnerGraphicsItems();
        if(!list.empty()){
            if(list.size() == 1) {
                currentMenu.addAction(createPortOwnerCameraLinksAction_);
                currentMenu.addAction(createPortOwnerPortSizeLinksAction_);
                currentMenu.addAction(deletePortOwnerLinksAction_);
                currentMenu.addAction(deletePortOwnerCameraLinksAction_);
                currentMenu.addAction(deletePortOwnerPortSizeLinksAction_);
            } else
            if(list.size() > 1) {
                currentMenu.addAction(createInnerCameraLinksAction_);
                currentMenu.addAction(createInnerPortSizeLinksAction_);
                currentMenu.addAction(deleteInnerLinksAction_);
                currentMenu.addAction(deleteInnerCameraLinksAction_);
                currentMenu.addAction(deleteInnerPortSizeLinksAction_);
                currentMenu.addSeparator();
                currentMenu.addAction(sortSubNetworkAction_);
                currentMenu.addSeparator();
            }
        }

        currentMenu.addSeparator();
        currentMenu.addAction(deleteAction_);
        currentMenu.addAction(copyAction_);

    }
    // no items have been hit
    else {
        scene()->clearSelection();
        currentMenu.addAction(deleteAllLinksAction_);
        currentMenu.addAction(deleteAllCameraLinksAction_);
        currentMenu.addAction(deleteAllPortSizeLinksAction_);
        currentMenu.addAction(createAllCameraLinksAction_);
        currentMenu.addAction(createAllPortSizeLinksAction_);

        if(clipboardHasValidContent()){
            currentMenu.addSeparator();
            rightClickPosition_ = mapToScene(mapFromGlobal(QCursor::pos()));
            currentMenu.addAction(pasteAction_);
        }
        currentMenu.addSeparator();
        currentMenu.addAction(createNewTextNoteAction_);
        currentMenu.addAction(createNewTextFrameAction_);
    }
    currentMenu.exec(event->globalPos());
}

} // namespace voreen
