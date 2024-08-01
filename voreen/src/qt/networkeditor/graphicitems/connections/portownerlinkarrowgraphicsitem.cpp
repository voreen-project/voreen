/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2024 University of Muenster, Germany,                        *
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

#include "voreen/qt/networkeditor/graphicitems/connections/portownerlinkarrowgraphicsitem.h"
#include "voreen/qt/networkeditor/styles/nwestyle_base.h"

//core
#include "voreen/core/processors/processor.h"
#include "voreen/core/properties/cameraproperty.h"
#include "voreen/core/ports/renderport.h"
//gi
#include "voreen/qt/networkeditor/graphicitems/core/portownergraphicsitem.h"
#include "voreen/qt/networkeditor/graphicitems/utils/propertylistbuttongraphicsitem.h"
#include "voreen/qt/networkeditor/graphicitems/core/propertygraphicsitem.h"


namespace voreen {

PortOwnerLinkArrowGraphicsItem::PortOwnerLinkArrowGraphicsItem(PropertyListButtonGraphicsItem* sourceItem, PropertyListButtonGraphicsItem* destinationItem)
    : ConnectionBaseGraphicsItem(sourceItem, destinationItem)
{
    setParent(sourceItem);
    sourceItem->scene()->addItem(this);
    setZValue(ZValuesPortOwnerLinkArrowGraphicsItem);
    updateNWELayerAndCursor();
}

PortOwnerLinkArrowGraphicsItem::~PortOwnerLinkArrowGraphicsItem() {}

//---------------------------------------------------------------------------------------------------------------
//                  getter and setter
//---------------------------------------------------------------------------------------------------------------
PropertyListButtonGraphicsItem* PortOwnerLinkArrowGraphicsItem::getSourceItem() const {
PropertyListButtonGraphicsItem* i = dynamic_cast<PropertyListButtonGraphicsItem*>(sourceItem_);
    tgtAssert(i, "sourceItem was of wrong type");
    return i;
}


PropertyListButtonGraphicsItem* PortOwnerLinkArrowGraphicsItem::getDestinationItem() const {
    if (destinationItem_) {
        PropertyListButtonGraphicsItem* i = dynamic_cast<PropertyListButtonGraphicsItem*>(destinationItem_);
        tgtAssert(i, "destinationItem was of wrong type");
        return i;
    }
    else
        return 0;
}

QPointF PortOwnerLinkArrowGraphicsItem::getSourcePoint() const{
    return getSourceItem()->dockingPoint();
}

QPointF PortOwnerLinkArrowGraphicsItem::getDestinationPoint() const{
    if(getDestinationItem())
        return getDestinationItem()->dockingPoint();
    else
        return destinationPoint_;
}

void PortOwnerLinkArrowGraphicsItem::setVisible(bool visible) {
    if(!visible){
        QGraphicsItem::setVisible(false);
    } else {
        switch (currentLayer()) {
        case NetworkEditorLayerDataFlow:
        case NetworkEditorLayerPortSizeLinking:
            QGraphicsItem::setVisible(false);
            break;
        case NetworkEditorLayerGeneralLinking: {
            QGraphicsItem::setVisible(false);
            if(getDestinationItem()){
                std::vector<const Property*> srcVec;
                std::vector<const Property*> dstVec;
                foreach (PropertyGraphicsItem* propItem, dynamic_cast<PortOwnerGraphicsItem*>(getSourceItem()->parent())->getPropertyList()->getAllPropertyItems()){
                    if(getNetworkEditor()->cameraLinksHidden() && dynamic_cast<const CameraProperty*>(propItem->getProperty())){
                        continue;
                    }
                    if(getNetworkEditor()->portSizeLinksHidden() && (dynamic_cast<const RenderSizeOriginProperty*>(propItem->getProperty())
                                                                 ||  dynamic_cast<const RenderSizeReceiveProperty*>(propItem->getProperty()))){
                        continue;
                    }
                    srcVec.push_back(propItem->getProperty());
                }
                foreach (PropertyGraphicsItem* propItem, dynamic_cast<PortOwnerGraphicsItem*>(getDestinationItem()->parent())->getPropertyList()->getAllPropertyItems()){
                    if(getNetworkEditor()->cameraLinksHidden() && dynamic_cast<const CameraProperty*>(propItem->getProperty())){
                        continue;
                    }
                    if(getNetworkEditor()->portSizeLinksHidden() && (dynamic_cast<const RenderSizeOriginProperty*>(propItem->getProperty())
                                                                 ||  dynamic_cast<const RenderSizeReceiveProperty*>(propItem->getProperty()))){
                        continue;
                    }
                    dstVec.push_back(propItem->getProperty());
                }
                for(size_t i = 0; i < srcVec.size(); ++i) {
                    for(size_t j = 0; j < dstVec.size(); ++j) {
                        if(srcVec[i]->isLinkedWith(dstVec[j],false) || dstVec[j]->isLinkedWith(srcVec[i])) { //was true, needed?
                            QGraphicsItem::setVisible(true);
                            goto EndOfLoop;
                        }
                    }
                }
                EndOfLoop:;
            }
            if(getIsHovered())
                QGraphicsItem::setVisible(true);
            } break;
        case NetworkEditorLayerCameraLinking: {
            QGraphicsItem::setVisible(false);
            if(getDestinationItem()) {
                std::vector<CameraProperty*> srcVec;
                std::vector<CameraProperty*> dstVec;

                foreach (Processor* processor, dynamic_cast<PortOwnerGraphicsItem*>(getSourceItem()->parent())->getProcessors()){
                    std::vector<CameraProperty*> hVec = processor->getPropertiesByType<CameraProperty>();
                    srcVec.insert(srcVec.end(), hVec.begin(), hVec.end());
                }
                foreach (Processor* processor, dynamic_cast<PortOwnerGraphicsItem*>(getDestinationItem()->parent())->getProcessors()){
                    std::vector<CameraProperty*> hVec = processor->getPropertiesByType<CameraProperty>();
                    dstVec.insert(dstVec.end(), hVec.begin(), hVec.end());
                }

                for(size_t i = 0; i < srcVec.size(); ++i)
                    for(size_t j = 0; j < dstVec.size(); ++j)
                        if(srcVec[i]->isLinkedWith(dstVec[j], true)){
                            QGraphicsItem::setVisible(true);
                            break;
                        }
            }
            if(getIsHovered())
                QGraphicsItem::setVisible(true);
            } break;
        default:
            tgtAssert(false, "shouldn't get here");
        }
        if((getSourceItem() && !getSourceItem()->isVisible()) ||
           (getDestinationItem() && !getDestinationItem()->isVisible()))
            QGraphicsItem::setVisible(false);
    }
}

//---------------------------------------------------------------------------------------------------------------
//                  nwebasegraphicsitem functions
//---------------------------------------------------------------------------------------------------------------
QRectF PortOwnerLinkArrowGraphicsItem::boundingRect() const {
    return currentStyle()->PortOwnerLinkArrowGI_boundingRect(this);
}

QPainterPath PortOwnerLinkArrowGraphicsItem::shape() const {
    return currentStyle()->PortOwnerLinkArrowGI_shape(this);
}

void PortOwnerLinkArrowGraphicsItem::initializePaintSettings() {
    currentStyle()->PortOwnerLinkArrowGI_initializePaintSettings(this);
}


void PortOwnerLinkArrowGraphicsItem::mainPaint(QPainter* painter, const QStyleOptionGraphicsItem* option, QWidget* widget, NWEItemSettings& setting) {
    currentStyle()->PortOwnerLinkArrowGI_paint(this, painter, option, widget, setting);
}

void PortOwnerLinkArrowGraphicsItem::updateNWELayerAndCursor() {
    switch(currentCursorMode()) {
    case NetworkEditorCursorSelectMode:
        switch (currentLayer()) {
        case NetworkEditorLayerDataFlow:
        case NetworkEditorLayerPortSizeLinking:
            setVisible(false);
            break;
        case NetworkEditorLayerGeneralLinking:
            setVisible(true);
            break;
        case NetworkEditorLayerCameraLinking:
            setVisible(true);
            break;
        default:
            tgtAssert(false, "shouldn't get here");
        }
        break;
    case NetworkEditorCursorMoveMode:
        switch (currentLayer()) {
        case NetworkEditorLayerDataFlow:
        case NetworkEditorLayerPortSizeLinking:
            setVisible(false);
            break;
        case NetworkEditorLayerGeneralLinking:
            setVisible(true);
            break;
        case NetworkEditorLayerCameraLinking:
            setVisible(true);
            break;
        default:
            tgtAssert(false, "shouldn't get here");
        }
        break;
    default:
        tgtAssert(false,"should not get here");
        break;
    }
}

} // namespace */
