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

#include "voreen/qt/networkeditor/graphicitems/connections/connectionbasegraphicsitem.h"
#include "tgt/assert.h"

namespace voreen {

ConnectionBaseGraphicsItem::ConnectionBaseGraphicsItem(NWEBaseGraphicsItem* sourceItem, NWEBaseGraphicsItem* destinationItem)
    : NWEBaseGraphicsItem(sourceItem->getNetworkEditor())
    , sourceArrowHeadDirection_(ArrowHeadDirectionNone)
    , destinationArrowHeadDirection_(ArrowHeadDirectionNone)
    , sourceItem_(sourceItem)
    , destinationItem_(destinationItem)
    , sourcePoint_(QPointF())
    , destinationPoint_(QPointF())
    , colorConnectableMode_(CCM_DEFAULT)
    , isHovered_(false)
{
    tgtAssert(sourceItem,"null pointer passed");
    if(!destinationItem)
        setDestinationPoint(getSourcePoint());
    setParent(sourceItem);
}

ConnectionBaseGraphicsItem::~ConnectionBaseGraphicsItem() {
}

//---------------------------------------------------------------------------------------------------------------
//                  getter and setter
//---------------------------------------------------------------------------------------------------------------
NWEBaseGraphicsItem* ConnectionBaseGraphicsItem::getSourceItem() const {
    return sourceItem_;
}

void ConnectionBaseGraphicsItem::setSourceItem(NWEBaseGraphicsItem* item) {
    tgtAssert(item, "passed null pointer");
    setParent(item);
    prepareGeometryChange();
    sourceItem_ = item;
    update();
}

NWEBaseGraphicsItem* ConnectionBaseGraphicsItem::getDestinationItem() const {
    return destinationItem_;
}

void ConnectionBaseGraphicsItem::setDestinationItem(NWEBaseGraphicsItem* item) {
    prepareGeometryChange();
    destinationItem_ = item;
    update();
}

void ConnectionBaseGraphicsItem::setSourcePoint(const QPointF& sourcePoint) {
    prepareGeometryChange();
    sourceItem_ = 0;
    sourcePoint_ = sourcePoint;
    update();
}

void ConnectionBaseGraphicsItem::setDestinationPoint(const QPointF& destinationPoint) {
    prepareGeometryChange();
    destinationItem_ = 0;
    destinationPoint_ = destinationPoint;
    update();
}

QPointF ConnectionBaseGraphicsItem::getSourcePoint() const {
    if (sourceItem_)
        return mapFromItem(sourceItem_, 0, 0);
    else
        return mapFromScene(sourcePoint_);
}

QPointF ConnectionBaseGraphicsItem::getDestinationPoint() const {
    if (destinationItem_)
        return mapFromItem(destinationItem_, 0, 0);
    else
        return mapFromScene(destinationPoint_);
}

ConnectionBaseGraphicsItem::ArrowHeadDirection ConnectionBaseGraphicsItem::getSourceArrowHeadDirection() const{
    return sourceArrowHeadDirection_;
}

ConnectionBaseGraphicsItem::ArrowHeadDirection ConnectionBaseGraphicsItem::getDestinationArrowHeadDirection() const{
    return destinationArrowHeadDirection_;
}

void ConnectionBaseGraphicsItem::setSourceArrowHeadDirection(ArrowHeadDirection ahd) {
    sourceArrowHeadDirection_ = ahd;
}

void ConnectionBaseGraphicsItem::setDestinationArrowHeadDirection(ArrowHeadDirection ahd) {
    destinationArrowHeadDirection_ = ahd;
}

bool ConnectionBaseGraphicsItem::getIsHovered() const{
    return isHovered_;
}

void ConnectionBaseGraphicsItem::setIsHovered(bool b) {
    isHovered_ = b;
}

} // namespace voreen
