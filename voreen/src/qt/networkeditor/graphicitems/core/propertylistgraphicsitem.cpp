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

#include "voreen/qt/networkeditor/graphicitems/core/propertylistgraphicsitem.h"
#include "voreen/qt/networkeditor/styles/nwestyle_classic.h"

//core
#include "voreen/core/properties/cameraproperty.h"
#include "voreen/core/ports/renderport.h"
//gi
#include "voreen/qt/networkeditor/graphicitems/core/propertyownergraphicsitem.h"
#include "voreen/qt/networkeditor/graphicitems/core/propertygraphicsitem.h"
#include "voreen/qt/networkeditor/graphicitems/connections/propertylinkarrowgraphicsitem.h"

#include <QBrush>



namespace voreen {

PropertyListGraphicsItem::PropertyListGraphicsItem(PropertyOwnerGraphicsItem* parent)
    : NWEBaseGraphicsItem(parent->getNetworkEditor())
    , parentPropertyOwnerItem_(parent)
    , headerLabelItem_(0)
    , superPropertyListItem_(0)
    , friendPropertyListItem_(0)
    , currentHeaderMode_(NO_HEADER)
    , currentPropertyMode_(SHOW_LINKED_PROPERTIES)
    , currentPropertyModifier_(HIDE_NO_PROPERTIES)
    , currentLinkArrowMode_(SHOW_NO_ARROWS)
    , currentHeight_(0.f)
    , currentWidth_(0.f)
    , isVisibleInEditor_(false)
{
    tgtAssert(parent != 0, "passed null pointer");
    //set z value
    setZValue(ZValuesPropertyListGraphicsItem);
    //set Qt parent
    setParent(parent);
    //add this to scene of parent
    parent->getNetworkEditor()->scene()->addItem(this);
    //generate label for header
    headerLabelItem_ = new QGraphicsTextItem(this);
    headerLabelItem_->setParentItem(this);
    //because of the dynamic behavior of the PropertyList, we have to set label text now
    initializePaintSettings();
    paintHasBeenInitialized_ = true;
    //create child items
    createChildItems();
}

PropertyListGraphicsItem::~PropertyListGraphicsItem() {
    //delete child items
    //deleteChildItems(); //is down automaticaly because the PropertyGraphicItems are children of this
    //set list parent of sub lists to zero
    foreach(PropertyListGraphicsItem* item,subPropertyListItems_) {
        removeSubPropertyListItem(item);
    }
}

//---------------------------------------------------------------------------------------------------------------
//                  nwebasegraphicsitem functions
//---------------------------------------------------------------------------------------------------------------
QRectF PropertyListGraphicsItem::boundingRect() const {
    return currentStyle()->PropertyListGI_boundingRect(this);
}

QPainterPath PropertyListGraphicsItem::shape() const{
    return currentStyle()->PropertyListGI_shape(this);
}

void PropertyListGraphicsItem::initializePaintSettings() {
    currentStyle()->PropertyListGI_initializePaintSettings(this);
    updateHeaderLabelItem();
}

void PropertyListGraphicsItem::mainPaint(QPainter* painter, const QStyleOptionGraphicsItem* option, QWidget* widget, NWEItemSettings& setting) {
    currentStyle()->PropertyListGI_paint(this,painter, option, widget, setting);
}


void PropertyListGraphicsItem::createChildItems(){
        addProperties(stdVectorToQList(parentPropertyOwnerItem_->getPropertyOwner()->getProperties()));
}

void PropertyListGraphicsItem::deleteChildItems(){
    foreach (PropertyGraphicsItem* item, propertyItems_)
        delete item;

}

void PropertyListGraphicsItem::layoutChildItems(){
    //prepare geometry change
    prepareGeometryChange();
    //set label position
    qreal x = (boundingRect().width() - headerLabelItem_->boundingRect().width()) / 2.0;
    qreal y = (currentHeaderMode_  - headerLabelItem_->boundingRect().height()) / 2.0;
    headerLabelItem_->setPos(x, y);
    //set current height depending on mode and number of visible properties
    if(getVisiblePropertyItems().empty()){
        headerLabelItem_->setVisible(false);
        currentHeight_ = 0.f;
    } else {
        if(currentHeaderMode_ != NO_HEADER)
            headerLabelItem_->setVisible(true);
        currentHeight_ = currentHeaderMode_;
    }
    //set position of all properties
    foreach(PropertyGraphicsItem* propItem, propertyItems_){
        if(propItem->isVisible()){
            propItem->setPos(0,currentHeight_);
            currentHeight_ += propItem->boundingRect().height();
        }
    }
    //set position of all sub lists
    foreach(PropertyListGraphicsItem* listItem, subPropertyListItems_){
            currentHeight_ += listItem->setChildListPosition(mapToItem(listItem->parentItem(),QPointF(0,currentHeight_)));
    }
}

//---------------------------------------------------------------------------------------------------------------
//                  getter and setter
//---------------------------------------------------------------------------------------------------------------
PropertyOwnerGraphicsItem* PropertyListGraphicsItem::getParentPropertyOwnerItem() const {
    return parentPropertyOwnerItem_;
}

void PropertyListGraphicsItem::setSuperPropertyListItem(PropertyListGraphicsItem* plgi) {
    superPropertyListItem_ = plgi;
}

PropertyListGraphicsItem* PropertyListGraphicsItem::getSuperPropertyListItem() const {
    return superPropertyListItem_;
}

const QGraphicsTextItem* PropertyListGraphicsItem::getHeaderLabelItem() const {
    return headerLabelItem_;
}

QGraphicsTextItem* PropertyListGraphicsItem::getHeaderLabelItem() {
    return headerLabelItem_;
}

void PropertyListGraphicsItem::updateHeaderLabelItem() {
    QString labelText = getParentPropertyOwnerItem()->getGuiName();
    getHeaderLabelItem()->setHtml(labelText);
    findNewListWidth();
}

PropertyGraphicsItem* PropertyListGraphicsItem::getPropertyItem(const Property* prop) {
    foreach(PropertyGraphicsItem* item, propertyItems_){
        if(item->getProperty() == prop)
            return item;
    }
    return 0;
}

QList<PropertyGraphicsItem*> PropertyListGraphicsItem::getPropertyItems() {
    return propertyItems_;
}

const QList<PropertyGraphicsItem*> PropertyListGraphicsItem::getPropertyItems() const{
    return propertyItems_;
}

QList<PropertyGraphicsItem*> PropertyListGraphicsItem::getVisiblePropertyItems() {
    QList<PropertyGraphicsItem*> list;
    foreach(PropertyGraphicsItem* item, propertyItems_)
        if(item->isVisible())
            list.append(item);
    return list;
}

const QList<PropertyGraphicsItem*> PropertyListGraphicsItem::getVisiblePropertyItems() const{
    QList<PropertyGraphicsItem*> list;
    foreach(PropertyGraphicsItem* item, propertyItems_)
        if(item->isVisible())
            list.append(item);
    return list;
}

QList<PropertyGraphicsItem*> PropertyListGraphicsItem::getAllPropertyItems() {
    QList<PropertyGraphicsItem*> list(propertyItems_);
    foreach(PropertyListGraphicsItem* item, getSubPropertyListItems())
        list.append(item->getPropertyItems());
    return list;
}

const QList<PropertyGraphicsItem*> PropertyListGraphicsItem::getAllPropertyItems() const{
    QList<PropertyGraphicsItem*> list(propertyItems_);
    foreach(PropertyListGraphicsItem* item, getSubPropertyListItems())
        list.append(item->getPropertyItems());
    return list;
}

QList<PropertyGraphicsItem*> PropertyListGraphicsItem::getAllVisiblePropertyItems() {
    QList<PropertyGraphicsItem*> list = getVisiblePropertyItems();
    foreach(PropertyListGraphicsItem* item, getSubPropertyListItems())
        list.append(item->getVisiblePropertyItems());
    return list;
}

const QList<PropertyGraphicsItem*> PropertyListGraphicsItem::getAllVisiblePropertyItems() const {
    QList<PropertyGraphicsItem*> list = getVisiblePropertyItems();
    foreach(PropertyListGraphicsItem* item, getSubPropertyListItems())
        list.append(item->getVisiblePropertyItems());
    return list;
}

PropertyListGraphicsItem::LinkArrowVisibleMode PropertyListGraphicsItem::getLinkArrowVisibleMode() const {
    return currentLinkArrowMode_;
}

PropertyListGraphicsItem::PropertyVisibleMode PropertyListGraphicsItem::getPropertyVisibleMode() const {
    return currentPropertyMode_;
}

PropertyListGraphicsItem::PropertyVisibleModeModifier PropertyListGraphicsItem::getPropertyVisibleModeModifier() const {
    return currentPropertyModifier_;
}
PropertyListGraphicsItem::HeaderMode PropertyListGraphicsItem::getHeaderMode() const {
    return currentHeaderMode_;
}

void PropertyListGraphicsItem::setLinkArrowMode(LinkArrowVisibleMode lavm, PropertyListGraphicsItem* friendList){
    currentLinkArrowMode_= lavm;
    friendPropertyListItem_ = friendList;
    //switch properties into right mode
    foreach(PropertyGraphicsItem* item, propertyItems_) {
        switch(lavm) {
        case SHOW_NO_ARROWS:
            item->setShowLinkArrows(false);
            break;
        case SHOW_ALL_ARROWS:
        case SHOW_FRIEND_LIST_ARROWS:
            item->setShowLinkArrows(true);
            break;
        }
        item->setPropertyListLinkedTo(friendList);
    }
    //switch sub lists into right mode
    foreach(PropertyListGraphicsItem* listItem, subPropertyListItems_){
        listItem->setLinkArrowMode(lavm, friendList);
    }
    //layout child items
    layoutChildItems();
}

void PropertyListGraphicsItem::setPropertyVisibleMode(PropertyVisibleMode pvm){
    currentPropertyMode_= pvm;
    //switch sub lists
    foreach(PropertyListGraphicsItem* listItem, subPropertyListItems_){
        listItem->setPropertyVisibleMode(pvm);
    }
    //layout child items
    layoutChildItems();
}

void PropertyListGraphicsItem::setPropertyVisibleModifier(PropertyVisibleModeModifier pvmm){
    currentPropertyModifier_= pvmm;
    //switch sub lists
    foreach(PropertyListGraphicsItem* listItem, subPropertyListItems_){
        listItem->setPropertyVisibleModifier(pvmm);
    }
    //layout child items
    layoutChildItems();
}

void PropertyListGraphicsItem::setHeaderMode(HeaderMode hm){
    if(hm == currentHeaderMode_)
        return;
    //set header mode
    currentHeaderMode_= hm;
    //layoutChildItems
    layoutChildItems();
}

bool PropertyListGraphicsItem::getIsVisibleInEditor() const {
    return isVisibleInEditor_;
}

void PropertyListGraphicsItem::setIsVisibleInEditor(bool b) {
    isVisibleInEditor_ = b;
}

QList<PropertyListGraphicsItem*> PropertyListGraphicsItem::getSubPropertyListItems() {
    return subPropertyListItems_;
}

const QList<PropertyListGraphicsItem*> PropertyListGraphicsItem::getSubPropertyListItems() const {
    return subPropertyListItems_;
}

void PropertyListGraphicsItem::setPosition(const QPointF& pos) {
    QGraphicsItem::setPos(pos);
    layoutChildItems();
}

void PropertyListGraphicsItem::setPosition(qreal x, qreal y) {
    QGraphicsItem::setPos(x,y);
    layoutChildItems();
}

qreal PropertyListGraphicsItem::getCurrentWidth() const {
    return currentWidth_;
}

void PropertyListGraphicsItem::setScene(QGraphicsScene* scene) {
    //change scene of this
    this->scene()->removeItem(this);
    scene->addItem(this);
    //change scene of all arrow items
    foreach(PropertyGraphicsItem* propItem, getPropertyItems()){
        foreach(PropertyLinkArrowGraphicsItem* arrow, propItem->getSourceLinkList()) {
            arrow->scene()->removeItem(arrow);
            scene->addItem(arrow);
        }
        foreach(PropertyLinkArrowGraphicsItem* arrow, propItem->getDestinationLinkList()) {
            arrow->scene()->removeItem(arrow);
            scene->addItem(arrow);
        }
    }
    //change scene of all sub lists
    foreach(PropertyListGraphicsItem* item, subPropertyListItems_){
        item->setScene(scene);
    }
};

//---------------------------------------------------------------------------------------------------------------
//                  utils
//---------------------------------------------------------------------------------------------------------------
void PropertyListGraphicsItem::addProperties(const QList<Property*>& items) {
    foreach(Property* prop, items) {
        PropertyGraphicsItem* item = new PropertyGraphicsItem(prop,this);
        propertyItems_.append(item);
        connect(item,SIGNAL(linkProperties(PropertyGraphicsItem*,PropertyGraphicsItem*)),this,SIGNAL(linkProperties(PropertyGraphicsItem*,PropertyGraphicsItem*)));
    }
    findNewListWidth();
}

void PropertyListGraphicsItem::addSubPropertyListItem(PropertyListGraphicsItem* list) {
    subPropertyListItems_.append(list);
    list->setSuperPropertyListItem(this);
    list->setHeaderMode(SMALL_HEADER);
    connect(list,SIGNAL(linkProperties(PropertyGraphicsItem*,PropertyGraphicsItem*)),this,SIGNAL(linkProperties(PropertyGraphicsItem*,PropertyGraphicsItem*)));
    findNewListWidth();
}

void PropertyListGraphicsItem::removeSubPropertyListItem(PropertyListGraphicsItem* list) {
    subPropertyListItems_.removeOne(list);
    list->setSuperPropertyListItem(0);
    list->setHeaderMode(BIG_HEADER);
    //since it is 'free' the width must be calculated new
    list->findNewListWidth();
    //same for this
    findNewListWidth();
}

float PropertyListGraphicsItem::setChildListPosition(QPointF point){
    setPos(point);
    layoutChildItems();
    return currentHeight_;
}

void PropertyListGraphicsItem::setVisible(bool visible) {
    switch(currentPropertyMode_){
    case SHOW_ALL_PROPERTIES:
        foreach (PropertyGraphicsItem* propertyItem, propertyItems_) {
            propertyItem->setVisible(visible);
        }
        break;
    case SHOW_LINKED_PROPERTIES:
        foreach (PropertyGraphicsItem* propertyItem, propertyItems_) {
            if(visible == false)
                propertyItem->setVisible(false);
            else {
                if(!propertyItem->getSourceLinkList().empty() || !propertyItem->getDestinationLinkList().empty())
                    propertyItem->setVisible(true);
                else
                    propertyItem->setVisible(false);
            }
        }
        break;
    case SHOW_CAMERA_PROPERTIES:
        foreach (PropertyGraphicsItem* propertyItem, propertyItems_) {
            if(visible == false)
                propertyItem->setVisible(false);
            else {
                if(dynamic_cast<const CameraProperty*>(propertyItem->getProperty()))
                    propertyItem->setVisible(true);
                else
                    propertyItem->setVisible(false);
            }
        }
        break;
    default:
        tgtAssert(false,"should not get here");
        break;
    }

    switch(currentPropertyModifier_) {
    case HIDE_CAMERA_AND_SIZE_PROPERTIES:
    case HIDE_CAMERA_PROPERTIES:
        foreach (PropertyGraphicsItem* propertyItem, propertyItems_) {
            if(dynamic_cast<CameraProperty*>(const_cast<Property*>(propertyItem->getProperty())))
                propertyItem->setVisible(false);
        }
        if(currentPropertyModifier_ == HIDE_CAMERA_PROPERTIES)
            break; //stop, if just camera porperties are hidden
    case HIDE_SIZE_PROPERTIES:
        foreach (PropertyGraphicsItem* propertyItem, propertyItems_) {
            if(dynamic_cast<RenderSizeOriginProperty*>(const_cast<Property*>(propertyItem->getProperty())) ||
               dynamic_cast<RenderSizeReceiveProperty*>(const_cast<Property*>(propertyItem->getProperty())))
                propertyItem->setVisible(false);
        }
        break;
    default:
        break;
    }

    if(currentHeaderMode_ == NO_HEADER || getVisiblePropertyItems().empty())
        headerLabelItem_->setVisible(false);
    else
        headerLabelItem_->setVisible(visible);

    //go into sub lists
    foreach(PropertyListGraphicsItem* listItem, subPropertyListItems_){
        listItem->setVisible(visible);
    }

    QGraphicsItem::setVisible(visible);
    if(visible == true)
        layoutChildItems();
}



void PropertyListGraphicsItem::prepareGeometryChange() {
    //prepare this
    NWEBaseGraphicsItem::prepareGeometryChange();
    //prepare properties
    foreach(PropertyGraphicsItem* item, getPropertyItems())
        item->prepareGeometryChange();
    //prepare sub lists
    foreach(PropertyListGraphicsItem* item, getSubPropertyListItems())
        item->prepareGeometryChange();
}

//---------------------------------------------------------------------------------------------------------------
//                  width management
//---------------------------------------------------------------------------------------------------------------
void PropertyListGraphicsItem::findNewListWidth() {
    qreal width = getHeaderLabelItem()->boundingRect().width();
    foreach (PropertyGraphicsItem* propertyItem, getPropertyItems()) {
        width = std::max(width, propertyItem->boundingRect().width());
    }
    foreach(PropertyListGraphicsItem* listItem, getSubPropertyListItems()){
        width = std::max(width,listItem->getCurrentWidth());
    }
    currentWidth_ = width;

    if(superPropertyListItem_ && width != superPropertyListItem_->getCurrentWidth())
        superPropertyListItem_->findNewListWidth();
    else {
        resizeListWidth(width);
        layoutChildItems();
    }
}

void PropertyListGraphicsItem::resizeListWidth(qreal width){
    currentWidth_ = width;
    foreach(PropertyGraphicsItem* item, propertyItems_){
        item->prepareGeometryChange();
        item->resetSize(width);
    }
    foreach(PropertyListGraphicsItem* item, subPropertyListItems_){
        item->prepareGeometryChange();
        item->resizeListWidth(width);
    }
}

} // namespace
