/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany,                        *
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

#include "bodyparts3dmodel.h"
#include "tgt/filesystem.h"
#include "voreen/core/utils/stringutils.h"
#include <QtGui>

namespace voreen{

BodyParts3DModel::BodyParts3DModel(const std::string& path, bool color, BodyPartsVersion version, QObject *parent)
    : QAbstractItemModel(parent)
{
    rootData_ << "Description" << "Render" << "Color";
    rootItem_ = new BodyParts3DItem();
    if(version == Version4_0 ){
        catalog_ = new BodyParts3DCatalog40(path, color);
        for(int i = 0; i < catalog_->getSize() ; i++){
            BodyPart* part = catalog_->getBodyPart(i);
            if(part->partOf == 0){
                rootItem_->appendChildren(1);
                rootItem_->getChild(rootItem_->getChildCount()-1)->setData(part);
            }
        }
    }else{
        catalog_ = new BodyParts3DCatalog30(path, color);
        BodyPart* head = catalog_->getHead();
        rootItem_->appendChildren(1);
        rootItem_->getChild(rootItem_->getChildCount()-1)->setData(head);
    }
}

BodyParts3DModel::BodyParts3DModel(QObject *parent)
    : QAbstractItemModel(parent)
{
    rootData_ << "Description" << "Render" << "Color";
    rootItem_ = new BodyParts3DItem();
}

BodyParts3DModel::~BodyParts3DModel()
{
    delete rootItem_;
}


QVariant BodyParts3DModel::data(const QModelIndex &index, int role) const
{
    if (!index.isValid())
        return QVariant();

    int col = index.column();
    BodyParts3DItem* item = static_cast<BodyParts3DItem*>(index.internalPointer());
    switch(role){
    case Qt::DisplayRole:
        switch(col){
        case 0:
            return QString(item->getDescription().c_str());
            break;
        case 1:
            return "";
            break;
        case 2:
            return "       ";
            break;
        }
        break;
    case Qt::BackgroundRole:
        if(col == 2){
            tgt::vec4 c = item->getColor();
            c*= 255.0f;
            return QColor(c.x,c.y,c.z,c.w);
        }
        break;
    case Qt::CheckStateRole:
        if(col == 1){
            if(item->getCheckState())
                return Qt::Checked;
            else
                return Qt::Unchecked;
        }
        break;
    }
    return QVariant();
}

QModelIndex BodyParts3DModel::index(int row, int column, const QModelIndex &parent) const
{
    if (parent.isValid() && parent.column() != 0)
        return QModelIndex();

    BodyParts3DItem *parentItem = getItem(parent);

    BodyParts3DItem *childItem = parentItem->getChild(row);
    if (childItem)
        return createIndex(row, column, childItem);
    else
        return QModelIndex();
}

QModelIndex BodyParts3DModel::parent(const QModelIndex &index) const
{
    if (!index.isValid())
        return QModelIndex();

    BodyParts3DItem *childItem = getItem(index);
    BodyParts3DItem *parentItem = childItem->getParent();

    if (parentItem == rootItem_)
        return QModelIndex();

    return createIndex(parentItem->getChildNumber(), 0, parentItem);
}

int BodyParts3DModel::rowCount(const QModelIndex &parent) const
{
    BodyParts3DItem *parentItem = getItem(parent);

    return parentItem->getChildCount();
}

int BodyParts3DModel::columnCount(const QModelIndex & ) const
{
    return 3;
}


QVariant BodyParts3DModel::headerData(int section, Qt::Orientation orientation,
                               int role) const
{
    if (orientation == Qt::Horizontal && role == Qt::DisplayRole)
        return rootData_[section];

    return QVariant();
}

Qt::ItemFlags BodyParts3DModel::flags(const QModelIndex &index) const
{
    if (!index.isValid())
        return 0;

    int col = index.column();

    switch(col){
    case 0:
        return Qt::ItemFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled);
    case 1:
        return Qt::ItemFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled | Qt::ItemIsUserCheckable);
    case 2:
        return Qt::ItemIsEnabled;
    }
    return Qt::ItemIsEnabled;
}

bool BodyParts3DModel::setData(const QModelIndex &index, const QVariant &value,
                        int role)
{
    if (!index.isValid())
        return false;

    int col = index.column();

    BodyParts3DItem *item = getItem(index);
    if(role == Qt::CheckStateRole && col == 1){
        checkAllChildren(value.toBool(),item);
        emit dataChanged(index, index.sibling(index.row()+1,col));
    }
    return false;
}


void BodyParts3DModel::checkAllChildren(bool chk, BodyParts3DItem* parent){
    parent->setCheckState(chk);
    for(int i = 0 ; i < parent->getChildCount() ; ++i){
        checkAllChildren(chk,parent->getChild(i));
    }
}

void BodyParts3DModel::exportCatalog(const std::string& fileName){
    catalog_->exportCatalog(fileName);
}

std::vector<std::string> BodyParts3DModel::getFileNames(){
    return catalog_->fileNameOutput();
}
std::vector<tgt::vec4> BodyParts3DModel::getColors(){
    return catalog_->colorOutput();
}
std::vector<std::string> BodyParts3DModel::getDescriptions(){
    return catalog_->descriptionOutput();
}


BodyParts3DItem *BodyParts3DModel::getItem(const QModelIndex &index) const
{
    if (index.isValid()) {
        BodyParts3DItem *item = static_cast<BodyParts3DItem*>(index.internalPointer());
        if (item)
            return item;
    }
    return rootItem_;
}

}//namespace voren
