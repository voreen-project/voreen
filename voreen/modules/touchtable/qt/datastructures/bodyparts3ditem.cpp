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

#include "bodyparts3ditem.h"
#include <QStringList>
#include <sstream>

namespace voreen{


BodyParts3DItem::BodyParts3DItem(BodyParts3DItem *parent)
{
    parentItem_ = parent;
}

BodyParts3DItem::~BodyParts3DItem()
{

}

int BodyParts3DItem::getID()const{
    return data_->id;
}

void BodyParts3DItem::setID(int iD){
    data_->id = iD;
}

tgt::vec4 BodyParts3DItem::getColor()const{
    return data_->color;
}

void BodyParts3DItem::setColor(tgt::vec4 c){
    data_->color = c;
}

bool BodyParts3DItem::getCheckState()const{
    return data_->checkState;
}

void BodyParts3DItem::setCheckState(bool c){
    data_->checkState = c;
}

std::string BodyParts3DItem::getDescription()const{
    return data_->description;
}

void BodyParts3DItem::setDescription(std::string desc){
    data_->description = desc;
}

BodyPart* BodyParts3DItem::getData()const{
    return data_;
}

void BodyParts3DItem::setData(BodyPart* dataNode){
    data_ = dataNode;
    if(data_->parts.size() != 0){
        for(int i = 0; i < data_->parts.size() ; ++i){
            childItems_.push_back(new BodyParts3DItem(this));
            childItems_.value(i)->setData(data_->parts[i]);
        }
    }
}

BodyParts3DItem *BodyParts3DItem::getChild(int number)
{
    return childItems_.value(number);
}

BodyParts3DItem *BodyParts3DItem::getParent()
{
    return parentItem_;
}

int BodyParts3DItem::getChildNumber() const
{
    if (parentItem_)
        return parentItem_->childItems_.indexOf(const_cast<BodyParts3DItem*>(this));

    return 0;
}

int BodyParts3DItem::getChildCount() const
{
    return childItems_.count();
}

void BodyParts3DItem::appendChildren(int position)
{
    childItems_.push_back(new BodyParts3DItem(this));
}

}//namespace voren
