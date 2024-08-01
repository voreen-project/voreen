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

#ifndef VRN_BODYPARTS3DITEM_H
#define VRN_BODYPARTS3DITEM_H

#include "voreen/core/voreenapplication.h"
#include "bodypartscatalogbase.h"
#include <QVector>
#include <QVariant>
#include <QString>
#include <QStandardItemModel>
#include <QList>

namespace voreen{

/**
 * This Class Represents a BodyPart in the BodyParts3DModel
 */
class BodyParts3DItem : public QStandardItem
{
public:
    BodyParts3DItem(BodyParts3DItem *parent = 0);
    ~BodyParts3DItem();

    int getID()const;
    void setID(int iD);
    tgt::vec4 getColor()const;
    void setColor(tgt::vec4 c);
    bool getCheckState()const;
    void setCheckState(bool c);
    std::string getDescription()const;
    void setDescription(std::string desc);
    BodyPart* getData()const;
    /**
     * This function calls itself recursively after the catalog has been build to construct the Item Model
     */
    void setData(BodyPart* dataNode);
    BodyParts3DItem* getChild(int number);
    BodyParts3DItem *getParent();
    /**
     * returns the index of current BodyPart
     */
    int getChildNumber() const;

    int getChildCount() const;
    /**
     * inserts new BodyParts3DItems in childitems
     */
    void appendChildren(int count);

private:
    BodyPart* data_;
    QList<BodyParts3DItem*> childItems_;
    BodyParts3DItem *parentItem_;
};

}//namespace voren

#endif //VRN_BODYPARTS3DITEM_H
