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

#ifndef VRN_BODYPARTS3DMODEL_H
#define VRN_BODYPARTS3DMODEL_H

#include <QAbstractItemModel>
#include <QModelIndex>
#include <QVariant>
#include "bodyparts3ditem.h"
#include "bodypartscatalogbase.h"
#include "bodypartscatalog30.h"
#include "bodypartscatalog40.h"
#include "modules/touchtable/processors/bodyparts3d/bodyparts3dsource.h"

namespace voreen{

class BodyParts3DItem;

class BodyParts3DModel : public QAbstractItemModel
{
    Q_OBJECT

public:
    BodyParts3DModel(const std::string& path, bool color, BodyPartsVersion version, QObject *parent = 0);
    BodyParts3DModel(QObject *parent = 0);
    ~BodyParts3DModel();

    QVariant data(const QModelIndex &index, int role) const; // = 0
    QModelIndex index(int row, int column, const QModelIndex &parent = QModelIndex()) const; // = 0
    QModelIndex parent(const QModelIndex &index) const; // = 0
    int rowCount(const QModelIndex &parent = QModelIndex()) const; // = 0
    int columnCount(const QModelIndex &parent = QModelIndex()) const; // = 0

    QVariant headerData(int section, Qt::Orientation orientation, int role = Qt::DisplayRole) const;

    Qt::ItemFlags flags(const QModelIndex &index) const;

    bool setData(const QModelIndex &index, const QVariant &value, int role = Qt::EditRole);

    /**
     * Sets the checkstate of all childitems to chk. This function calls itself recursively.
     *
     * @param chk value to set the check state
     * @param parent all children of the parent get checked
     */
    void checkAllChildren(bool chk, BodyParts3DItem* parent);

    /**
     * passes a save request to the BodyParts model
     *
     * @param fileName name and path of file to be saved
     */
    void exportCatalog(const std::string& fileName);

    /**
     * Returns all relevant filenames of the Bodyparts that are checked
     */
    std::vector<std::string> getFileNames();

    /**
     * Returns the Colors of checked BodyParts
     */
    std::vector<tgt::vec4> getColors();

    /**
     * Returns the Descritions (names) of checked BodyParts
     */
    std::vector<std::string> getDescriptions();
private:

    /**
     * Returns the Item at index index
     */
    BodyParts3DItem *getItem(const QModelIndex &index) const;

    QVector<QVariant> rootData_;
    BodyParts3DItem *rootItem_;
    BodyParts3DCatalogBase* catalog_;
};

}//namespace voren


#endif // VRN_BODYPARTS3DMODEL_H
