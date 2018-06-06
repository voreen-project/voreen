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

#include "dicomhierarchymodel.h"

namespace voreen {

/* DicomHierarchyModel
   ------------------- */

const std::string DicomHierarchyModel::loggerCat_ = "voreen.gdcm.DicomHierarchyModel";

DicomHierarchyModel::DicomHierarchyModel(QObject * parent) :
    QAbstractItemModel(parent)
{
    rootNode_ = 0;
}

DicomHierarchyModel::~DicomHierarchyModel() {
    if (rootNode_)
        delete rootNode_;
}

AbstractDicomHierarchyNode* DicomHierarchyModel::nodeFromIndex(const QModelIndex &index) const
{
    if (index.isValid()) {
        return static_cast<AbstractDicomHierarchyNode*>(index.internalPointer());
    } else {
        return rootNode_;
    }
}

void DicomHierarchyModel::clear() {
    beginResetModel();
    if (rootNode_)
        delete rootNode_;

    PatientInfo p;
    rootNode_ = new PatientHierarchyNode(AbstractDicomHierarchyNode::ROOT, p);

    endResetModel();
}

AbstractDicomHierarchyNode* DicomHierarchyModel::getInvisibleRoot() const {
    return rootNode_;
}

bool DicomHierarchyModel::canFetchMore ( const QModelIndex & parent ) const {

    if (!parent.isValid())
        return false;

    AbstractDicomHierarchyNode* node = nodeFromIndex(parent);

    if (!node || node->type_ == AbstractDicomHierarchyNode::SERIES)
        return false;

    return true;
}

bool DicomHierarchyModel::hasChildren ( const QModelIndex & parent) const {

    AbstractDicomHierarchyNode* node = nodeFromIndex(parent);

    if (!node)
        return false;

    if (node->type_ == AbstractDicomHierarchyNode::SERIES)
        return false;
    else
        return true;
}

int DicomHierarchyModel::rowCount ( const QModelIndex & parent) const {

    AbstractDicomHierarchyNode *node = nodeFromIndex(parent);

    if (!node)
        return 0;

    return node->children_.count();
}

int DicomHierarchyModel::columnCount(const QModelIndex &/*parent*/) const {
    return 1;
}

void DicomHierarchyModel::fetchMore ( const QModelIndex & parent ) {

    //eventually: check for Volume type first to instead use C-MOVE ;)
    if (!canFetchMore(parent))
        return;

    AbstractDicomHierarchyNode *node = nodeFromIndex(parent);

    //if node already has children: return without contacting the server
    //as data has already been fetched
    if (node->children_.size() > 0)
        return;
    /*
    //remove old studies
    if (node->children_.size() > 0) {
        LERROR("Already found children ;)");
        beginRemoveRows(parent,0,node->children_.size()-1);
        node->children_.clear();
        endRemoveRows();
    }*/
    switch (node->type_) {
    case AbstractDicomHierarchyNode::PATIENT: {
        std::vector<StudyInfo> studies;
        try {
            //find all studies
            studies = gdcmReader_->findNetworkStudies(url_,scpAet_,scpPort_,node->getDicomHierarchyID());
        }
        catch (tgt::FileException e) {
            LERROR(e.what());
            QMessageBox::warning(0,"Could not find studies", e.what());
            return;
        }

        //insert all studies
        beginInsertRows(parent, 0, static_cast<int>(studies.size()-1));

        for (std::vector<StudyInfo>::iterator it = studies.begin(); it != studies.end(); ++it) {
            StudyHierarchyNode* child = new StudyHierarchyNode(AbstractDicomHierarchyNode::STUDY,*it);
            child->parent_ = node;
            node->children_.append(child);
        }
        endInsertRows();
        break;
    }
    case AbstractDicomHierarchyNode::STUDY: {
        std::vector<SeriesInfo> series;
        try {
            //find all series
            series = gdcmReader_->findNetworkSeries(url_,scpAet_,scpPort_,node->parent_->getDicomHierarchyID(),node->getDicomHierarchyID());
        }
        catch (tgt::FileException e) {
            LERROR(e.what());
            QMessageBox::warning(0,"Could not find series", e.what());
            return;
        }

        //insert all studies
        beginInsertRows(parent,0, static_cast<int>(series.size()-1));

        for (std::vector<SeriesInfo>::iterator it = series.begin(); it != series.end(); ++it) {
            SeriesHierarchyNode* child = new SeriesHierarchyNode(AbstractDicomHierarchyNode::SERIES,*it);
            child->parent_ = node;
            node->children_.append(child);
        }
        endInsertRows();
        break;
    }
    default:
        //Do nothing ;)
        break;
    }
}

QModelIndex DicomHierarchyModel::index(int row, int column, const QModelIndex &parent) const {
    if (!rootNode_)
            return QModelIndex();

    AbstractDicomHierarchyNode *parentNode = nodeFromIndex(parent);

    return createIndex(row, column, parentNode->children_[row]);
}

QModelIndex DicomHierarchyModel::parent(const QModelIndex& index) const {
    AbstractDicomHierarchyNode *node = nodeFromIndex(index);

    if (!node)
        return QModelIndex();

    AbstractDicomHierarchyNode *parentNode = node->parent_;

    if (!parentNode)
        return QModelIndex();

    AbstractDicomHierarchyNode *grandparentNode = parentNode->parent_;

    if (!grandparentNode)
        return QModelIndex();

    int row = grandparentNode->children_.indexOf(parentNode);

    return createIndex(row, 0, parentNode);
}

QVariant DicomHierarchyModel::data(const QModelIndex &index, int role) const
{
    if (role != Qt::DisplayRole)
        return QVariant();

    AbstractDicomHierarchyNode *node = nodeFromIndex(index);

    if (!node)
        return QVariant();


    return node->getDataString();
}

QVariant DicomHierarchyModel::headerData(int section, Qt::Orientation orientation, int role) const {
    if (role == Qt::DisplayRole)
        return QString("PatientID (Name) -> StudyInstanceUID (Description) -> SeriesInstanceUID (Description)");

    return QAbstractItemModel::headerData(section,orientation,role);
}

} //end namespace
