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

#include <QTreeView>
#include <QStandardItemModel>
#include <QStandardItem>
#include <QMessageBox>

#include <QAbstractItemModel>

#include "../io/gdcmvolumereader.h"
#include "voreen/qt/voreenqtapi.h"
#include "voreen/core/datastructures/volume/volume.h"
#include "tgt/types.h"

#ifndef DICOMHIERARCHYMODEL_H
#define DICOMHIERARCHYMODEL_H

namespace voreen {

/*
  HierarchyNode classes
                        */

class VRN_QT_API AbstractDicomHierarchyNode {

public:

    enum DicomType {ROOT, PATIENT, STUDY, SERIES};

    virtual ~AbstractDicomHierarchyNode() {
        qDeleteAll(children_);
    }

    /// Returns the data of the item to be displayed within the GUI
    virtual QString getDataString() const = 0;

    /// Returns the ID within the DICOM hierarchy
    virtual std::string getDicomHierarchyID() const = 0;

    DicomType type_; ///< type of the node
    AbstractDicomHierarchyNode *parent_; ///< parent node
    QList<AbstractDicomHierarchyNode *> children_; ///< list of child nodes

};

template <class DICOMTYPE> class VRN_QT_API DicomHierarchyNode :
    public AbstractDicomHierarchyNode
{

public:
    DicomHierarchyNode() {}

    DICOMTYPE data_; ///< the actual data within the node
};

class VRN_QT_API PatientHierarchyNode :
        public DicomHierarchyNode<PatientInfo>
{
public:
    PatientHierarchyNode(AbstractDicomHierarchyNode::DicomType type, PatientInfo patient) {
        type_ = type;
        data_ = patient;
        parent_ = 0;
    }

    /// implements abstract method from AbstractDicomHierarchyNode
    virtual QString getDataString() const {
        std::stringstream s;
        s << data_.patientID_;
        if ((!data_.patientName_.empty()) && !(data_.patientName_ == " "))
            s << " (" << data_.patientName_ << ")";
        return QString(s.str().c_str());
    }

    /// returns the PatientID
    virtual std::string getDicomHierarchyID() const {
        return data_.patientID_;
    }
};

class VRN_QT_API StudyHierarchyNode :
        public DicomHierarchyNode<StudyInfo>
{
public:
    StudyHierarchyNode(AbstractDicomHierarchyNode::DicomType type, StudyInfo study) {
        type_ = type;
        data_ = study;
        parent_ = 0;
    }

    /// implements abstract method from AbstractDicomHierarchyNode
    virtual QString getDataString() const {
        std::stringstream s;
        s << data_.studyInstanceUID_;
        if ((!data_.studyDescription_.empty()) && !(data_.studyDescription_ == " "))
            s << " (" << data_.studyDescription_ << ")";
        return QString(s.str().c_str());
    }

    /// returns the StudyInstanceUID
    virtual std::string getDicomHierarchyID() const {
        return data_.studyInstanceUID_;
    }
};

class VRN_QT_API SeriesHierarchyNode :
        public DicomHierarchyNode<SeriesInfo>
{
public:
    SeriesHierarchyNode(AbstractDicomHierarchyNode::DicomType type, SeriesInfo series) {
        type_ = type;
        data_ = series;
        parent_ = 0;
    }

    /// implements abstract method from AbstractDicomHierarchyNode
    virtual QString getDataString() const {
        std::stringstream s;
        s << data_.seriesInstanceUID_;
        if ((!data_.seriesDescription_.empty()) && !(data_.seriesDescription_== " "))
            s << " (" << data_.seriesDescription_ << ")";

        s << " - " << itos(data_.numberOfImages_) << " image(s)";

        return QString(s.str().c_str());
    }

    /// returns the SeriesInstanceUID
    virtual std::string getDicomHierarchyID() const {
        return data_.seriesInstanceUID_;
    }
};


/*
  DicomHierarchyModel
                        */

class VRN_QT_API DicomHierarchyModel :
       public QAbstractItemModel
{
    Q_OBJECT
    public:
        DicomHierarchyModel(QObject * parent = 0);
        ~DicomHierarchyModel();

        /// implements abstract method from QAbstractItemModel
        virtual int rowCount ( const QModelIndex & parent = QModelIndex() ) const;

        /// implements abstract method from QAbstractItemModel
        virtual int columnCount ( const QModelIndex & parent = QModelIndex() ) const;

        /// implements abstract method from QAbstractItemModel
        virtual QModelIndex index(int row, int column, const QModelIndex &parent) const;

        /// implements abstract method from QAbstractItemModel
        virtual QModelIndex parent ( const QModelIndex & index ) const;

        /// implements abstract method from QAbstractItemModel
        virtual QVariant data ( const QModelIndex & index, int role = Qt::DisplayRole ) const;

        /// reimplemented from QAbstractItemModel
        virtual bool canFetchMore ( const QModelIndex & parent ) const;

        /// reimplemented from QAbstractItemModel
        virtual void fetchMore ( const QModelIndex & parent );

        /// reimplemented from QAbstractItemModel
        virtual bool hasChildren ( const QModelIndex & parent = QModelIndex() ) const;

        /// reimplemented from QAbstractItemModel
        virtual QVariant headerData(int section, Qt::Orientation orientation, int role = Qt::DisplayRole ) const;

        /**
         * Clears the whole model, i.e. removes all entries and data
         */
        void clear();

        /**
         * Returns a pointer to the root node that is not visible in the view (i.e. GUI)
         */
        AbstractDicomHierarchyNode* getInvisibleRoot() const;

        /**
         * Encapsulates casting an index within the model into an actual node
         */
        AbstractDicomHierarchyNode* nodeFromIndex(const QModelIndex &index) const;

        void setUrl(const std::string& url) {
            url_ = url;
        }

        void setScpAet(const std::string& aet) {
            scpAet_ = aet;
        }

        void setScpPort(const uint16_t& port) {
            scpPort_ = port;
        }

        std::string getUrl() const {
            return url_;
        }

        std::string getScpAet() const {
            return scpAet_;
        }

        uint16_t getScpPort() const {
            return scpPort_;
        }

        void setGdcmVolumeReader(GdcmVolumeReader* reader) {
            gdcmReader_ = reader;
        }

    private:
        static const std::string loggerCat_;

        AbstractDicomHierarchyNode* rootNode_; ///< root node of the model (invisible in view)

        std::string url_; ///< url for connecting to a DICOM SCP
        std::string scpAet_; ///< AET of a DICOM SCP to connect to it
        uint16_t scpPort_; ///< port where a DICOM SCP should be contacted

        GdcmVolumeReader* gdcmReader_; ///< GdcmVolumeReader that is used to dynamically contact the server

};

} //namespace

#endif // DICOMHIERARCHYMODEL_H
