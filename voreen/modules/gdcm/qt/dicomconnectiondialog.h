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

#ifndef VRN_DICOMCONNECTIONDIALOG_H
#define VRN_DICOMCONNECTIONDIALOG_H

#include <QDialog>
#include <QPushButton>
#include <QLineEdit>
#include <QTreeView>
#include <QStandardItemModel>
#include <QStandardItem>
#include <QTableWidget>
#include <QHeaderView>
#include <QScrollBar>
#include <QMessageBox>

#include <QHBoxLayout>
#include <QVBoxLayout>
#include <QLabel>
#include <QWidget>

#include <QAbstractItemModel>

#include "tgt/types.h"

#include "../io/gdcmvolumereader.h"
#include "voreen/qt/voreenqtapi.h"
#include "voreen/core/datastructures/volume/volume.h"
#include "./dicomhierarchymodel.h"
#include "voreen/core/voreenapplication.h"

#include "voreen/core/utils/stringutils.h"

namespace voreen {

class VRN_QT_API ConnectionInformationWidget : public QWidget {
    //Q_OBJECT
public:
    ConnectionInformationWidget(QWidget* parent = 0);
    ~ConnectionInformationWidget();

    void setConnectionInformation(const std::string& url, const std::string& scpaet, uint16_t port);

private:

    QLabel* currentConnectionTitle_;
    QLabel* currentUrlLabel_; ///< displays the url of the SCP currently connected to
    QLabel* currentAetLabel_; ///< displays the AET of the SCP currently connected to
    QLabel* currentPortLabel_; ///< displays the port of the SCP currently connected to
};

class VolumeReader;

class VRN_QT_API DicomConnectionDialog : public QDialog {
    Q_OBJECT
public:
    DicomConnectionDialog(QWidget* parent);
    ~DicomConnectionDialog();

signals:
    void loadFromURL(const std::string& url, VolumeReader* reader) const;

private:
    static const std::string loggerCat_;

    QLineEdit* urlField_; ///< contains the DICOM SCP url
    QLineEdit* scpAetField_; ///< contains the AE title of the DICOM SCP
    QLineEdit* scpPortField_; ///< port to contact the DICOM SCP

    QPushButton* submitButton_; ///< button for contacting the DICOM SCP
    QPushButton* cancelButton_; ///< button for closing the dialog
    QPushButton* refreshButton_; ///< button for refreshing the connection to the DICOM SCP

    ConnectionInformationWidget* connectInfo_; ///< widget for displaying information about current connection

    GdcmVolumeReader* gdcmReader_; ///< GdcmVolumeReader (with network support: provides the interface to contact the DICOM SCP)

    QTreeView* tree_; ///< the tree view for the hierarchical DICOM data of the SCP
    DicomHierarchyModel* model_; ///< the (hierarchical) data model

    QTableWidget* detailTable_; ///< show details about the selected data set

private slots:
    /// when called with network support: contacts the DICOM SCP and gets a list of available patients
    void submitFindPatients();

    /// when called with network support: contacts the DICOM SCP slready connected to and refreshes the list of available patients
    void refreshPatients();

    /// is called when an item of the tree is selected so that details can be displayed in the table
    void treeItemSelected(const QModelIndex & /*current*/, const QModelIndex & /*previous*/);

    /// emits a signal to the VolumeReader to open a specific data set
    void emitLoading(const QModelIndex& index) /*const*/;
};

} // namespace

#endif // VRN_DICOMCONNECTIONDIALOG_H
