/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2019 University of Muenster, Germany,                        *
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

#include "dicomconnectiondialog.h"

#include "../io/dicomnetworkconnector.h"
#include "voreen/core/voreenapplication.h"

namespace voreen {

/*
  DicomConnectionInformationWidget
                                    */

ConnectionInformationWidget::ConnectionInformationWidget(QWidget* parent) : QWidget(parent)
{

    QVBoxLayout* mainLayout = new QVBoxLayout();
    setLayout(mainLayout);
    currentConnectionTitle_ = new QLabel("Current Connection:");
    mainLayout->addWidget(currentConnectionTitle_);
    QHBoxLayout* hLayout = new QHBoxLayout();
    currentUrlLabel_ = new QLabel("-");
    hLayout->addWidget(currentUrlLabel_);
    currentAetLabel_ = new QLabel("-");
    //currentAetLabel_->setAlignment(Qt::AlignCenter);
    hLayout->addWidget(currentAetLabel_);
    currentPortLabel_ = new QLabel("-");
    //currentPortLabel_->setAlignment(Qt::AlignRight);
    hLayout->addWidget(currentPortLabel_);
    mainLayout->addLayout(hLayout);
}

ConnectionInformationWidget::~ConnectionInformationWidget(){}

void ConnectionInformationWidget::setConnectionInformation(const std::string& url, const std::string& scpaet, uint16_t port) {
    //display information about current connection
    std::string urlString = "URL: " + url;
    currentUrlLabel_->setText(urlString.c_str());
    std::string aetString = "SCP AE Title: " + scpaet;
    currentAetLabel_->setText(aetString.c_str());
    std::string portString = "Port: " + itos(static_cast<int>(port));
    currentPortLabel_->setText(portString.c_str());
}

/*
  DicomConnectionDialog
                        */

const std::string DicomConnectionDialog::loggerCat_ = "voreen.gdcm.DicomConnectionDialog";

DicomConnectionDialog::DicomConnectionDialog(QWidget* parent)
    : QDialog(parent)
    , gdcmReader_(new GdcmVolumeReader(VoreenApplication::app()->createProgressDialog()))
{
    setWindowTitle("Connect to DICOM Server");

    QVBoxLayout* mainLayout = new QVBoxLayout();
    setLayout(mainLayout);

    QHBoxLayout* scpLayout = new QHBoxLayout();
    mainLayout->addLayout(scpLayout);

    connectInfo_ = new ConnectionInformationWidget();
    mainLayout->addWidget(connectInfo_);
    connectInfo_->setVisible(false);

    //add widgets for SCP connection data
    scpLayout->addWidget(new QLabel("URL: "));
    urlField_ = new QLineEdit();
    //get default SCP url
    std::string scpUrl = dynamic_cast<const GdcmModule*>(VoreenApplication::app()->getModule("gdcm"))->getDefaultScpUrl();
    urlField_->setText(scpUrl.c_str());
    scpLayout->addWidget(urlField_);

    scpLayout->addWidget(new QLabel("SCP AE Title: "));
    scpAetField_ = new QLineEdit();
    //get default SCP AE Title
    std::string scpAet = dynamic_cast<const GdcmModule*>(VoreenApplication::app()->getModule("gdcm"))->getDefaultScpAet();
    scpAetField_->setText(scpAet.c_str());
    scpLayout->addWidget(scpAetField_);

    scpLayout->addWidget(new QLabel("Port: "));
    scpPortField_ = new QLineEdit();
    //get default SCP port
    int scpPort = dynamic_cast<const GdcmModule*>(VoreenApplication::app()->getModule("gdcm"))->getDefaultScpPortNumber();
    scpPortField_->setText(itos(scpPort).c_str());
    scpLayout->addWidget(scpPortField_);

    QHBoxLayout* buttonLayout = new QHBoxLayout();
    submitButton_ = new QPushButton("Open");
    buttonLayout->addWidget(submitButton_);
    refreshButton_ = new QPushButton("Refresh");
    buttonLayout->addWidget(refreshButton_);
    cancelButton_ = new QPushButton("Cancel");
    buttonLayout->addWidget(cancelButton_);
    mainLayout->addLayout(buttonLayout);

    connect(submitButton_, SIGNAL(clicked()),
        this, SLOT(submitFindPatients()));
    connect(refreshButton_, SIGNAL(clicked()),
        this, SLOT(refreshPatients()));
    connect(cancelButton_, SIGNAL(clicked()),
        this, SLOT(close()));

    //Model
    model_ = new DicomHierarchyModel(this);
    model_->setGdcmVolumeReader(gdcmReader_);

    //View
    tree_ = new QTreeView(this);
    tree_->setModel(model_);
    //tree_->setExpandsOnDoubleClick(false);
    mainLayout->addWidget(tree_,1);
    tree_->setVisible(false);

    // ...
    connect(tree_, SIGNAL(doubleClicked(const QModelIndex& /*index*/)),
                  this,   SLOT(emitLoading(const QModelIndex& /*index*/)));

    connect(tree_->selectionModel(), SIGNAL(currentChanged(const QModelIndex & /*current*/, const QModelIndex & /*previous*/)),
            this, SLOT(treeItemSelected(const QModelIndex & /*current*/, const QModelIndex & /*previous*/)));

    detailTable_ = new QTableWidget(1,1,this);
    detailTable_->setMaximumHeight(detailTable_->rowHeight(0) + detailTable_->horizontalHeader()->height() + detailTable_->horizontalScrollBar()->height());
    detailTable_->setEditTriggers(QAbstractItemView::NoEditTriggers);
    detailTable_->removeColumn(0); detailTable_->removeRow(0);
    detailTable_->setDisabled(true);
    detailTable_->setVisible(false);
    mainLayout->addWidget(detailTable_,0);
}

DicomConnectionDialog::~DicomConnectionDialog() {

    model_->clear();
    detailTable_->clear();

    delete gdcmReader_;
    gdcmReader_ = 0;
}

void DicomConnectionDialog::emitLoading(const QModelIndex& index) /*const*/ {

    //blockSignals(true);
    setDisabled(true);
    update();

    AbstractDicomHierarchyNode* node = model_->nodeFromIndex(index);
    if (node->type_ != AbstractDicomHierarchyNode::SERIES) {
        setDisabled(false);
        return;
    }

    std::string patientID = node->parent_->parent_->getDicomHierarchyID();
    std::string studyID = node->parent_->getDicomHierarchyID();
    std::string seriesID = node->getDicomHierarchyID();

    //construct url
    VolumeURL origin("dicom-scp", model_->getUrl());
    origin.addSearchParameter("port", itos(static_cast<int>(model_->getScpPort())));
    origin.addSearchParameter("SCP-AET", model_->getScpAet());
    origin.addSearchParameter("PatientID",patientID);
    origin.addSearchParameter("StudyInstanceUID",studyID);
    origin.addSearchParameter("SeriesInstanceUID",seriesID);

    emit(loadFromURL(origin.getURL(), gdcmReader_));

    //blockSignals(false);
    setDisabled(false);
    close();
}

void DicomConnectionDialog::submitFindPatients() {
    tgtAssert(gdcmReader_, "no gdcm reader");

    /*submitButton_->setDisabled(true);
    refreshButton_->setDisabled(true);
    cancelButton_->setDisabled(true);*/
    setDisabled(true);
    update();

    //set connection widget visibility false
    connectInfo_->setVisible(false);

    //clear detail table and set visibility false
    detailTable_->setVisible(false);
    detailTable_->setDisabled(true);
    detailTable_->clear();
    if (detailTable_->rowCount() > 0)
        detailTable_->removeRow(0);

    while (detailTable_->columnCount() > 0)
        detailTable_->removeColumn(0);

    std::vector<PatientInfo> patients;
    tree_->setVisible(false);
    model_->clear();

    try {
        //find all patients
        patients = gdcmReader_->findNetworkPatients(urlField_->text().toStdString(),scpAetField_->text().toStdString(),static_cast<uint16_t>(stoi(scpPortField_->text().toStdString())));
    }
    catch (tgt::FileException e) {
        LERROR(e.what());
        QMessageBox::warning(this,"Could not establish connection", e.what());

        /*submitButton_->setDisabled(false);
        refreshButton_->setDisabled(false);
        cancelButton_->setDisabled(false);*/
        adjustSize();
        setDisabled(false);
        return;
    }

    //set server data
    model_->setScpAet(scpAetField_->text().toStdString());
    model_->setUrl(urlField_->text().toStdString());
    model_->setScpPort(static_cast<uint16_t>(stoi(scpPortField_->text().toStdString())));

    //insert all patients and enable submit button
    AbstractDicomHierarchyNode* root = model_->getInvisibleRoot();

    for (std::vector<PatientInfo>::iterator it = patients.begin(); it != patients.end(); ++it) {
        PatientHierarchyNode *patient = new PatientHierarchyNode(AbstractDicomHierarchyNode::PATIENT, *it);
        patient->parent_ = root;
        root->children_.append(patient);
    }
    tree_->setVisible(true);
    detailTable_->setVisible(true);
    /*submitButton_->setDisabled(false);
    refreshButton_->setDisabled(false);
    cancelButton_->setDisabled(false);*/
    setDisabled(false);

    //set connection widget information and set visible true
    connectInfo_->setConnectionInformation(model_->getUrl(),model_->getScpAet(),model_->getScpPort());
    connectInfo_->setVisible(true);
    adjustSize();
}

void DicomConnectionDialog::treeItemSelected(const QModelIndex & current, const QModelIndex & /*previous*/) {

    if (!current.isValid())
        detailTable_->setDisabled(true);

    AbstractDicomHierarchyNode* node = model_->nodeFromIndex(current);

    if (!node)
        return;

    //clear detail table
    detailTable_->clear();
    if (detailTable_->rowCount() > 0)
        detailTable_->removeRow(0);
    while (detailTable_->columnCount() > 0)
        detailTable_->removeColumn(0);

    QStringList columnLabels;
    QStringList rowLabel;

    //fill table
    switch(node->type_) {
    case AbstractDicomHierarchyNode::PATIENT: {
        //add columns
        detailTable_->insertColumn(0);
        detailTable_->insertColumn(1);
        columnLabels << QString("PatientID") << QString("PatientName");
        detailTable_->setHorizontalHeaderLabels(columnLabels);

        //add rows
        detailTable_->insertRow(0);
        rowLabel << QString("Patient");
        detailTable_->setVerticalHeaderLabels(rowLabel);

        //insert data
        QTableWidgetItem* id = new QTableWidgetItem(node->getDicomHierarchyID().c_str());
        //id->setFlags(!Qt::ItemIsEditable);
        detailTable_->setItem(0,0,id);

        QTableWidgetItem* name;
        PatientHierarchyNode* patientNode = dynamic_cast<PatientHierarchyNode*>(node);
        if (patientNode)
                name = new QTableWidgetItem(patientNode->data_.patientName_.c_str());
        else
            name = new QTableWidgetItem("");

        //name->setFlags(!Qt::ItemIsEditable);

        detailTable_->setItem(0,1,name);

        //resize and set visible
        detailTable_->resizeColumnsToContents();
        detailTable_->setVisible(true);
        break;
    }
    case AbstractDicomHierarchyNode::STUDY: {
        //add columns
        detailTable_->insertColumn(0);
        detailTable_->insertColumn(1);
        detailTable_->insertColumn(2);
        detailTable_->insertColumn(3);
        columnLabels << QString("StudyInstanceUID") << QString("StudyDescription") << QString("StudyDate") << QString("StudyTime");
        detailTable_->setHorizontalHeaderLabels(columnLabels);

        //add rows
        detailTable_->insertRow(0);
        rowLabel << QString("Study");
        detailTable_->setVerticalHeaderLabels(rowLabel);

        //insert data
        QTableWidgetItem* id = new QTableWidgetItem(node->getDicomHierarchyID().c_str());
        //id->setFlags(!Qt::ItemIsEditable);
        detailTable_->setItem(0,0,id);

        QTableWidgetItem* desc;
        QTableWidgetItem* date;
        QTableWidgetItem* time;
        StudyHierarchyNode* studyNode = dynamic_cast<StudyHierarchyNode*>(node);
        if (studyNode) {
                desc = new QTableWidgetItem(studyNode->data_.studyDescription_.c_str());
                date = new QTableWidgetItem(studyNode->data_.studyDate_.c_str());
                time = new QTableWidgetItem(studyNode->data_.studyTime_.c_str());
        }
        else {
            desc = new QTableWidgetItem("");
            date = new QTableWidgetItem("");
            time = new QTableWidgetItem("");
        }
        //desc->setFlags(!Qt::ItemIsEditable);
        //date->setFlags(!Qt::ItemIsEditable);
        //time->setFlags(!Qt::ItemIsEditable);
        detailTable_->setItem(0,1,desc);
        detailTable_->setItem(0,2,date);
        detailTable_->setItem(0,3,time);

        //resize and set visible
        detailTable_->resizeColumnsToContents();
        detailTable_->setVisible(true);
        break;
    }
    case AbstractDicomHierarchyNode::SERIES: {
        //add columns
        detailTable_->insertColumn(0);
        detailTable_->insertColumn(1);
        detailTable_->insertColumn(2);
        detailTable_->insertColumn(3);
        detailTable_->insertColumn(4);
        columnLabels << QString("SeriesInstanceUID") << QString("SeriesDescription") << QString("SeriesDate") << QString("SeriesTime") << QString("Number of Images");
        detailTable_->setHorizontalHeaderLabels(columnLabels);

        //add rows
        detailTable_->insertRow(0);
        rowLabel << QString("Series");
        detailTable_->setVerticalHeaderLabels(rowLabel);

        //insert data
        QTableWidgetItem* id = new QTableWidgetItem(node->getDicomHierarchyID().c_str());
        //id->setFlags(!Qt::ItemIsEditable);
        detailTable_->setItem(0,0,id);

        QTableWidgetItem* desc;
        QTableWidgetItem* date;
        QTableWidgetItem* time;
        QTableWidgetItem* images;
        SeriesHierarchyNode* seriesNode = dynamic_cast<SeriesHierarchyNode*>(node);
        if (seriesNode) {
                desc = new QTableWidgetItem(seriesNode->data_.seriesDescription_.c_str());
                date = new QTableWidgetItem(seriesNode->data_.seriesDate_.c_str());
                time = new QTableWidgetItem(seriesNode->data_.seriesTime_.c_str());
                images = new QTableWidgetItem(itos(seriesNode->data_.numberOfImages_).c_str());
        }
        else {
            desc = new QTableWidgetItem("");
            date = new QTableWidgetItem("");
            time = new QTableWidgetItem("");
            images = new QTableWidgetItem("");
        }
        //desc->setFlags(!Qt::ItemIsEditable);
        //date->setFlags(!Qt::ItemIsEditable);
        //time->setFlags(!Qt::ItemIsEditable);
        //images->setFlags(!Qt::ItemIsEditable);
        detailTable_->setItem(0,1,desc);
        detailTable_->setItem(0,2,date);
        detailTable_->setItem(0,3,time);
        detailTable_->setItem(0,4,images);

        //resize and set visible
        detailTable_->resizeColumnsToContents();
        /*detailTable_->setMaximumHeight(detailTable_->horizontalHeader()->height()+
                                       detailTable_->rowHeight(0) +
                                       detailTable_->horizontalScrollBar()->height());*/
        detailTable_->setVisible(true);
        break;
    }
    default:
        ;
    }

    detailTable_->adjustSize();
    detailTable_->setDisabled(false);
}

void DicomConnectionDialog::refreshPatients() {
    tgtAssert(gdcmReader_, "no gdcm reader");

    //nothing to refresh: return
    if (!connectInfo_->isVisible())
        return;

    setDisabled(true);
    submitButton_->setDisabled(true);
    refreshButton_->setDisabled(true);
    cancelButton_->setDisabled(true);

    update();

    //clear detail table and set visibility false
    detailTable_->clear();
    detailTable_->setDisabled(true);
    detailTable_->setVisible(false);
    if (detailTable_->rowCount() > 0)
        detailTable_->removeRow(0);

    while (detailTable_->columnCount() > 0)
        detailTable_->removeColumn(0);

    std::vector<PatientInfo> patients;
    tree_->setVisible(false);
    model_->clear();

    try {
        //find all patients
        patients = gdcmReader_->findNetworkPatients(model_->getUrl(),model_->getScpAet(),model_->getScpPort());
    }
    catch (tgt::FileException e) {
        LERROR(e.what());
        QMessageBox::warning(this,"Could not establish connection", e.what());
        submitButton_->setDisabled(false);
        refreshButton_->setDisabled(false);
        cancelButton_->setDisabled(false);
        setDisabled(false);
        return;
    }

    //insert all patients and enable submit button
    AbstractDicomHierarchyNode* root = model_->getInvisibleRoot();

    for (std::vector<PatientInfo>::iterator it = patients.begin(); it != patients.end(); ++it) {
        PatientHierarchyNode *patient = new PatientHierarchyNode(AbstractDicomHierarchyNode::PATIENT, *it);
        patient->parent_ = root;
        root->children_.append(patient);
    }
    tree_->setVisible(true);
    detailTable_->setVisible(true);
    submitButton_->setDisabled(false);
    refreshButton_->setDisabled(false);
    cancelButton_->setDisabled(false);
    setDisabled(false);
}

} // namespace
