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

#include "voreen/qt/widgets/volumeiohelper.h"

#include "voreen/core/io/volumeserializerpopulator.h"
#include "voreen/core/io/volumeserializer.h"
#include "voreen/core/io/volumereader.h"
#include "voreen/core/io/volumewriter.h"
#include "voreen/core/voreenapplication.h"
#include "voreen/core/utils/stringutils.h"
#include "voreen/core/datastructures/volume/volume.h"
#include "modules/core/io/rawvolumereader.h"

#include "voreen/qt/widgets/rawvolumewidget.h"
#include "voreen/qt/widgets/volumereaderselectiondialog.h"
#include "voreen/qt/widgets/volumelistingdialog.h"
#include "voreen/qt/voreenapplicationqt.h"
#include "voreen/qt/progressdialog.h"

#include "tgt/filesystem.h"

#include <QAction>
#include <QApplication>
#include <QDesktopServices>
#include <QErrorMessage>
#include <QFile>
#include <QFileDialog>
#include <QInputDialog>
#include <QMainWindow>
#include <QMenu>
#include <QMessageBox>
#include <QUrl>
#include <QSettings>

namespace voreen {

const std::string VolumeIOHelper::loggerCat_("voreen.qt.VolumeIOHelper");

VolumeIOHelper::VolumeIOHelper(QWidget* parent, FileMode fileMode)
    : parent_(parent)
    , fileMode_(fileMode)
    , volumeSerializerPopulator_(0)
    , readerSelectionDialog_(0)
    , volumeListingDialog_(0)
{
    progressBar_ = VoreenApplication::app()->createProgressDialog();
    if (progressBar_) {
        progressBar_->setTitle("Loading volume");
        progressBar_->setProgressMessage("Loading volume ...");
    }
    volumeSerializerPopulator_ = new VolumeSerializerPopulator(progressBar_);

    readerSelectionDialog_ = new VolumeReaderSelectionDialog(VoreenApplicationQt::qtApp()->getMainWindow());
    connect(readerSelectionDialog_, SIGNAL(readerSelected(const std::string&, VolumeReader*)),
        this, SLOT(loadURL(const std::string&, VolumeReader*)));

    volumeListingDialog_ = new VolumeListingDialog(VoreenApplicationQt::qtApp()->getMainWindow(),
        fileMode == SINGLE_FILE ? VolumeListingDialog::SINGLE_SELECTION : VolumeListingDialog::MULTI_SELECTION);
    connect(volumeListingDialog_, SIGNAL(originsSelected(const std::vector<VolumeURL>&, VolumeReader*)),
        this, SLOT(loadOrigins(const std::vector<VolumeURL>&, VolumeReader*)));
}

VolumeIOHelper::~VolumeIOHelper() {
    delete progressBar_;
    progressBar_ = 0;

    delete volumeSerializerPopulator_;
    volumeSerializerPopulator_ = 0;
}

void VolumeIOHelper::showFileOpenDialog() {

    tgtAssert(volumeSerializerPopulator_, "no volumeserializerpopulator");
    tgtAssert(volumeSerializerPopulator_->getVolumeSerializer(), "no volume serializer");

    //
    // 1. Setup and execute file open dialog
    //

    // open path
    std::string openPath;
    QSettings settings;
    settings.beginGroup("VolumeIOHelper");
    if (settings.contains("lastVolumePath"))
        openPath = settings.value("lastVolumePath").toString().toStdString();
    else
        openPath = VoreenApplication::app()->getUserDataPath("volumes");
    // we want absolute path name to prevent problems when the file dialog perform chdir()
    QDir openDir(QString::fromStdString(openPath));

    // filters
    QStringList filters;
    QString allVolumesFilter = QString::fromStdString(getVolumeReaderFilterString());
    QString allFilesFilter = "All Files (*)";
    filters << allVolumesFilter << allFilesFilter;
    // additionally, create one filter per volume reader
    std::map<std::string, VolumeReader*> filterToReaderMap;
    std::vector<VolumeReader*> volumeReaders =
        volumeSerializerPopulator_->getVolumeSerializer()->getReaders();
    for (size_t i=0; i<volumeReaders.size(); i++) {
        VolumeReader* curReader = volumeReaders.at(i);
        std::vector<std::string> extensionVec = curReader->getSupportedExtensions();
        std::vector<std::string> filenameVec = curReader->getSupportedFilenames();
        if (extensionVec.empty() && filenameVec.empty())
            continue;

        std::string filterStr = curReader->getFormatDescription() + " (";
        if (!extensionVec.empty())
            filterStr += "*." + strJoin(extensionVec, " *.");
        if (!filenameVec.empty())
            filterStr += " " + strJoin(filenameVec, " ");
        filterStr += ")";

        if (filterStr.find("vvd") != std::string::npos)
            filters.insert(1, QString::fromStdString(filterStr)); //< put VVD format at front
        else
            filters << QString::fromStdString(filterStr);
        filterToReaderMap.insert(std::pair<std::string, VolumeReader*>(filterStr, curReader));
    }

    // sidebar URLs
    QList<QUrl> urls;
    urls << QUrl::fromLocalFile(VoreenApplication::app()->getApplicationResourcePath("volumes").c_str());
    urls << QUrl::fromLocalFile(VoreenApplication::app()->getApplicationResourcePath("datasets").c_str());
    urls << QUrl::fromLocalFile(VoreenApplication::app()->getUserDataPath().c_str());
    urls << QUrl::fromLocalFile(VoreenApplication::app()->getBasePath("modules").c_str());
    if (QDir(VoreenApplication::app()->getBasePath("custommodules").c_str()).exists())
        urls << QUrl::fromLocalFile(VoreenApplication::app()->getBasePath("custommodules").c_str());
    if (VoreenApplication::app()->getTestDataPath() != "")
        urls << QUrl::fromLocalFile(VoreenApplication::app()->getTestDataPath().c_str());
    for (auto f : QStandardPaths::standardLocations(QStandardPaths::DesktopLocation))
        urls << QUrl::fromLocalFile(f);
    for (auto f : QStandardPaths::standardLocations(QStandardPaths::HomeLocation))
        urls << QUrl::fromLocalFile(f);

    // create dialog
    QFileDialog fileDialog(parent_, tr("Load Volume..."), openDir.absolutePath(), "");
    fileDialog.setNameFilters(filters);
    fileDialog.setOption(QFileDialog::DontUseNativeDialog);
    fileDialog.setSidebarUrls(urls);
    fileDialog.setViewMode(QFileDialog::Detail);
    fileDialog.setAcceptMode(QFileDialog::AcceptOpen);
    if (fileMode_ == SINGLE_FILE)
        fileDialog.setFileMode(QFileDialog::ExistingFile);
    else
        fileDialog.setFileMode(QFileDialog::ExistingFiles);

    // execute dialog
    if (fileDialog.exec() != QDialog::Accepted)
        return;

    // retrieve selected filenames from dialog
    std::vector<std::string> filenames;
    const QStringList& lst = fileDialog.selectedFiles();
    QStringList::const_iterator it = lst.begin();
    for (; it != lst.end(); ++it)
        filenames.push_back(it->toStdString());
    if (filenames.empty()) {
        LWARNING("no files selected");
        return;
    }

    // retrieve the user selected file filter
    QString selectedFilter = fileDialog.selectedNameFilter();

    // store dialog path
    settings.setValue("lastVolumePath", fileDialog.directory().absolutePath());
    settings.endGroup();


    //
    // 2. Load selected volume files (depending on selected filter)
    //
    tgtAssert(volumeSerializerPopulator_, "no volumeserializerpopulator");
    if (selectedFilter == allVolumesFilter || selectedFilter == allFilesFilter) {
        // combined filter for all supported volume files => pass each filename to matching reader
        const VolumeSerializer* volumeSerializer = volumeSerializerPopulator_->getVolumeSerializer();
        for (size_t i = 0; i < filenames.size(); i++) {
            std::vector<VolumeReader*> volumeReaders = volumeSerializer->getReaders(filenames.at(i));
            tgtAssert(!volumeReaders.empty(), "no volume readers"); //< otherwise exception expected
            if (volumeReaders.size() == 1) { // only one compatible reader
                loadURL(filenames.at(i), volumeReaders.front());
            }
            else { // multiple possible readers => query user
                tgtAssert(readerSelectionDialog_, "no volumereaderselectiondialog");
                readerSelectionDialog_->setURL(filenames.at(i));
                readerSelectionDialog_->setReaders(volumeReaders);
                readerSelectionDialog_->show();
            }
        }
    }
    else {
        // one of the single volume reader filters directly selected
        if (filterToReaderMap.find(selectedFilter.toStdString()) == filterToReaderMap.end()) {
            LWARNING("unknown filter selected");
            return;
        }
        VolumeReader* reader = filterToReaderMap.find(selectedFilter.toStdString())->second;
        for (size_t i = 0; i < filenames.size(); i++)
            loadURL(filenames[i], reader);
    }
}

void VolumeIOHelper::showFileSaveDialog(const VolumeBase* volume) {
    tgtAssert(volume, "null pointer passed");

    tgtAssert(volumeSerializerPopulator_, "no volumeserializerpopulator");
    tgtAssert(volumeSerializerPopulator_->getVolumeSerializer(), "no volume serializer");

    //
    // 1. Setup and execute file save dialog
    //

    // open path
    std::string savePath;
    QSettings settings;
    settings.beginGroup("VolumeIOHelper");
    if (settings.contains("lastVolumePath"))
        savePath = settings.value("lastVolumePath").toString().toStdString();
    else
        savePath = VoreenApplication::app()->getUserDataPath("volumes");
    // we want absolute path name to prevent problems when the file dialog perform chdir()
    QDir openDir(QString::fromStdString(savePath));

    // create one filter per volume writer and extension
    QStringList filters;
    std::map<std::string, std::pair<VolumeWriter*, std::string> > filterToWriterMap;
    std::vector<VolumeWriter*> volumeWriters =
        volumeSerializerPopulator_->getVolumeSerializer()->getWriters();
    for (size_t i=0; i<volumeWriters.size(); i++) {
        VolumeWriter* curWriter = volumeWriters.at(i);

        // extensions
        std::vector<std::string> extensionVec = curWriter->getSupportedExtensions();
        for (size_t j=0; j<extensionVec.size(); j++) {
            std::string extension = extensionVec.at(j);
            std::string filterStr = curWriter->getFormatDescription() + " (*." + extension + ")";
            if (extension == "vvd")
                filters.insert(0, QString::fromStdString(filterStr)); //< put VVD format at front
            else
                filters << QString::fromStdString(filterStr);
            filterToWriterMap.insert(
                std::pair<std::string, std::pair<VolumeWriter*, std::string> >(filterStr, std::pair<VolumeWriter*, std::string>(curWriter, extension)));
        }

        // filename
        std::vector<std::string> filenamesVec = curWriter->getSupportedFilenames();
        for (size_t j=0; j<filenamesVec.size(); j++) {
            std::string filename = filenamesVec.at(j);
            std::string filterStr = curWriter->getFormatDescription() + " (" + filename + ")";
            filters << QString::fromStdString(filterStr);
            filterToWriterMap.insert(
                std::pair<std::string, std::pair<VolumeWriter*, std::string> >(filterStr, std::pair<VolumeWriter*, std::string>(curWriter, filename)));
        }

    }

    // sidebar URLs
    QList<QUrl> urls;
    urls << QUrl::fromLocalFile(VoreenApplication::app()->getApplicationResourcePath("volumes").c_str());
    urls << QUrl::fromLocalFile(VoreenApplication::app()->getApplicationResourcePath("datasets").c_str());
    urls << QUrl::fromLocalFile(VoreenApplication::app()->getUserDataPath().c_str());
    urls << QUrl::fromLocalFile(VoreenApplication::app()->getBasePath("modules").c_str());
    if (QDir(VoreenApplication::app()->getBasePath("custommodules").c_str()).exists())
        urls << QUrl::fromLocalFile(VoreenApplication::app()->getBasePath("custommodules").c_str());
    if (VoreenApplication::app()->getTestDataPath() != "")
        urls << QUrl::fromLocalFile(VoreenApplication::app()->getTestDataPath().c_str());
    for (auto f : QStandardPaths::standardLocations(QStandardPaths::DesktopLocation))
        urls << QUrl::fromLocalFile(f);
    for (auto f : QStandardPaths::standardLocations(QStandardPaths::HomeLocation))
        urls << QUrl::fromLocalFile(f);

    // create dialog
    QFileDialog saveAsDialog(parent_, tr("Save Volume..."), openDir.absolutePath(), "");
    saveAsDialog.setNameFilters(filters);
    saveAsDialog.setOption(QFileDialog::DontUseNativeDialog);
    saveAsDialog.setSidebarUrls(urls);
    saveAsDialog.setViewMode(QFileDialog::Detail);
    saveAsDialog.setAcceptMode(QFileDialog::AcceptSave);
    saveAsDialog.setFileMode(QFileDialog::AnyFile);
    saveAsDialog.setConfirmOverwrite(true);

    // execute dialog
    if (saveAsDialog.exec() != QDialog::Accepted)
        return;

    // retrieve selected filename from dialog
    const QStringList& lst = saveAsDialog.selectedFiles();
    if (lst.empty()) {
        LWARNING("no file selected");
        return;
    }
    else if (lst.size() > 1) {
        LWARNING("more than one file selected");
        return;
    }
    std::string filepath = lst.first().toStdString();

    // retrieve the user selected file filter
    QString selectedFilter = saveAsDialog.selectedNameFilter();

    // store dialog path
    settings.setValue("lastVolumePath", saveAsDialog.directory().absolutePath());
    settings.endGroup();


    //
    // 2. Save file to specified path (depending on selected filter)
    //
    tgtAssert(volumeSerializerPopulator_, "no volumeserializerpopulator");
    if (filterToWriterMap.find(selectedFilter.toStdString()) == filterToWriterMap.end()) {
        LWARNING("unknown filter selected");
        return;
    }
    VolumeWriter* writer = filterToWriterMap.find(selectedFilter.toStdString())->second.first;

    // if file path does not end with extension, add extension to path
    std::string selectedExtension = filterToWriterMap.find(selectedFilter.toStdString())->second.second;
    std::string filename = tgt::FileSystem::fileName(filepath);
    if (filename != selectedExtension && !QString::fromStdString(filename).endsWith(QString::fromStdString("." + selectedExtension)))
        filepath += "." + selectedExtension;

    LDEBUG("saving volume: " << filepath);
    saveVolumeToPath(volume, writer, filepath);
}

void VolumeIOHelper::loadURL(const std::string& url, VolumeReader* reader) {
    tgtAssert(reader, "null pointer passed");

    // In case we have selected a raw volume reader, open raw volume file dialog.
    if (dynamic_cast<RawVolumeReader*>(reader)) {
        loadRawVolume(url);
        return;
    }

    try {
        std::vector<VolumeURL> origins = reader->listVolumes(url);
        if (origins.size() > 1) { // if more than one volume listed at URL, show listing dialog
            tgtAssert(volumeListingDialog_, "no volumelistingdialog");
            volumeListingDialog_->setOrigins(origins, reader);
            volumeListingDialog_->show();
        }
        else if (origins.size() == 1) { // load single volume directly
            loadOrigin(origins.front(), reader);
        }
        else {
            LWARNING("No volumes found at URL '" << url << "'");
            QErrorMessage* errorMessageDialog = new QErrorMessage(VoreenApplicationQt::qtApp()->getMainWindow());
            errorMessageDialog->showMessage(QString::fromStdString("No volumes found at '" + url + "'"));
        }
    }
    catch (const tgt::FileException& e) {
        LWARNING(e.what());
        QErrorMessage* errorMessageDialog = new QErrorMessage(VoreenApplicationQt::qtApp()->getMainWindow());
        errorMessageDialog->showMessage(e.what());
    }
}

void VolumeIOHelper::loadOrigin(const VolumeURL& origin, VolumeReader* reader) {
    tgtAssert(reader, "null pointer passed");

    try {
        VolumeBase* handle = reader->read(origin);
        if (handle) {
            handle->getOrigin().addSearchParameter("preferredReader", reader->getClassName());
            emit(volumeLoaded(handle));
        }
        else {
            LERROR("Reader '" << reader->getClassName() << "' returned null pointer as volume (exception expected)");
        }
    }
    catch (const tgt::FileException& e) {
        LERROR(e.what());
        QErrorMessage* errorMessageDialog = new QErrorMessage(VoreenApplicationQt::qtApp()->getMainWindow());
        errorMessageDialog->showMessage(e.what());
    }
    catch (std::bad_alloc&) {
        LERROR("bad allocation while reading file: " << origin.getURL());
        QErrorMessage* errorMessageDialog = new QErrorMessage(VoreenApplicationQt::qtApp()->getMainWindow());
        errorMessageDialog->showMessage("Bad allocation while reading file: " + QString::fromStdString(origin.getURL()));
    }
    catch (std::exception& e) {
        LERROR("unknown exception while reading file '" << origin.getURL() << "':" << e.what());
        QErrorMessage* errorMessageDialog = new QErrorMessage(VoreenApplicationQt::qtApp()->getMainWindow());
        errorMessageDialog->showMessage("Unknown exception while reading file '" +
            QString::fromStdString(origin.getURL()) + "': " + e.what());
    }

    if (progressBar_)
        progressBar_->hide();
}

void VolumeIOHelper::loadOrigins(const std::vector<VolumeURL>& origins, VolumeReader* reader) {
    tgtAssert(reader, "null pointer passed");
    for (size_t i=0; i<origins.size(); i++) {
        loadOrigin(origins.at(i), reader);
    }
}

void VolumeIOHelper::loadRawVolume(const std::string& filenameStd) {

    QString filename = QString::fromStdString(filenameStd);
    if (filename.isEmpty())
        return;

    // query raw parameters via dialog
    std::string objectModel;
    std::string format;
        int numFrames;
    tgt::ivec3 dim;
    tgt::vec3 spacing;
    int headerSkip;
    bool bigEndian;
    tgt::mat4 trafoMat = tgt::mat4::identity;
    RawVolumeWidget* rawVW = new RawVolumeWidget(parent_, tr("Please enter the properties for <br><strong>") + filename + "</strong>",
        objectModel, format, numFrames, dim, spacing, headerSkip, bigEndian, trafoMat);
    if (rawVW->exec() != QDialog::Accepted)
        return;

    // derive expected file size from provided properties
    uint formatBytes = 1;;
    if (format == "USHORT" || format == "USHORT_12" || format == "SHORT")
        formatBytes = 2;
    else if (format == "FLOAT")
        formatBytes = 4;
    else if (format == "UINT" || format == "INT")
        formatBytes = 4;

    int numChannels = 1;
    if (objectModel == "RGB")
        numChannels = 3;
    else if (objectModel == "RGBA")
        numChannels = 4;
    else if (objectModel.find("TENSOR_") == 0)
        numChannels = 6;

    uint rawSize = headerSkip + formatBytes * numChannels * (dim.x * dim.y * dim.z) * numFrames;

    // inform/query user, if file size does not match
    if (QFile(filename).size() != rawSize) {
        QMessageBox::StandardButton retButton = QMessageBox::Yes;
        if(QFile(filename).size() > rawSize) {
            QString msg = tr("The provided properties result in a size smaller\nthan the actual file size. Do you want to continue?");
            retButton = QMessageBox::question(parent_, tr("Size mismatch"), msg,
                QMessageBox::Yes | QMessageBox::Cancel, QMessageBox::Yes);
        }
        else if (QFile(filename).size() < rawSize) {
            QString msg = tr("The provided properties result in a size\ngreater than the actual file size.");
            retButton = QMessageBox::warning(parent_, tr("Size mismatch"), msg,
                QMessageBox::Cancel);
        }
        if (retButton != QMessageBox::Yes && retButton != QMessageBox::Ok)
            return;
    }

    qApp->processEvents();

    // load raw volume
    try {
        for (int frame=0; frame < numFrames; ++frame) {
            RawVolumeReader rawReader(progressBar_);
            rawReader.setReadHints(dim, spacing, objectModel, format, frame, headerSkip, bigEndian);
            VolumeList* collection = rawReader.read(filename.toStdString());
            if (collection && !collection->empty()) {
                tgtAssert(collection->size() == 1, "More than one raw volume returned");
                Volume* volumeHandle = static_cast<Volume*>(collection->first());
                oldVolumePosition(volumeHandle);
                volumeHandle->setPhysicalToWorldMatrix(trafoMat);
                volumeHandle->setTimestep(static_cast<float>(frame));
                emit(volumeLoaded(volumeHandle));
            }
            delete collection;
        }
    }
    catch (const tgt::FileException& e) {
        LERROR(e.what());
        QErrorMessage* errorMessageDialog = new QErrorMessage(VoreenApplicationQt::qtApp()->getMainWindow());
        errorMessageDialog->showMessage(e.what());
    }
    catch (std::bad_alloc&) {
        LERROR("bad allocation while reading file: " << filename.toStdString());
        QErrorMessage* errorMessageDialog = new QErrorMessage(VoreenApplicationQt::qtApp()->getMainWindow());
        errorMessageDialog->showMessage("Bad allocation while reading file: " + filename);
    }
}

void VolumeIOHelper::saveVolumeToPath(const VolumeBase* volume, VolumeWriter* writer, const std::string& filepath) {
    tgtAssert(volume, "null pointer passed");
    tgtAssert(writer, "null pointer passed");

    try {
        writer->write(filepath, volume);
        emit(volumeSaved(volume, filepath));
    }
    catch (const tgt::FileException& e) {
        LERROR(e.what());
        QErrorMessage* errorMessageDialog = new QErrorMessage(VoreenApplicationQt::qtApp()->getMainWindow());
        errorMessageDialog->showMessage(e.what());
    }
    catch (std::exception& e) {
        LERROR("unknown exception while writing file '" << filepath << "':" << e.what());
        QErrorMessage* errorMessageDialog = new QErrorMessage(VoreenApplicationQt::qtApp()->getMainWindow());
        errorMessageDialog->showMessage("Unknown exception while writing file '" +
            QString::fromStdString(filepath) + "': " + e.what());
    }
}

std::string VolumeIOHelper::getVolumeReaderFilterString() const {
    tgtAssert(volumeSerializerPopulator_, "no volumeserializerpopulator");

    std::string filterStr;

    std::vector<std::string> extensionVec = volumeSerializerPopulator_->getSupportedReadExtensions();
    if (!extensionVec.empty()) {
        filterStr += "*." + strJoin(extensionVec, " *.");
    }

    std::vector<std::string> filenameVec = volumeSerializerPopulator_->getSupportedReadFilenames();
    if (!filenameVec.empty()) {
        filterStr += " " + strJoin(filenameVec, " ");
    }

    return "Volume Files (" + filterStr + ")";
}

VolumeIOHelper::FileMode VolumeIOHelper::getFileMode() const {
    return fileMode_;
}

} // namespace voreen
