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

#include "voreen/qt/utils/voreenqtworkspacehandler.h"

#include "tgt/gpucapabilities.h"

#include "voreen/core/network/workspace.h"
#include "voreen/core/network/networkevaluator.h"
#include "voreen/core/properties/link/linkevaluatorid.h"

#include "voreen/qt/voreenapplicationqt.h"

#include "voreen/qt/mainwindow/voreenqtmainwindow.h"
#include "voreen/qt/widgets/propertylistwidget.h"
#include "voreen/qt/networkeditor/networkeditor.h"

#include <QApplication>
#include <QFileDialog>
#include <QMessageBox>
#include <QString>
#include <QErrorMessage>
#include <QStandardPaths>

#include "tgt/glcontextmanager.h"

namespace voreen {

const std::string VoreenQtWorkspaceHandler::loggerCat_ = "voreenqt.WorkspaceHandler";

VoreenQtWorkspaceHandler::VoreenQtWorkspaceHandler()
    : QObject()
    , workspace_(new Workspace())
    , evaluator_(0)
    , mainWindow_(0)
    , networkEditorWidget_(0)
    //, propertyListWidget_(0)
    , readOnlyWorkspace_(false)
{
    workspace_->addObserver(static_cast<WorkspaceObserver*>(this));
}

VoreenQtWorkspaceHandler::~VoreenQtWorkspaceHandler() {
    // Ensure network evaluator is unlocked for deinitialization.
    evaluator_->unlock();
    propagateWorkspace(0);
    delete workspace_;
    workspace_ = 0;
}

Workspace* VoreenQtWorkspaceHandler::getWorkspace() const {
    return workspace_;
}

void VoreenQtWorkspaceHandler::registerNetworkEvaluator(NetworkEvaluator* evaluator) {
    tgtAssert(evaluator, "null pointer passed");
    evaluator_ = evaluator;
    evaluator_->addObserver(static_cast<NetworkEvaluatorObserver*>(this));
}

void VoreenQtWorkspaceHandler::registerMainWindow(VoreenQtMainWindow* mainWindow) {
    tgtAssert(mainWindow, "null pointer passed");
    mainWindow_ = mainWindow;
}

void VoreenQtWorkspaceHandler::registerNetworkEditorWidget(NetworkEditor* networkEditorWidget) {
    tgtAssert(networkEditorWidget, "null pointer passed");
    networkEditorWidget_ = networkEditorWidget;
}

/*void VoreenQtWorkspaceHandler::registerPropertyListWidget(PropertyListWidget* propertyListWidget) {
    tgtAssert(propertyListWidget, "null pointer passed");
    propertyListWidget_ = propertyListWidget;
}*/

void VoreenQtWorkspaceHandler::registerWorkspaceUsingWidget(UsesWorkspace* workspaceUser) {
    tgtAssert(workspaceUser, "null pointer passed");
    std::pair<std::set<UsesWorkspace*>::iterator, bool> res = workspaceUsers_.insert(workspaceUser);
    tgtAssert(res.second, "Widget has been registered already!");
}

void VoreenQtWorkspaceHandler::unregisterWorkspaceUsingWidget(UsesWorkspace* workspaceUser) {
    tgtAssert(workspaceUser, "null pointer passed");
    size_t res = workspaceUsers_.erase(workspaceUser);
    tgtAssert(res == 1, "Widget was not registered!");
}


void VoreenQtWorkspaceHandler::newWorkspace() {
    tgtAssert(workspace_, "no workspace (not initialized)");
    tgtAssert(evaluator_, "no network evaluator");

    if (!askSave())
        return;

    blockSignals(true);

    // Unlock for proper deinitialization.
    bool locked = evaluator_->isLocked();
    if (locked)
        evaluator_->unlock();

    // remove workspace/network from registered entities
    propagateWorkspace(0);

    // clear workspace resources
    delete workspace_;
    workspace_ = new Workspace();
    workspace_->addObserver(static_cast<WorkspaceObserver*>(this));
    // generate new resources
    workspace_->setProcessorNetwork(new ProcessorNetwork());

    // Restore lock (forced by toggleNetworkEvaluator).
    if (locked)
        evaluator_->lock();

    blockSignals(false);

    readOnlyWorkspace_ = false;

    // propagate new workspace/network to registered entities
    propagateWorkspace(workspace_);
}

void VoreenQtWorkspaceHandler::openWorkspace(const QString& filename) {
    tgtAssert(workspace_, "no workspace (not initialized)");

    if (!askSave())
        return;

    // check if workspace to load is a test workspace (i.e. located inside a */test/* subdirectory)
    // if so, query user whether the test data path should be used as workDir
    QString workDir;
    std::vector<std::string> pathComponents = tgt::FileSystem::splitPath(filename.toStdString());
    if (std::find(pathComponents.begin(), pathComponents.end(), "test") != pathComponents.end()) {
        tgtAssert(VoreenApplication::app(), "Voreen app not instantiated");
        if (VoreenApplication::app()->getTestDataPath() != "") {
            QMessageBox::StandardButton button =
                QMessageBox::question(mainWindow_, "Regression Test Workspace",
                    "You are apparently opening a regression test workspace. "
                    "Should the test data directory be used as working directory for this workspace?",
                    static_cast<QMessageBox::StandardButtons>(QMessageBox::Yes | QMessageBox::No), QMessageBox::Yes);
            if (button == QMessageBox::Yes)
                workDir = QString::fromStdString(VoreenApplication::app()->getTestDataPath());
        }
        else {
            QMessageBox::warning(mainWindow_, "Regression Test Workspace",
                "You are apparently opening a regression test workspace, "
                "but the test data directory has not been specified in the application settings. "
                "Therefore, the file paths inside the workspace will most likely be invalid.");
        }
    }

    LINFO("Loading workspace " << tgt::FileSystem::absolutePath(filename.toStdString()));
    if (!workDir.isEmpty())
        LINFO("Workspace working path: " << tgt::FileSystem::cleanupPath(workDir.toStdString()));

    QApplication::setOverrideCursor(QCursor(Qt::WaitCursor));
    qApp->processEvents();
    blockSignals(true);

    // Unlock for proper deinitialization.
    bool locked = evaluator_->isLocked();
    if (locked)
        evaluator_->unlock();

    // remove workspace/network from registered entities
    propagateWorkspace(0);
    // clear workspace resources
    delete workspace_;

    workspace_ = new Workspace();
    workspace_->addObserver(static_cast<WorkspaceObserver*>(this));

    // generate new resources
    workspace_->setProcessorNetwork(new ProcessorNetwork());

    // Restore lock (forced by toggleNetworkEvaluator).
    if (locked)
        evaluator_->lock();

    blockSignals(false);

    // load workspace
    try {
        workspace_->load(filename.toStdString(), workDir.toStdString());
    }
    catch (SerializationException& e) {
        QApplication::restoreOverrideCursor();
        LERROR("Could not open workspace: " << e.what());

        QErrorMessage* errorMessageDialog = new QErrorMessage(mainWindow_);
        errorMessageDialog->showMessage(tr("Could not open workspace:\n") + e.what());

        newWorkspace();
        return;
    }

    readOnlyWorkspace_ = workspace_->readOnly();
    lastUsedWorkspacePath_ = QString::fromStdString(tgt::FileSystem::dirName(workspace_->getFilename()));

    // propagate workspace/network to registered entities
    propagateWorkspace(workspace_);

    QApplication::restoreOverrideCursor();

    // Emit warning if version is beyond specific threshold defined by the Processor Network.
    if(workspace_->getProcessorNetwork()->getVersion() < ProcessorNetwork::WARNING_VERSION) {
        VoreenApplication::app()->showMessageBox("Outdated Workspace",
                                                 "Your workspace was created using a version of Voreen prior to " + ProcessorNetwork::WARNING_VERSION_STRING + ".\n"
                                                 "Therefore this workspace is not guaranteed to work as expected.\n\n"
                                                 "To disable this warning, please save the workspace again.\n"
        );
    }

    // show deserialization errors
    showWorkspaceErrors();
    showNetworkErrors();
}

void VoreenQtWorkspaceHandler::openWorkspace() {
    tgtAssert(workspace_, "no workspace (not initialized)");

    /*if (!askSave())
        return; */

    QFileDialog fileDialog(mainWindow_, tr("Open Workspace..."), QDir(lastUsedWorkspacePath_).absolutePath());
    fileDialog.setOption(QFileDialog::DontUseNativeDialog);

    QStringList filters;
    filters << "Voreen workspaces (*.vws)";
    fileDialog.setNameFilters(filters);

    QList<QUrl> urls;
    urls << QUrl::fromLocalFile(VoreenApplication::app()->getApplicationResourcePath("workspaces").c_str());
    urls << QUrl::fromLocalFile(VoreenApplication::app()->getUserDataPath().c_str());
    urls << QUrl::fromLocalFile(VoreenApplication::app()->getBasePath("modules").c_str());
    if (QDir(VoreenApplication::app()->getBasePath("custommodules").c_str()).exists())
        urls << QUrl::fromLocalFile(VoreenApplication::app()->getBasePath("custommodules").c_str());
    for (auto f : QStandardPaths::standardLocations(QStandardPaths::DesktopLocation))
        urls << QUrl::fromLocalFile(f);
    for (auto f : QStandardPaths::standardLocations(QStandardPaths::HomeLocation))
        urls << QUrl::fromLocalFile(f);
    fileDialog.setSidebarUrls(urls);

    if (!fileDialog.exec())
        return;

    openWorkspace(fileDialog.selectedFiles().at(0));
}

bool VoreenQtWorkspaceHandler::saveWorkspace() {
    tgtAssert(workspace_, "no workspace (not initialized)");

    QString filename = QString::fromStdString(workspace_->getFilename());

    if (filename.isEmpty() || readOnlyWorkspace_) {
        return saveWorkspaceAs();
    }

    // check if workspace is saved as a test workspace (i.e. located inside a test subdirectory)
    // if so, query user whether the test data path should be used as workDir
    QString workDir;
    std::vector<std::string> pathComponents = tgt::FileSystem::splitPath(filename.toStdString());
    if (std::find(pathComponents.begin(), pathComponents.end(), "test") != pathComponents.end()) {
        tgtAssert(VoreenApplication::app(), "Voreen app not instantiated");
        if (VoreenApplication::app()->getTestDataPath() != "") {
            QMessageBox::StandardButton button =
                QMessageBox::question(mainWindow_, "Regression Test Workspace",
                    "You are apparently saving a regression test workspace. "
                    "Should the test data directory be used as working directory for this workspace?",
                    static_cast<QMessageBox::StandardButtons>(QMessageBox::Yes | QMessageBox::No), QMessageBox::Yes);
            if (button == QMessageBox::Yes)
                workDir = QString::fromStdString(VoreenApplication::app()->getTestDataPath());
        }
        else {
            QMessageBox::warning(mainWindow_, "Regression Test Workspace",
                "You are apparently saving a regression test workspace, "
                "but the test data directory has not been specified in the application settings. "
                "Therefore, the file paths inside the workspace will most likely be invalid.");
        }
    }

    if (workDir.isEmpty())
        LINFO("Saving workspace to " << tgt::FileSystem::cleanupPath(filename.toStdString()));
    else
        LINFO("Saving workspace to " << tgt::FileSystem::cleanupPath(filename.toStdString()) << " with working path " << tgt::FileSystem::cleanupPath(workDir.toStdString()));

    try {
        readOnlyWorkspace_ = false;
        workspace_->save(filename.toStdString(), true, workDir.toStdString());
    }
    catch (SerializationException& e) {
        LERROR("Could not save workspace: " << e.what());
        QErrorMessage* errorMessageDialog = new QErrorMessage(mainWindow_);
        errorMessageDialog->showMessage(tr("Could not save workspace:\n") + e.what());
        return false;
    }

    lastUsedWorkspacePath_ = QString::fromStdString(tgt::FileSystem::dirName(filename.toStdString()));
    if (mainWindow_)
        mainWindow_->setWorkspace(workspace_);

    return true;
}

bool VoreenQtWorkspaceHandler::saveWorkspaceAs() {
    tgtAssert(workspace_, "no workspace (not initialized)");

    QFileDialog fileDialog(mainWindow_, tr("Save Workspace As..."), QDir(lastUsedWorkspacePath_).absolutePath());
    fileDialog.setFileMode(QFileDialog::AnyFile);
    fileDialog.setAcceptMode(QFileDialog::AcceptSave);
    fileDialog.setConfirmOverwrite(true);
    fileDialog.setOption(QFileDialog::DontUseNativeDialog);

    QStringList filters;
    filters << "Voreen workspaces (*.vws)";
    fileDialog.setNameFilters(filters);
    QList<QUrl> urls;
    urls << QUrl::fromLocalFile(VoreenApplication::app()->getApplicationResourcePath("workspaces").c_str());
    urls << QUrl::fromLocalFile(VoreenApplication::app()->getUserDataPath().c_str());
    urls << QUrl::fromLocalFile(VoreenApplication::app()->getBasePath("modules").c_str());
    if (QDir(VoreenApplication::app()->getBasePath("custommodules").c_str()).exists())
        urls << QUrl::fromLocalFile(VoreenApplication::app()->getBasePath("custommodules").c_str());
    for (auto f : QStandardPaths::standardLocations(QStandardPaths::DesktopLocation))
        urls << QUrl::fromLocalFile(f);
    for (auto f : QStandardPaths::standardLocations(QStandardPaths::HomeLocation))
        urls << QUrl::fromLocalFile(f);
    fileDialog.setSidebarUrls(urls);

    if (fileDialog.exec()) {
        // set workspace filename and call standard save method
        QString name = fileDialog.selectedFiles().at(0);
        std::string newpath = (!name.endsWith(".vws") ? name.toStdString() + ".vws" : name.toStdString());
        bool readOnly = false;
        try {
            readOnly = checkForReadOnly(newpath);
        }
        catch (...) {
            readOnly = false;
        }
        if(!readOnly) {
            workspace_->setFilename(newpath);
            readOnlyWorkspace_ = false;
            return saveWorkspace();
        }
        else {
            QMessageBox::warning(mainWindow_, "ReadOnly Workspace",
                "You are trying to override a read only workspace.\n"
                "Please try another filename or save under a new filename.");
            return false;
        }
    }
    else {
        return false;
    }
}

bool VoreenQtWorkspaceHandler::askSave() {
    //Workspace must be modified and the application property must be true
    if (workspace_->isModified() && VoreenApplication::app()->getAskForSave()) {
        switch (QMessageBox::question(mainWindow_, tr("Modified Workspace"), tr("Save the current workspace?"),
            QMessageBox::Yes | QMessageBox::No | QMessageBox::Cancel, QMessageBox::Yes))
        {
        case QMessageBox::Yes:
            return saveWorkspace();
        case QMessageBox::No:
            return true;
        default:
            return false;
        }
    }
    return true;
}

bool VoreenQtWorkspaceHandler::checkForReadOnly(const std::string& wsFilepath) {
     // open file for reading
    std::fstream fileStream(wsFilepath.c_str(), std::ios_base::in);
    if (fileStream.fail()) {
        //LERROR("Failed to open file '" << tgt::FileSystem::absolutePath(filename) << "' for reading.");
        throw SerializationException("Failed to open workspace file '" + tgt::FileSystem::absolutePath(wsFilepath) + "' for reading.");
    }

    // read data stream into deserializer
    std::stringbuf buffer;
    do
    {
        // Use 0 character instead of '\n' to minimize the number of get-calls...
        fileStream.get(buffer, 0);
    } while (fileStream.good() && !fileStream.eof()
        && (buffer.sputc(fileStream.get()) != std::stringbuf::traits_type::eof()));

    // Parse input...
    TiXmlDocument document;
    document.Parse(buffer.str().c_str());
    TiXmlElement* root = document.RootElement();
    // Is there no root element?
    if (!root)
        throw SerializationFormatException(std::string("No root node found."));

    // Has root node incorrect name?
    if (root->ValueStr() != XmlSerializationConstants::ROOTNODE) {
        throw SerializationFormatException("XML root node name is '" + root->ValueStr()
                                              + "' instead of '" + XmlSerializationConstants::ROOTNODE + "'.");
    }
    const std::string* version = root->Attribute(XmlSerializationConstants::VERSIONATTRIBUTE);
    // Is serialization version not set?
    if (!version)
        throw SerializationFormatException("XML root node has no version attribute.");
    // Does XmlSerializer and XmlDeserializer version not match the XML document version?
    if (*version != XmlSerializationConstants::VERSION) {
        throw SerializationVersionMismatchException("XML document has version " + *version
                                                        + " instead of " + XmlSerializationConstants::VERSION + ".");
    }
    // deserialize readonly from workspace
    bool readOnly = false;
    if(TiXmlElement* wsElem = root->FirstChildElement("Workspace")) {
        const char* ro = wsElem->Attribute("readonly");
        if(ro && std::strcmp(ro, "true") == 0)
            return true;
        else
            return false;
    } else
        return false;
}

void VoreenQtWorkspaceHandler::showNetworkErrors() {
    tgtAssert(workspace_, "no workspace (not initialized)");
    tgtAssert(workspace_->getProcessorNetwork(), "no network");

    // alert about errors in the Network
    std::vector<std::string> errors = workspace_->getProcessorNetwork()->getErrors();
    if (!errors.empty()) {
        QString msg;
        for (size_t i=0; i < errors.size(); i++) {
            msg += "<li>" + QString(errors[i].c_str()) + "</li>\n";
            LWARNING(errors[i]);
        }

        QErrorMessage* errorMessageDialog = new QErrorMessage(mainWindow_);
        errorMessageDialog->resize(600, 300);
        errorMessageDialog->setWindowTitle(tr("Network Deserialization"));
        errorMessageDialog->showMessage(tr("There were <b>%1 errors</b> loading the network:\n<ul>").arg(errors.size())
            + msg + "\n</ul>");

        qApp->processEvents();
    }
}

void VoreenQtWorkspaceHandler::showWorkspaceErrors() {
    tgtAssert(workspace_, "no workspace (not initialized)");

    // alert about errors in the Network
    std::vector<std::string> errors = workspace_->getErrors();
    if (!errors.empty()) {
        QString msg;
        for (size_t i=0; i < errors.size(); i++) {
            msg += "<li>" + QString(errors[i].c_str()) + "</li>\n";
            LWARNING(errors[i]);
        }

        QErrorMessage* errorMessageDialog = new QErrorMessage(mainWindow_);
        errorMessageDialog->resize(600, 300);
        errorMessageDialog->setWindowTitle(tr("Workspace Deserialization"));
        errorMessageDialog->showMessage(tr("There were <b>%1 errors</b> loading the workspace %2:\n<ul>").arg(
            errors.size()).arg(QString::fromStdString(workspace_->getFilename()))
            + msg + "\n</ul>");

        qApp->processEvents();
    }
}
void VoreenQtWorkspaceHandler::networkAssigned(ProcessorNetwork* network, ProcessorNetwork* /*previousNetwork*/) {
    if (!network)
        return;

    // check whether the assigned network has a different workspace than the current one
    // (happens when the network has not been assigned by the WorkspaceHandler but a module for example)
    if (network->getWorkspace() != workspace_) {
        workspace_ = network->getWorkspace();

        // propagate workspace
        propagateWorkspace(workspace_);

        // emit gui mode change
        if (network->getMetaDataContainer().hasMetaData("uiMode")) {
            qApp->processEvents();
            std::string guiMode = network->getMetaDataContainer().getMetaData("uiMode")->toString();
            if (guiMode == "network")
                emit(guiModeChanged(VoreenQtMainWindow::MODE_NETWORK));
            else if (guiMode == "application")
                emit(guiModeChanged(VoreenQtMainWindow::MODE_APPLICATION));
            else {
                LWARNING("unknown uiMode: " << guiMode);
            }
        }
    }

}

void VoreenQtWorkspaceHandler::propagateWorkspace(Workspace* workspace) {
    tgtAssert(evaluator_, "no network evaluator");
    tgtAssert(mainWindow_, "no mainwindow");
    tgtAssert(networkEditorWidget_, "no network editor");

    ProcessorNetwork* network = workspace ? workspace->getProcessorNetwork() : 0;

    // register as network observer (no longer needed)
    if (workspace && !workspace->isObservedBy(static_cast<WorkspaceObserver*>(this))) {
        workspace->addObserver(static_cast<WorkspaceObserver*>(this));
    }
    // propagate to network editor first to convey to the user immediately that the workspace has been loaded
    networkEditorWidget_->setWorkspace(workspace);
    qApp->processEvents();

    // assign network to evaluator, also initializes the network
    evaluator_->setProcessorNetwork(network);
    // propagate workspace to further widgets
    for (std::set<UsesWorkspace*>::iterator it = workspaceUsers_.begin(); it != workspaceUsers_.end(); it++) {
        (*it)->setWorkspace(workspace);
    }
    // notify mainwindow //main window must be notified after propertywidget to set prodessor selection right
    mainWindow_->setWorkspace(workspace);
}

bool VoreenQtWorkspaceHandler::rebuildShaders() {
    bool allSuccessful = true;

    std::vector<Processor*> procs = workspace_->getProcessorNetwork()->getProcessors();
    for(size_t i = 0; i < procs.size(); i++) {
        std::vector<ShaderProperty*> props = procs.at(i)->getPropertiesByType<ShaderProperty>();
        for(size_t j = 0; j < props.size(); j++) {
            if(!props.at(j)->rebuild())
                allSuccessful = false;
        }
    }

    if (!ShdrMgr.rebuildAllShadersFromFile())
        allSuccessful = false;

    if (allSuccessful) {
        evaluator_->invalidateProcessors();
        return true;
    }

    return false;
}

} // namespace
