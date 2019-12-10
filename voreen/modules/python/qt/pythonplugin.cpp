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

// include Python stuff first
#include "modules/python/pythonmodule.h"
#include "pythonplugin.h"
#include "pythonhighlighter.h"
#include "../properties/pythonproperty.h"

#include "voreen/core/voreenapplication.h"
#include "voreen/qt/widgets/codeedit.h"

#include <QObject>
#include <QApplication>
#include <QDesktopServices>
#include <QFileDialog>
#include <QFrame>
#include <QHBoxLayout>
#include <QMessageBox>
#include <QToolButton>

namespace voreen {

const std::string PythonPlugin::loggerCat_ = "voreen.Python.PythonPlugin";

PythonPlugin::PythonPlugin(PythonProperty* property, QWidget* parent)
    : QWidget(parent)
    , property_(property)
    , script_(nullptr)
    , highlighter_(nullptr)
    , fontSize_(9)
    , scriptOwner_(false)
{

    QHBoxLayout* hbox = new QHBoxLayout();
    hbox->setContentsMargins(0,0,0,0);

    if(!property_) {
        runBt_ = new QToolButton();
        runBt_->setIcon(QIcon(":/modules/python/python.png"));
        runBt_->setIconSize(QSize(24, 24));
        runBt_->setShortcut(QKeySequence("Ctrl+R"));
        runBt_->setToolTip("Run Script (Ctrl+R)");
        hbox->addWidget(runBt_);
        QObject::connect(runBt_, SIGNAL(clicked()), this, SLOT(runScript()));

        QFrame* sep = new QFrame();
        sep->setFrameShape(QFrame::VLine);
        hbox->addWidget(sep);

        newBt_ = new QToolButton();
        newBt_->setIcon(QIcon(":/modules/python/python_2.png"));
        newBt_->setIconSize(QSize(24, 24));
        newBt_->setToolTip("New Script");
        hbox->addWidget(newBt_);

        QObject::connect(newBt_, SIGNAL(clicked()), this, SLOT(newScript()));
    }

    openBt_ = new QToolButton();
    openBt_->setIcon(QIcon(":/qt/icons/open.png"));
    openBt_->setIconSize(QSize(24, 24));
    openBt_->setToolTip("Load Script");
    hbox->addWidget(openBt_);
    QObject::connect(openBt_, SIGNAL(clicked()), this, SLOT(openScript()));

    if(!property_) {
        reloadBt_ = new QToolButton();
        reloadBt_->setIcon(QIcon(":/qt/icons/refresh.png"));
        reloadBt_->setIconSize(QSize(24, 24));
        reloadBt_->setToolTip("Reload Script from Disk");
        hbox->addWidget(reloadBt_);
        QObject::connect(reloadBt_, SIGNAL(clicked()), this, SLOT(reloadScript()));
    }

    saveBt_ = new QToolButton();
    saveBt_->setIcon(QIcon(":/qt/icons/save.png"));
    saveBt_->setIconSize(QSize(24, 24));
    saveBt_->setToolTip("Save Script");
    hbox->addWidget(saveBt_);
    QObject::connect(saveBt_, SIGNAL(clicked()), this, SLOT(saveScript()));

    if(!property_) {
        saveAsBt_ = new QToolButton();
        saveAsBt_->setIcon(QIcon(":/qt/icons/saveas.png"));
        saveAsBt_->setIconSize(QSize(24, 24));
        saveAsBt_->setToolTip("Save Script As");
        hbox->addWidget(saveAsBt_);

        QObject::connect(saveAsBt_, SIGNAL(clicked()), this, SLOT(saveScriptAs()));
    }

    hbox->addStretch();

    QToolButton* increaseFontSizeBt = new QToolButton();
    increaseFontSizeBt->setIcon(QIcon(":/qt/icons/viewmag+.png"));
    increaseFontSizeBt->setIconSize(QSize(24, 24));
    increaseFontSizeBt->setToolTip("Increase Font Size");
    hbox->addWidget(increaseFontSizeBt);
    QObject::connect(increaseFontSizeBt, SIGNAL(clicked()), this, SLOT(increaseFontSize()));

    QToolButton* decreaseFontSizeBt = new QToolButton();
    decreaseFontSizeBt->setIcon(QIcon(":/qt/icons/viewmag_.png"));
    decreaseFontSizeBt->setIconSize(QSize(24, 24));
    decreaseFontSizeBt->setToolTip("Decrease Font Size");
    hbox->addWidget(decreaseFontSizeBt);
    QObject::connect(decreaseFontSizeBt, SIGNAL(clicked()), this, SLOT(decreaseFontSize()));

    QWidget* toolButtonBar = new QWidget();
    toolButtonBar->setLayout(hbox);

    QFont font;
    font.setFamily("Courier");
    font.setFixedPitch(true);
    font.setPointSize(fontSize_);

    codeEdit_ = new CodeEdit(false, fontSize_);
    codeEdit_->setFont(font);
    highlighter_ = new PythonHighlighter(codeEdit_->document());

    compilerLogWidget_ = new QTextEdit();
    compilerLogWidget_->setFont(font);
    compilerLogWidget_->setReadOnly(true);
    compilerLogWidget_->setFixedHeight(150);
    QVBoxLayout* vbox = new QVBoxLayout();
    vbox->addWidget(toolButtonBar);
    vbox->addWidget(codeEdit_);
    vbox->addWidget(compilerLogWidget_);
    setLayout(vbox);

    QObject::connect(codeEdit_, SIGNAL(textChanged()), this, SIGNAL(modified()));

    setMinimumSize(300, 400);

    // Register as output listener.
    if(PythonModule::getInstance()) {
        PythonModule::getInstance()->addOutputListener(this);
    }
}

PythonPlugin::~PythonPlugin() {
    if(PythonModule::getInstance()) {
        PythonModule::getInstance()->removeOutputListener(this);
    }

    clearScript();
    delete highlighter_;
}

void PythonPlugin::runScript() {
    tgtAssert(!property_, "Invalid function call");

    if (!script_)
        return;

    compilerLogWidget_->setTextColor(Qt::black);
    compilerLogWidget_->clear();

    // retrieve script source from code editor and assign it to script object
    QString source = codeEdit_->toPlainText();
    source = source.replace("\t", "    ");   //< replace tabs with four spaces
    script_->setSource(source.toStdString());

    if (script_->compile(false)) {
        compilerLogWidget_->setPlainText("Running script...");
        setEnabled(false);
        codeEdit_->updateHighlight();
        qApp->setOverrideCursor(Qt::WaitCursor);
        qApp->processEvents();

        bool success = script_->run(false);

        qApp->restoreOverrideCursor();
        setEnabled(true);
        codeEdit_->updateHighlight();
        codeEdit_->setFocus();
        qApp->processEvents();

        if (success) {
            compilerLogWidget_->append("\nfinished.");
            //compilerLogWidget_->setPlainText(compilerLogWidget_->toPlainText() + QString("finished."));
        }
        else {
            compilerLogWidget_->setTextColor(Qt::red);
            compilerLogWidget_->setPlainText("Runtime Error:\n" + QString::fromStdString(script_->getLog()));
            codeEdit_->moveCursorToPosition(script_->getErrorLine() - 1);
        }
    }
    else { // compilation failed
        compilerLogWidget_->setTextColor(Qt::red);
        compilerLogWidget_->setPlainText("Compile Error:\n" + QString::fromStdString(script_->getLog()));
        codeEdit_->moveCursorToPosition(script_->getErrorLine() - 1, script_->getErrorCol() - 1);
    }
}

void PythonPlugin::newScript() {
    clearScript();
    if (PythonModule::getInstance() && PythonModule::getInstance()->isInitialized()) {
        script_ = PythonModule::getInstance()->loadScript("template.py", false);
        if (script_) {
            QString source = QString::fromStdString(script_->getSource());
            source = source.replace("\t", "    ");  //< replace tabs with four spaces
            PythonModule::getInstance()->dispose(script_);

            script_ = new PythonScript();
            script_->setSource(source.toStdString());
            codeEdit_->setPlainText(source);
            scriptOwner_ = true;
        }
    }
    else {
        QMessageBox::critical(this, tr("Python Error"), tr("Failed to create Python script: module not initialized"));
    }
    updateGuiState();
}
void PythonPlugin::tryLoadScript(QString scriptPath) {
    if (!scriptPath.isEmpty()) {
        clearScript();

        script_ = PythonModule::getInstance()->loadScript(scriptPath.toStdString(), false);
        if (!script_) {
            QMessageBox::critical(this, tr("Python Error"), tr("Python script '%1' could not be loaded.").arg(scriptPath));
        }
        else {
            // retrieve script source from script object and assign it to code editor
            QString source = QString::fromStdString(script_->getSource());
            source = source.replace("\t", "    ");  //< replace tabs with four spaces
            codeEdit_->setPlainText(source);
            scriptOwner_ = false;
        }
        updateGuiState();
    }
}

void PythonPlugin::openScript() {
    if (!PythonModule::getInstance()) {
        LERROR("PythonModule not instantiated");
        return;
    }

    //create filter with supported file formats
    QString filter = "Python Script (*.py)";
    QString scriptPath = selectOpenFileName(filter);
    tryLoadScript(scriptPath);
}

void PythonPlugin::reloadScript() {
    if (!PythonModule::getInstance()) {
        LERROR("PythonModule not instantiated");
        return;
    }

    if(script_) {
        tryLoadScript(QString::fromStdString(script_->getFilename()));
    }
}

void PythonPlugin::saveScript() {
    tgtAssert(codeEdit_, "No code editor");

    if (property_) {

        // In case we are connected to a property, we only care about the current content of the code editor.
        // The script based on an actual file from disk, if even available, can be deleted after having stored
        // the code in the property's script.
        PythonScript script;
        script.setSource(codeEdit_->toPlainText().toStdString());
        property_->set(script);
        // Invalidate, even if the value didn't change. This forces script execution.
        property_->invalidate();

        // Now, we no longer need the script!
        if(script_) {
            clearScript();
            codeEdit_->setPlainText(QString::fromStdString(script.getSource()));
        }
    }
    else {

        if (!script_) {
            // This should not happen since the save button should be disabled.
            LERROR("No script");
            return;
        }

        if (script_->getFilename().empty()) {
            LERROR("Script has no filename");
        } else {
            if (!saveScriptInternal(QString::fromStdString(script_->getFilename()), codeEdit_->toPlainText())) {
                QString message = tr("Failed to save script to file '%1'").arg(
                        QString::fromStdString(script_->getFilename()));
                QMessageBox::critical(this, tr("Python Error"), message);
            }
        }
    }

    updateGuiState();
}

void PythonPlugin::saveScriptAs() {
    if (!PythonModule::getInstance()) {
        LERROR("PythonModule not instantiated");
        return;
    }

    if (!script_) {
        LERROR("No script");
        return;
    }

    //create filter with supported file formats
    QStringList filter;
    filter << "Python Script (*.py)";

    QString fileName = selectSaveFileName(filter);
    if (!fileName.isEmpty()) {
        if (saveScriptInternal(fileName, codeEdit_->toPlainText())) {
            clearScript();
            script_ = PythonModule::getInstance()->loadScript(fileName.toStdString(), false);
            scriptOwner_ = false;
            if (script_) {
                codeEdit_->setPlainText(QString::fromStdString(script_->getSource()));
            }
            else {
                QMessageBox::critical(this, tr("Python Error"), tr("Saved script '%1' could not be loaded.").arg(fileName));
            }
        }
        else {
            QString message = tr("Failed to save script to file '%1'").arg(fileName);
            QMessageBox::critical(this, tr("Python Error"), message);
        }
    }

    updateGuiState();
}

const QString PythonPlugin::selectOpenFileName(QString filter) {
    QFileDialog fileDialog(this);
    fileDialog.setWindowTitle(tr("Choose a Python script to open"));
    fileDialog.setNameFilter(filter);
    fileDialog.setOption(QFileDialog::DontUseNativeDialog);

    QList<QUrl> urls;
    urls << QUrl::fromLocalFile(VoreenApplication::app()->getUserDataPath().c_str());
    urls << QUrl::fromLocalFile(VoreenApplication::app()->getBasePath("modules").c_str());
    urls << QUrl::fromLocalFile(VoreenApplication::app()->getBasePath("custommodules").c_str());
    urls << QUrl::fromLocalFile(VoreenApplication::app()->getCoreResourcePath("scripts").c_str());
    for (const auto& f : QStandardPaths::standardLocations(QStandardPaths::DesktopLocation))
        urls << QUrl::fromLocalFile(f);
    for (const auto& f : QStandardPaths::standardLocations(QStandardPaths::HomeLocation))
        urls << QUrl::fromLocalFile(f);
    fileDialog.setSidebarUrls(urls);

    if (fileDialog.exec() && !fileDialog.selectedFiles().empty()) {
        return fileDialog.selectedFiles()[0];
    }

    return QString();
}

const QString PythonPlugin::selectSaveFileName(QStringList filters) {
    QFileDialog fileDialog(this);
    fileDialog.setWindowTitle(tr("Choose a filename to save script"));
    fileDialog.setDirectory(VoreenApplication::app()->getUserDataPath("scripts").c_str());
    fileDialog.setNameFilters(filters);
    fileDialog.setAcceptMode(QFileDialog::AcceptSave);
    fileDialog.setOption(QFileDialog::DontUseNativeDialog);

    QList<QUrl> urls;
    urls << QUrl::fromLocalFile(VoreenApplication::app()->getUserDataPath().c_str());
    urls << QUrl::fromLocalFile(VoreenApplication::app()->getBasePath("modules").c_str());
    urls << QUrl::fromLocalFile(VoreenApplication::app()->getCoreResourcePath("scripts").c_str());
    for (const auto& f : QStandardPaths::standardLocations(QStandardPaths::DesktopLocation))
        urls << QUrl::fromLocalFile(f);
    for (const auto& f : QStandardPaths::standardLocations(QStandardPaths::HomeLocation))
        urls << QUrl::fromLocalFile(f);
    fileDialog.setSidebarUrls(urls);

    QStringList fileList;
    if (fileDialog.exec() && !fileDialog.selectedFiles().empty()) {
        QString endingFilter = fileDialog.selectedNameFilter();
        int pos = endingFilter.lastIndexOf(".");
        //removes closing bracket
        endingFilter.chop(1);
        endingFilter = endingFilter.mid(pos);

        //look whether the user specified an ending
        std::string fileExtension;
        fileList = fileDialog.selectedFiles();
        size_t dotPosition = fileList[0].toStdString().rfind(".");
        if (dotPosition == std::string::npos) {
            // no ending given -> add ending of selected filter
            fileList[0].append(endingFilter);
        }
        else {
            // an ending was given -> test whether it matches the selected filter
            if (fileList[0].mid(static_cast<int>(dotPosition)) != endingFilter)
                fileList[0].append(endingFilter);
        }
        return fileList[0];
    }
    return QString();
}

bool PythonPlugin::saveScriptInternal(QString filename, QString source) {

    // retrieve script source from code editor and replace tabs with four spaces
    source = source.replace("\t", "    ");

    // save script to disk
    QFile outputFile(filename);
    if (outputFile.open(QIODevice::WriteOnly)) {
        outputFile.write(source.toStdString().c_str(), source.size());
        outputFile.close();
        LINFO("Saved script to file '" << filename.toStdString() << "'");
        return true;
    }
    else {
        LERROR("Failed to open file '" << filename.toStdString() << "' for writing");
        return false;
    }

    //updateGuiState(); //cannot be reached
}

void PythonPlugin::clearScript() {
    if (!script_)
        return;

    if (!PythonModule::getInstance()) {
        LERROR("PythonModule not instantiated");
        return;
    }

    if (!scriptOwner_ && PythonModule::getInstance())
        PythonModule::getInstance()->dispose(script_);
    else
        delete script_;
    script_ = 0;
    scriptOwner_ = false;

    codeEdit_->setPlainText("");
    compilerLogWidget_->setPlainText("");
}

void PythonPlugin::updateFromProperty() {
    tgtAssert(property_, "property was null");

    QString source = QString(property_->get().getSource().c_str());
    if(codeEdit_->toPlainText() != source)
        codeEdit_->setPlainText(source);

    updateGuiState();
}

void PythonPlugin::increaseFontSize() {
    fontSize_++;
    updateFont();
}

void PythonPlugin::decreaseFontSize() {
    fontSize_ = std::max(fontSize_-1, 7);
    updateFont();
}

void PythonPlugin::updateFont() {
    QFont font;
    font.setFamily("Courier");
    font.setFixedPitch(true);
    font.setPointSize(fontSize_);
    codeEdit_->setFont(font);
    compilerLogWidget_->setFont(font);
}

void PythonPlugin::updateGuiState() {

    compilerLogWidget_->setTextColor(Qt::black);
    compilerLogWidget_->clear();

    std::string windowTitle;

    if(property_) {
        // No need to modify button states.
        if (property_->getOwner())
            windowTitle += property_->getOwner()->getGuiName() + " - ";
        windowTitle += property_->getGuiName();
        if(!property_->get().getFilename().empty()) {
            windowTitle += " (" + property_->get().getFilename() + ")";
        }
    }
    else {
        tgtAssert(script_, "Either property or script must be available");

        windowTitle = "Python Script Editor";
        saveAsBt_->setEnabled(true);
        saveBt_->setEnabled(!script_->getFilename().empty());
        if(!script_->getFilename().empty()) {
            windowTitle += " (" + script_->getFilename() + ")";
        }
    }

    if(window()) {
        window()->setWindowTitle(QString::fromStdString(windowTitle));
    }
}

void PythonPlugin::pyStdout(const std::string& out, const std::string& id) {
    // Check if this actual script has been executed.
    if(((property_ && property_->get().getId() == id) || (script_ && script_->getId() == id))) {
        compilerLogWidget_->setTextColor(Qt::black);
        compilerLogWidget_->append(QString::fromStdString(out));
    }
}

void PythonPlugin::pyStderr(const std::string& err, const std::string& id) {
    // Check if this actual script has been executed.
    if(((property_ && property_->get().getId() == id) || (script_ && script_->getId() == id))) {
        compilerLogWidget_->setTextColor(Qt::red);
        compilerLogWidget_->append(QString::fromStdString(err));
    }
}


} // namespace voreen
