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

        QFrame* sep = new QFrame();
        sep->setFrameShape(QFrame::VLine);
        hbox->addWidget(sep);

        newBt_ = new QToolButton();
        newBt_->setIcon(QIcon(":/modules/python/python_2.png"));
        newBt_->setIconSize(QSize(24, 24));
        newBt_->setToolTip("New Script");
        hbox->addWidget(newBt_);

        QObject::connect(runBt_, SIGNAL(clicked()), this, SLOT(runScript()));
        QObject::connect(newBt_, SIGNAL(clicked()), this, SLOT(newScript()));
    }
    openBt_ = new QToolButton();
    openBt_->setIcon(QIcon(":/qt/icons/open.png"));
    openBt_->setIconSize(QSize(24, 24));
    openBt_->setToolTip("Load Script");
    reloadBt_ = new QToolButton();
    reloadBt_->setIcon(QIcon(":/qt/icons/refresh.png"));
    reloadBt_->setIconSize(QSize(24, 24));
    reloadBt_->setToolTip("Reload Script from Disk");
    saveBt_ = new QToolButton();
    saveBt_->setIcon(QIcon(":/qt/icons/save.png"));
    saveBt_->setIconSize(QSize(24, 24));
    saveBt_->setToolTip("Save Script");
    saveAsBt_ = new QToolButton();
    saveAsBt_->setIcon(QIcon(":/qt/icons/saveas.png"));
    saveAsBt_->setIconSize(QSize(24, 24));
    saveAsBt_->setToolTip("Save Script As");

    increaseFontSizeBt_ = new QToolButton();
    increaseFontSizeBt_->setIcon(QIcon(":/qt/icons/viewmag+.png"));
    increaseFontSizeBt_->setIconSize(QSize(24, 24));
    increaseFontSizeBt_->setToolTip("Increase Font Size");
    decreaseFontSizeBt_ = new QToolButton();
    decreaseFontSizeBt_->setIcon(QIcon(":/qt/icons/viewmag_.png"));
    decreaseFontSizeBt_->setIconSize(QSize(24, 24));
    decreaseFontSizeBt_->setToolTip("Decrease Font Size");

    hbox->addWidget(openBt_);
    hbox->addWidget(reloadBt_);
    hbox->addWidget(saveBt_);
    hbox->addWidget(saveAsBt_);
    hbox->addStretch();
    hbox->addWidget(increaseFontSizeBt_);
    hbox->addWidget(decreaseFontSizeBt_);
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

    QObject::connect(openBt_, SIGNAL(clicked()), this, SLOT(openScript()));
    QObject::connect(reloadBt_, SIGNAL(clicked()), this, SLOT(reloadScript()));
    QObject::connect(saveBt_, SIGNAL(clicked()), this, SLOT(saveScript()));
    QObject::connect(saveAsBt_, SIGNAL(clicked()), this, SLOT(saveScriptAs()));
    QObject::connect(increaseFontSizeBt_, SIGNAL(clicked()), this, SLOT(increaseFontSize()));
    QObject::connect(decreaseFontSizeBt_, SIGNAL(clicked()), this, SLOT(decreaseFontSize()));
    QObject::connect(codeEdit_, SIGNAL(textChanged()), this, SIGNAL(modified()));

    setMinimumSize(300, 400);
}

PythonPlugin::~PythonPlugin() {
    delete highlighter_;
}

void PythonPlugin::runScript() {
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
        PythonScript script;
        script.setSource(codeEdit_->toPlainText().toStdString());
        property_->set(script);
        return;
    }

    if (!script_) {
        LERROR("No script");
        return;
    }

    if (script_->getFilename().empty()) {
        LERROR("Script has no filename");
    }
    else {
        if (!saveScriptInternal(QString::fromStdString(script_->getFilename()), codeEdit_->toPlainText())) {
            QString message = tr("Failed to save script to file '%1'").arg(QString::fromStdString(script_->getFilename()));
            QMessageBox::critical(this, tr("Python Error"), message);
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
    for (auto f : QStandardPaths::standardLocations(QStandardPaths::DesktopLocation))
        urls << QUrl::fromLocalFile(f);
    for (auto f : QStandardPaths::standardLocations(QStandardPaths::HomeLocation))
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
    for (auto f : QStandardPaths::standardLocations(QStandardPaths::DesktopLocation))
        urls << QUrl::fromLocalFile(f);
    for (auto f : QStandardPaths::standardLocations(QStandardPaths::HomeLocation))
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

    std::string windowTitle = "";
    if (property_->getOwner())
        windowTitle += property_->getOwner()->getGuiName() + " - ";
    windowTitle += property_->getGuiName() + " (" + property_->get().getFilename() + ")";
    window()->setWindowTitle(QString::fromStdString(windowTitle));
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
    if(!property_) {
        runBt_->setEnabled(script_);
        saveAsBt_->setEnabled(script_);
    }

    QString title;
    if (property_ || (script_ && !script_->getFilename().empty())) {
        saveBt_->setEnabled(true);
        title = "Python Script Editor - " + QString::fromStdString(script_->getFilename());
    }
    else {
        saveBt_->setEnabled(false);
        title = "Python Script Editor";
    }

    //toolWindow_->setWindowTitle(title); // TODO: implement
    /*if (parentWidget() && parentWidget()->parentWidget())
        parentWidget()->parentWidget()->setWindowTitle(title);
    else
        setWindowTitle(title);*/
}

void PythonPlugin::pyStdout(const std::string& out) {
    if (compilerLogWidget_)
        compilerLogWidget_->append(QString::fromStdString(out));
}

void PythonPlugin::pyStderr(const std::string& err){
    if (compilerLogWidget_)
        compilerLogWidget_->append(QString::fromStdString(err));
}


} // namespace voreen
