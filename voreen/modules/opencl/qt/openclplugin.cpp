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

#include "modules/opencl/qt/openclplugin.h"
#include "modules/opencl/openclmodule.h"

#include <QFileDialog>
#include <QTextEdit>
#include <QTextStream>
#include <QVBoxLayout>

namespace voreen {

OpenCLPlugin::OpenCLPlugin(OpenCLProperty* prop, QWidget* parent)
    : QWidget(parent)
    , property_(prop)
    , highlighter_(nullptr)
{
}

OpenCLPlugin::~OpenCLPlugin() {
    delete highlighter_;
}

void OpenCLPlugin::createWidgets() {
    updateBt_ = new QToolButton();
    updateBt_->setIcon(QIcon(":/qt/icons/rebuild.png"));
    updateBt_->setIconSize(QSize(24, 24));
    updateBt_->setShortcut(QKeySequence("Ctrl+R"));
    updateBt_->setToolTip("Update program (Ctrl+R)");
    undoBt_ = new QToolButton();
    undoBt_->setIcon(QIcon(":/qt/icons/revert.png"));
    undoBt_->setIconSize(QSize(24, 24));
    undoBt_->setShortcut(QKeySequence("Ctrl+U"));
    undoBt_->setToolTip("Revert to last loaded source");
    fullUndoBt_ = new QToolButton();
    fullUndoBt_->setIcon(QIcon(":/qt/icons/revert.png"));
    fullUndoBt_->setIconSize(QSize(24, 24));
    fullUndoBt_->setToolTip("Revert to original source");
    openBt_ = new QToolButton();
    openBt_->setIcon(QIcon(":/qt/icons/open.png"));
    openBt_->setIconSize(QSize(24, 24));
    openBt_->setToolTip("Load program");
    saveBt_ = new QToolButton();
    saveBt_->setIcon(QIcon(":/qt/icons/save.png"));
    saveBt_->setIconSize(QSize(24, 24));
    saveBt_->setToolTip("Export program");
    fontSizeBox_ = new QSpinBox();
    fontSizeBox_->setMinimum(6);
    fontSizeBox_->setMaximum(24);
    fontSizeBox_->setValue(9);
    fontSizeBox_->setToolTip("Choose font size");

    QHBoxLayout* hbox = new QHBoxLayout();
    hbox->setContentsMargins(0,0,0,0);
    hbox->addWidget(updateBt_);
    hbox->addWidget(undoBt_);
    hbox->addWidget(fullUndoBt_);
    hbox->addWidget(openBt_);
    hbox->addWidget(saveBt_);
    hbox->addWidget(fontSizeBox_);
    hbox->addStretch();
    QWidget* toolButtonBar = new QWidget();
    toolButtonBar->setLayout(hbox);

    codeEdit_ = new CodeEdit(false, fontSizeBox_->value());
    highlighter_ = new OpenCLHighlighter(codeEdit_->document());

    compilerLogWidget_ = new QTextEdit();
    QFont font;
    font.setFamily("Courier");
    font.setFixedPitch(true);
    font.setPointSize(10);
    compilerLogWidget_->setFont(font);
    compilerLogWidget_->setReadOnly(true);
    compilerLogWidget_->setFixedHeight(150);

    QVBoxLayout* vbox = new QVBoxLayout();
    vbox->addWidget(toolButtonBar);
    vbox->addWidget(codeEdit_);
    vbox->addWidget(compilerLogWidget_);
    setLayout(vbox);

    updateFromProperty();
}

void OpenCLPlugin::createConnections() {
    connect(undoBt_, SIGNAL(clicked()), this, SLOT(undoProgram()));
    connect(fullUndoBt_, SIGNAL(clicked()), this, SLOT(fullUndoProgram()));
    connect(openBt_, SIGNAL(clicked()), this, SLOT(openProgram()));
    connect(saveBt_, SIGNAL(clicked()), this, SLOT(saveProgram()));
    connect(updateBt_, SIGNAL(clicked()), this, SLOT(setProperty()));
    //connect(codeEdit_, SIGNAL(textChanged()), this, SIGNAL(modified()));
    connect(fontSizeBox_, SIGNAL(valueChanged(int)), this, SLOT(changeFontSize()));
}

void OpenCLPlugin::changeFontSize() {
    codeEdit_->updateFontSize(fontSizeBox_->value());
}

const QString OpenCLPlugin::getOpenFileName(QString filter) {
    QFileDialog fileDialog(this);
    fileDialog.setWindowTitle(tr("Choose a program to open"));
    fileDialog.setDirectory(VoreenApplication::app()->getBasePath().c_str());
    fileDialog.setNameFilter(filter);
    fileDialog.setOption(QFileDialog::DontUseNativeDialog);

    QList<QUrl> urls;
    urls << QUrl::fromLocalFile(VoreenApplication::app()->getBasePath().c_str());
    fileDialog.setSidebarUrls(urls);

    if (fileDialog.exec() && !fileDialog.selectedFiles().empty()) {
        return fileDialog.selectedFiles()[0];
    }

    return QString();
}

const QString OpenCLPlugin::getSaveFileName(QStringList filters) {
    QFileDialog fileDialog(this);
    fileDialog.setWindowTitle(tr("Choose a filename to save program"));
    fileDialog.setDirectory(VoreenApplication::app()->getBasePath().c_str());
    fileDialog.setNameFilters(filters);
    fileDialog.setAcceptMode(QFileDialog::AcceptSave);
    fileDialog.setOption(QFileDialog::DontUseNativeDialog);

    QList<QUrl> urls;
    urls << QUrl::fromLocalFile(VoreenApplication::app()->getBasePath().c_str());
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

void OpenCLPlugin::undoProgram() {
    property_->resetProgramSource();
}

void OpenCLPlugin::fullUndoProgram() {
    property_->resetProgramFilename();
}

void OpenCLPlugin::openProgram() {
    //create filter with supported file formats
    QString filter = "OpenCL Program (*.cl)";
    QString fileName = getOpenFileName(filter);
    if (!fileName.isEmpty()) {
        property_->setProgramFilename(fileName.toStdString());
    }
}

void OpenCLPlugin::saveProgram() {
    //create filter with supported file formats
    QStringList filter;
    filter << "OpenCL Program (*.cl)";

    QString fileName = getSaveFileName(filter);
    if (!fileName.isEmpty()) {
        //save program to disk
        QFile outputFile(fileName);
        outputFile.open(QIODevice::WriteOnly);
        outputFile.write(codeEdit_->toPlainText().toStdString().c_str(), codeEdit_->toPlainText().size());
        outputFile.close();
    }
}

void OpenCLPlugin::setProperty() {
    property_->setProgramSource(codeEdit_->toPlainText().toStdString());
}

void OpenCLPlugin::updateFromProperty() {
    if(codeEdit_->toPlainText() != QString(property_->get().programSource_.c_str()))
        codeEdit_->setPlainText(property_->get().programSource_.c_str());

    const cl::Program* prog = property_->getProgram();
    if(prog)
        compilerLogWidget_->setText(prog->getBuildLog(OpenCLModule::getInstance()->getCLDevice()).c_str());

    std::string mod = "";
    if(property_->get().programModified_)
        mod = "original source: ";

    std::string windowTitle = "";
    if (property_->getOwner())
        windowTitle += property_->getOwner()->getGuiName() + " - ";
    windowTitle += property_->getGuiName() + " (" + mod + property_->get().programFilename_ + ")";
    window()->setWindowTitle(QString::fromStdString(windowTitle));
}

} // namespace voreen
