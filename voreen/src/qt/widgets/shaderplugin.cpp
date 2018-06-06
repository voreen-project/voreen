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

#include "voreen/qt/widgets/shaderplugin.h"
#include "voreen/core/utils/stringutils.h"

#include "tgt/shadermanager.h"

#include <QFileDialog>
#include <QTextEdit>
#include <QTextStream>
#include <QVBoxLayout>
#include <QLabel>
#include <QShortcut>
#include <QDesktopServices>

namespace voreen {

ShaderTab::ShaderTab(CodeEdit* codeEdit, int sourceComponentIndex)
    : codeEdit_(codeEdit)
    , highlighter_(new GLSLHighlighter(codeEdit_->document()))
    , sourceComponentIndex_(sourceComponentIndex)
{
}

ShaderTab::ShaderTab(ShaderTab&& other)
    : codeEdit_(other.codeEdit_)
    , highlighter_(std::move(other.highlighter_))
    , sourceComponentIndex_(other.sourceComponentIndex_)
{
    other.codeEdit_ = nullptr;
    other.highlighter_ = nullptr;
}

const int ShaderTab::NO_SOURCE_COMPONENT_INDEX = -1;

bool ShaderTab::showsImported() const {
    return sourceComponentIndex_ == NO_SOURCE_COMPONENT_INDEX;
}

ShaderPlugin::ShaderPlugin(ShaderProperty* prop, QWidget* parent)
    : QWidget(parent)
    , property_(prop)
{
}

ShaderPlugin::~ShaderPlugin() {
}

const std::string ShaderPlugin::loggerCat_("voreen.shaderplugin");

void ShaderPlugin::createWidgets() {
    updateBt_ = new QToolButton();
    updateBt_->setIcon(QIcon(":/qt/icons/save.png"));
    updateBt_->setIconSize(QSize(24, 24));
    updateBt_->setShortcut(QKeySequence("Ctrl+R"));
    updateBt_->setToolTip("Update shader (Ctrl+R)");
    updateAllBt_ = new QToolButton();
    updateAllBt_->setIcon(QIcon(":/qt/icons/saveall.png"));
    updateAllBt_->setIconSize(QSize(24, 24));
    //updateAllBt_->setShortcut(QShortcut(QKeySequence("Ctrl+Shift+R"), updateAllBt_, 0, 0, Qt::WidgetWithChildrenShortcut));
    updateAllBt_->setShortcut(QKeySequence("Ctrl+Shift+R"));
    updateAllBt_->setToolTip("Update shaders (Ctrl+Shift+R)");
    undoBt_ = new QToolButton();
    undoBt_->setIcon(QIcon(":/qt/icons/delete.png"));
    undoBt_->setIconSize(QSize(24, 24));
    undoBt_->setToolTip("Revert to last loaded source");
    fullUndoBt_ = new QToolButton();
    fullUndoBt_->setIcon(QIcon(":/qt/icons/revert.png"));
    fullUndoBt_->setIconSize(QSize(24, 24));
    fullUndoBt_->setToolTip("Revert to original source");
    openBt_ = new QToolButton();
    openBt_->setIcon(QIcon(":/qt/icons/open.png"));
    openBt_->setIconSize(QSize(24, 24));
    openBt_->setToolTip("Load shader");
    exportBt_ = new QToolButton();
    exportBt_->setIcon(QIcon(":/qt/icons/saveas.png"));
    exportBt_->setIconSize(QSize(24, 24));
    exportBt_->setToolTip("Export shader");
    fontSizeBox_ = new QSpinBox();
    fontSizeBox_->setMinimum(6);
    fontSizeBox_->setMaximum(24);
    fontSizeBox_->setValue(9);
    fontSizeBox_->setToolTip("Choose font size");
    QLabel* fontSizeLabel = new QLabel("Font Size:");

    QHBoxLayout* hbox = new QHBoxLayout();
    hbox->setContentsMargins(0,0,0,0);
    hbox->addWidget(openBt_);
    hbox->addWidget(updateBt_);
    hbox->addWidget(updateAllBt_);
    hbox->addWidget(exportBt_);
    hbox->addWidget(undoBt_);
    hbox->addWidget(fullUndoBt_);
    hbox->addWidget(fontSizeLabel);
    hbox->addWidget(fontSizeBox_);
    hbox->addStretch();
    QWidget* toolButtonBar = new QWidget();
    toolButtonBar->setLayout(hbox);

    compilerLogWidget_ = new QTextEdit();
    QFont font;
    font.setFamily("Courier");
    font.setFixedPitch(true);
    font.setPointSize(10);
    compilerLogWidget_->setFont(font);
    compilerLogWidget_->setReadOnly(false);
    compilerLogWidget_->setFixedHeight(150);

    QVBoxLayout* vbox = new QVBoxLayout();
    vbox->addWidget(toolButtonBar);

    tabWidget_ = new QTabWidget(this);
    tabWidget_->setMovable(false); // Important! Otherwise the internal shaderTabs_ model would not match anymore
    tabWidget_->setTabsClosable(true);
    vbox->addWidget(tabWidget_);

    vbox->addWidget(compilerLogWidget_);
    setLayout(vbox);

    for(size_t i = 0; i < property_->get().numComponents(); ++i) {
        openSourceComponentTab(i);
    }
}

void ShaderPlugin::createConnections() {
    connect(undoBt_, SIGNAL(clicked()), this, SLOT(undoShader()));
    connect(fullUndoBt_, SIGNAL(clicked()), this, SLOT(fullUndoShader()));
    connect(openBt_, SIGNAL(clicked()), this, SLOT(openShader()));
    connect(exportBt_, SIGNAL(clicked()), this, SLOT(exportShader()));
    connect(updateBt_, SIGNAL(clicked()), this, SLOT(saveShader()));
    connect(updateAllBt_, SIGNAL(clicked()), this, SLOT(saveAllShaders()));
    connect(fontSizeBox_, SIGNAL(valueChanged(int)), this, SLOT(changeFontSize()));
    connect(tabWidget_, SIGNAL(tabCloseRequested(int)), this, SLOT(closeTab(int)));
    connect(tabWidget_, SIGNAL(currentChanged(int)), this, SLOT(adaptButtonsToTab(int)));
}

void ShaderPlugin::urlClicked(const QUrl& link) {
    std::string filename = link.path().toStdString();
    openImportTab(filename);
}
void ShaderPlugin::changeFontSize() {
    for(auto& tab : shaderTabs_) {
        tab.codeEdit_->updateFontSize(fontSizeBox_->value());
    }
}

const QString ShaderPlugin::getOpenFileName(QString filter) {
    QFileDialog fileDialog(this);
    fileDialog.setWindowTitle(tr("Choose a shader to open"));
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

const QString ShaderPlugin::getSaveFileName(QStringList filters) {
    QFileDialog fileDialog(this);
    fileDialog.setWindowTitle(tr("Choose a filename to save shader"));
    fileDialog.setDirectory(VoreenApplication::app()->getUserDataPath("shaders").c_str());
    fileDialog.setNameFilters(filters);
    fileDialog.setAcceptMode(QFileDialog::AcceptSave);
    fileDialog.setOption(QFileDialog::DontUseNativeDialog);

    QList<QUrl> urls;
    urls << QUrl::fromLocalFile(VoreenApplication::app()->getUserDataPath().c_str());
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
QString ShaderPlugin::getFilter(tgt::ShaderObject::ShaderType type) const {
    switch(type) {
        case tgt::ShaderObject::ShaderType::VERTEX_SHADER:
            return "Vertex Shader (*.vert)";
        case tgt::ShaderObject::ShaderType::FRAGMENT_SHADER:
            return "Fragment Shader (*.frag)";
        case tgt::ShaderObject::ShaderType::GEOMETRY_SHADER:
            return "Geometry Shader (*.geom)";
        case tgt::ShaderObject::ShaderType::COMPUTE_SHADER:
            return "Compute Shader (*.comp)";
        default:
            tgtAssert(false, "Unimplemented shader type");
            return "";
    }
}
ShaderTab& ShaderPlugin::getCurrentShaderTab() {
    int index = tabWidget_->currentIndex();
    tgtAssert(index >= 0 && index < static_cast<int>(shaderTabs_.size()), "Invalid tab index");
    auto it = shaderTabs_.begin();
    std::advance(it, index);
    return *it;
}

void ShaderPlugin::undoShader() {
    compilerLogWidget_->clear();
    compilerLogWidget_->insertPlainText("Reverting to last loaded source...\n");

    ShaderTab& tab = getCurrentShaderTab();
    tgtAssert(!tab.showsImported(), "Undoing imported shader");
    property_->getMutator()->at(tab.sourceComponentIndex_).reloadSource();

    updateTabTexts();
}

void ShaderPlugin::fullUndoShader() {
    compilerLogWidget_->clear();
    compilerLogWidget_->insertPlainText("Reverting to original loaded source...\n");

    ShaderTab& tab = getCurrentShaderTab();
    tgtAssert(!tab.showsImported(), "Undoing imported shader");
    property_->getMutator()->at(tab.sourceComponentIndex_).reset();

    updateTabTexts();
}

void ShaderPlugin::openShader() {
    ShaderTab& tab = getCurrentShaderTab();
    CodeEdit* curEdit = tab.codeEdit_;
    tgtAssert(!tab.showsImported(), "Opening imported shader");

    QString fileName = getOpenFileName(getFilter(property_->get()[tab.sourceComponentIndex_].getType()));
    if (!fileName.isEmpty()) {
        property_->getMutator()->at(tab.sourceComponentIndex_).setExternalFilename(fileName.toStdString());
    }

    updateTabTexts();
}


void ShaderPlugin::exportShader() {
    ShaderTab& tab = getCurrentShaderTab();
    CodeEdit* curEdit = tab.codeEdit_;
    tgtAssert(!tab.showsImported(), "Saving imported shader");

    // create filter with supported file formats
    QStringList filter;
    filter << getFilter(property_->get()[tab.sourceComponentIndex_].getType());

    QString fileName = getSaveFileName(filter);
    if (!fileName.isEmpty()) {
        // save shader to disk
        QFile outputFile(fileName);
        outputFile.open(QIODevice::WriteOnly);
#ifdef _WIN32
        std::string outputStr = convertNewlinesUnixToWindows(curEdit->toPlainText().toStdString());
#else
        std::string outputStr = curEdit->toPlainText().toStdString();
#endif
        outputFile.write(outputStr.c_str(), outputStr.size());
        outputFile.close();

        // The new saved shader shall be the new external filename
        property_->getMutator()->at(tab.sourceComponentIndex_).setExternalFilename(fileName.toStdString());
    }

    updateTabTexts();
}

void ShaderPlugin::saveShader() {
    compilerLogWidget_->clear();
    compilerLogWidget_->insertPlainText("Saving current shader...\n");

    const ShaderTab& tab = getCurrentShaderTab();
    if(!tab.showsImported()) {
        property_->getMutator()->at(tab.sourceComponentIndex_).setSource(tab.codeEdit_->toPlainText().toStdString());
    }

    updateTabTexts();
}

void ShaderPlugin::saveAllShaders() {
    compilerLogWidget_->clear();
    compilerLogWidget_->insertPlainText("Saving all shaders...\n");

    auto mut = property_->getMutator();
    for(auto& tab : shaderTabs_) {
        if(!tab.showsImported()) {
            mut->at(tab.sourceComponentIndex_).setSource(tab.codeEdit_->toPlainText().toStdString());
        }
    }

    updateTabTexts();
}

void ShaderPlugin::updateFromProperty() {
    const ShaderSource& source = property_->get();

    // We only want to show compile messages if the shader has been rebuilt already.
    if(!property_->requiresRebuild()) {
        compilerLogWidget_->insertPlainText("Compiling...\n");
        for(auto& tab : shaderTabs_) {
            if(!tab.showsImported()) {
                const std::string currentSrc = source[tab.sourceComponentIndex_].getSource();
                if(tab.codeEdit_->toPlainText() != QString(currentSrc.c_str())) {
                    tab.codeEdit_->setShaderSource(currentSrc);
                }
            }
        }
        tgt::Shader* shader = property_->getShader();
        for(tgt::ShaderObject* shaderObj : shader->getObjects()) {
            std::string compilerLog = shaderObj->getCompilerLog();
            compilerLogWidget_->insertPlainText(QString(compilerLog.c_str()));
        }
        if(shader) {
            std::string linkerLog = shader->getLinkerLog();
            compilerLogWidget_->insertPlainText(QString(linkerLog.c_str()));
        }
        if(property_->hasValidShader()) {
            compilerLogWidget_->insertPlainText("Compilation successful!\n");
        } else {
            compilerLogWidget_->insertPlainText("Compilation failed!\n");
        }
    }

    tgtAssert(property_->get().numComponents() > 0, "No shader source components");

    const ShaderSourceComponent& firstComponent = property_->get()[0];
    std::string mod = "";
    if(firstComponent.isModified())
        mod = "original source: ";

    std::string windowTitle = "";
    if (property_->getOwner())
        windowTitle += property_->getOwner()->getGuiName() + " - ";
    windowTitle += property_->getGuiName() + " (" + mod + firstComponent.getCurrentFileName() + ")";
    window()->setWindowTitle(QString::fromStdString(windowTitle));

    updateTabTexts();
}

void ShaderPlugin::openImportTab(const std::string& filename) {
    CodeEdit* codeEdit = new CodeEdit(true, fontSizeBox_->value());

    tgtAssert(!filename.empty(), "filename empty");
    tgtAssert(tgt::Singleton<tgt::ShaderManager>::isInited(), "ShaderManager not instantiated");
    tgtAssert(tgt::Singleton<tgt::FileSystem>::isInited(), "FileSystem not instantiated");

    std::string completeFilename = ShdrMgr.completePath(filename);

    std::unique_ptr<tgt::File> file = std::unique_ptr<tgt::File>(FileSys.open(completeFilename));

    if (!file || !file->isOpen()) {
        LERROR("File not found: " << filename);
        return;
    }

    codeEdit->setShaderSource(file->getAsString());

    openTab(filename, ShaderTab(codeEdit));
}

void ShaderPlugin::openSourceComponentTab(int sourceComponentIndex) {
    tgtAssert(sourceComponentIndex != ShaderTab::NO_SOURCE_COMPONENT_INDEX, "invalid sourceComponentIndex");
    CodeEdit* codeEdit = new CodeEdit(false, fontSizeBox_->value());
    const ShaderSourceComponent& component = property_->get()[sourceComponentIndex];
    codeEdit->setShaderSource(component.getSource());

    openTab(component.getCurrentFileName(), ShaderTab(codeEdit, sourceComponentIndex));
}

void ShaderPlugin::openTab(const std::string& filename, ShaderTab&& newTab) {
    QWidget* tab = new QWidget();
    tabWidget_->addTab(tab, QString((filename +
                (newTab.showsImported() ? " [readonly]" : "")).c_str()));

    QVBoxLayout* layout = new QVBoxLayout(tab);

    tab->setLayout(layout);

    layout->addWidget(newTab.codeEdit_);

    connect(newTab.codeEdit_, SIGNAL(anchorClicked(const QUrl&)), this, SLOT(urlClicked(const QUrl&)));
    connect(newTab.codeEdit_, SIGNAL(textChanged()), this, SIGNAL(modified()));

    shaderTabs_.emplace_back(std::move(newTab));

    tabWidget_->setCurrentIndex(shaderTabs_.size()-1);
}
void ShaderPlugin::updateTabTexts() {
    int i = 0;
    for(auto& shaderTab : shaderTabs_) {
        tgtAssert(i < tabWidget_->count(), "tabWidget_ <-> shaderTabs_ size mismatch");
        if(!shaderTab.showsImported()) {
            const ShaderSourceComponent& ssc = property_->get()[shaderTab.sourceComponentIndex_];
            std::string text = ssc.getCurrentFileName();
            if(ssc.isModified()) {
                text += " [modified]";
            }
            if(ssc.isExternal()) {
                text += " [external]";
            }
            tabWidget_->setTabText(i, QString(text.c_str()));
        }
        ++i;
    }
}

void ShaderPlugin::closeTab(int index) {
    auto it = shaderTabs_.begin();
    std::advance(it, index);
    if(it->showsImported()) {
        shaderTabs_.erase(it);
        tabWidget_->removeTab(index);
    } else {
        LWARNING("Non imported Tabs cannot be closed.");
    }
}

void ShaderPlugin::adaptButtonsToTab(int index) {
    bool imported = getCurrentShaderTab().showsImported();
    undoBt_->setEnabled(!imported);
    fullUndoBt_->setEnabled(!imported);
    openBt_->setEnabled(!imported);
    updateBt_->setEnabled(!imported);
    updateAllBt_->setEnabled(!imported);
    exportBt_->setEnabled(!imported);
}

} // namespace voreen
