/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2024 University of Muenster, Germany,                        *
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

#include "voreen/qt/mainwindow/menuentities/screenshotmenuentity/screenshotpluginbase.h"

#include "voreen/core/network/processornetwork.h"
#include "voreen/core/datastructures/meta/templatemetadata.h"

#include "tgt/gpucapabilities.h"
#include "voreen/core/voreenapplication.h"

#include <QComboBox>
#include <QGridLayout>
#include <QLabel>
#include <QSpinBox>
#include <QLineEdit>
#include <QToolButton>
#include <QApplication>
#include <QFileDialog>
#include <QUrl>
#include <QVBoxLayout>
#include <QDesktopServices>
#include <QMessageBox>

namespace voreen {


ScreenshotPluginBase::ScreenshotPluginBase(QWidget* parent, Workspace* workspace)
    : QWidget(parent)
    , path_("")
    , currentWorkspace_(workspace)
    , metaDataPrefix_("base")
{
    //set flags
    setWindowTitle(tr("Screenshot"));
    setWindowFlags(Qt::Tool);
    setAttribute(Qt::WA_DeleteOnClose,false);

    //load supported resolutions
    resolutions_.push_back("native");
    resolutions_.push_back("user-defined");
    resolutions_.push_back("512x512");
    resolutions_.push_back("800x600");
    resolutions_.push_back("1024x768");
    resolutions_.push_back("1024x1024");
    resolutions_.push_back("1280x1024");
    resolutions_.push_back("1600x1200");
    resolutions_.push_back("1920x1080");
    resolutions_.push_back("1920x1200");
    resolutions_.push_back("2048x2048");

    //Qt layout here
    createLayout();

    //set to fix size
    adjustSize();
    setFixedSize(sizeHint());

    //make connections
    connect(screenshotButton_, SIGNAL(clicked()), this, SLOT(screenshotButtonPressedSlot()));
    connect(prefixLineEdit_, SIGNAL(textChanged(QString)), this, SLOT(prefixLineEditChangedSlot(QString)));
    connect(userDefinedWidthSpinBox_, SIGNAL(valueChanged(int)), this, SLOT(userDefinedWidthSpinBoxChangedSlot(int)));
    connect(userDefinedHeightSpinBox_, SIGNAL(valueChanged(int)), this, SLOT(userDefinedHeightSpinBoxChangedSlot(int)));
    connect(resolutionComboBox_, SIGNAL(currentIndexChanged(int)), this, SLOT(resolutionComboBoxChangedSlot(int)));
}

ScreenshotPluginBase::~ScreenshotPluginBase() {
}


//--------------------------------------------------------
//      Layout
//--------------------------------------------------------)
void ScreenshotPluginBase::createLayout() {
    QVBoxLayout* vboxLayout = new QVBoxLayout();
    QGridLayout* gridLayout = new QGridLayout();

    // set labels
    gridLayout->addWidget(new QLabel(tr("Preset:")), 0, 0);
    gridLayout->addWidget(new QLabel(tr("Resolution:")), 1, 0);
    gridLayout->addWidget(new QLabel(tr("Prefix:")), 2, 0);

    //add resolution box
    resolutionComboBox_ = new QComboBox();
    resolutionComboBox_->addItems(resolutions_);
    gridLayout->addWidget(resolutionComboBox_, 0, 1, 1, 3);

    //add user-defines resolution
    userDefinedWidthSpinBox_ = new QSpinBox();
    gridLayout->addWidget(userDefinedWidthSpinBox_, 1, 1);
    userDefinedWidthSpinBox_->setRange(2, GpuCaps.getMaxTextureSize());
    userDefinedWidthSpinBox_->setValue(800);
    gridLayout->addWidget(new QLabel("x"), 1, 2, Qt::AlignCenter);
    userDefinedHeightSpinBox_ = new QSpinBox();
    gridLayout->addWidget(userDefinedHeightSpinBox_, 1, 3);
    userDefinedHeightSpinBox_->setRange(2, GpuCaps.getMaxTextureSize());
    userDefinedHeightSpinBox_->setValue(600);

    //add prefix
    prefixLineEdit_ = new QLineEdit();
    prefixLineEdit_->setText("snapshot");
    gridLayout->addWidget(prefixLineEdit_, 2, 1, 1, 3);

    //set stretches (why? to get x small in user defined settings?)
    gridLayout->setColumnStretch(1, 5);
    gridLayout->setColumnStretch(2, 1);
    gridLayout->setColumnStretch(3, 5);
    vboxLayout->addLayout(gridLayout);

    //add screenshot button
    screenshotButton_ = new QToolButton();
    screenshotButton_->setToolButtonStyle(Qt::ToolButtonTextBesideIcon);
    screenshotButton_->setText(tr(" Save snapshot as..."));
    screenshotButton_->setIcon(QIcon(":/qt/icons/saveas.png"));
    screenshotButton_->setToolTip(tr("Save snapshot as..."));
    screenshotButton_->setSizePolicy(QSizePolicy(QSizePolicy::MinimumExpanding, QSizePolicy::Preferred));
    vboxLayout->addWidget(screenshotButton_);
    vboxLayout->addStretch();

    //set layout
    setLayout(vboxLayout);

    //update from workspace
    loadMetaData();
}


//--------------------------------------------------------
//      Callbacks
//--------------------------------------------------------
void ScreenshotPluginBase::resolutionComboBoxChangedSlot(int index) {
    // if user-defined
    if (index == 1) { //user-defined
        userDefinedWidthSpinBox_->setEnabled(true);
        userDefinedHeightSpinBox_->setEnabled(true);
        return;
    }
    else {
        userDefinedWidthSpinBox_->setEnabled(false);
        userDefinedHeightSpinBox_->setEnabled(false);
    }

    userDefinedWidthSpinBox_->blockSignals(true);
    userDefinedHeightSpinBox_->blockSignals(true);

     if(resolutionComboBox_->currentIndex() > 1) { //"bla x blub"
            QString curText = resolutionComboBox_->currentText();
            int xIndex = curText.indexOf("x");
            int width = curText.left(xIndex).toInt();
            int height = curText.right(curText.size()-xIndex-1).toInt();
            userDefinedWidthSpinBox_->setValue(width);
            userDefinedHeightSpinBox_->setValue(height);
        }

    userDefinedWidthSpinBox_->blockSignals(false);
    userDefinedHeightSpinBox_->blockSignals(false);

    //update native settings
    nativeSizeHasChangedSlot();
    //save meta
    saveMetaData();
}

void ScreenshotPluginBase::userDefinedWidthSpinBoxChangedSlot(int /*value*/) {
    if(resolutionComboBox_->currentIndex() == 1) //user-defined
        saveMetaData();
}

void ScreenshotPluginBase::userDefinedHeightSpinBoxChangedSlot(int /*value*/) {
    if(resolutionComboBox_->currentIndex() == 1) //user-defined
        saveMetaData();
}

void ScreenshotPluginBase::prefixLineEditChangedSlot(QString prefix) {
    //save meta in workspace
    saveMetaData();
}

void ScreenshotPluginBase::screenshotButtonPressedSlot() {
    QFileDialog filedialog(this);
    filedialog.setWindowTitle(tr("Save Screenshot"));
    //set to screenshots or last used
    if(path_.isEmpty())
        filedialog.setDirectory(VoreenApplication::app()->getUserDataPath("screenshots").c_str());
    else
        filedialog.setDirectory(path_);

    //set filter
    filedialog.setDefaultSuffix(tr("png"));
    QStringList filter;
    filter << tr("PNG image (*.png)");
    filter << tr("JPEG image (*.jpg)");
    filter << tr("Windows Bitmap (*.bmp)");
    filter << tr("TIFF image (*.tif)");
    filedialog.setNameFilters(filter);
    filedialog.setAcceptMode(QFileDialog::AcceptSave);
    filedialog.setOption(QFileDialog::DontUseNativeDialog);

    //set sidebar
    QList<QUrl> urls;
    urls << QUrl::fromLocalFile(VoreenApplication::app()->getUserDataPath("screenshots").c_str());
    for (auto f : QStandardPaths::standardLocations(QStandardPaths::DesktopLocation))
        urls << QUrl::fromLocalFile(f);
    for (auto f : QStandardPaths::standardLocations(QStandardPaths::HomeLocation))
        urls << QUrl::fromLocalFile(f);
    filedialog.setSidebarUrls(urls);

    //get time and create name
    struct tm* Tm;
    time_t currentTime = time(NULL);
    Tm = localtime(&currentTime);
    std::stringstream timestamp;
    timestamp << prefixLineEdit_->text().toStdString() << "_" << (Tm->tm_year+1900) << "-" << (Tm->tm_mon+1) << "-" << Tm->tm_mday << "-" << Tm->tm_hour << "-" << Tm->tm_min << "-" << Tm->tm_sec;
    timestamp << ".png";
    filedialog.selectFile(tr(timestamp.str().c_str()));

    //open dialog
    QStringList fileList;
    if (filedialog.exec())
        fileList = filedialog.selectedFiles();
    if (fileList.empty())
        return;

    //save lastused path
    QString file = fileList.at(0);
    path_ = filedialog.directory().absolutePath();

    //check, if right file has been specified
    if (!file.endsWith(".jpg", Qt::CaseInsensitive) && !file.endsWith(".png", Qt::CaseInsensitive) &&
        !file.endsWith(".tif", Qt::CaseInsensitive) && !file.endsWith(".bmp", Qt::CaseInsensitive)) {
        QString text = tr("Image file could not be saved.\n");
        int index = file.lastIndexOf(".");
        if ((index == -1) || (index+1 == fileList[0].size()))
            text += tr("No file extension specified.");
        else
            text += tr("Invalid file extension: ") + file.right(file.size() - index - 1);

        QMessageBox::critical(this, tr("Error saving snapshot"), text);
        return;
    }

    //call pure virtual function to actuelly save the image
    saveScreenshot(file, userDefinedWidthSpinBox_->value(), userDefinedHeightSpinBox_->value());
}

//--------------------------------------------------------
//      Meta Handling
//--------------------------------------------------------
void ScreenshotPluginBase::updateWorkspace(Workspace* workspace) {
    currentWorkspace_ = workspace;
    //laod meta
    loadMetaData();
    //update native
    resolutionComboBoxChangedSlot(resolutionComboBox_->currentIndex());
}

void ScreenshotPluginBase::loadMetaData() {
    std::string newLineEditValue = "snapshot";
    int resolutionIndex = 0;
    tgt::ivec2 lastUsedSize = tgt::ivec2(256);
    if(currentWorkspace_) {
        if(ProcessorNetwork* network = currentWorkspace_->getProcessorNetwork()) {
            if(MetaDataBase* meta = network->getMetaDataContainer().getMetaData(metaDataPrefix_ + "_Prefix")) {
                newLineEditValue = static_cast<StringMetaData*>(meta)->getValue();
            }
            if(MetaDataBase* meta = network->getMetaDataContainer().getMetaData(metaDataPrefix_ + "_Index")) {
                resolutionIndex = static_cast<IntMetaData*>(meta)->getValue();
            }
            if(MetaDataBase* meta = network->getMetaDataContainer().getMetaData(metaDataPrefix_ + "_Resolution")) {
                lastUsedSize = static_cast<IVec2MetaData*>(meta)->getValue();
            }
        }
    }
    //set value
    prefixLineEdit_->setText(newLineEditValue.c_str());
    resolutionComboBox_->setCurrentIndex(resolutionIndex);
    if(resolutionIndex == 1) { //user-defined
        userDefinedWidthSpinBox_->setValue(lastUsedSize.x);
        userDefinedHeightSpinBox_->setValue(lastUsedSize.y);
    }

}

void ScreenshotPluginBase::saveMetaData() {
    if(currentWorkspace_) {
        if(ProcessorNetwork* network = currentWorkspace_->getProcessorNetwork()) {
            network->getMetaDataContainer().addMetaData(metaDataPrefix_ + "_Prefix",new StringMetaData(prefixLineEdit_->text().toStdString()));
            network->getMetaDataContainer().addMetaData(metaDataPrefix_ + "_Index",new IntMetaData(resolutionComboBox_->currentIndex()));
            if(resolutionComboBox_->currentIndex() == 1) // user-defined
                network->getMetaDataContainer().addMetaData(metaDataPrefix_ + "_Resolution",new IVec2MetaData(tgt::ivec2(userDefinedWidthSpinBox_->value(),
                                                                                                                         userDefinedHeightSpinBox_->value())));
        }
    }
}

void ScreenshotPluginBase::removeMetaData() {
    if(currentWorkspace_) {
        if(ProcessorNetwork* network = currentWorkspace_->getProcessorNetwork()) {
            network->getMetaDataContainer().removeMetaData(metaDataPrefix_ + "_Prefix");
            network->getMetaDataContainer().removeMetaData(metaDataPrefix_ + "_Index");
            network->getMetaDataContainer().removeMetaData(metaDataPrefix_ + "_Resolution");
        }
    }
}

} // namespace voreen
