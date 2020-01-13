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

#include "voreen/qt/widgets/property/volumeurlpropertywidget.h"

#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/voreenapplication.h"
#include "voreen/core/io/volumereader.h"
#include "voreen/core/io/volumeserializer.h"
#include "voreen/core/io/volumeserializerpopulator.h"

#include "voreen/qt/widgets/volumeviewhelper.h"
#include "voreen/qt/voreenapplicationqt.h"

#include <QMessageBox>
#include <QMenu>
#include <QVBoxLayout>
#include <QCheckBox>

#include "tgt/vector.h"
#include "tgt/matrix.h"

#include <math.h>

#ifdef VRN_MODULE_GDCM
#include "modules/gdcm/qt/dicomconnectiondialog.h"
#endif

namespace voreen {


const std::string VolumeURLPropertyWidget::loggerCat_("voreen.qt.VolumeHandlePropertyWidget");

VolumeURLPropertyWidget::VolumeURLPropertyWidget(VolumeURLProperty* volumeHandleProp, QWidget* parent)
    : QPropertyWidget(volumeHandleProp, parent, false)
    , volumeIOHelper_(parent, VolumeIOHelper::SINGLE_FILE)
    , fileWatchCheckBox_(nullptr)
#ifdef VRN_MODULE_GDCM
    , dicomConnectionDialog_(new DicomConnectionDialog(this))
#else
    , dicomConnectionDialog_(0)
#endif
{

    if (!volumeHandleProp) {
        tgtAssert(false, "No volume property");
        LERROR("No volume property");
        return;
    }
    QVBoxLayout* mainLayout = new QVBoxLayout();
    //mainLayout->setAlignment(Qt::AlignLeft);
    //mainLayout->setContentsMargins(0, 2, 0, 0);

    layout_->addLayout(mainLayout);

    //loadButton_ = new QToolButton(tr("Load Volume..."));
    loadButton_ = new QToolButton();
    loadButton_->setIcon(QPixmap(":/qt/icons/open-volume.png"));
    loadButton_->setText(tr("Load Volume..."));
    loadButton_->setToolButtonStyle(Qt::ToolButtonTextBesideIcon);
    loadButton_->setPopupMode(QToolButton::MenuButtonPopup);

    QMenu* loadMenu = new QMenu(loadButton_);
    loadButton_->setMenu(loadMenu);
    QAction* loadVolumeAction = new QAction(tr("Load Volume..."), loadButton_);
    loadVolumeAction->setIcon(QPixmap(":/qt/icons/open-volume.png"));
    loadMenu->addAction(loadVolumeAction);
#ifdef VRN_MODULE_GDCM
    QAction* dicomServerAction = new QAction(tr("Load from DICOM Server..."), loadButton_);
    //dicomServerAction->setIcon(QPixmap(":/qt/icons/open-volume.png"));
    loadMenu->addAction(dicomServerAction);
#endif
    clearButton_ = new QPushButton(tr("Clear Volume"));

    QHBoxLayout* volumeLayout = new QHBoxLayout();
    volumeLayout->setContentsMargins(0, 0, 4, 0);
    volumeLayout->setSpacing(4);
    volumeLayout->setMargin(0);
    volumeLayout->addWidget(loadButton_, 1);
    volumeLayout->addWidget(clearButton_, 1);

    loadButton_->setSizePolicy(QSizePolicy::MinimumExpanding, QSizePolicy::Minimum);
    loadButton_->setMinimumWidth(105);
    clearButton_->setSizePolicy(QSizePolicy::MinimumExpanding, QSizePolicy::Minimum);
    //clearButton_->setMinimumWidth(105);

    mainLayout->addLayout(volumeLayout);

    if (volumeHandleProp->getWatchMode() != VoreenFileWatchListener::ALWAYS_OFF) {
        fileWatchCheckBox_ = new QCheckBox(tr("Watch File Change?"), this);
        fileWatchCheckBox_->setChecked(volumeHandleProp->isFileWatchEnabled());
        if(volumeHandleProp->isFileWatchEditable())
            connect(fileWatchCheckBox_, SIGNAL(clicked(bool)), this, SLOT(changeFileWatchOption(bool)));
        else
            fileWatchCheckBox_->setEnabled(false);
        mainLayout->addWidget(fileWatchCheckBox_);
    }

    connect(loadButton_, SIGNAL(clicked()),
        &volumeIOHelper_, SLOT(showFileOpenDialog()));
    connect(loadVolumeAction, SIGNAL(triggered()),
        &volumeIOHelper_, SLOT(showFileOpenDialog()));
    connect(&volumeIOHelper_, SIGNAL(volumeLoaded(const VolumeBase*)),
        this, SLOT(volumeLoaded(const VolumeBase*)));
    connect(clearButton_, SIGNAL(clicked()),
        this, SLOT(clearVolume()));

#ifdef VRN_MODULE_GDCM
    tgtAssert(dicomConnectionDialog_, "DicomConnectionDialog is null");
    connect(dicomServerAction, SIGNAL(triggered()),
        this, SLOT(showDicomConnectionDialog()));
    connect(dicomConnectionDialog_, SIGNAL(loadFromURL(const std::string&, VolumeReader*)),
        &volumeIOHelper_, SLOT(loadURL(const std::string&, VolumeReader*)));
#endif

    updateFromPropertySlot();
}

void VolumeURLPropertyWidget::updateFromPropertySlot() {
    VolumeURLProperty* handleProp = dynamic_cast<VolumeURLProperty*>(prop_);
    if (!handleProp) {
        LWARNING("No volume property");
        return;
    }

    bool supportFileWatching = fileWatchCheckBox_ != nullptr && handleProp->isFileWatchEditable() && !handleProp->getURL().empty();
    if (supportFileWatching) {
        std::string preferredReader = VolumeURL(handleProp->getURL()).getSearchParameter("preferredReader");
        if (preferredReader.empty()) {
            if (fileWatchCheckBox_->isEnabled())
                LWARNING("No preferred reader set. Disabling file watching.");
            supportFileWatching = false;
        }
        else {
            VolumeSerializerPopulator populator;
            VolumeReader* reader = populator.getVolumeSerializer()->getReaderByName(preferredReader);
            supportFileWatching = reader->canSupportFileWatching();
        }
    }

    // Handle file watching ability.
    if (supportFileWatching) {
        fileWatchCheckBox_->setEnabled(true);
    }
    else {
        handleProp->setFileWatchEnabled(false);
        fileWatchCheckBox_->setChecked(false);
        fileWatchCheckBox_->setEnabled(false);
    }

    VolumeBase* handle = handleProp->getVolume();
    if (handle) {
        clearButton_->setEnabled(true);
    }
    else {
        clearButton_->setEnabled(false);
    }
}

void VolumeURLPropertyWidget::changeFileWatchOption(bool checked) {
    dynamic_cast<VolumeURLProperty*>(prop_)->setFileWatchEnabled(checked);
}

void VolumeURLPropertyWidget::showNameLabel(bool) {
    if (nameLabel_)
        nameLabel_->hide();
}

CustomLabel* VolumeURLPropertyWidget::getOrCreateNameLabel() const {
    return 0;
}

// private slot
void VolumeURLPropertyWidget::volumeLoaded(const VolumeBase* handle) {
    tgtAssert(handle, "null pointer passed");

    VolumeURLProperty* handleProp = dynamic_cast<VolumeURLProperty*>(prop_);
    if (!handleProp) {
        tgtAssert(false, "No volume property");
        LERROR("No volume property");
        return;
    }

    handleProp->setVolume(const_cast<VolumeBase*>(handle), true);
}

void VolumeURLPropertyWidget::clearVolume() {
    VolumeURLProperty* handleProp = dynamic_cast<VolumeURLProperty*>(prop_);
    if (!handleProp) {
        tgtAssert(false, "No volume property");
        LERROR("No volume property");
        return;
    }

    handleProp->setVolume(0);
}

void VolumeURLPropertyWidget::showDicomConnectionDialog() {
#ifdef VRN_MODULE_GDCM
    if (!dicomConnectionDialog_)
        dicomConnectionDialog_ = new DicomConnectionDialog(this);
    dicomConnectionDialog_->show();
#else
    QMessageBox::information(this, "DICOM Server Connection",
        "Loading DICOM data from a server requires the GDCM module.");
#endif
}

} //namespace voreen
