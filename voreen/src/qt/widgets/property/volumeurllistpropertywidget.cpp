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

#include "voreen/core/datastructures/volume/volumelist.h"

#include "voreen/core/properties/volumeurllistproperty.h"
#include "voreen/core/io/volumereader.h"
#include "voreen/core/io/volumeserializer.h"
#include "voreen/core/io/volumeserializerpopulator.h"

#include "voreen/qt/widgets/volumeviewhelper.h"
#include "voreen/qt/widgets/property/volumeurllistpropertywidget.h"

#include <QLabel>
#include <QPushButton>
#include <QVBoxLayout>

namespace {
#ifdef __APPLE__
    int fontSize = 13;
#else
    int fontSize = 8;
#endif
}

namespace voreen {

const std::string VolumeURLListPropertyWidget::loggerCat_("voreenqt.VolumeURLListPropertyWidget");

VolumeURLListPropertyWidget::VolumeURLListPropertyWidget(VolumeURLListProperty* volumeListProp, QWidget* parent)
    : QPropertyWidget(volumeListProp, parent, false)
    , volumeIOHelper_(parent, VolumeIOHelper::MULTI_FILE)
    , fileWatchCheckBox_(nullptr)
{
    urlListProperty_ = volumeListProp;
    tgtAssert(urlListProperty_, "No volume collection property");

    setFocusPolicy(Qt::StrongFocus);
    QVBoxLayout* mainLayout = new QVBoxLayout();
    layout_->addLayout(mainLayout);

    QHBoxLayout* buttonLayout = new QHBoxLayout();
    loadButton_ = new QPushButton(tr("Load Volumes..."));
    loadButton_->setIcon(QPixmap(":/qt/icons/open-volume.png"));
    loadButton_->setMinimumWidth(110);
    clearButton_ = new QPushButton(tr("Clear Volumes"));
    buttonLayout->addWidget(loadButton_);
    buttonLayout->addWidget(clearButton_);
    mainLayout->addLayout(buttonLayout);

    if (urlListProperty_->getWatchMode() != VoreenFileWatchListener::ALWAYS_OFF) {
        fileWatchCheckBox_ = new QCheckBox(tr("Watch File Change?"), this);
        fileWatchCheckBox_->setChecked(volumeListProp->isFileWatchEnabled());
        if(urlListProperty_->isFileWatchEditable())
            connect(fileWatchCheckBox_, SIGNAL(clicked(bool)), this, SLOT(changeFileWatchOption(bool)));
        else
            fileWatchCheckBox_->setEnabled(false);
        mainLayout->addWidget(fileWatchCheckBox_);
    }

    selectAll_ = new QCheckBox(tr("Select All"), this);
    //selectAll_->move(8, 0);
    mainLayout->addWidget(selectAll_);

    volumeTreeWidget_ = new QTreeWidget(this);
    QTreeWidgetItem* header = volumeTreeWidget_->headerItem();
    header->setText(0, tr(""));
    volumeTreeWidget_->setColumnCount(1);
    volumeTreeWidget_->show();
    volumeTreeWidget_->setIconSize(QSize(50,50));
    volumeTreeWidget_->setSizePolicy(QSizePolicy::Minimum, QSizePolicy::Minimum);
    setSizePolicy(QSizePolicy::Minimum, QSizePolicy::Minimum);

    connect(loadButton_, SIGNAL(clicked()), &volumeIOHelper_, SLOT(showFileOpenDialog()));
    connect(&volumeIOHelper_, SIGNAL(volumeLoaded(const VolumeBase*)), this, SLOT(volumeLoaded(const VolumeBase*)));
    connect(clearButton_, SIGNAL(clicked()), this, SLOT(clearVolumes()));

    connect(volumeTreeWidget_, SIGNAL(itemClicked(QTreeWidgetItem*, int)), this, SLOT(itemSelected(QTreeWidgetItem*, int)));
    connect(volumeTreeWidget_, SIGNAL(itemClicked(QTreeWidgetItem*, int)), this, SIGNAL(widgetChanged()));
    mainLayout->addWidget(volumeTreeWidget_);

    connect(selectAll_, SIGNAL(toggled(bool)), this, SLOT(selectAll(bool)));

    updateFromPropertySlot();
}

CustomLabel* VolumeURLListPropertyWidget::getOrCreateNameLabel() const {
    return 0;
}

void VolumeURLListPropertyWidget::updateFromPropertySlot() {
    if (!urlListProperty_)
        return;
    if (!volumeTreeWidget_->updatesEnabled())
        return;

    VolumeList* collection = urlListProperty_->getVolumes(false);
    tgtAssert(collection, "null pointer returned");

    volumeTreeWidget_->clear();

    bool supportFileWatching = fileWatchCheckBox_ != nullptr && urlListProperty_->isFileWatchEditable() && !collection->empty();
    int numSelected = 0;
    for (size_t i = 0 ; i< collection->size(); i++) {
        VolumeBase* handle = collection->at(i);
        std::string url = handle->getOrigin().getURL();

        // Only allow file watching, if all volume formats are watchable.
        if (supportFileWatching) {
            std::string preferredReader = handle->getOrigin().getSearchParameter("preferredReader");
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

        QTreeWidgetItem* qtwi = new QTreeWidgetItem(volumeTreeWidget_);
        qtwi->setFont(0, QFont(QString("Arial"), fontSize));
        QString info = QString::fromStdString(VolumeViewHelper::getStrippedVolumeName(handle)
            + "\n"
            + VolumeViewHelper::getVolumePath(handle));
        qtwi->setText(0, info);
        if(urlListProperty_->getPreviewsVisible())
            qtwi->setIcon(0, QIcon(VolumeViewHelper::generateBorderedPreview(handle, 27, 0)));
        qtwi->setSizeHint(0,QSize(27,27));
        qtwi->setFlags(static_cast<Qt::ItemFlags>(Qt::ItemIsUserCheckable | Qt::ItemIsSelectable | Qt::ItemIsEnabled));

        // set tree widget to checked, if the corresponding volume is contained by the property's collection
        bool selected = urlListProperty_->isSelected(url);
        qtwi->setCheckState(0, selected ? Qt::Checked : Qt::Unchecked);
        if (selected)
            numSelected++;

        volumeTreeWidget_->addTopLevelItem(qtwi);
    }

    // Handle file watching ability.
    if (supportFileWatching) {
        fileWatchCheckBox_->setEnabled(true);
    }
    else {
        urlListProperty_->setFileWatchEnabled(false);
        fileWatchCheckBox_->setChecked(false);
        fileWatchCheckBox_->setEnabled(false);
    }

    clearButton_->setEnabled(!collection->empty());
    selectAll_->setEnabled(!collection->empty());
    if (numSelected == 0)
        selectAll_->setChecked(false);
    else if (numSelected == (int)collection->size())
        selectAll_->setChecked(true);

    delete collection;
}

void VolumeURLListPropertyWidget::updateSelection() {
    if (!prop_)
        return;

    VolumeList* collection = urlListProperty_->getVolumes(false);
    tgtAssert(collection, "null pointer returned");

    volumeTreeWidget_->setUpdatesEnabled(false);

    QList<QTreeWidgetItem*> items = volumeTreeWidget_->findItems("", Qt::MatchContains);
    for(size_t i = 0; i < collection->size(); i++) {
        bool selected = items.at(static_cast<int>(i))->checkState(0) == Qt::Checked;
        urlListProperty_->setSelected(collection->at(i)->getOrigin().getURL(), selected);
    }
    volumeTreeWidget_->setUpdatesEnabled(true);

    delete collection;
}

void VolumeURLListPropertyWidget::changeFileWatchOption(bool checked) {
    urlListProperty_->setFileWatchEnabled(checked);
}

void VolumeURLListPropertyWidget::volumeLoaded(const VolumeBase* handle) {
    urlListProperty_->addVolume(const_cast<VolumeBase*>(handle), true, true);
}

void VolumeURLListPropertyWidget::clearVolumes() {
    urlListProperty_->clear();
}

void VolumeURLListPropertyWidget::itemSelected(QTreeWidgetItem*, int) {
    updateSelection();
}

void VolumeURLListPropertyWidget::selectAll(bool toggle) {
    urlListProperty_->setAllSelected(toggle);
}

} //namespace
