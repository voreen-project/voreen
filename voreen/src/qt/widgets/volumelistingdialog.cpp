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

#include "voreen/qt/widgets/volumelistingdialog.h"

#include "voreen/core/datastructures/meta/metadatabase.h"
#include "voreen/core/utils/stringutils.h"

#include <QTableWidget>
#include <QComboBox>
#include <QLineEdit>
#include <QPushButton>
#include <QLabel>
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QFontMetrics>

#include <set>

namespace voreen {

VolumeListingDialog::VolumeListingDialog(QWidget* parent, SelectionMode selectionMode)
    : QDialog(parent)
    , selectionMode_(selectionMode)
    , reader_(0)
{
    setWindowFlags(Qt::Dialog);
    setWindowModality(Qt::WindowModal);
    setWindowTitle(tr("Volume Selection"));

    // table configuration
    table_ = new QTableWidget();
    table_->setSelectionBehavior(QAbstractItemView::SelectRows);
    table_->setSelectionMode(selectionMode == SINGLE_SELECTION ? QAbstractItemView::SingleSelection : QAbstractItemView::MultiSelection);
    table_->setSortingEnabled(true);
    table_->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);

    // top label
    QVBoxLayout* mainLayout = new QVBoxLayout();
    titleLabel_ = new QLabel(tr("The selected file contains multiple volumes. Please select the volumes to load:"));
    QFont labelFont = titleLabel_->font();
    labelFont.setPointSize(titleLabel_->font().pointSize()+1);
    labelFont.setBold(true);
    titleLabel_->setFont(labelFont);
    mainLayout->addWidget(titleLabel_);

    // filter bar
    QHBoxLayout* searchLayout = new QHBoxLayout();
    filterLabel_ = new QLabel(tr("Filter: "));
    searchLayout->addWidget(filterLabel_);
    filterTextBox_ = new QLineEdit();
    searchLayout->addWidget(filterTextBox_);
    comboBoxFilterAttribute_ = new QComboBox();
    comboBoxFilterAttribute_->setMinimumWidth(100);
    searchLayout->addWidget(comboBoxFilterAttribute_);
    mainLayout->addLayout(searchLayout);

    // table
    mainLayout->addWidget(table_);

    // buttons
    QHBoxLayout* buttonLayout = new QHBoxLayout();
    buttonLayout->addStretch();
    cancelButton_ = new QPushButton(tr("Cancel"));
    cancelButton_->setFixedWidth(100);
    buttonLayout->addWidget(cancelButton_);
    if (selectionMode_ == MULTI_SELECTION) {
        selectAllButton_ = new QPushButton(tr("Select All"));
        selectAllButton_->setFixedWidth(100);
        buttonLayout->addWidget(selectAllButton_);
    }
    else
        selectAllButton_ = 0;
    loadButton_ = new QPushButton(tr("Load"));
    loadButton_->setFixedWidth(100);
    QFont loadFont(loadButton_->font());
    loadFont.setBold(true);
    loadButton_->setFont(loadFont);
    buttonLayout->addWidget(loadButton_);
    mainLayout->addLayout(buttonLayout);

    setLayout(mainLayout);

    connect(table_, SIGNAL(itemSelectionChanged()), this, SLOT(updateGuiState()));
    connect(table_, SIGNAL(cellDoubleClicked(int, int)), this, SLOT(cellDoubleClicked()));

    connect(filterTextBox_, SIGNAL(textChanged(const QString&)), this, SLOT(updateTableRows()));
    connect(filterTextBox_, SIGNAL(editingFinished()), this, SLOT(updateTableRows()));
    connect(comboBoxFilterAttribute_, SIGNAL(currentIndexChanged(int)), this, SLOT(updateTableRows()));

    connect(loadButton_, SIGNAL(clicked()), this, SLOT(loadClicked()));
    if (selectionMode_ == MULTI_SELECTION)
        connect(selectAllButton_, SIGNAL(clicked()), this, SLOT(selectAllClicked()));
    connect(cancelButton_, SIGNAL(clicked()), this, SLOT(cancelClicked()));

    resize(600, 300);
    updateGuiState();
}

void VolumeListingDialog::setOrigins(const std::vector<VolumeURL>& origins, VolumeReader* reader) {
    origins_ = origins;
    reader_ = reader;

    QString labelText;
    if (selectionMode_ == SINGLE_SELECTION)
        labelText = tr("The selected file contains %1 volumes. Please select the volume to load:").arg(origins.size());
    else {
        tgtAssert(selectionMode_ == MULTI_SELECTION, "unknown selection mode");
        labelText = tr("The selected file contains %1 volumes. Please select one or multiple volumes to load:").arg(origins.size());
    }
    titleLabel_->setText(labelText);

    // clear table and search widgets
    table_->blockSignals(true);
    filterTextBox_->blockSignals(true);
    comboBoxFilterAttribute_->blockSignals(true);
    table_->clear();
    filterTextBox_->clear();
    comboBoxFilterAttribute_->clear();

    metaDataKeys_.clear();
    filteredOrigins_.clear();

    // determine key set of origin meta data
    std::set<std::string> keySet;
    for (size_t i=0; i<origins_.size(); i++) {
        std::vector<std::string> keys = origins_.at(i).getMetaDataContainer().getKeys();
        keySet.insert(keys.begin(), keys.end());
    }
    for (std::set<std::string>::iterator it = keySet.begin(); it != keySet.end(); ++it)
        metaDataKeys_.push_back(QString::fromStdString(*it));
    metaDataKeys_.sort();

    // add keys to attribute combobox
    comboBoxFilterAttribute_->addItem(tr("All attributes"));
    comboBoxFilterAttribute_->addItems(metaDataKeys_);
    comboBoxFilterAttribute_->addItem("URL");

    // initialize table
    QStringList headers(metaDataKeys_);
    headers.push_back("URL");
    table_->setColumnCount(headers.size());
    table_->setHorizontalHeaderLabels(headers);

    table_->blockSignals(false);
    filterTextBox_->blockSignals(false);
    comboBoxFilterAttribute_->blockSignals(false);

    updateTableRows();
    adjustColumnWidths();
}

VolumeListingDialog::SelectionMode VolumeListingDialog::getSelectionMode() const {
    return selectionMode_;
}

void VolumeListingDialog::showEvent(QShowEvent* event) {
    QDialog::showEvent(event);
    adjustColumnWidths();
}

void VolumeListingDialog::updateTableRows() {
    tgtAssert(table_->columnCount() == metaDataKeys_.size() + 1 /*+1 for URL*/, "column count mismatch");

    QString filterAttribute;
    if (comboBoxFilterAttribute_->currentIndex() > 0)
        filterAttribute = comboBoxFilterAttribute_->currentText();

    std::vector<VolumeURL> filteredOriginsNew = getFilteredOrigins(filterTextBox_->text().trimmed(), filterAttribute);
    if (filteredOrigins_ == filteredOriginsNew)
        return;
    else
        filteredOrigins_ = filteredOriginsNew;

    table_->clearContents();
    table_->setRowCount(static_cast<int>(filteredOrigins_.size()));

    // disable sorting during insertion
    table_->setSortingEnabled(false);

    // iterate over filter origins create a row for each
    for (size_t i=0; i<filteredOrigins_.size(); i++) {
        VolumeURL& origin = filteredOrigins_.at(i);
        // fill row cells with meta data entries
        for (int m=0; m<metaDataKeys_.size(); m++) {
            std::string key = metaDataKeys_.at(m).toStdString();
            std::string cellText;
            const MetaDataBase* metaDate = origin.getMetaDataContainer().getMetaData(key);
            if (metaDate)
                cellText = metaDate->toString();

            // set cell content type according to type of meta data item
            QTableWidgetItem* cellItem = new QTableWidgetItem();
            if (dynamic_cast<const IntMetaData*>(metaDate))
                cellItem->setData(Qt::DisplayRole, static_cast<const IntMetaData*>(metaDate)->getValue());
            else if (dynamic_cast<const SizeTMetaData*>(metaDate))
                cellItem->setData(Qt::DisplayRole, static_cast<qulonglong>(static_cast<const SizeTMetaData*>(metaDate)->getValue()));
            else if (dynamic_cast<const FloatMetaData*>(metaDate))
                cellItem->setData(Qt::DisplayRole, static_cast<const FloatMetaData*>(metaDate)->getValue());
            else if (dynamic_cast<const DoubleMetaData*>(metaDate))
                cellItem->setData(Qt::DisplayRole, static_cast<const DoubleMetaData*>(metaDate)->getValue());
            else
                cellItem->setData(Qt::DisplayRole, QString::fromStdString(cellText));

            cellItem->setToolTip(QString::fromStdString(key) + ": " + QString::fromStdString(cellText));
            //cellItem->setTextAlignment(Qt::AlignLeft);
            cellItem->setFlags(static_cast<Qt::ItemFlags>(Qt::ItemIsSelectable|Qt::ItemIsEnabled));
            table_->setItem(static_cast<int>(i), m, cellItem);
        }

        // put URL into last cell
        QTableWidgetItem* urlItem = new QTableWidgetItem(QString::fromStdString(origin.getURL()));
        urlItem->setToolTip("URL: " + QString::fromStdString(origin.getURL()));
        //urlItem->setTextAlignment(Qt::AlignLeft);
        urlItem->setFlags(static_cast<Qt::ItemFlags>(Qt::ItemIsSelectable|Qt::ItemIsEnabled));
        table_->setItem(static_cast<int>(i), table_->columnCount()-1, urlItem);
    }

    table_->setSortingEnabled(true);

    updateGuiState();
}

void VolumeListingDialog::adjustColumnWidths() {
    if (table_->columnCount() == 0)
        return;

    // let columns fill the table, but use minimum width
    int colWidth = std::max((table_->width()-20) / table_->columnCount(), 100);
    for (int i=0; i<table_->columnCount(); i++)
        table_->setColumnWidth(i, colWidth);
}

std::vector<VolumeURL> VolumeListingDialog::getFilteredOrigins(QString filterStr, QString filterAttribute /*= ""*/) const {
    if (filterStr == "")
        return origins_;

    filterStr = filterStr.toLower();
    std::vector<VolumeURL> result;

    // iterate over all origins and check for match with filter string
    for (size_t i=0; i<origins_.size(); i++) {
        const VolumeURL& origin = origins_.at(i);
        bool match = false;
        if (filterAttribute == "") {
            // no filter attribute specified => check all meta data entries + URL for filter string
            std::vector<std::string> keys = origin.getMetaDataContainer().getKeys();
            for (size_t k=0; k<keys.size() && !match; k++) {
                const MetaDataBase* metaDate =  origin.getMetaDataContainer().getMetaData(keys.at(k));
                if (metaDate)
                    match |= (QString::fromStdString(metaDate->toString()).toLower().contains(filterStr) == true);
            }
            match |= (QString::fromStdString(origin.getURL()).toLower().contains(filterStr) == true);
        }
        else if (filterAttribute == "URL") {
            // check only URL
            match |= (QString::fromStdString(origin.getURL()).toLower().contains(filterStr) == true);
        }
        else {
            // check meta data entry whose key matches the filterAttribute
            const MetaDataBase* metaDate = origin.getMetaDataContainer().getMetaData(filterAttribute.toStdString());
            if (metaDate)
                match |= (QString::fromStdString(metaDate->toString()).toLower().contains(filterStr) == true);
        }

        if (match)
            result.push_back(origin);
    }

    return result;
}

void VolumeListingDialog::cellDoubleClicked() {
    if (!table_->selectedItems().empty())
        loadClicked();
}

void VolumeListingDialog::updateGuiState() {
    loadButton_->setEnabled(!table_->selectedItems().empty());
    if (selectAllButton_)
        selectAllButton_->setEnabled(table_->rowCount() > 0);
}

void VolumeListingDialog::loadClicked() {
    std::vector<VolumeURL> selectedOrigins;

    // map selected rows to origins by comparing url column to origin url
    QList<QTableWidgetItem*> selectedItems = table_->selectedItems();
    int urlColumn = table_->columnCount()-1; // last columns contains URL
    for (int i=0; i<selectedItems.size(); i++) {
        QTableWidgetItem* item = selectedItems.at(i);
        if (item->column() == urlColumn) {
            std::string url = item->text().toStdString();
            for (size_t j=0; j<filteredOrigins_.size(); j++) {
                VolumeURL& origin = filteredOrigins_.at(j);
                if (origin.getURL() == url) {
                    selectedOrigins.push_back(origin);
                    break;
                }
            }
        }
    }

    if (!selectedOrigins.empty()) {
        close();
        emit(originsSelected(selectedOrigins, reader_));
    }
    else {
        LWARNINGC("voreen.qt.VolumeListingDialog", "no origins selected");
    }
}

void VolumeListingDialog::selectAllClicked() {
    table_->selectAll();
    table_->setFocus();
}

void VolumeListingDialog::cancelClicked() {
    close();
    emit(selectionCancelled());
}

} // namespace
