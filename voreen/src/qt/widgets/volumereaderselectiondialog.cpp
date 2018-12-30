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

#include "voreen/qt/widgets/volumereaderselectiondialog.h"

#include "voreen/core/utils/stringutils.h"

#include <QString>
#include <QStringList>
#include <QLabel>
#include <QPushButton>
#include <QTableWidget>
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QFontMetrics>

namespace {
bool compareReaderClassName(voreen::VolumeReader* i, voreen::VolumeReader* j) {
    return i->getClassName() < j->getClassName();
}
}

namespace voreen {

VolumeReaderSelectionDialog::VolumeReaderSelectionDialog(QWidget* parent)
    : QWidget(parent)
{

    setWindowFlags(Qt::Dialog);
    setWindowModality(Qt::WindowModal);
    setWindowTitle("Volume Reader Selection");

    // table configuration
    table_ = new QTableWidget(1, 3);
    QString headers = tr("Volume Reader,Extensions,Format Description");
    table_->setHorizontalHeaderLabels(headers.split(","));
    table_->setColumnWidth(0, 150);
    table_->setColumnWidth(1, 100);
    table_->setColumnWidth(2, 260);
    table_->setSelectionBehavior(QAbstractItemView::SelectRows);
    table_->setSelectionMode(QAbstractItemView::SingleSelection);
    table_->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);

    // labels
    QVBoxLayout* mainLayout = new QVBoxLayout();
    QLabel* topLabel = new QLabel(tr("Please select the volume reader that should be used for loading the file:"));
    QFont labelFont = topLabel->font();
    labelFont.setPointSize(topLabel->font().pointSize()+1);
    labelFont.setBold(true);
    topLabel->setFont(labelFont);
    mainLayout->addWidget(topLabel);
    fileLabel_ = new QLabel();
    fileLabel_->setFixedWidth(530);
    mainLayout->addWidget(fileLabel_);
    mainLayout->addWidget(table_);

    // buttons
    QHBoxLayout* buttonLayout = new QHBoxLayout();
    buttonLayout->addStretch();
    cancelButton_ = new QPushButton(tr("Cancel"));
    cancelButton_->setFixedWidth(100);
    buttonLayout->addWidget(cancelButton_);
    loadButton_ = new QPushButton(tr("Load"));
    loadButton_->setFixedWidth(100);
    QFont loadFont(loadButton_->font());
    loadFont.setBold(true);
    loadButton_->setFont(loadFont);
    buttonLayout->addWidget(loadButton_);
    mainLayout->addLayout(buttonLayout);

    setLayout(mainLayout);
    setFixedSize(550, 250);

    connect(table_, SIGNAL(itemSelectionChanged()), this, SLOT(updateGuiState()));
    connect(table_, SIGNAL(cellDoubleClicked(int, int)), this, SLOT(cellDoubleClicked()));
    connect(loadButton_, SIGNAL(clicked()), this, SLOT(loadClicked()));
    connect(cancelButton_, SIGNAL(clicked()), this, SLOT(cancelClicked()));

    updateGuiState();
}

void VolumeReaderSelectionDialog::setURL(const std::string& url) {
    url_ = url;

    // if the string overflows the label, add points in the middle (ElideMiddle) and set tooltip
    QString str = QString::fromStdString(url_);
    QFontMetrics fm(fileLabel_->font());
    if (fileLabel_->width() < fm.width(str) ){
        fileLabel_->setText(fileLabel_->fontMetrics().elidedText(str, Qt::ElideMiddle, fileLabel_->width()));
        fileLabel_->setToolTip(str);
    }
    else {
        fileLabel_->setText(str);
    }

    updateGuiState();
}

void VolumeReaderSelectionDialog::setReaders(const std::vector<VolumeReader*>& readers) {
    readers_ = readers;
    std::sort(readers_.begin(), readers_.end(), compareReaderClassName);

    table_->clearContents();
    table_->setRowCount(static_cast<int>(readers_.size()));

    // add one row per reader
    for (size_t i=0; i<readers_.size(); i++) {
        const VolumeReader* reader = readers_.at(i);

        // name cell
        QTableWidgetItem* nameItem =
            new QTableWidgetItem(QString::fromStdString(reader->getClassName()));
        nameItem->setFlags(static_cast<Qt::ItemFlags>(Qt::ItemIsSelectable|Qt::ItemIsEnabled));
        table_->setItem(static_cast<int>(i), 0, nameItem);

        // extension cell
        std::string extensionStr = strJoin(reader->getSupportedExtensions(), " .");
        if (!extensionStr.empty())
            extensionStr = "." + extensionStr;
        std::string filenameStr = strJoin(reader->getSupportedFilenames(), " ");
        if (!filenameStr.empty())
            filenameStr = " " + filenameStr;
        extensionStr += filenameStr;
        QTableWidgetItem* extensionItem = new QTableWidgetItem(QString::fromStdString(extensionStr));
        extensionItem->setToolTip(tr("Supported extensions: ") + QString::fromStdString(extensionStr));
        extensionItem->setTextAlignment(Qt::AlignCenter);
        extensionItem->setFlags(static_cast<Qt::ItemFlags>(Qt::ItemIsSelectable|Qt::ItemIsEnabled));
        table_->setItem(static_cast<int>(i), 1, extensionItem);

        // description cell
        QTableWidgetItem* descItem =
            new QTableWidgetItem(QString::fromStdString(reader->getFormatDescription()));
        descItem->setFlags(static_cast<Qt::ItemFlags>(Qt::ItemIsSelectable|Qt::ItemIsEnabled));
        table_->setItem(static_cast<int>(i), 2, descItem);
    }

    updateGuiState();
}

void VolumeReaderSelectionDialog::updateGuiState() {
    loadButton_->setEnabled(!table_->selectedItems().empty());
}

void VolumeReaderSelectionDialog::loadClicked() {
    int row = table_->currentRow();
    if (row < 0 || row >= static_cast<int>(readers_.size())) {
        LWARNINGC("voreen.qt.VolumeReaderSelectionDialog", "invalid row: " << row);
        return;
    }
    close();
    emit(readerSelected(url_, readers_.at(row)));
}

void VolumeReaderSelectionDialog::cellDoubleClicked() {
    if (table_->currentRow() >= 0)
        loadClicked();
}

void VolumeReaderSelectionDialog::cancelClicked() {
    close();
    emit(selectionCancelled());
}

} // namespace
