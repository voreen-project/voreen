/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2021 University of Muenster, Germany,                        *
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

#include "voreen/qt/widgets/property/filedialogpropertywidget.h"

#include "voreen/core/properties/filedialogproperty.h"

#include <QDir>
#include <QFileDialog>
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QPushButton>
#include <QCheckBox>
#include <QLabel>

#include "tgt/filesystem.h"

namespace voreen {

FileDialogPropertyWidget::FileDialogPropertyWidget(FileDialogProperty* prop, QWidget* parent)
    : QPropertyWidget(prop, parent)
    , property_(prop)
    , openFileDialogBtn_(new QPushButton(this))
    , fileWatchCheckBox_(nullptr)
{
    updateButtonText(prop->get());

    mainLayout_ = new QVBoxLayout();
    layout_->addLayout(mainLayout_);

    QHBoxLayout* buttonLayout = new QHBoxLayout();
    mainLayout_->addLayout(buttonLayout);
    buttonLayout->addWidget(openFileDialogBtn_);

    // Add button for file dialog.
    connect(openFileDialogBtn_, SIGNAL(clicked(void)), this, SLOT(setProperty(void)));
    QFontInfo fontInfo(font());
    openFileDialogBtn_->setFont(QFont(fontInfo.family(), QPropertyWidget::fontSize_));

    // Add file watch checkbox, if not safe file - dialog.
    if(prop->getFileMode() != FileDialogProperty::SAVE_FILE && prop->getWatchMode() != VoreenFileWatchListener::ALWAYS_OFF) {
        fileWatchCheckBox_ = new QCheckBox(tr("Watch File Change?"), this);
        fileWatchCheckBox_->setChecked(prop->isFileWatchEnabled());

        if(property_->isFileWatchEditable())
            connect(fileWatchCheckBox_, SIGNAL(clicked(bool)), this, SLOT(changeFileWatchOption(bool)));
        else
            fileWatchCheckBox_->setEnabled(false);
        mainLayout_->addWidget(fileWatchCheckBox_);
    }
}

FileDialogPropertyWidget::~FileDialogPropertyWidget() {
}

void FileDialogPropertyWidget::setProperty() {
    if (!disconnected_) {
        QString dialogCaption = QString::fromStdString(property_->getDialogCaption());
        QString directory;
        // use directory of current property value if any, default directory otherwise
        if (!property_->get().empty())
            directory = QString::fromStdString(property_->get());
        else
            directory = QString::fromStdString(property_->getDirectory());

        QString fileFilter;
        if (property_->getFileFilter().empty())
            fileFilter = tr("All files (*)");
        else
            fileFilter = QString::fromStdString(property_->getFileFilter()) + ";;" + tr("All files (*)");

        QString filename;
        if (property_->getFileMode() == FileDialogProperty::OPEN_FILE) {
            filename = QFileDialog::getOpenFileName(QWidget::parentWidget(), dialogCaption, directory, fileFilter, nullptr, QFileDialog::DontUseNativeDialog);
        }
        else if (property_->getFileMode() == FileDialogProperty::SAVE_FILE) {
            QString selectedFilter;
            filename = QFileDialog::getSaveFileName(QWidget::parentWidget(), dialogCaption, directory, fileFilter, &selectedFilter, QFileDialog::DontUseNativeDialog);

            // Create regular expression to parse the file filter string for the extension.
            // This is necessary, since the QFileDialog will not add the file extension by default.
            QRegExp filter_regex(QLatin1String("(?:^\\*\\.(?!.*\\()|\\(\\*\\.)(\\w+(.\\w+)*)"));

            // add the first extension of the selected file filter to the filename, if it is not present in the file name
            QFileInfo info(filename);
            if (info.suffix().isEmpty() && !selectedFilter.isEmpty()) {
                if (filter_regex.indexIn(selectedFilter) != -1) {
                    QString extension = filter_regex.cap(1);
                    filename += QLatin1String(".") + extension;
                }
            }
        }
        else if (property_->getFileMode() == FileDialogProperty::DIRECTORY) {
            filename = QFileDialog::getExistingDirectory(QWidget::parentWidget(), dialogCaption, QString::fromStdString(property_->get()), QFileDialog::DontUseNativeDialog);
        }

        if (!filename.isEmpty()) {
            property_->set(filename.toStdString());
            emit valueModifiedByUser();
        }
    }
}

void FileDialogPropertyWidget::changeFileWatchOption(bool checked) {
    property_->setFileWatchEnabled(checked);
}

void FileDialogPropertyWidget::updateButtonText(const std::string& filename) {

    if (!filename.empty()) {
        if ((property_->getFileMode() == FileDialogProperty::OPEN_FILE) || (property_->getFileMode() == FileDialogProperty::SAVE_FILE)) {
            size_t index = filename.find_last_of('/');
            if (index == filename.npos)
                index = filename.find_last_of('\\');
            std::string endFilename = filename;
            if (index != filename.npos)
                endFilename = filename.substr(index + 1, filename.length());
            openFileDialogBtn_->setText(QString::fromStdString(endFilename));
        }
        else if (property_->getFileMode() == FileDialogProperty::DIRECTORY) {
            std::string directory = filename;
            if (directory.length() >= 20)
                directory = "..." + directory.substr(directory.length()-20);
            openFileDialogBtn_->setText(QString::fromStdString(directory));
        }
    }
    else {
        if (property_->getFileMode() == FileDialogProperty::OPEN_FILE)
             openFileDialogBtn_->setText(tr("Select File"));
        else if (property_->getFileMode() == FileDialogProperty::SAVE_FILE)
             openFileDialogBtn_->setText(tr("Select File"));
        else if (property_->getFileMode() == FileDialogProperty::DIRECTORY)
            openFileDialogBtn_->setText(tr("Select Directory"));
    }

    openFileDialogBtn_->update();
}

void FileDialogPropertyWidget::updateFromPropertySlot() {
    updateButtonText(property_->get());

    if (fileWatchCheckBox_)
        fileWatchCheckBox_->setChecked(property_->isFileWatchEnabled());
}

} // namespace
