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

#ifndef VRN_VOLUMEREADERSELECTIONDIALOG_H
#define VRN_VOLUMEREADERSELECTIONDIALOG_H

#include "voreen/core/io/volumereader.h"
#include "voreen/qt/voreenqtapi.h"

#include <string>
#include <vector>

#include <QWidget>

class QTableWidget;
class QPushButton;
class QLabel;

namespace voreen {

/**
 * Presents a list of VolumeReaders and allows the user to select one of them.
 * The selected reader is passed through a signal.
 */
class VRN_QT_API VolumeReaderSelectionDialog : public QWidget {
    Q_OBJECT
public:
    VolumeReaderSelectionDialog(QWidget* parent = 0);

    void setURL(const std::string& url);
    void setReaders(const std::vector<VolumeReader*>& readers);

signals:
    void readerSelected(const std::string& url, VolumeReader* reader);
    void selectionCancelled();

private:
    std::string url_; ///< url of the file to load, shown in the dialog
    std::vector<VolumeReader*> readers_; ///< volume readers to present to the user

    QTableWidget* table_;
    QLabel* fileLabel_;
    QPushButton* loadButton_;
    QPushButton* cancelButton_;

private slots:
    void updateGuiState();
    void loadClicked();
    void cellDoubleClicked();
    void cancelClicked();
};

} // namespace

#endif
