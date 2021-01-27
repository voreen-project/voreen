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

#ifndef VRN_VOLUMELISTINGDIALOG_H
#define VRN_VOLUMELISTINGDIALOG_H

#include "voreen/core/io/volumereader.h"
#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/qt/voreenqtapi.h"

#include <string>
#include <vector>

#include <QWidget>
#include <QString>
#include <QStringList>
#include <QShowEvent>
#include <QDialog>

class QTableWidget;
class QComboBox;
class QLineEdit;
class QPushButton;
class QLabel;

namespace voreen {

/**
 * Presents a list of VolumeURL items and allows the user to select one or multiple
 * for loading. The selected origins are emitted through a signal.
 */
class VRN_QT_API VolumeListingDialog : public QDialog {
    Q_OBJECT
public:

    /// Specifies whether one or multiple files can be selected
    enum SelectionMode {
        SINGLE_SELECTION,
        MULTI_SELECTION
    };

    VolumeListingDialog(QWidget* parent = 0, SelectionMode selectionMode = SINGLE_SELECTION);

    /**
     * Assigns the origins that should be presented to the user.
     * The \p reader is not used internally, but it is passed to the originsSelected() signal.
     */
    void setOrigins(const std::vector<VolumeURL>& origins, VolumeReader* reader = 0);

    /// Returns the dialog's SelectionMode.
    SelectionMode getSelectionMode() const;

signals:
    void originsSelected(const std::vector<VolumeURL>& origins, VolumeReader* reader);
    void selectionCancelled();

protected:
    void showEvent(QShowEvent* event);

private:
    /**
     * Returns all origins whose meta data entries contain the filter string.
     *
     * @param filterStr the string to search for
     * @param filterAttribute the key of the meta data entry that should be considered.
     *  If an empty string is passed, all meta data items are searched.
     */
    std::vector<VolumeURL> getFilteredOrigins(QString filterStr, QString filterAttribute) const;

    SelectionMode selectionMode_;

    std::vector<VolumeURL> origins_; ///< the origins to present to the user
    VolumeReader* reader_;              ///< only needed for passing it to originsSelected() signal

    /// meta data keys of the currently assigned origins, extracted by setOrigins()
    QStringList metaDataKeys_;
    /// currently visible origins, set by updateTableRows
    std::vector<VolumeURL> filteredOrigins_;

    QTableWidget* table_;

    QLabel* titleLabel_;
    QLabel* filterLabel_;
    QLineEdit* filterTextBox_;
    QComboBox* comboBoxFilterAttribute_;

    QPushButton* loadButton_;
    QPushButton* selectAllButton_;
    QPushButton* cancelButton_;

private slots:
    /// Updates the table content. Called after origin assignment and after each filter change.
    void updateTableRows();

    void adjustColumnWidths();
    void updateGuiState();

    void cellDoubleClicked();
    void loadClicked();
    void selectAllClicked();
    void cancelClicked();
};

} // namespace

#endif
