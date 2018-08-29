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

#ifndef VRN_VOLUMEIOHELPER_H
#define VRN_VOLUMEIOHELPER_H

#include "voreen/qt/voreenqtapi.h"

#include "voreen/core/datastructures/volume/volume.h"

#include <vector>
#include <QToolButton>

namespace voreen {

class VolumeBase;
class VolumeReader;
class VolumeWriter;
class ProgressBar;
class VolumeSerializerPopulator;
class VolumeReaderSelectionDialog;
class VolumeListingDialog;
class VolumeList;

/**
 * Helper class for loading and saving volumes.
 * All loaded volumes are emitted via a signal.
 */
class VRN_QT_API VolumeIOHelper : public QObject {
    Q_OBJECT
public:

    /// Specifies whether one or multiple files can be loaded during a single load operation
    enum FileMode {
        SINGLE_FILE,   ///< only a single file can be loaded
        MULTI_FILE     ///< multiple files can be loaded
    };

    VolumeIOHelper(QWidget* parent = 0, FileMode fileMode = SINGLE_FILE);
    ~VolumeIOHelper();

signals:
    /**
     * Is emitted after a volume has been loaded.
     *
     * @param handle the loaded volume
     */
    void volumeLoaded(const VolumeBase* handle);

    /**
     * Is emitted after a volume has been written to disk.
     *
     * @param handle the save volume
     * @param exportPath the file path the volume has been written to
     */
    void volumeSaved(const VolumeBase* handle, const std::string& exportPath);

public slots:
    /**
     * Presents a file open dialog to the user
     * with filters for all supported volume formats.
     * Additionally, filters for loading raw volumes
     * and raw slice stacks are added.
     * The user selected files are loaded then
     * using the resp. volume readers.
     */
    void showFileOpenDialog();

    /**
     * Presents a file save dialog to the user
     * with filters for all supported volume formats.
     * The passed volume is saved using the user selected
     * volume writer under the user specified path.
     */
    void showFileSaveDialog(const VolumeBase* volume);

    /**
     * Loads volume(s) from the passed URL using the passed reader.
     * If the URL references a container file with more than one volume,
     * a VolumeListingDialog is shown for selecting the volume(s) to load.
     * For each loaded volume, a volumeLoaded() signal is emitted.
     */
    void loadURL(const std::string& url, VolumeReader* reader);

    /**
     * Loads the volume from the passed origin using the passed reader.
     * The loaded volume is emitted via a volumeLoaded() signal.
     */
    void loadOrigin(const VolumeURL& origin, VolumeReader* reader);

    /**
     * Helper slot: passes all origins and the reader to loadOrigin().
     */
    void loadOrigins(const std::vector<VolumeURL>& origins, VolumeReader* reader);

    /**
     * Loads the raw volume with the passed filename.
     * The required parameters are queried from the user
     * via a RawVolumeWidget. The loaded volume is emitted
     * via a volumeLoaded() signal.
     */
    void loadRawVolume(const std::string& filename);

    /**
     * Saves the passed volume using the passed writer to the passed path.
     */
    void saveVolumeToPath(const VolumeBase* volume, VolumeWriter* writer, const std::string& filepath);

    /// Returns the helper's FileMode.
    FileMode getFileMode() const;

private:
    /// Returns the filter string used in the file open dialog.
    std::string getVolumeReaderFilterString() const;

    QWidget* parent_;
    FileMode fileMode_;

    VolumeSerializerPopulator* volumeSerializerPopulator_;
    ProgressBar* progressBar_;

    VolumeReaderSelectionDialog* readerSelectionDialog_;
    VolumeListingDialog* volumeListingDialog_;

    static const std::string loggerCat_;
};

} // namespace

#endif // VRN_VOLUMEIOHELPER_H
