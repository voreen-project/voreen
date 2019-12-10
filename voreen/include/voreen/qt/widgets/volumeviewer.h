/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2019 University of Muenster, Germany,                        *
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

#ifndef VRN_VOLUMEVIEWER_H
#define VRN_VOLUMEVIEWER_H

#include "voreen/core/datastructures/volume/volumeram.h"
#include "voreen/core/ports/port.h"
#include "voreen/core/network/networkevaluator.h"

#include "voreen/qt/voreenqtapi.h"
#include "voreen/qt/widgets/volumeiohelper.h"

#include <QDialog>
#include <QLabel>
#include <QPushButton>
#include <QTreeWidget>
#include <QTreeWidgetItem>
#include <QWidget>
#include <QPushButton>
#include <QCheckBox>
#include <QLabel>
#include <QShowEvent>
#include <QPaintEvent>

class QPixmap;

namespace voreen {

class QRCTreeWidget : public QTreeWidget {
    Q_OBJECT
public:
    QRCTreeWidget(QWidget* parent = 0);

signals:
    void refresh(QTreeWidgetItem*);
    void exportVolumes();

protected:
    void keyPressEvent(QKeyEvent* event);
    void contextMenuEvent(QContextMenuEvent*);
};

// ---------------------------------------------------------------------------

/**
 * Widget that traverses the network and displays all contained volumes in a tree.
 */
class VRN_QT_API VolumeViewer : public QWidget, public VolumeObserver, public NetworkEvaluatorObserver {
    Q_OBJECT
public:

    VolumeViewer(QWidget* parent = 0);
    ~VolumeViewer();

    void setNetworkEvaluator(NetworkEvaluator* evaluator);

    /// @see VolumeObserver
    virtual void volumeDelete(const VolumeBase* source);

    /// @see VolumeObserver
    virtual void volumeChange(const VolumeBase* source);

    /// @see NetworkEvaluator::ProcessWrapper
    void afterNetworkProcess();

    /// @see NetworkEvaluator::ProcessWrapper
    void afterNetworkInitialize();

    /// @see NetworkEvaluator::ProcessWrapper
    void beforeNetworkDeinitialize();

    virtual QSize sizeHint() const;

    /// Clears the widget.
    void clear();

protected:
    /// triggers an update from the network, if the updateRequires_ flag is true
    virtual void paintEvent(QPaintEvent* event);

    virtual void showEvent(QShowEvent* event);

private:
    QTreeWidgetItem* createTreeWidgetItem(const VolumeBase* handle, const Port* port);
    std::string calculateVolumeSizeString(const std::vector<const VolumeBase*>& volumes);

    NetworkEvaluator* evaluator_;

    QPushButton* updateButton_;
    QCheckBox* autoUpdateCheckbox_;
    QLabel* containerInfo_;

    VolumeIOHelper volumeIOHelper_;

    /// QTreeWidget with custom context Menu
    QRCTreeWidget* volumeInfos_;

    /// stores the current displayed handles
    std::vector<const VolumeBase*> volumeHandles_;

    /// maps from a displayed handle to its hash (currently not in use)
    std::map<const VolumeBase*, std::string> handleToHashMap_;

    /// maps from a displayed handle to the outport it originates from
    std::map<const VolumeBase*, Port*> handleToPortMap_;

    // if true, the widget needs to be updated from the network
    bool updateRequired_;

    static const std::string loggerCat_;

private slots:
    /// updates the widget from the current network state
    void updateFromNetwork();

    /// exports currently selected volumes
    void exportVolumes();

    /// issued when the auto-update checkbox changes
    void updateStateChanged(int state);

    /// shows a message dialog to the user
    void volumeSaved(const VolumeBase* handle, const std::string& exportPath);

    void resizeOnCollapse(bool);

};

} // namespace

#endif
