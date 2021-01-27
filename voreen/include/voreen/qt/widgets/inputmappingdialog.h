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

#ifndef VRN_INPUTMAPPINGDIALOG_H
#define VRN_INPUTMAPPINGDIALOG_H

#include "voreen/core/network/workspace.h"

#include "voreen/core/network/processornetworkobserver.h"
#include "voreen/qt/voreenqtapi.h"

#include <QWidget>
#include <QScrollArea>
#include <QMap>

class QVBoxLayout;
class QLabel;
class QGroupBox;
class QSpacerItem;
class QShowEvent;

namespace voreen {

class EventPropertyBase;
class EventPropertyWidget;
class ProcessorNetwork;
class Processor;

class VRN_QT_API InputMappingDialog : public QWidget, public ProcessorNetworkObserver, public UsesWorkspace {
Q_OBJECT
public:
    InputMappingDialog(QWidget* parent = 0);

public slots:
    /**
     * This method passes the Workspace whose network's event properties
     * are to be displayed by the InputMappingDialog.
     */
    void setWorkspace(Workspace* workspace);
    void rebuildWidgets();

    // Implementation of the ProcessorNetworkObserver interface
    virtual void networkChanged();
    virtual void processorAdded(const Processor* processor);
    virtual void processorRemoved(const Processor* processor);

protected:
    ProcessorNetwork* getProcessorNetwork();

    /// We need to create the widgets on show (lazy instantiation).
    void showEvent(QShowEvent* event);

    void createWidgets();

    void addProcessorToLayout(const Processor* processor);

    Workspace* workspace_;

    QMap<Processor*, QGroupBox*> processorBoxMap_;
    QScrollArea* scrollArea_;
    QVBoxLayout* scrollLayout_;
    QSpacerItem* scrollStretchItem_;

    bool widgetsValid_;     ///< Determines, whether widgets have to be rebuild on next show
};

} // namespace

#endif // VRN_INPUTMAPPINGDIALOG_H
