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

#ifndef VRN_PROPERTYLISTWIDGET_H
#define VRN_PROPERTYLISTWIDGET_H

#include "voreen/core/network/processornetwork.h"
#include "voreen/core/network/workspace.h"
#include "voreen/core/properties/property.h"
#include "voreen/core/properties/optionproperty.h"

#include "voreen/qt/voreenqtapi.h"

#include <map>
#include <QScrollArea>
#include <QToolButton>
#include <QVBoxLayout>

#include <QStackedWidget>

class QLabel;

namespace voreen {

class ProcessorPropertiesWidget;
class AppModePropertyGroupWidget;

/**
 * This class displays a list of ProcessorPropertiesWidgets/AppModePropertyGroupWidget representing the processors' properties in a network.
 *
 * PropertyListWidget objects register themselves as observers at the assigned workspace's ProcessorNetwork.
 */
class VRN_QT_API PropertyListWidget : public QScrollArea, public ProcessorNetworkObserver, public WorkspaceObserver, public UsesWorkspace {
    Q_OBJECT

public:
    /**
     * Determines the widget's GUI mode.
     */
    enum WidgetMode {
        NETWORK,        ///< All properties of the selected network processors are shown in an expandable list.
        APPLICATION     ///< Property visibility and grouping is determined by the workspace's ApplicationModeConfiguration.
    };

    /**
     * Constructor - creates an empty layout
     *
     * @param parent parent widget
     * @param mode the widget mode
     */
    PropertyListWidget(QWidget* parent = 0, WidgetMode mode = NETWORK);
    ~PropertyListWidget();

    /**
     * This method passes the Workspace whose processors' properties
     * are to be displayed by the PropertyListWidget.
     */
    void setWorkspace(Workspace* workspace);

    /// Changes the WidgetMode.
    void setWidgetMode(WidgetMode mode);

    /// Returns the WidgetMode.
    WidgetMode getWidgetMode() const;

    /// Returns the property list's current level of detail.
    Property::LevelOfDetail getLevelOfDetail() const;

    /**
     * Convenience function for setting the WidgetMode and LevelOfDetail
     * settings in combination.
     */
    void setState(WidgetMode mode, Property::LevelOfDetail lod);

    /// Deletes all property-dependent widgets.
    void clear();

    // Implementation of the ProcessorNetworkObserver interface
    void networkChanged();
    void processorAdded(const Processor* processor);
    void processorRemoved(const Processor* processor);
    void processorRenamed(const Processor* p, const std::string& prevName);
    void connectionsChanged();
    void propertyLinkAdded(const PropertyLink* link);
    void propertyLinkRemoved(const PropertyLink* link);
    void portConnectionAdded(const Port* outport, const Port* inport);
    void portConnectionRemoved(const Port* outport, const Port* inport);

    // implementation of WorkspaceObserver interface
    void applicationModeConfigurationChanged();

    virtual QSize sizeHint() const;

public slots:
    /**
     * This method is called by the NetworkEditor when a processor is selected or unselected.
     * It will unfold the header widget of the selected processors.
     *
     * @param p processors that are selected
     */
    void processorsSelected(const QList<Processor*>& processors);

    /// Adjusts the level of detail of all property widgets.
    void setLevelOfDetail(Property::LevelOfDetail lod);
signals:
    /// Emitted when a processor (i.e. its properties) was modified.
    void modified();

    /**
     * Signal that is emitted when a property was been changed by the user. It is fired
     * when a repaint of the rendering is necessary.
     */
    void repaintSignal();

protected:
    void levelOfDetailOnChange();
    /// Generates the stack widget and adds the container widgets of the two widget modes.
    void createCommonWidgets();

    /// Generates the property widgets for the network mode (called lazy on first activation of the mode).
    void createNetworkModeWidgets();
    /// Deleted all network mode widgets.
    void clearNetworkModeWidgets();

    /// Generates the property widgets for the application mode (called lazy on first activation of the mode).
    void createAppModeWidgets();
    /// Deletes all application mode widgets.
    void clearAppModeWidgets();

    /// The assigned workspace whose network's processor properties are displayed.
    Workspace* workspace_;

    WidgetMode widgetMode_;
     //property storing the current lod. Will be serialized by the main window.
    OptionProperty<Property::LevelOfDetail> currentLevelOfDetail_;

    /// main widget (created on construction)
    QStackedWidget* stackWidget_;

    // network mode widgets/resources (created lazy)
    QWidget* networkModeContainer_;
    QVBoxLayout* networkModeContainerLayout_;
    std::map<const Processor*, ProcessorPropertiesWidget*> networkModeProcessorWidgetMap_; ///< mapping from processors to property widgets
    QList<Processor*> previouslySelectedProcessors_;

    // application mode widgets (created lazy)
    QWidget* appModeContainer_;
    QVBoxLayout* appModeContainerLayout_;
    QLabel* appModeStatusLabel_;
    std::vector<AppModePropertyGroupWidget*> appModeGroupWidgets_; /// one widget per group
    bool appModeWidgetsValid_;

};

} // namespace voreen

#endif //VRN_PROPERTYLISTWIDGET_H
