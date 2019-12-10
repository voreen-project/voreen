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

#ifndef VRN_APPLICATIONMODECONFIGDIALOG_H
#define VRN_APPLICATIONMODECONFIGDIALOG_H

#include "voreen/core/network/workspace.h"
#include "voreen/core/network/processornetworkobserver.h"

#include "voreen/qt/mainwindow/voreenqtmainwindow.h"

#include <QWidget>
#include <QCheckBox>
#include <QComboBox>
#include <QTreeWidgetItem>

class QVBoxLayout;
class QLabel;
class QGroupBox;
class QSpacerItem;
class QShowEvent;

class QTreeWidget;
class QTreeWidgetItem;
class QListWidget;
class QPushButton;

namespace voreen {

class ProcessorNetwork;
class Processor;
class Property;
class ApplicationModeConfigProcessorItem;

/**
 * Dialog to configurate the property and menu entity visiblility in the applicatuon mode.
 */
class ApplicationModeConfigDialog : public QWidget, public ProcessorNetworkObserver, public UsesWorkspace {
    Q_OBJECT
public:
    /**
     * Constructor
     *
     * @param parent normal parent
     * @param workspace current workspace
     */
    ApplicationModeConfigDialog(QWidget* parent = 0, Workspace* workspace = 0);

    /** Sets the main window. */
    void setMainWindow(VoreenQtMainWindow* window);
private:

    //-----------------------
    //    layout functions
    //-----------------------
private:
    /** Creates the group control elements. */
    QWidget* layoutGroups();
    /** Creates the menu entities control elements. */
    QWidget* layoutMenuEntities();
    /** Creates the combo box for main canvas selection. */
    QWidget* layoutMainCanvas();
    /** Creates the tree widget. */
    QWidget* layoutProperties();

    //-----------------------
    //  initialize functions
    //-----------------------
private:
    /**
     * Re-initializes all widgets.
     * Should be called after a new workspace has been loaded.
     */
    void reinitAllWidgets();
    /** Updates the group widget to the current configuration */
    void initGroupListWidget();
    /** Updates the entity widget to the current configuration */
    void initEntityTreeWidget();
    /** Updates the main canvas widget to the current configuration */
    void initMainCanvasWidget();
    /** Updates the property widget to the current configuration */
    void initPropertyTreeWidget();

    //-----------------------
    //    observer functions
    //-----------------------
public slots:
    /** @override ProcessorNetworkObserver */
    virtual void networkChanged();
    /** @override ProcessorNetworkObserver */
    virtual void processorAdded(const Processor* processor);
    /** @override ProcessorNetworkObserver */
    virtual void processorRemoved(const Processor* processor);
    /** @override ProcessorNetworkObserver */
    virtual void processorRenamed(const Processor* processor, const std::string& prevName);
    /** @override UsesWorkspace */
    virtual void setWorkspace(Workspace* workspace);

    //-----------------------
    //    callback functions
    //-----------------------
private slots:
    void handleAddGroup();
    void handleRemoveGroup();
    void handleGroupRename();
    void handleGroupMoveUp();
    void handleGroupMoveDown();

    void handlePropertyPriorityChange(QTreeWidgetItem* propertyItem, Property* property, std::string groupName, int priority);
    void handleEntityStateChange(std::string entityName, bool visible);
    void handleMainCanvasChange(const QString & text);

    //-----------------------
    //    helper functions
    //-----------------------
private slots:
    /**
     * Updates all widget states, e.i., enable/disable
     */
    void updateWidgetState();
    void expandProcessor(QTreeWidgetItem*);
protected:
    /**
     * Initializes widgets on show (lazy instantiation).
     * @see QWidget::showEvent(QShowEvent* event)
     */
    virtual void showEvent(QShowEvent* event);

    //-----------------------
    //    TWidget updates
    //-----------------------
private:
    void updateTreeHeaders();
    void addProcessorToTreeWidget(const Processor* processor);
    void addProcessorPropertiesToWidget(ApplicationModeConfigProcessorItem* processorItem);
    void removeProcessorFromTreeWidget(const Processor* processor);
    void renameProcessorInTreeWidget(const std::string& prevName, const std::string& newName);

    //-----------------------
    //    member
    //-----------------------
private:
    Workspace* workspace_;
    VoreenQtMainWindow* mainWindow_;

    QListWidget* groupListWidget_;
    QPushButton* addGroupButton_;
    QPushButton* deleteGroupButton_;
    QPushButton* renameGroupButton_;
    QPushButton* moveGroupUpButton_;
    QPushButton* moveGroupDownButton_;

    QTreeWidget* propertyTreeWidget_;
    QTreeWidget* entityTreeWidget_;
    QComboBox* mainCanvasRendererCB_;

    bool widgetsValid_;     ///< determines, whether widgets have to be rebuild on next show

    static const std::string loggerCat_;
};



/**
 * CheckBox of MenuEntity Visibility.
 * Used to get SINGNALS with correct parameters.
 */
class MenuEntityVisibilityBox : public QCheckBox {
    Q_OBJECT
public:
    /**
     * Constructor
     */
    MenuEntityVisibilityBox(std::string name);
private slots:
    void handleStateChanged(int);
signals:
    void visibilityChanged(std::string, bool);
private:
    std::string menuEntityName_;
};

class ApplicationModeConfigProcessorItem: public QTreeWidgetItem{
public:
    ApplicationModeConfigProcessorItem(const Processor* processor);
    const Processor* getProcessor() const;
    void setInitialized(bool b);
    bool getInitialized() const;
private:
    const Processor* processor_;
    bool initialized_;
};

/**
 * ComboBox of fixed size containing numbers for priority selection.
 * This class is used in the ApplicationModeConfigDialog
 */
class PrioritySelectionBox : public QComboBox {
    Q_OBJECT
public:
    /**
     * Constructor
     *
     * @param treeWidgetItem
     * @param property
     * @param groupName
     * @param numProperties
     */
    PrioritySelectionBox(QTreeWidgetItem* treeWidgetItem, Property* property, const std::string& groupName, int numPriorities);

    void setPriority(int priority);

private:
    int priorityCount_;               ///< number of stored priorities
    Property* property_;              ///< associated property
    std::string groupName_;           ///< associated group name
    QTreeWidgetItem* treeWidgetItem_; ///< associated tree widget item (i.e. property row)

private slots:
    /**
     * Slot called on selection change.
     * Emits priorityChanged signal.
     * @see priorityChanged(...)
     */
    void handleSelectionChange();

signals:
    /**
     * Emited by selection change.
     * Is connected with die ApplicationModeConfigDialog.
     *
     * @param treeWidgetItem
     * @param property
     * @param groupName
     * @param priority
     */
    void priorityChanged(QTreeWidgetItem* treeWidgetItem, Property* property, std::string groupName, int priority);
};

} // namespace

#endif // VRN_APPLICATIONMODECONFIGDIALOG_H
