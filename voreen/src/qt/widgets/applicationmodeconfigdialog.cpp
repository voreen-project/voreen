/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2020 University of Muenster, Germany,                        *
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

#include "voreen/qt/widgets/applicationmodeconfigdialog.h"

#include "voreen/qt/mainwindow/menuentities/voreenqtmenuentity.h"

#include "voreen/core/voreenapplication.h"
#include "voreen/core/processors/processor.h"
#include "voreen/core/network/processornetwork.h"
#include "modules/core/processors/output/canvasrenderer.h"
#include "voreen/core/network/workspace.h"

#include <QApplication>
#include <QPushButton>
#include <QLineEdit>
#include <QTreeWidgetItem>
#include <QLabel>
#include <QVBoxLayout>
#include <QGroupBox>
#include <QScrollArea>
#include <QShowEvent>
#include <QTreeWidget>
#include <QHeaderView>
#include <QListWidget>
#include <QInputDialog>
#include <QMessageBox>
#include <QCheckBox>

namespace voreen {
/**
 * Magic number of weight count
 */
const int WEIGHT_COUNT = 50;

const std::string ApplicationModeConfigDialog::loggerCat_("voreen.qt.ApplicationModeConfigDialog");

ApplicationModeConfigDialog::ApplicationModeConfigDialog(QWidget* parent, Workspace* workspace)
    : QWidget(parent)
    , workspace_(workspace)
    , mainWindow_(0)
    , groupListWidget_(0)
    , addGroupButton_(0)
    , deleteGroupButton_(0)
    , renameGroupButton_(0)
    , moveGroupUpButton_(0)
    , moveGroupDownButton_(0)
    , propertyTreeWidget_(0)
    , entityTreeWidget_(0)
    , mainCanvasRendererCB_(0)
    , widgetsValid_(false)
{
    //main layout
    QVBoxLayout* mainLayout = new QVBoxLayout();
    setLayout(mainLayout);
    // create components
    QWidget* groups = layoutGroups();
    QWidget* entities = layoutMenuEntities();
    QWidget* canvas = layoutMainCanvas();
    QWidget* properties = layoutProperties();
    //layout components
    QVBoxLayout* tmp1Layout = new QVBoxLayout();
    QHBoxLayout* tmp2Layout = new QHBoxLayout();
    tmp1Layout->addWidget(entities);
    tmp1Layout->addWidget(canvas);
    tmp2Layout->addWidget(groups);
    tmp2Layout->addLayout(tmp1Layout);
    tmp2Layout->addStretch();
    mainLayout->addLayout(tmp2Layout);
    mainLayout->addWidget(properties);
    //resize
    resize(625, 600);
    setMinimumSize(625, 500);
    // add this as observer
    if (workspace_ && workspace_->getProcessorNetwork())
        workspace_->getProcessorNetwork()->addObserver(this);
}

void ApplicationModeConfigDialog::setMainWindow(VoreenQtMainWindow* window) {
    mainWindow_ = window;
}

    //-----------------------
    //    layout functions
    //-----------------------
QWidget* ApplicationModeConfigDialog::layoutGroups() {
    // group box
    QGroupBox* groupsGroupBox = new QGroupBox("Property Groups");
    QHBoxLayout* groupLayout = new QHBoxLayout();
    groupsGroupBox->setLayout(groupLayout);
    groupsGroupBox->setFixedSize(300, 230);


    groupListWidget_ = new QListWidget();
    connect(groupListWidget_, SIGNAL(itemSelectionChanged()), this, SLOT(updateWidgetState()));
    groupLayout->addWidget(groupListWidget_);

    QVBoxLayout* buttonLayout = new QVBoxLayout();
    groupLayout->addLayout(buttonLayout);
    addGroupButton_ = new QPushButton(" Add Group ");
    buttonLayout->addWidget(addGroupButton_);
    connect(addGroupButton_, SIGNAL(clicked()), this, SLOT(handleAddGroup()));
    deleteGroupButton_ = new QPushButton(" Remove Group ");
    buttonLayout->addWidget(deleteGroupButton_);
    connect(deleteGroupButton_, SIGNAL(clicked()), this, SLOT(handleRemoveGroup()));

    QFrame* separator = new QFrame();
    separator->setFrameShape(QFrame::HLine);
    buttonLayout->addWidget(separator);

    renameGroupButton_ = new QPushButton(" Rename Group ");
    buttonLayout->addWidget(renameGroupButton_);
    connect(renameGroupButton_, SIGNAL(clicked()), this, SLOT(handleGroupRename()));
    moveGroupUpButton_ = new QPushButton(" Move Up ");
    buttonLayout->addWidget(moveGroupUpButton_);
    connect(moveGroupUpButton_, SIGNAL(clicked()), this, SLOT(handleGroupMoveUp()));
    moveGroupDownButton_ = new QPushButton(" Move Down ");
    buttonLayout->addWidget(moveGroupDownButton_);
    connect(moveGroupDownButton_, SIGNAL(clicked()), this, SLOT(handleGroupMoveDown()));

    buttonLayout->addStretch();

    groupLayout->addStretch();
    return groupsGroupBox;
}

QWidget* ApplicationModeConfigDialog::layoutMenuEntities() {
    // group box
    QGroupBox* entityGroupBox = new QGroupBox("Menu Entities");
    QHBoxLayout* groupLayout = new QHBoxLayout();
    entityGroupBox->setLayout(groupLayout);
    entityGroupBox->setFixedSize(300, 165);

    // internal tree widget
    entityTreeWidget_ = new QTreeWidget();
        //configure header
    entityTreeWidget_->header()->setSectionsMovable(false);
    QStringList headers;
    headers.push_back("Name");
    headers.push_back("Value");
    entityTreeWidget_->setHeaderLabels(headers);
    entityTreeWidget_->setAlternatingRowColors(true);
    entityTreeWidget_->setColumnWidth(0, 210);
    entityTreeWidget_->setColumnWidth(1, 40);
    groupLayout->addWidget(entityTreeWidget_);

    return entityGroupBox;
}

QWidget* ApplicationModeConfigDialog::layoutMainCanvas() {
    // group box
    QGroupBox* mainCanvasGroupBox = new QGroupBox("Main Canvas");
    QHBoxLayout* groupLayout = new QHBoxLayout();
    mainCanvasGroupBox->setLayout(groupLayout);
    mainCanvasGroupBox->setFixedSize(300, 55);

    // internal widget
    mainCanvasRendererCB_ = new QComboBox(this);
    groupLayout->addWidget(mainCanvasRendererCB_);
    connect(mainCanvasRendererCB_,SIGNAL(currentIndexChanged(QString)),this,SLOT(handleMainCanvasChange(QString)));
    return mainCanvasGroupBox;
}

QWidget* ApplicationModeConfigDialog::layoutProperties() {
    propertyTreeWidget_ = new QTreeWidget();
    propertyTreeWidget_->header()->setSectionsMovable(false);
    connect(propertyTreeWidget_, SIGNAL(itemExpanded(QTreeWidgetItem*)), this, SLOT(expandProcessor(QTreeWidgetItem*)));


    return propertyTreeWidget_;
}

    //-----------------------
    //  initialize functions
    //-----------------------
void ApplicationModeConfigDialog::reinitAllWidgets() {
    initGroupListWidget();
    initEntityTreeWidget();
    initMainCanvasWidget();
    initPropertyTreeWidget();

    widgetsValid_ = true;
}

void ApplicationModeConfigDialog::initGroupListWidget() {
    tgtAssert(groupListWidget_, "no group list widget");

    groupListWidget_->clear();

    if (!workspace_)
        return;
    const ApplicationModeConfiguration& config = workspace_->getApplicationModeConfig();

    std::vector<std::string> propertyGroups = config.getPropertyGroups();
    for (size_t i=0; i<propertyGroups.size(); i++) {
        groupListWidget_->addItem(QString::fromStdString(propertyGroups.at(i)));
    }
}

void ApplicationModeConfigDialog::initEntityTreeWidget() {
    tgtAssert(entityTreeWidget_, "no entity tree widget");

    entityTreeWidget_->clear();

    if (!workspace_ || !workspace_->getProcessorNetwork())
        return;

    //get all menu entities
    tgtAssert(mainWindow_, "no main window!");
    std::vector<VoreenQtMenuEntity*> entities = mainWindow_->getAllMenuEntities();

    //set entities checked according to the current config
    const ApplicationModeConfiguration& config = workspace_->getApplicationModeConfig();

    //add and update
    for(size_t i = 0; i < entities.size(); i++) {
        QStringList entityName;
        entityName.push_back(QString::fromStdString(entities[i]->getName()));
        QTreeWidgetItem* entityItem = new QTreeWidgetItem(entityName);
        entityTreeWidget_->addTopLevelItem(entityItem);
        entityItem->setIcon(0, entities[i]->getIcon());
        MenuEntityVisibilityBox* cBox = new MenuEntityVisibilityBox(entities[i]->getName());
            //update state
        std::pair<bool, bool> entityConfig = config.isMenuEntityVisible(entities[i]->getName());
        if(entityConfig.first) // if entity is present
            cBox->setChecked(entityConfig.second);
        else //entity not present, use default setting
            cBox->setChecked(entities[i]->getDefaultVisibilityInApplicationMode());
        entityTreeWidget_->setItemWidget(entityItem,1,cBox);
            //make connection
        connect(cBox,SIGNAL(visibilityChanged(std::string,bool)),this,SLOT(handleEntityStateChange(std::string,bool)));
    }

}

void ApplicationModeConfigDialog::initMainCanvasWidget() {
    tgtAssert(mainCanvasRendererCB_, "no canvas cb widget");

    mainCanvasRendererCB_->clear();

    if (!workspace_ || !workspace_->getProcessorNetwork())
        return;

    const ApplicationModeConfiguration& config = workspace_->getApplicationModeConfig();
    ProcessorNetwork* network = workspace_->getProcessorNetwork();
    std::vector<CanvasRenderer*> renderers(network->getProcessorsByType<CanvasRenderer>());

    //save old selection
    std::string selectedCanvas = config.getMainCanvas();

    for(size_t i = 0; i < renderers.size(); i++) {
        mainCanvasRendererCB_->addItem((renderers[i]->getClassName() + "::" + renderers[i]->getGuiName()).c_str());
    }

    int index = mainCanvasRendererCB_->findText(selectedCanvas.c_str());
    if(index >= 0) {
        mainCanvasRendererCB_->setCurrentIndex(index);
    } else {
        mainCanvasRendererCB_->setCurrentIndex(std::min(0, mainCanvasRendererCB_->count() -1));
    }
}

void ApplicationModeConfigDialog::initPropertyTreeWidget() {
    tgtAssert(propertyTreeWidget_, "no tree widget");

    propertyTreeWidget_->clear();

    if (!workspace_ || !workspace_->getProcessorNetwork())
        return;

    const ApplicationModeConfiguration& config = workspace_->getApplicationModeConfig();
    ProcessorNetwork* network = workspace_->getProcessorNetwork();
    std::vector<Processor*> processors(network->getProcessors());
    const std::vector<std::string> propertyGroups = config.getPropertyGroups();

    // order processors by name
    for (size_t i=0; i<processors.size(); i++) {
        for (size_t j=i+1; j<processors.size(); j++) {
            if (processors.at(i)->getGuiName() > processors.at(j)->getGuiName())
                std::swap(processors.at(i), processors.at(j));
        }
    }

    // adjust column count and headers to property groups
    propertyTreeWidget_->setColumnCount((int)propertyGroups.size() + 1);
    updateTreeHeaders();

    QApplication::setOverrideCursor(Qt::WaitCursor);
    qApp->processEvents();

    // add processors/properties
    for (size_t i=0; i<processors.size(); i++) {
        addProcessorToTreeWidget(processors.at(i));
    }

    propertyTreeWidget_->setColumnWidth(0, 200);
    for (size_t i=1; i<static_cast<size_t>(propertyTreeWidget_->columnCount()); i++)
        propertyTreeWidget_->setColumnWidth((int)i, 75);

    QApplication::restoreOverrideCursor();
}

    //-----------------------
    //    observer functions
    //-----------------------
void ApplicationModeConfigDialog::networkChanged() {
    // nothing
}

void ApplicationModeConfigDialog::processorAdded(const Processor* processor) {
    tgtAssert(processor, "null pointer passed");
    //update config
    if(workspace_) {
        const CanvasRenderer* cproc = 0;
        if((cproc = dynamic_cast<const CanvasRenderer*>(processor)) &&
                                                    workspace_->getApplicationModeConfig().getMainCanvas().empty()) {
           workspace_->getApplicationModeConfig().setMainCanvas(cproc->getClassName() + "::" + cproc->getGuiName());
        }
    }

    //update gui
    if (isVisible()) {
        if(const CanvasRenderer* cproc = dynamic_cast<const CanvasRenderer*>(processor)) {
            mainCanvasRendererCB_->addItem((cproc->getClassName() + "::" + cproc->getGuiName()).c_str());
        }

        addProcessorToTreeWidget(processor);
        updateWidgetState();
    }
    else
        widgetsValid_ = false;
}

void ApplicationModeConfigDialog::processorRemoved(const Processor* processor) {
    tgtAssert(processor, "null pointer passed");

    // remove config
    if (workspace_) {
        workspace_->getApplicationModeConfig().removeProcessorProperties(processor);
        const CanvasRenderer* cproc = 0;
        if((cproc = dynamic_cast<const CanvasRenderer*>(processor)) &&
                        !workspace_->getApplicationModeConfig().getMainCanvas().compare(cproc->getClassName() + "::" + cproc->getGuiName())) {
            workspace_->getApplicationModeConfig().setMainCanvas("");
        }
    }

    //remove gui
    if (isVisible()) {
        if(const CanvasRenderer* cproc = dynamic_cast<const CanvasRenderer*>(processor)) {
            mainCanvasRendererCB_->removeItem(mainCanvasRendererCB_->findText((cproc->getClassName() + "::" + cproc->getGuiName()).c_str()));
            workspace_->getApplicationModeConfig().setMainCanvas(mainCanvasRendererCB_->currentText().toStdString());
        }
        removeProcessorFromTreeWidget(processor);
        updateWidgetState();
    }
    else
        widgetsValid_ = false;
}

void ApplicationModeConfigDialog::processorRenamed(const Processor* processor, const std::string& prevName) {
    tgtAssert(processor, "null pointer passed");

    // rename config
    if (workspace_) {
        const CanvasRenderer* cproc = 0;
        if((cproc = dynamic_cast<const CanvasRenderer*>(processor)) &&
                                !workspace_->getApplicationModeConfig().getMainCanvas().compare(cproc->getClassName() + "::" + prevName)) {
           workspace_->getApplicationModeConfig().setMainCanvas(cproc->getClassName() + "::" + cproc->getGuiName());
        }
    }

    //rename gui
    if (isVisible()) {
        if(const CanvasRenderer* cproc = dynamic_cast<const CanvasRenderer*>(processor)) {
            mainCanvasRendererCB_->setItemText(mainCanvasRendererCB_->findText((cproc->getClassName() + "::" + prevName).c_str()),(cproc->getClassName() + "::" + cproc->getGuiName()).c_str());
        }
        renameProcessorInTreeWidget(prevName, processor->getGuiName());
        updateWidgetState();
    }
    else
        widgetsValid_ = false;
}

void ApplicationModeConfigDialog::setWorkspace(Workspace* workspace) {
    // stop observation of previously assigned workspace/network
    if (workspace_)
        stopObservation(workspace_->getProcessorNetwork());

    workspace_ = workspace;

    if (workspace_) {
        tgtAssert(workspace_->getProcessorNetwork(), "workspace has no processor network");
        workspace_->getProcessorNetwork()->addObserver(this);
    }

    updateWidgetState();

    // only rebuild widgets immediately, if widget is visible or config is empty
    if (isVisible() || (workspace_ && workspace_->getApplicationModeConfig().getPropertyGroups().empty()))
        reinitAllWidgets();
    else
        widgetsValid_ = false;
}

    //-----------------------
    //    callback functions
    //-----------------------
void ApplicationModeConfigDialog::handleAddGroup() {
    tgtAssert(groupListWidget_, "no group list widget");
    if (!workspace_)
        return;
    ApplicationModeConfiguration& config = workspace_->getApplicationModeConfig();
    tgtAssert(groupListWidget_->count() == static_cast<int>(config.getPropertyGroups().size()), "group count mismatch");

    // query group name from user
    QString newGroup;
    bool inputOk = false;
    do {
        newGroup = QInputDialog::getText(this, "Add Property Group", "Please enter the name of the new property group: ", QLineEdit::Normal, newGroup, &inputOk);
    }
    while (inputOk && config.hasPropertyGroup(newGroup.toStdString()));
    if (!inputOk)
        return;

    // add group
    config.addPropertyGroup(newGroup.toStdString());

    // update widgets
    groupListWidget_->addItem(newGroup);
    initPropertyTreeWidget();
}

void ApplicationModeConfigDialog::handleRemoveGroup() {
    tgtAssert(groupListWidget_, "no group list widget");
    if (!workspace_)
        return;
    ApplicationModeConfiguration& config = workspace_->getApplicationModeConfig();
    tgtAssert(groupListWidget_->count() == static_cast<int>(config.getPropertyGroups().size()), "group count mismatch");

    // determine selected group
    int selectedGroupID = groupListWidget_->currentRow();
    if (selectedGroupID < 0)
        return;
    QString groupName = QString::fromStdString(config.getPropertyGroups().at(selectedGroupID));

    // user confirmation
    if (QMessageBox::question(this, "Remove Property Group", "Remove property group '" + groupName + "'?", QMessageBox::Yes | QMessageBox::Cancel) != QMessageBox::Yes)
        return;

    // remove group
    config.removePropertyGroup(selectedGroupID);

    // update widgets
    delete groupListWidget_->takeItem(selectedGroupID);
    initPropertyTreeWidget();
}

void ApplicationModeConfigDialog::handleGroupMoveUp() {
    tgtAssert(groupListWidget_, "no group list widget");
    if (!workspace_)
        return;
    ApplicationModeConfiguration& config = workspace_->getApplicationModeConfig();

    // determine selected group
    int selectedGroupID = groupListWidget_->currentRow();
    if (selectedGroupID <= 0)
        return;

    // swap groups
    size_t swapID = selectedGroupID-1;
    config.swapPropertyGroups(selectedGroupID, swapID);

    // update widgets
    initPropertyTreeWidget();
    initGroupListWidget();
    groupListWidget_->setCurrentRow((int)swapID);
}

void ApplicationModeConfigDialog::handleGroupMoveDown() {
    tgtAssert(groupListWidget_, "no group list widget");
    if (!workspace_)
        return;
    ApplicationModeConfiguration& config = workspace_->getApplicationModeConfig();

    // determine selected group
    int selectedGroupID = groupListWidget_->currentRow();
    if (selectedGroupID < 0 || selectedGroupID >= static_cast<int>(config.getPropertyGroups().size())-1)
        return;

    // swap groups
    size_t swapID = selectedGroupID+1;
    config.swapPropertyGroups(selectedGroupID, swapID);

    // update widgets
    initPropertyTreeWidget();
    initGroupListWidget();
    groupListWidget_->setCurrentRow((int)swapID);
}

void ApplicationModeConfigDialog::handleGroupRename() {
    tgtAssert(groupListWidget_, "no group list widget");
    if (!workspace_)
        return;
    ApplicationModeConfiguration& config = workspace_->getApplicationModeConfig();
    tgtAssert(groupListWidget_->count() == static_cast<int>(config.getPropertyGroups().size()), "group count mismatch");

    // determine selected group
    int selectedGroupID = groupListWidget_->currentRow();
    if (selectedGroupID < 0)
        return;

    std::string prevName = config.getPropertyGroups().at(selectedGroupID);

    // query group name from user
    QString newName;
    bool inputOk = false;
    do {
        newName = QInputDialog::getText(this, "Rename Property Group", "Enter new name for property group '" + QString::fromStdString(prevName) + "': ", QLineEdit::Normal, newName, &inputOk);
    }
    while (inputOk && config.hasPropertyGroup(newName.toStdString()));
    if (!inputOk)
        return;

    // rename group
    config.renamePropertyGroup(selectedGroupID, newName.toStdString());

    // update widgets
    tgtAssert(groupListWidget_->item(selectedGroupID), "no list item at index");
    groupListWidget_->item(selectedGroupID)->setText(newName);
    initPropertyTreeWidget();
    //updateTreeHeaders();
}

void ApplicationModeConfigDialog::handlePropertyPriorityChange(QTreeWidgetItem* propertyItem, Property* property, std::string groupName, int priority) {
    tgtAssert(propertyItem, "null pointer passed");
    tgtAssert(property, "null pointer passed");
    tgtAssert(propertyTreeWidget_, "no tree widget");

    if (!workspace_)
        return;
    ApplicationModeConfiguration& config = workspace_->getApplicationModeConfig();

    // update group membership
    config.setPropertyGroupMembership(property, groupName, priority);

    // tree widget: remove property from all other groups
    std::vector<std::string> propertyGroups = config.getPropertyGroups();
    tgtAssert(propertyTreeWidget_->columnCount() == static_cast<int>(config.getPropertyGroups().size() + 1), "invalid tree widget column count");
    for (size_t col=1; col<static_cast<size_t>(propertyTreeWidget_->columnCount()); col++) {
        PrioritySelectionBox* colPriorityBox = dynamic_cast<PrioritySelectionBox*>(propertyTreeWidget_->itemWidget(propertyItem, (int)col));
        if (colPriorityBox) {
            if (propertyGroups.at(col-1) != groupName)
                colPriorityBox->setPriority(-1);
        }
        else {
            tgtAssert(false, "item widget is not of type PrioritySelectionBox");
        }
    }

}

void ApplicationModeConfigDialog::handleEntityStateChange(std::string entityName, bool visible) {
    tgtAssert(workspace_,"no workspace");
    ApplicationModeConfiguration& config = workspace_->getApplicationModeConfig();
    config.setMenuEntityVisibility(entityName,visible);
}

void ApplicationModeConfigDialog::handleMainCanvasChange(const QString & text ) {
    if(!workspace_) return;
    ApplicationModeConfiguration& config = workspace_->getApplicationModeConfig();
    config.setMainCanvas(text.toStdString());
}

    //-----------------------
    //    help functions
    //-----------------------
void ApplicationModeConfigDialog::updateWidgetState() {
    tgtAssert(groupListWidget_, "no group list widget");
    tgtAssert(propertyTreeWidget_, "no tree widget");
    tgtAssert(addGroupButton_ && deleteGroupButton_ && moveGroupUpButton_ && moveGroupDownButton_ && renameGroupButton_, "no buttons");

    if (!workspace_ || !workspace_->getProcessorNetwork() || workspace_->getProcessorNetwork()->empty()) {
        addGroupButton_->setEnabled(false);
        deleteGroupButton_->setEnabled(false);
        renameGroupButton_->setEnabled(false);
        moveGroupUpButton_->setEnabled(false);
        moveGroupDownButton_->setEnabled(false);

        groupListWidget_->setEnabled(false);
        propertyTreeWidget_->setEnabled(false);
        mainCanvasRendererCB_->setEnabled(false);
        entityTreeWidget_->setEnabled(false);
    }
    else {
        groupListWidget_->setEnabled(true);
        propertyTreeWidget_->setEnabled(true);
        mainCanvasRendererCB_->setEnabled(true);
        entityTreeWidget_->setEnabled(true);

        addGroupButton_->setEnabled(true);
        bool groupSelected = groupListWidget_->count() > 0 && groupListWidget_->currentRow() >= 0;
        deleteGroupButton_->setEnabled(groupSelected);
        renameGroupButton_->setEnabled(groupSelected);
        moveGroupUpButton_->setEnabled(groupSelected);
        moveGroupDownButton_->setEnabled(groupSelected);
    }
}


void ApplicationModeConfigDialog::expandProcessor(QTreeWidgetItem* qitem){
    ApplicationModeConfigProcessorItem *item = dynamic_cast<ApplicationModeConfigProcessorItem*>(qitem);
    if (!item->getInitialized()){
        item->removeChild(item->child(0));
        addProcessorPropertiesToWidget(item);
        item->setInitialized(true);
    }

}

void ApplicationModeConfigDialog::showEvent(QShowEvent* event) {
    if (!widgetsValid_) {
        qApp->processEvents();
        reinitAllWidgets();
    }

    updateWidgetState();
}

//-----------------------
//    TWidget updates
//-----------------------
void ApplicationModeConfigDialog::addProcessorToTreeWidget(const Processor* processor) {
    tgtAssert(processor, "null pointer passed");
    tgtAssert(workspace_, "no workspace");

    ApplicationModeConfigProcessorItem* processorItem = new ApplicationModeConfigProcessorItem(processor);
    propertyTreeWidget_->addTopLevelItem(processorItem);

    processorItem->addChild(new QTreeWidgetItem());
    //addProcessorPropertiesToWidget(processorItem);
}

void ApplicationModeConfigDialog::addProcessorPropertiesToWidget(ApplicationModeConfigProcessorItem* processorItem) {
    const Processor* processor = processorItem->getProcessor();
    const ApplicationModeConfiguration& config = workspace_->getApplicationModeConfig();
    const std::vector<std::string> propertyGroups = config.getPropertyGroups();
    bool processorHasVisibleProps = false;
    for (size_t j=0; j<processor->getProperties().size(); j++) {
        Property* prop = processor->getProperties().at(j);
        QStringList propName;
        propName.push_back(QString::fromStdString(prop->getGuiName()));
        QTreeWidgetItem* propItem = new QTreeWidgetItem(propName);
        processorItem->addChild(propItem);

        for (size_t k=0; k<propertyGroups.size(); k++) {
            PrioritySelectionBox* priorityBox = new
                PrioritySelectionBox(propItem, prop, propertyGroups.at(k),
                        WEIGHT_COUNT);
            int priority = config.getPropertyGroupPriority(prop,
                    propertyGroups.at(k));
            //int priority = 0;
            priorityBox->setPriority(priority);
            processorHasVisibleProps |= priority >= 0;
            propertyTreeWidget_->setItemWidget(propItem, (int)k+1, priorityBox);
            connect(priorityBox,
                    SIGNAL(priorityChanged(QTreeWidgetItem*, Property*, std::string, int)),
                    this,
                    SLOT(handlePropertyPriorityChange(QTreeWidgetItem*, Property*, std::string, int)));
        }
    }

}

void ApplicationModeConfigDialog::removeProcessorFromTreeWidget(const Processor* processor) {
    tgtAssert(processor, "null pointer passed");
    tgtAssert(propertyTreeWidget_, "no tree widget");

    // find top level item corresponding to removed processor
    QTreeWidgetItem* procItem = 0;
    for (size_t i=0; i<static_cast<size_t>(propertyTreeWidget_->topLevelItemCount()); i++) {
        if (propertyTreeWidget_->topLevelItem((int)i)->text(0).toStdString() == processor->getGuiName()) {
            procItem = propertyTreeWidget_->topLevelItem((int)i);
            break;
        }
    }
    if (!procItem) {
        LWARNINGC("voreen.qt.ApplicationModeConfigDialog", "removeProcessorFromTreeWidget(): processor item not found");
        return;
    }

    // remove item
    propertyTreeWidget_->invisibleRootItem()->removeChild(procItem);
}

void ApplicationModeConfigDialog::renameProcessorInTreeWidget(const std::string& prevName, const std::string& newName) {
    // find processor item (top level)
    QTreeWidgetItem* procItem = 0;
    for (size_t i=0; i<static_cast<size_t>(propertyTreeWidget_->topLevelItemCount()); i++) {
        if (propertyTreeWidget_->topLevelItem((int)i)->text(0).toStdString() == prevName) {
            procItem = propertyTreeWidget_->topLevelItem((int)i);
            break;
        }
    }
    if (!procItem) {
        LWARNINGC("voreen.qt.ApplicationModeConfigDialog", "renameProcessorInTreeWidget(): processor item not found");
        return;
    }

    procItem->setText(0, QString::fromStdString(newName));
}

void ApplicationModeConfigDialog::updateTreeHeaders() {
    if (!workspace_)
        return;
    const ApplicationModeConfiguration& config = workspace_->getApplicationModeConfig();
    const std::vector<std::string> propertyGroups = config.getPropertyGroups();

    QStringList headers;
    headers.push_back("Property");
    for (size_t i=0; i<propertyGroups.size(); i++)
        headers.push_back(QString::fromStdString(propertyGroups.at(i)));
    propertyTreeWidget_->setHeaderLabels(headers);
}


//-----------------------------------------------------------------------------
//                 MenuEntityVisibilityBox
//-----------------------------------------------------------------------------
MenuEntityVisibilityBox::MenuEntityVisibilityBox(std::string name)
    : menuEntityName_(name)
{
    connect(this,SIGNAL(stateChanged(int)),this,SLOT(handleStateChanged(int)));
}

void MenuEntityVisibilityBox::handleStateChanged(int checkState) {
    emit visibilityChanged(menuEntityName_, checkState); //0 = unchecked, 2 = checked
}
//-----------------------------------------------------------------------------
//                 ApplicationModeConfigProcessorItem
//-----------------------------------------------------------------------------
ApplicationModeConfigProcessorItem::ApplicationModeConfigProcessorItem(const Processor* processor)
    :QTreeWidgetItem(QStringList(QString::fromStdString(processor->getGuiName())))
    , processor_(processor)
    , initialized_(false){

}

const Processor* ApplicationModeConfigProcessorItem::getProcessor() const{
    return processor_;
}

bool ApplicationModeConfigProcessorItem::getInitialized() const{
    return initialized_;
}

void ApplicationModeConfigProcessorItem::setInitialized(bool b){
    initialized_ = b;
}
//-----------------------------------------------------------------------------
//                 PrioritySelectionBox
//-----------------------------------------------------------------------------
PrioritySelectionBox::PrioritySelectionBox(QTreeWidgetItem* treeWidgetItem, Property* property, const std::string& groupName, int numPriorities)
    : property_(property)
    , groupName_(groupName)
    , priorityCount_(numPriorities)
    , treeWidgetItem_(treeWidgetItem)
{
    tgtAssert(property_, "null pointer passed");
    tgtAssert(priorityCount_ > 0, "invalid priority count");

    addItem(""); //< is mapped to priority -1
    for (int i=0; i<priorityCount_; i++) {
        addItem(QString::number(i));
    }

    setFixedWidth(40);

    connect(this, SIGNAL(currentIndexChanged(int)), this, SLOT(handleSelectionChange()));
}

void PrioritySelectionBox::setPriority(int priority) {
    blockSignals(true);

    if (priority < 0)
        setCurrentIndex(0); // index 0 is equal to ""
    else if (priority < priorityCount_)
        setCurrentIndex(priority+1); // +1 because priority 0 has index 1 etc.
    else {
        tgtAssert(false, "invalid priority");
    }

    blockSignals(false);
}

void PrioritySelectionBox::handleSelectionChange() {
    emit(priorityChanged(treeWidgetItem_, property_, groupName_, currentIndex()-1));
}

} // namespace
