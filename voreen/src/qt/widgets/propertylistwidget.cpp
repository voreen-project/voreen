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

#include "voreen/qt/widgets/propertylistwidget.h"

#include "voreen/core/voreenapplication.h"

#include "voreen/qt/widgets/property/qpropertywidget.h"
#include "voreen/qt/widgets/property/processorpropertieswidget.h"
#include "voreen/qt/widgets/property/appmodepropertygroupwidget.h"
#include "voreen/qt/widgets/customlabel.h"

#include "voreen/core/network/workspace.h"
#include "voreen/core/datastructures/callback/memberfunctioncallback.h"

#include <QVBoxLayout>
#include <QGridLayout>
#include <QToolButton>
#include <QLabel>

#include "tgt/glcontextmanager.h"

namespace voreen {

PropertyListWidget::PropertyListWidget(QWidget* parent, PropertyListWidget::WidgetMode mode)
    : QScrollArea(parent)
    , workspace_(0)
    , widgetMode_(mode)
    , currentLevelOfDetail_("currentLevelOfDetail","Level of Detail")
    , stackWidget_(0)
    , networkModeContainer_(0)
    , networkModeContainerLayout_(0)
    , appModeContainer_(0)
    , appModeContainerLayout_(0)
    , appModeStatusLabel_(0)
    , appModeWidgetsValid_(false)
{
    setWidgetResizable(true);
    setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
    setFrameStyle(QFrame::NoFrame);
    setContentsMargins(0, 0, 0, 0);

    setMinimumWidth(250);

    currentLevelOfDetail_.addOption("default","Standard",Property::LOD_DEFAULT);
    currentLevelOfDetail_.addOption("advanced","Advanced",Property::LOD_ADVANCED);
    currentLevelOfDetail_.addOption("debug","Debug",Property::LOD_DEBUG);

    currentLevelOfDetail_.onChange(MemberFunctionCallback<PropertyListWidget>(this,&PropertyListWidget::levelOfDetailOnChange));

    createCommonWidgets();
}

PropertyListWidget::~PropertyListWidget() {
    //clear();
}

void PropertyListWidget::levelOfDetailOnChange() {
    setLevelOfDetail(currentLevelOfDetail_.getValue());
}

void PropertyListWidget::setWorkspace(Workspace* workspace) {
    // stop observation of previously assigned workspace and network, if new workspace is assigned
    if (workspace_ && workspace_ != workspace) {
        WorkspaceObserver::stopObservation(workspace_);
        if (workspace_->getProcessorNetwork())
            ProcessorNetworkObserver::stopObservation(workspace_->getProcessorNetwork());
    }

    workspace_ = workspace;

    // update network and register as observer
    if (workspace_) {
        workspace_->addObserver(static_cast<WorkspaceObserver*>(this));
        if (workspace_->getProcessorNetwork())
            workspace_->getProcessorNetwork()->addObserver(static_cast<ProcessorNetworkObserver*>(this));
    }

    // update the widget's state
    networkChanged();
}

void PropertyListWidget::setWidgetMode(WidgetMode mode) {
    //set mode
    widgetMode_ = mode;

    //disable inactive size policy
    if (stackWidget_->currentWidget() !=0) {
            stackWidget_->currentWidget()->setSizePolicy(QSizePolicy::Ignored,
                                                         QSizePolicy::Ignored);
    }
    //load/set right widget
    switch(widgetMode_)  {
    case NETWORK:
        stackWidget_->setCurrentIndex(0);

        if (networkModeProcessorWidgetMap_.empty())
            createNetworkModeWidgets();
        processorsSelected(previouslySelectedProcessors_);
        break;
    case APPLICATION:
        stackWidget_->setCurrentIndex(1);

        if (!appModeWidgetsValid_) {
            clearAppModeWidgets();
            createAppModeWidgets();
        }
        break;
    default:
        tgtAssert(false,"unknown widget mode!");
    }
    //resizes the property list
    stackWidget_->currentWidget()->setSizePolicy(QSizePolicy::Expanding,
                                                 QSizePolicy::Expanding);
    adjustSize();
}

PropertyListWidget::WidgetMode PropertyListWidget::getWidgetMode() const {
    return widgetMode_;
}

void PropertyListWidget::setLevelOfDetail(Property::LevelOfDetail lod) {
    currentLevelOfDetail_.selectByValue(lod);

    std::map<const Processor*, ProcessorPropertiesWidget*>::iterator it;
    for (it = networkModeProcessorWidgetMap_.begin(); it != networkModeProcessorWidgetMap_.end(); ++it) {
        if(std::find(previouslySelectedProcessors_.begin(), previouslySelectedProcessors_.end(), it->first) != previouslySelectedProcessors_.end())
            it->second->setLevelOfDetail(lod);
    }
}

Property::LevelOfDetail PropertyListWidget::getLevelOfDetail() const {
    return currentLevelOfDetail_.getValue();
}

void PropertyListWidget::setState(WidgetMode mode, Property::LevelOfDetail lod) {
    setWidgetMode(mode);
    setLevelOfDetail(lod);
}

void PropertyListWidget::clear() {
    clearNetworkModeWidgets();
    clearAppModeWidgets();
}

void PropertyListWidget::networkChanged() {
    clear();
    setWidgetMode(widgetMode_);
}

void PropertyListWidget::processorAdded(const Processor* processor) {
    if (widgetMode_ == NETWORK) { // generate header widget for the new processor and insert it into the network mode widget map
        tgtAssert(networkModeContainerLayout_, "no network mode container layout");

        ProcessorPropertiesWidget* headerWidget = new ProcessorPropertiesWidget(const_cast<Processor*>(processor), this);
        networkModeProcessorWidgetMap_.insert(std::make_pair(processor, headerWidget));
        networkModeContainerLayout_->addWidget(headerWidget);

        // Assign volume container to headerWidgets for further propagation to volumehandleproperty- and volumecollectionwidgets
        headerWidget->instantiateWidgets();
        headerWidget->setLevelOfDetail(currentLevelOfDetail_.getValue());
        headerWidget->setVisible(true);

        appModeWidgetsValid_ = false;
    }
    else if (widgetMode_ == APPLICATION) {
        appModeWidgetsValid_ = false;
        //clearAppModeWidgets();
        //createAppModeWidgets();
    }
    else {
        tgtAssert(false, "unknown widget mode");
    }
}

void PropertyListWidget::processorRemoved(const Processor* processor){
    if (widgetMode_ == NETWORK) {
        // delete processor's header widget and remove it from the map
        std::map<const Processor*, ProcessorPropertiesWidget*>::iterator it = networkModeProcessorWidgetMap_.find(processor);
        if (it != networkModeProcessorWidgetMap_.end()) {
            if (widget() && widget()->layout())
                widget()->layout()->removeWidget(it->second);
            delete it->second;
            networkModeProcessorWidgetMap_.erase(it);
        }

        appModeWidgetsValid_ = false;
    }
    else if (widgetMode_ == APPLICATION) {
        appModeWidgetsValid_ = false;
        //clearAppModeWidgets();
        //createAppModeWidgets();
    }
    else {
        tgtAssert(false, "unknown widget mode");
    }
}

void PropertyListWidget::processorRenamed(const Processor* processor, const std::string& /*prevName*/){
    std::map<const Processor*, ProcessorPropertiesWidget*>::iterator it = networkModeProcessorWidgetMap_.find(processor);
    if (it != networkModeProcessorWidgetMap_.end()) {
        it->second->updateHeaderTitle();
    }
}

void PropertyListWidget::connectionsChanged() {}

void PropertyListWidget::propertyLinkAdded(const PropertyLink*) {}

void PropertyListWidget::propertyLinkRemoved(const PropertyLink*) {}

void PropertyListWidget::portConnectionAdded(const Port*, const Port*) {}

void PropertyListWidget::portConnectionRemoved(const Port*, const Port*) {}

void PropertyListWidget::applicationModeConfigurationChanged() {
    if (widgetMode_ == NETWORK) {
        appModeWidgetsValid_ = false;
    }
    else if (widgetMode_ == APPLICATION) {
        clearAppModeWidgets();
        createAppModeWidgets();
    }
    else {
        tgtAssert(false, "unknown widget mode");
    }
}

QSize PropertyListWidget::sizeHint() const {
    //return size of current widget(i.e. application or network widget)
    if(stackWidget_ && stackWidget_->currentWidget())
        return stackWidget_->currentWidget()->size();
    //return default value, if widget is not initialized
    return QSize(300, 0);
}

void PropertyListWidget::processorsSelected(const QList<Processor*>& selectedProcessors) {
    if (widgetMode_ == NETWORK) {
        std::map<const Processor*, ProcessorPropertiesWidget*>::iterator it;
        for (it = networkModeProcessorWidgetMap_.begin(); it != networkModeProcessorWidgetMap_.end(); ++it) {
            bool isSelected = std::find(selectedProcessors.begin(), selectedProcessors.end(), it->first) != selectedProcessors.end();
            it->second->setLevelOfDetail(currentLevelOfDetail_.getValue());
            it->second->setVisible(isSelected);
            it->second->setExpanded(isSelected);
        }
    }

    previouslySelectedProcessors_ = selectedProcessors;
}

void PropertyListWidget::createCommonWidgets() {
    tgtAssert(!stackWidget_, "stack widget already created");

    if (!stackWidget_) {
        stackWidget_ = new QStackedWidget();
        setWidget(stackWidget_);

        // network mode container widget
        networkModeContainer_ = new QWidget;
        QVBoxLayout* outerLayout = new QVBoxLayout(networkModeContainer_);
        outerLayout->setContentsMargins(5, 5, 5, 0);
        networkModeContainerLayout_ = new QVBoxLayout();
        networkModeContainerLayout_->setContentsMargins(0, 5, 0, 0);//5,5,5,0
        networkModeContainerLayout_->setAlignment(Qt::AlignTop);
        outerLayout->addLayout(networkModeContainerLayout_);
        outerLayout->addStretch(10);

        //create and add lod property
        PropertyWidget* currentLevelOfDetailWidget = VoreenApplication::app()->createPropertyWidget(&currentLevelOfDetail_);
        if (currentLevelOfDetailWidget) {
            currentLevelOfDetail_.addWidget(currentLevelOfDetailWidget);

            QWidget* lodWidget = new QWidget();
            lodWidget->setObjectName("PLW_LODWIDGET");
            lodWidget->setStyleSheet("QWidget#PLW_LODWIDGET {background-color:white; border: 1px solid black }");
            QGridLayout* gridLayout = new QGridLayout();
            gridLayout->setContentsMargins(4, 4, 4, 4); // 4 4 4 4
            gridLayout->setSpacing(2);
            gridLayout->setColumnStretch(0, 1);
            gridLayout->setColumnStretch(1, 2);

            gridLayout->addWidget(dynamic_cast<QPropertyWidget*>(currentLevelOfDetailWidget)->getOrCreateNameLabel(), 0, 0, 1, 1);
            gridLayout->addWidget(dynamic_cast<QPropertyWidget*>(currentLevelOfDetailWidget), 0, 1, 1, 1);

            lodWidget->setLayout(gridLayout);
            networkModeContainerLayout_->addWidget(lodWidget);
        }

        stackWidget_->addWidget(networkModeContainer_);

        // application mode container widget
        appModeContainer_ = new QWidget();
        appModeContainerLayout_ = new QVBoxLayout();
        appModeContainerLayout_->setContentsMargins(5, 5, 5, 0);
        appModeContainer_->setLayout(appModeContainerLayout_);
        stackWidget_->addWidget(appModeContainer_);
    }

}

void PropertyListWidget::createNetworkModeWidgets() {
    tgtAssert(networkModeContainerLayout_, "no network mode container layout");

    if (!workspace_ || !workspace_->getProcessorNetwork())
        return;

    ProcessorNetwork* network = workspace_->getProcessorNetwork();

    // generate processor property widgets: each processor property contains a processor's properties
    // along with a expansion header
    for (size_t i=0; i<network->numProcessors(); ++i) {
        Processor* proc = network->getProcessors().at(i);
        ProcessorPropertiesWidget* headerWidget = new ProcessorPropertiesWidget(proc, this, true, false);
        networkModeContainerLayout_->addWidget(headerWidget);
        networkModeProcessorWidgetMap_.insert(std::make_pair(proc, headerWidget));
    }
}

void PropertyListWidget::createAppModeWidgets() {
    tgtAssert(appModeContainer_, "no app mode container");
    tgtAssert(appModeContainerLayout_, "no app mode container layout");

    tgtAssert(!appModeStatusLabel_, "app mode status label already created");

    // create group widgets
    if (workspace_ && workspace_->getProcessorNetwork() && !workspace_->getProcessorNetwork()->empty()) {
        ApplicationModeConfiguration config = workspace_->getApplicationModeConfig();
        const std::vector<std::string> propertyGroups = config.getPropertyGroups();
        if (!propertyGroups.empty()) {
            for (size_t i=0; i<propertyGroups.size(); i++) {
                AppModePropertyGroupWidget* groupWidget = new AppModePropertyGroupWidget(propertyGroups.at(i), appModeContainer_, true);
                groupWidget->setProperties(config.getGroupProperties(propertyGroups.at(i)));
                appModeGroupWidgets_.push_back(groupWidget);
                appModeContainerLayout_->addWidget(groupWidget);
            }

            for (size_t i=0; i<appModeGroupWidgets_.size(); i++) {
                appModeGroupWidgets_.at(i)->instantiateWidgets();
                appModeGroupWidgets_.at(i)->setExpanded(false);
            }
        }
        else {
            appModeStatusLabel_ = new QLabel();
            appModeStatusLabel_->setText(
                "No visible properties in application mode. \n\n"
                "Please configure the property visibility \nusing the 'Application Mode Config' tool.");
            appModeContainerLayout_->addWidget(appModeStatusLabel_);
        }

        appModeContainerLayout_->addStretch(10);

        appModeWidgetsValid_ = true;
    }
}

void PropertyListWidget::clearNetworkModeWidgets() {
    tgtAssert(networkModeContainerLayout_, "no network mode container widget");

    setUpdatesEnabled(false);

    // delete network mode widgets and layout item
    std::map<const Processor*, ProcessorPropertiesWidget*>::iterator it;
    for (it = networkModeProcessorWidgetMap_.begin(); it != networkModeProcessorWidgetMap_.end(); ++it) {
        delete it->second;
    }
    // remove remaining layout items
    while (networkModeContainerLayout_->count() > 1) //at pos 0 is the lod proprty
        networkModeContainerLayout_->removeItem(networkModeContainerLayout_->itemAt(1));

    setUpdatesEnabled(true);

    networkModeProcessorWidgetMap_.clear();
    previouslySelectedProcessors_.clear();
}

void PropertyListWidget::clearAppModeWidgets() {
    tgtAssert(appModeContainerLayout_, "no app mode container layout");

    //setUpdatesEnabled(false);

    // delete application mode group widgets and layout items
    delete appModeStatusLabel_;
    appModeStatusLabel_ = 0;

    for (size_t i = 0; i<appModeGroupWidgets_.size(); i++) {
        //appModeContainerLayout_->removeWidget(appModeGroupWidgets_.at(i));
        delete appModeGroupWidgets_.at(i);
    }
    // remove remaining layout items
    while (appModeContainerLayout_->count() > 0)
        appModeContainerLayout_->removeItem(appModeContainerLayout_->itemAt(0));

    appModeGroupWidgets_.clear();
    appModeWidgetsValid_ = false;

    //setUpdatesEnabled(true);
}

} // namespace voreen
