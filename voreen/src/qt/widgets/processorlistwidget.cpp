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

#include "voreen/qt/widgets/processorlistwidget.h"

#include "voreen/core/voreenmodule.h"
#include "voreen/core/ports/port.h"
#include "voreen/core/properties/property.h"

#include "voreen/qt/widgets/lineeditresetwidget.h"
#include "voreen/qt/voreenapplicationqt.h"

#include "voreen/core/utils/stringutils.h"

#include <QMimeData>
#include <QDrag>
#include <QAction>
#include <QApplication>
#include <QContextMenuEvent>
#include <QHeaderView>
#include <QMenu>
#include <QMouseEvent>
#include <QPushButton>
#include <QSplitter>
#include <QStringList>
#include <QTextBrowser>
#include <QVBoxLayout>
#include <QMap>
#include <QList>

namespace voreen {

ProcessorListWidget::ProcessorListWidget(QWidget* parent)
    : QWidget(parent)
    , tree_(new ProcessorListTreeWidget(this))
    , resetIcon_()
    , info_(new QTextBrowser())
    , splitter_(new QSplitter(Qt::Vertical))
    , recentlyUsedModule_(0)
    , recentlyUsedProcessor_(0)
    , resetSettings_(false)
{
    QPushButton* sortButton = new QPushButton(this);
    sortButton->setIcon(QIcon(":/qt/icons/configure.png"));
    sortButton->setGeometry(0, 0, 32, 32);
    sortButton->setFlat(true);
    sortButton->setToolTip(tr("Sort processor list"));

    edit_ = new LineEditResetWidget(this);
    edit_->setToolTip(tr("Filter processor list"));
    QAction* searchAction = new QAction(this);
    searchAction->setShortcut(tr("Ctrl+F"));
    connect(searchAction, SIGNAL(triggered()), edit_, SLOT(setFocus()));
    edit_->addAction(searchAction);

    QVBoxLayout* mainLayout = new QVBoxLayout(this);
    QHBoxLayout* filterSortLayout = new QHBoxLayout();
    info_->resize(120, 100);
    info_->setContextMenuPolicy(Qt::CustomContextMenu);
    info_->setReadOnly(true);
    info_->setOpenExternalLinks(true);

    connect(edit_, SIGNAL(textChanged(const QString&)), tree_, SLOT(filter(const QString&)));
    connect(tree_, SIGNAL(itemSelectionChanged()), this, SLOT(setInfo()));
    connect(sortButton, SIGNAL(clicked()), this, SLOT(sortMenu()));
    connect(this, SIGNAL(sort(int)), tree_, SLOT(sort(int)));
    connect(this, SIGNAL(hideStatus(bool)), tree_, SLOT(hideStatus(bool)));

    connect(this, SIGNAL(showModule(QString, bool)), tree_, SLOT(setModuleNameVisibility(QString, bool)));

    filterSortLayout->addWidget(edit_);
    filterSortLayout->addWidget(sortButton);
    splitter_->addWidget(tree_);
    splitter_->addWidget(info_);
    QList<int> widgetSpacing;
    widgetSpacing << 500 << 100;
    splitter_->setSizes(widgetSpacing);
    mainLayout->addLayout(filterSortLayout);
    mainLayout->addWidget(splitter_);

    showCodeState_ = new QAction(tr("Show Code State"), this);
    showCodeState_->setCheckable(true);
    showCodeState_->setChecked(true);

    const std::vector<const Processor*> processors = VoreenApplication::app()->getSerializableTypes<Processor>();

    //get all moduleNames and set bool for visibility filtration
    std::string moduleName;
    for (size_t i = 0; i < processors.size(); ++i) {
        moduleName = processors.at(i)->getModuleName();
        moduleVisibility_[moduleName] = true;
        tree_->moduleVisibility_[moduleName] = true;
    }

    sortByCategory_ = new QAction("Sort by Category", this);
    sortByModule_ = new QAction("Sort by Module", this);
    sortByModuleThenCategory_ = new QAction("Sort by Module, then Category", this);
    sortByCategory_->setCheckable(true);
    sortByModule_->setCheckable(true);
    sortByModuleThenCategory_->setCheckable(true);

    showCodeStateStable_ = new QAction("Show Code State 'stable'", this);
    showCodeStateExperimental_ = new QAction("Show Code State 'testing'", this);
    showCodeStateBroken_ = new QAction("Show Code State 'experimental'", this);
    showCodeStateStable_->setCheckable(true);
    showCodeStateExperimental_->setCheckable(true);
    showCodeStateBroken_->setCheckable(true);
    moduleVisibilityMenu_ = new QMenu("Show Modules", this);
    std::map<std::string, bool>::iterator it = moduleVisibility_.begin();
    moduleVisibilityMenu_->addAction("Check All");
    moduleVisibilityMenu_->addAction("Uncheck All");
    moduleVisibilityMenu_->addSeparator();
    loadSettings();
    while (it != moduleVisibility_.end()) {
        QAction* tempAction = new QAction(QString::fromStdString((*it).first), this);
        moduleVisibilityActions_[(*it).first] = tempAction;
        moduleVisibilityMenu_->addAction(tempAction);
        tempAction->setCheckable(true);
        tempAction->setChecked((*it).second);
        ++it;
    }
}

ProcessorListWidget::~ProcessorListWidget() {
    if (!resetSettings_)
        saveSettings();
}

void ProcessorListWidget::loadSettings()
{
    QSettings settings;
    settings.beginGroup("ProcessorListWidget");
    if (!settings.contains("core")) {      // first time load without any data written in the settings
        sortByModuleThenCategory_->setChecked(true);         // first time and no data is stored in qsettings
        std::map<std::string, bool>::iterator it = moduleVisibility_.begin();
        while (it != moduleVisibility_.end()) {
            moduleVisibility_[(*it).first] = true;
            tree_->setModuleNameVisibility(QString::fromStdString((*it).first), true);
            it++;
        }
        showCodeState_->setChecked(true);
        showCodeStateStable_->setChecked(true);
        showCodeStateExperimental_->setChecked(true);
        showCodeStateBroken_->setChecked(true);
        tree_->codeState_.insert(Processor::CODE_STATE_STABLE);
        tree_->codeState_.insert(Processor::CODE_STATE_TESTING);
        tree_->codeState_.insert(Processor::CODE_STATE_BROKEN);
        tree_->codeState_.insert(Processor::CODE_STATE_OBSOLETE);
        tree_->codeState_.insert(Processor::CODE_STATE_EXPERIMENTAL);
        setModuleNameVisibility("Check All", true);
    }
    else {
        if (settings.value("sbc").toBool()) {
            sortByCategory_->setChecked(true);
        }
        else if (settings.value("sbm").toBool()) {
            sortByModule_->setChecked(true);
        }
        else if (settings.value("sbmtc").toBool()) {
            sortByModuleThenCategory_->setChecked(true);
        }

        showCodeState_->setChecked(settings.value("showCodeState").toBool());
    }

    std::map<std::string, bool>::iterator it = moduleVisibility_.begin();
    while (it != moduleVisibility_.end()) {
        if (settings.value(QString::fromStdString((*it).first)).toString() == "") {
            moduleVisibility_[(*it).first] = true;
            tree_->setModuleNameVisibility(QString::fromStdString((*it).first), true);
        }
        else {
            moduleVisibility_[(*it).first] = settings.value(QString::fromStdString((*it).first)).toBool();
            tree_->setModuleNameVisibility(QString::fromStdString((*it).first), settings.value(QString::fromStdString((*it).first)).toBool());
        }
        it++;
    }

    if (settings.contains("cst")) {
        showCodeStateStable_->setChecked(settings.value("csw").toBool());
        showCodeStateExperimental_->setChecked(settings.value("cst").toBool());
        showCodeStateBroken_->setChecked(settings.value("csb").toBool());
    }
    if (showCodeStateStable_->isChecked()) {
        tree_->codeState_.insert(Processor::CODE_STATE_STABLE);
        showCodeStateStable_->setChecked(true);
    }
    else {
        tree_->codeState_.erase(Processor::CODE_STATE_STABLE);
        showCodeStateStable_->setChecked(false);
    }
    if (showCodeStateExperimental_->isChecked()) {
        tree_->codeState_.insert(Processor::CODE_STATE_TESTING);
        showCodeStateExperimental_->setChecked(true);
    }
    else {
        tree_->codeState_.erase(Processor::CODE_STATE_TESTING);
        showCodeStateExperimental_->setChecked(false);
    }

    if (showCodeStateBroken_->isChecked()) {
        tree_->codeState_.insert(Processor::CODE_STATE_BROKEN);
        tree_->codeState_.insert(Processor::CODE_STATE_OBSOLETE);
        tree_->codeState_.insert(Processor::CODE_STATE_EXPERIMENTAL);
        showCodeStateBroken_->setChecked(true);
    }
    else {
        tree_->codeState_.erase(Processor::CODE_STATE_BROKEN);
        tree_->codeState_.erase(Processor::CODE_STATE_OBSOLETE);
        tree_->codeState_.erase(Processor::CODE_STATE_EXPERIMENTAL);
        showCodeStateBroken_->setChecked(false);
    }

    if (sortByCategory_->isChecked()) {
        sort(ProcessorListTreeWidget::GroupTypeCategory);
    }
    else if (sortByModule_->isChecked()) {
        sort(ProcessorListTreeWidget::GroupTypeModule);
    }
    else if (sortByModuleThenCategory_->isChecked()) {
        sort(ProcessorListTreeWidget::GroupTypeModuleCategory);
    }

    if (settings.contains("infoBoxSize")) {
        QList<int> widgetSpacing;
        widgetSpacing << settings.value("processorListSize").toInt() << settings.value("infoBoxSize").toInt();
        splitter_->setSizes(widgetSpacing);
    }

    settings.endGroup();

    emit hideStatus(showCodeState_->isChecked());
}

void ProcessorListWidget::saveSettings() {
    QSettings settings;
    settings.beginGroup("ProcessorListWidget");
    std::map<std::string, bool>::iterator it = moduleVisibility_.begin();
    while (it != moduleVisibility_.end()) {
        QVariant val = QVariant(moduleVisibility_[(*it).first]);
        settings.setValue(QString::fromStdString((*it).first), val);
        it++;
    }
    settings.setValue("filtertext", QString::fromStdString(tree_->filterText_));

    settings.setValue("sbc", sortByCategory_->isChecked());
    settings.setValue("sbm", sortByModule_->isChecked());
    settings.setValue("sbmtc", sortByModuleThenCategory_->isChecked());

    settings.setValue("csw", showCodeStateStable_->isChecked());
    settings.setValue("cst", showCodeStateExperimental_->isChecked());
    settings.setValue("csb", showCodeStateBroken_->isChecked());

    settings.setValue("showCodeState", showCodeState_->isChecked());
    settings.setValue("processorListSize", tree_->height());
    settings.setValue("infoBoxSize", info_->height());
    settings.endGroup();
}

void ProcessorListWidget::setInfo() {
    clearInfo();

    QTreeWidgetItem* currentItem = tree_->currentItem();
    if (currentItem) {
        std::string selectionName = currentItem->text(0).toStdString();
        const Processor* processor = dynamic_cast<const Processor*>(VoreenApplication::app()->getSerializableType(selectionName));
        if (processor)
            setInfo(processor);
        else {
            // selected item is a module or a category
            const std::vector<VoreenModule*> modules = VoreenApplication::app()->getModules();
            foreach(VoreenModule* module, modules) {
                if (module->getGuiName() == selectionName)
                    setInfo(module);
            }
        }
    }
}

void ProcessorListWidget::reloadInfoText() {
    if (recentlyUsedProcessor_)
        setInfo(recentlyUsedProcessor_);
    else if (recentlyUsedModule_)
        setInfo(recentlyUsedModule_);
}

void ProcessorListWidget::clearInfo() {
    recentlyUsedModule_ = 0;
    recentlyUsedProcessor_ = 0;
    info_->setHtml("");
}

void ProcessorListWidget::setInfo(const Processor* processor) {
    tgtAssert(processor, "null pointer passed");
    std::string processorClass = processor->getClassName();
    std::string processorModule = processor->getModuleName();

    clearInfo();
    recentlyUsedProcessor_ = processor;

    std::string info;
    info = "<b>" + processorClass + "</b><br>";
    info += processor->getDescription();

    std::vector<std::pair<std::string, std::string> > ports;
    const std::vector<Port*>& tPorts = processor->getPorts();
    for (size_t i = 0; i<tPorts.size(); i++) {
        if (tPorts[i]->getDescription() != "")
            ports.push_back(std::pair<std::string, std::string>(tPorts[i]->getGuiName(), tPorts[i]->getDescription()));
    }
    if (!ports.empty()) {
        info += "<br><br>";
        info += "<i><b>Ports</b></i><br>";
        for (size_t i = 0; i < ports.size(); ++i) {
            std::pair<std::string, std::string> p = ports[i];
            info += "<b>" + p.first + "</b>: " + p.second;
            if (i != ports.size() - 1) {
                info += "<br>";
            }
        }
    }

    std::vector<std::pair<std::string, std::string> > properties;
    const std::vector<Property*>&  props = processor->getProperties();
    for (size_t i = 0; i<props.size(); i++) {
        if (props[i]->getDescription() != "")
            properties.push_back(std::pair<std::string, std::string>(props[i]->getGuiName(), props[i]->getDescription()));
    }
    if (!properties.empty()) {
        info += "<br><br>";
        info += "<i><b>Properties</b></i><br>";
        std::pair<std::string, std::string> p;
        for (size_t i = 0; i < properties.size(); ++i) {
            std::pair<std::string, std::string> p = properties[i];
            info += "<b>" + p.first + "</b>: " + p.second;
            if (i != properties.size() - 1) {
                info += "<br>";
            }
        }
    }

    info_->setHtml(QString::fromStdString(info));
}

void ProcessorListWidget::setInfo(VoreenModule* module) {
    clearInfo();
    recentlyUsedModule_ = module;
    info_->setHtml(QString::fromStdString("<b>" + module->getGuiName() + "</b><br>" + module->getDescription()));
}

void ProcessorListWidget::sortMenu() {
    QMenu* menu = new QMenu();
    QActionGroup* sortActions = new QActionGroup(this);
    sortActions->addAction(sortByCategory_);
    sortActions->addAction(sortByModule_);
    sortActions->addAction(sortByModuleThenCategory_);
    // sort by category
    menu->addAction(sortByCategory_);
    // sort by module
    menu->addAction(sortByModule_);
    //sort by module then category
    menu->addAction(sortByModuleThenCategory_);
    menu->addSeparator();

    menu->addAction(showCodeState_);
    menu->addSeparator();

    menu->addAction(showCodeStateStable_);
    menu->addAction(showCodeStateExperimental_);
    menu->addAction(showCodeStateBroken_);
    menu->addSeparator();

    menu->addMenu(moduleVisibilityMenu_);
    QAction* action = menu->exec(QCursor::pos());
    if (action) {
        if (action == sortByCategory_) {
            emit sort(ProcessorListTreeWidget::GroupTypeCategory);
        }
        else if (action == sortByModule_) {
            emit sort(ProcessorListTreeWidget::GroupTypeModule);
        }
        else if (action == sortByModuleThenCategory_) {
            emit sort(ProcessorListTreeWidget::GroupTypeModuleCategory);
        }
        if (action == showCodeStateStable_ || action == showCodeStateExperimental_ || action == showCodeStateBroken_) {
            if (action == showCodeStateStable_) {
                if (showCodeStateStable_->isChecked())
                    tree_->codeState_.insert(Processor::CODE_STATE_STABLE);
                else
                    tree_->codeState_.erase(Processor::CODE_STATE_STABLE);
            }
            else if (action == showCodeStateExperimental_) {
                if (showCodeStateExperimental_->isChecked())
                    tree_->codeState_.insert(Processor::CODE_STATE_TESTING);
                else
                    tree_->codeState_.erase(Processor::CODE_STATE_TESTING);
            }
            else if (action == showCodeStateBroken_) {
                if (showCodeStateBroken_->isChecked()) {
                    tree_->codeState_.insert(Processor::CODE_STATE_BROKEN);
                    tree_->codeState_.insert(Processor::CODE_STATE_OBSOLETE);
                    tree_->codeState_.insert(Processor::CODE_STATE_EXPERIMENTAL);
                }
                else {
                    tree_->codeState_.erase(Processor::CODE_STATE_BROKEN);
                    tree_->codeState_.erase(Processor::CODE_STATE_OBSOLETE);
                    tree_->codeState_.erase(Processor::CODE_STATE_EXPERIMENTAL);
                }
            }
            if (sortByCategory_->isChecked())
                emit sort(ProcessorListTreeWidget::GroupTypeCategory);
            if (sortByModule_->isChecked())
                emit sort(ProcessorListTreeWidget::GroupTypeModule);
            if (sortByModuleThenCategory_->isChecked())
                emit sort(ProcessorListTreeWidget::GroupTypeModuleCategory);
        }
        else if (action == showCodeState_) {
            emit hideStatus(showCodeState_->isChecked());
        }
        else {
            emit showModule(action->text(), action->isChecked());
            setModuleNameVisibility(action->text().toStdString(), action->isChecked());
        }
    }
}

void ProcessorListWidget::setModuleNameVisibility(std::string name, bool visibility) {
    if (name == "Uncheck All" || name == "Check All"){
        bool checked;
        if (name == "Uncheck All") {
            checked = false;
        }
        else
            checked = true;

        std::map<std::string, bool>::iterator it = moduleVisibility_.begin();
        while (it != moduleVisibility_.end()) {
            (*it).second = checked;
            it++;
        }
        std::map<std::string, QAction*>::iterator it2 = moduleVisibilityActions_.begin();
        while (it2 != moduleVisibilityActions_.end()) {
            (*it2).second->setChecked(checked);
            it2++;
        }

    }
    else
        moduleVisibility_[name] = visibility;
    //moduleVisibilityActions_[name]->setChecked(visibility);
    //saveSettings();

}

void ProcessorListWidget::processorsSelected(const QList<Processor*>& processors) {
    if (processors.size() == 1)
        setInfo(processors.front());
    else
        clearInfo();
}

void ProcessorListWidget::resetSettings() {
    QSettings settings;
    settings.remove("ProcessorListWidget");

    resetSettings_ = true;
}

void ProcessorListWidget::keyReleaseEvent(QKeyEvent* event) {

    if (event->key() != Qt::Key_Return) {
        QWidget::keyReleaseEvent(event);
        return;
    }

    ProcessorListItem* processor = nullptr;
    if (tree_->currentItem() && tree_->hasFocus())
        processor = dynamic_cast<ProcessorListItem*>(tree_->currentItem());
    else {

        // Don't add anything, in case no input exists.
        if (edit_->text().trimmed().isEmpty()) {
            QWidget::keyReleaseEvent(event);
            return;
        }

        // Search the first valid result.
        QTreeWidgetItem* item = tree_->topLevelItem(0);
        while (!processor && (item = tree_->itemBelow(item))) {
            processor = dynamic_cast<ProcessorListItem*>(item);
        }
    }

    // A valid result has been found.
    if (processor) {
        // Create a new processor in the Networkeditor.
        event->accept();
        emit processorAdded(QString(processor->getId().c_str()));
    }
    else
        QWidget::keyReleaseEvent(event);
}

// ----------------------------------------------------------------------------

ProcessorListTreeWidget::ProcessorListTreeWidget(ProcessorListWidget* processorListWidget, QWidget* parent)
    : QTreeWidget(parent)
    , processorListWidget_(processorListWidget)
    , sortType_(ProcessorListTreeWidget::GroupTypeCategory)
    , showCodeState_(true)
{
    setHeaderHidden(true);
    QStringList headeritems = QStringList();
    headeritems << "Processor" << "State";
    setHeaderLabels(headeritems);

    header()->setSectionResizeMode(0, QHeaderView::Stretch);
    header()->setSectionResizeMode(1, QHeaderView::ResizeToContents);
    header()->setStretchLastSection(false);

    setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);

    setColumnCount(2);

    // Setup Cache for Icons.
    codeStateIcons_[Processor::CODE_STATE_BROKEN]  = QIcon(":/qt/icons/processor-broken.png");
    codeStateIcons_[Processor::CODE_STATE_TESTING] = QIcon(":/qt/icons/processor-testing.png");
    codeStateIcons_[Processor::CODE_STATE_STABLE]  = QIcon(":/qt/icons/processor-stable.png");

    const std::vector<VoreenModule*> modules = VoreenApplication::app()->getModules();
    foreach(VoreenModule* module, modules) {
        std::string moduleName = module->getGuiName();
        moduleVisibility_[moduleName] = true;
    }
    sortByCategory();
}

void ProcessorListTreeWidget::filter(const QString& text) {
    filterText_ = text.toStdString();
    items_.clear();
    clear();
    if (sortType_ == ProcessorListTreeWidget::GroupTypeCategory)
        sortByCategory();
    else if (sortType_ == ProcessorListTreeWidget::GroupTypeModule)
        sortByModuleName();
    else if (sortType_ == ProcessorListTreeWidget::GroupTypeModuleCategory)
        sortByModuleThenCategory();
}

void ProcessorListTreeWidget::mousePressEvent(QMouseEvent *event) {
    QTreeWidget::mousePressEvent(event);

    // if not left button - return
    if (!(event->buttons() & Qt::LeftButton))
        return;

    // if no item selected, return
    if (!currentItem())
        return;

    ProcessorListItem* item = dynamic_cast<ProcessorListItem*>(itemAt(event->pos()));
    // if no ProcessorListItem selected, return
    if (!item) {
        const std::vector<VoreenModule*> modules = VoreenApplication::app()->getModules();
        foreach(VoreenModule* module, modules) {
            QTreeWidgetItem* item = itemAt(event->pos());

            if (item && (item->text(0).toStdString() == module->getGuiName()))
                processorListWidget_->setInfo(module);
        }
        return;

    }
    QString idStr(item->getId().c_str());

    QMimeData* mimeData = new QMimeData;
    mimeData->setText(idStr);

    QDrag* drag = new QDrag(processorListWidget_);
    drag->setMimeData(mimeData);
    //drag->setPixmap(pixmap);
    //drag->setHotSpot(event->pos()-item->pos());

    drag->start(Qt::CopyAction);

}

void ProcessorListTreeWidget::setModuleNameVisibility(QString module, bool visible) {
    if (module == "Uncheck All" || module == "Check All"){
        bool checked = module != "Unchecked All";

        std::map<std::string, bool>::iterator it = moduleVisibility_.begin();
        while (it != moduleVisibility_.end()) {
            (*it).second = checked;
            it++;
        }
    }
    else
        moduleVisibility_[module.toStdString()] = visible;
    //filter(QString::fromStdString(filterText_));
    sort(sortType_);

}

std::vector<const Processor*> ProcessorListTreeWidget::getVisibleProcessors() const {
    const std::vector<const Processor*>& knownProcessors = VoreenApplication::app()->getSerializableTypes<Processor>();
    std::vector<const Processor*> visibleProcessors;
    for (size_t i = 0; i < knownProcessors.size(); ++i) {
        const Processor* processor = knownProcessors[i];
        QString procClassName = QString::fromStdString(processor->getClassName());
        if (codeState_.find(processor->getCodeState()) != codeState_.end()) {
            std::string moduleName = processor->getModuleName();
            std::map<std::string, bool>::const_iterator it = moduleVisibility_.find(moduleName);
            if (it != moduleVisibility_.end() && (it->second == true)) {
                if (procClassName.indexOf(QString::fromStdString(filterText_), 0, Qt::CaseInsensitive) != -1)
                    visibleProcessors.push_back(processor);
            }
        }
    }

    return visibleProcessors;
}

void ProcessorListTreeWidget::sortByCategory() {
    sortType_ = ProcessorListTreeWidget::GroupTypeCategory;
    items_.clear();
    clear();

    QList<QTreeWidgetItem*> categoryItems = createCategoryHierarchy(getVisibleProcessors(), 0);
    addTopLevelItems(categoryItems);

    expandAll();
    if (!showCodeState_)
        header()->hide();
    else
        header()->show();
}
void ProcessorListTreeWidget::sortByModuleName() {
    sortType_ = ProcessorListTreeWidget::GroupTypeModule;
    items_.clear();

    const std::vector<const Processor*> visibleProcessors = getVisibleProcessors();
    clear();

    QString moduleIdentifier("");
    std::map<QString, QTreeWidgetItem*> modules;

    for (size_t i = 0; i < visibleProcessors.size(); ++i) {
        moduleIdentifier = QString::fromStdString(visibleProcessors[i]->getModuleName());
        if (modules.find(moduleIdentifier) == modules.end()) {
            QTreeWidgetItem* item = new QTreeWidgetItem(QStringList(moduleIdentifier));
            QFont font = item->font(0);
            font.setBold(true);
            item->setFont(0, font);
            modules[moduleIdentifier] = item;
        }
    }

    for (size_t i = 0; i < visibleProcessors.size(); ++i) {
        moduleIdentifier = QString::fromStdString(visibleProcessors[i]->getModuleName());
        ProcessorListItem* tempListItem = new ProcessorListItem(visibleProcessors[i]->getClassName());

        modules[moduleIdentifier]->addChild(tempListItem);

        if (showCodeState_) {
            setCodeStateIcon(visibleProcessors[i]->getClassName(), tempListItem);
        }
    }

    std::map<QString, QTreeWidgetItem*>::iterator it = modules.begin();
    while (it != modules.end()) {
        addTopLevelItem((*it).second);
        ++it;
    }

    expandAll();
    if (!showCodeState_)
        header()->hide();
    else
        header()->show();
}

void ProcessorListTreeWidget::sortByModuleThenCategory() {
    sortType_ = ProcessorListTreeWidget::GroupTypeModuleCategory;
    items_.clear();
    clear();

    const std::vector<const Processor*> visibleProcessors = getVisibleProcessors();

    // group visible processors by module
    QMap<QString, std::vector<const Processor*> > moduleToProcessorsMap;
    for (size_t i = 0; i < visibleProcessors.size(); ++i) {
        QString moduleIdentifier = QString::fromStdString(visibleProcessors[i]->getModuleName());
        if (!moduleToProcessorsMap.contains(moduleIdentifier))
            moduleToProcessorsMap.insert(moduleIdentifier, std::vector<const Processor*>());
        moduleToProcessorsMap[moduleIdentifier].push_back(visibleProcessors[i]);
    }

    // create tree item with sub-hierarchy for each visible module
    QList<QString> modules = moduleToProcessorsMap.keys();
    QTreeWidgetItem* coreModuleItem = 0;
    QTreeWidgetItem* baseModuleItem = 0;
    for (int i = 0; i<modules.size(); i++) {
        QString curModule = modules.at(i);
        std::vector<const Processor*> moduleProcessors = moduleToProcessorsMap[curModule];

        // module tree item
        QTreeWidgetItem* moduleItem = new QTreeWidgetItem(QStringList(curModule));
        QFont font = moduleItem->font(0);
        font.setBold(true);
        moduleItem->setFont(0, font);

        // category sub-hierarchy
        QList<QTreeWidgetItem*> categoryHierarchy = createCategoryHierarchy(moduleProcessors, 0);
        moduleItem->addChildren(categoryHierarchy);

        if (curModule == "Core")
            coreModuleItem = moduleItem; //< will be inserted later
        else if (curModule == "Base")
            baseModuleItem = moduleItem; //< will be inserted later
        else
            items_.append(moduleItem);
    }

    // add and sort top level items
    addTopLevelItems(items_);
    QTreeWidget::sortItems(0, Qt::AscendingOrder);

    // insert core and base module items at top
    if (baseModuleItem)
        insertTopLevelItem(0, baseModuleItem);
    if (coreModuleItem)
        insertTopLevelItem(0, coreModuleItem);
    /*else
    LWARNINGC("voreenqt.ProcessorListTreeWidget", "no core module item"); */

    expandAll();
    if (!showCodeState_)
        header()->hide();
    else
        header()->show();

}

QList<QTreeWidgetItem*> ProcessorListTreeWidget::createCategoryHierarchy(
    const std::vector<const Processor*>& processors, int categoryLevel) const {
    tgtAssert(categoryLevel >= 0, "invalid category level");

    // separate processors at current category level from those with more sub-categories
    std::map<std::string, QTreeWidgetItem*> categoryToTreeItemMap;
    std::map<std::string, std::vector<const Processor*> > categoryToProcessorsMap; //< for processors with more sub-cats
    for (size_t i = 0; i<processors.size(); i++) {
        const Processor* curProc = processors.at(i);

        // extract sub-category string for current level
        std::vector<std::string> subcategories = strSplit(curProc->getCategory(), '/');
        tgtAssert(subcategories.size() >(size_t)categoryLevel, "too few sub-categories for this category level");
        std::string subCategoryStr = subcategories.at(categoryLevel);

        // create tree item for sub-category, if not already present
        if (categoryToTreeItemMap.find(subCategoryStr) == categoryToTreeItemMap.end()) {
            categoryToTreeItemMap.insert(std::make_pair(subCategoryStr,
                new QTreeWidgetItem(QStringList(QString::fromStdString(subCategoryStr)))));
        }
        QTreeWidgetItem* curCategoryItem = categoryToTreeItemMap[subCategoryStr];

        // if processor has no more sub-categories => add processor item directly
        if (subcategories.size() == static_cast<size_t>(categoryLevel + 1)) {
            ProcessorListItem* processorItem = new ProcessorListItem(curProc->getClassName());
            curCategoryItem->addChild(processorItem);
            if (showCodeState_)
                setCodeStateIcon(curProc->getClassName(), processorItem);
        }
        // else queue processor for sub-category creation
        else {
            if (categoryToProcessorsMap.find(subCategoryStr) == categoryToProcessorsMap.end())
                categoryToProcessorsMap.insert(std::make_pair(subCategoryStr, std::vector<const Processor*>()));
            categoryToProcessorsMap[subCategoryStr].push_back(curProc);
        }
    }

    // recursively create sub-hierarchies for processors with more sub-categories than current category level
    for (std::map<std::string, std::vector<const Processor*> >::const_iterator it = categoryToProcessorsMap.begin();
        it != categoryToProcessorsMap.end(); it++) {
        std::string subCategory = it->first;
        QTreeWidgetItem* categoryItem = categoryToTreeItemMap[subCategory];
        tgtAssert(categoryItem, "not category tree item");

        QList<QTreeWidgetItem*> subHierarchy = createCategoryHierarchy(it->second, categoryLevel + 1);
        categoryItem->addChildren(subHierarchy);
    }

    // return created top-level tree items
    QList<QTreeWidgetItem*> resultList;
    for (std::map<std::string, QTreeWidgetItem*>::const_iterator it = categoryToTreeItemMap.begin();
        it != categoryToTreeItemMap.end(); it++) {
        resultList.append(it->second);
    }
    return resultList;
}

void ProcessorListTreeWidget::setCodeStateIcon(const std::string& classname, QTreeWidgetItem* item) const {
    tgtAssert(dynamic_cast<const Processor*>(VoreenApplication::app()->getSerializableType(classname)), "unknown class name");
    // Set Tool tip.
    Processor::CodeState codeState = dynamic_cast<const Processor*>(VoreenApplication::app()->getSerializableType(classname))->getCodeState();
    switch (codeState) {
    case Processor::CODE_STATE_BROKEN:
        item->setToolTip(1, "Code state: broken");
        break;
    case Processor::CODE_STATE_OBSOLETE:
        item->setToolTip(1, "Code state: obsolete");
        break;
    case Processor::CODE_STATE_EXPERIMENTAL:
        item->setToolTip(1, "Code state: experimental");
        break;
    case Processor::CODE_STATE_TESTING:
        item->setToolTip(1, "Code state: testing");
        break;
    case Processor::CODE_STATE_STABLE:
        item->setToolTip(1, "Code state: stable");
        break;
    default:
        item->setToolTip(1, "Code state: unknown");
        break;
    }

    // Set Icon. Use CODE_STATE_BROKEN as fallback.
    if (codeStateIcons_.find(codeState) != codeStateIcons_.end())
        item->setIcon(1, codeStateIcons_.at(codeState));
    else
        item->setIcon(1, codeStateIcons_.at(Processor::CODE_STATE_BROKEN));
}

void ProcessorListTreeWidget::hideStatus(bool hidden) {
    showCodeState_ = hidden;
    filter(QString::fromStdString(filterText_));
}

void ProcessorListTreeWidget::sort(int type) {
    if (type == ProcessorListTreeWidget::GroupTypeCategory)
        sortByCategory();
    if (type == ProcessorListTreeWidget::GroupTypeModule)
        sortByModuleName();
    if (type == ProcessorListTreeWidget::GroupTypeModuleCategory)
        sortByModuleThenCategory();
}


// ----------------------------------------------------------------------------

ProcessorListItem::ProcessorListItem(const std::string& id)
    : QTreeWidgetItem(QStringList(QString(id.c_str()))),
    id_(id)
{}


} //namespace voreen
