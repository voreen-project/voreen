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

#include "tgt/matrix.h"

#include "voreen/qt/widgets/volumeviewer.h"
#include "voreen/qt/widgets/volumeviewhelper.h"
#include "voreen/qt/widgets/rawvolumewidget.h"

#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/datastructures/volume/volumedecorator.h"
#include "voreen/core/datastructures/volume/volumedisk.h"
#include "voreen/core/datastructures/volume/volumegl.h"
#include "voreen/core/io/progressbar.h"
#include "voreen/core/io/volumeserializerpopulator.h"
#include "voreen/core/network/processornetwork.h"
#include "voreen/core/ports/volumeport.h"
#include "voreen/core/voreenapplication.h"
#include "voreen/core/datastructures/octree/volumeoctreebase.h"

#include "voreen/qt/voreenapplicationqt.h"
#include "voreen/qt/progressdialog.h"

#include <set>

#include <QCheckBox>
#include <QErrorMessage>
#include <QFile>
#include <QFileDialog>
#include <QHeaderView>
#include <QInputDialog>
#include <QKeyEvent>
#include <QLabel>
#include <QListWidget>
#include <QListWidgetItem>
#include <QMenu>
#include <QMessageBox>
#include <QPixmap>
#include <QToolButton>
#include <QUrl>
#include <QVBoxLayout>
#include <QMessageBox>

#include <QApplication>

namespace {
#ifdef __APPLE__
    int fontSize = 13;
#else
    int fontSize = 8;
#endif
}

namespace voreen {

QRCTreeWidget::QRCTreeWidget(QWidget* parent) : QTreeWidget(parent) {
    setSelectionMode(QAbstractItemView::ExtendedSelection);
    setTextElideMode(Qt::ElideMiddle);
}

void QRCTreeWidget::contextMenuEvent(QContextMenuEvent* cmevent) {
    if (itemAt(cmevent->pos())){
        QMenu* menu = new QMenu();
        menu->addAction(QPixmap(":/voreenve/icons/export.png"), "Export...");
        QAction* action = menu->exec(mapToGlobal(cmevent->pos()));
        if(action) {
            if(action->text() == "Export...")
                emit exportVolumes();
            }
    }
}

void QRCTreeWidget::keyPressEvent(QKeyEvent* event) {
    QTreeWidget::keyPressEvent(event);
}

// ---------------------------------------------------------------------------

const std::string VolumeViewer::loggerCat_("voreen.qt.VolumeContainerWidget");

VolumeViewer::VolumeViewer(QWidget* parent)
    : QWidget(parent, Qt::Tool)
    , evaluator_(0)
    , volumeIOHelper_(parent)
    , updateRequired_(true)
{
    setFocusPolicy(Qt::StrongFocus);

    QVBoxLayout* mainLayout = new QVBoxLayout(this);
    QHBoxLayout* buttonLayout = new QHBoxLayout();

    volumeInfos_ = new QRCTreeWidget();
    QTreeWidgetItem* header = volumeInfos_->headerItem();
    header->setText(0, tr("Volumes"));
    volumeInfos_->header()->hide();
    volumeInfos_->setColumnCount(1);
    volumeInfos_->show();
    volumeInfos_->setIconSize(QSize(63,63));

    updateButton_ = new QPushButton(tr("Update"), this);
    updateButton_->setIcon(QPixmap(":/icons/reload.png"));
    updateButton_->setToolTip(tr("Update view from network"));

    autoUpdateCheckbox_ = new QCheckBox(tr("auto"));
    autoUpdateCheckbox_->setToolTip(tr("Automatically update view after each network processing (may be slow!)"));
    autoUpdateCheckbox_->setChecked(true);

    containerInfo_ = new QLabel(this);

    buttonLayout->addWidget(updateButton_);
    buttonLayout->addWidget(autoUpdateCheckbox_);
    buttonLayout->addStretch();
    buttonLayout->addWidget(containerInfo_);

    mainLayout->addWidget(volumeInfos_);
    mainLayout->addLayout(buttonLayout);

    connect(updateButton_, SIGNAL(clicked()),
        this, SLOT(updateFromNetwork()));
    connect(autoUpdateCheckbox_, SIGNAL(stateChanged(int)),
        this, SLOT(updateStateChanged(int)));

    connect(volumeInfos_, SIGNAL(exportVolumes()),
        this, SLOT(exportVolumes()));
    connect(&volumeIOHelper_, SIGNAL(volumeSaved(const VolumeBase*, const std::string&)),
        this, SLOT(volumeSaved(const VolumeBase*, const std::string&)));

    setMinimumWidth(250);
    setMinimumHeight(125);
}

VolumeViewer::~VolumeViewer() {
}

void VolumeViewer::resizeOnCollapse(bool) {
    //adjustSize();
}

QSize VolumeViewer::sizeHint() const {
    return QSize(300, 300);
}

void VolumeViewer::setNetworkEvaluator(NetworkEvaluator* evaluator) {
    if (evaluator_)
        evaluator_->removeObserver(static_cast<NetworkEvaluatorObserver*>(this));

    evaluator_ = evaluator;

    if (evaluator_)
        evaluator_->addObserver(static_cast<NetworkEvaluatorObserver*>(this));

    updateRequired_ = true;
    update();
}

void VolumeViewer::volumeDelete(const VolumeBase* /*source*/) {
    // force update even with auto-update disabled,
    // since we do not want to store invalid handle pointers
    if (isVisible()) {
        updateRequired_ = true;
        update();
    }
}

void VolumeViewer::volumeChange(const VolumeBase* /*source*/) {
    // force update even with auto-update disabled,
    // since we do not want to store invalid handle pointers
    if (isVisible()) {
        updateRequired_ = true;
        update();
    }
}

void VolumeViewer::afterNetworkProcess() {
    if (isVisible() && autoUpdateCheckbox_->isChecked()) {
        updateRequired_ = true;
        update();
    }
}

void VolumeViewer::afterNetworkInitialize() {
    if (isVisible() && autoUpdateCheckbox_->isChecked()) {
        updateRequired_ = true;
        update();
    }
}

void VolumeViewer::beforeNetworkDeinitialize() {
    clear();
}

void VolumeViewer::clear() {
    volumeHandles_.clear();
    handleToHashMap_.clear();
    handleToPortMap_.clear();

    tgtAssert(volumeInfos_, "no volumeInfos_");
    tgtAssert(containerInfo_, "no containerInfo_");
    volumeInfos_->clear();
    containerInfo_->clear();

    if (parentWidget() && parentWidget()->parentWidget())
        parentWidget()->parentWidget()->setWindowTitle("Volume Viewer");
}

void VolumeViewer::paintEvent(QPaintEvent* event) {
    if (updateRequired_) {
        if (!evaluator_) {
            updateRequired_ = false;
            clear();
        }
        else if (!evaluator_->isLocked()) {
            updateRequired_ = false;
            updateFromNetwork();
        }
    }

    QWidget::paintEvent(event);
}

void VolumeViewer::showEvent(QShowEvent* /*event*/) {
    if (isVisible() && autoUpdateCheckbox_->isChecked()) {
        updateRequired_ = true;
        update();
    }
}

void VolumeViewer::updateFromNetwork() {

    if (!evaluator_ || !evaluator_->getProcessorNetwork()) {
        clear();
        return;
    }

    const ProcessorNetwork* network = evaluator_->getProcessorNetwork();
    tgtAssert(network, "no network");

    // retrieve volumes from network (outports)
    std::vector<const VolumeBase*> tempVolumes;
    std::map<const VolumeBase*, std::string> tempHashMap_;
    std::map<const VolumeBase*, Port*> tempPortMap_;
    for (size_t i=0; i<network->getProcessors().size(); i++) {
        Processor* processor = network->getProcessors().at(i);
        tgtAssert(processor, "no processor");
        const std::vector<Port*>& ports = processor->getOutports();
        for (size_t j=0; j<ports.size(); j++) {
            if (dynamic_cast<VolumePort*>(ports.at(j))) {
                const VolumeBase* handle = static_cast<VolumePort*>(ports.at(j))->getData();

                // a VolumeDecorator without a decorated volume is not a regular volume and is ignored here
                const VolumeDecoratorIdentity* vdp = dynamic_cast<const VolumeDecoratorIdentity*>(handle);
                bool isRegularVolume = (vdp == 0 || vdp->getDecorated() != 0);

                if (handle && isRegularVolume && !tempHashMap_.count(handle)) {
                    tempVolumes.push_back(handle);
                    std::string hashString = handle->getHash();
                    // add representations to hash string so representation changes are detected
                    if (handle->hasRepresentation<VolumeDisk>())
                        hashString += "-Disk";
                    if (handle->hasRepresentation<VolumeRAM>())
                        hashString += "-RAM";
                    if (handle->hasRepresentation<VolumeOctreeBase>())
                        hashString += "-Octree";
                    if (handle->hasRepresentation<VolumeGL>())
                        hashString += "-VolumeGL";

                    tempHashMap_.insert(std::make_pair(handle, hashString));
                    tempPortMap_.insert(std::make_pair(handle, ports.at(j)));
                }
            }
        }
    }

    // check wether network volumes have changed
    bool changed = false;
    if (tempVolumes.size() != volumeHandles_.size())
        changed = true;
    else {
        for (size_t i=0; i<tempVolumes.size() && !changed; i++) {
            if (tempVolumes.at(i) != volumeHandles_.at(i))
                changed = true;
            else {
                tgtAssert(tempHashMap_.find(tempVolumes.at(i)) != tempHashMap_.end(), "missing entry in tempHashMap");
                tgtAssert(handleToHashMap_.find(volumeHandles_.at(i)) != handleToHashMap_.end(), "missing entry in handleToHashMap_");
                if (tempHashMap_[tempVolumes.at(i)] != handleToHashMap_[volumeHandles_.at(i)])
                    changed = true;
            }
        }
    }

    if (!changed)
        return;

    clear();

    // update widget from new volumes
    volumeHandles_ = tempVolumes;
    handleToHashMap_ = tempHashMap_;
    handleToPortMap_ = tempPortMap_;

    for (size_t i=0; i<volumeHandles_.size(); i++) {
        const VolumeBase* handle = volumeHandles_.at(i);
        if (!handle->Observable<VolumeObserver>::isObservedBy(static_cast<VolumeObserver*>(this)))
            handle->Observable<VolumeObserver>::addObserver(static_cast<VolumeObserver*>(this));

        tgtAssert(handleToPortMap_.find(handle) != handleToPortMap_.end(),
            "missing entry in handleToPortMap");
        QTreeWidgetItem* treeWidgetItem = createTreeWidgetItem(handle, handleToPortMap_[handle]);
        volumeInfos_->addTopLevelItem(treeWidgetItem);
    }

    if (!volumeHandles_.empty()) {
        containerInfo_->setText(QString::fromStdString(calculateVolumeSizeString(volumeHandles_)));
        if (parentWidget() && parentWidget()->parentWidget()) {
            size_t numVolumes = volumeHandles_.size();
            parentWidget()->parentWidget()->setWindowTitle(QString::fromStdString("Volume Viewer (" + itos(volumeHandles_.size()) + (numVolumes > 1 ? " Volumes)" : " Volume)")));
        }
    }
}

QTreeWidgetItem* VolumeViewer::createTreeWidgetItem(const VolumeBase* handle, const Port* port) {
    tgtAssert(handle, "no handle");
    tgtAssert(port, "no port");

    QTreeWidgetItem* treeItem = new QTreeWidgetItem(volumeInfos_);
    std::string name = VolumeViewHelper::getStrippedVolumeName(handle);
    std::string path = VolumeViewHelper::getVolumePath(handle);
    /*if (name.empty())
        name = "-";
    if (path.empty())
        path += "-"; */

    //QLabel* infos = new QLabel(QString(name.c_str()) + " (" + QString(VolumeViewHelper::getVolumeType(volume).c_str()) + ")" + QString(QChar::LineSeparator)
    //              + QString(path.c_str()) + QString(QChar::LineSeparator)
    //              + "Dimension: " + QString(VolumeViewHelper::getVolumeDimension(volume).c_str()));
    //infos->setWordWrap(true);

    QFontInfo fontInfo(treeItem->font(0));
    treeItem->setFont(0, QFont(fontInfo.family(), fontSize));

    // description lines
    std::vector<std::string> descLines;
    descLines.push_back(port->getProcessor()->getGuiName() + "." + port->getGuiName());
    descLines.push_back("File: " + VolumeViewHelper::getStrippedVolumeName(handle));
    descLines.push_back("Data Type: " + VolumeViewHelper::getVolumeType(handle));
    descLines.push_back("Dimension: " + VolumeViewHelper::getVolumeDimension(handle));

    std::vector<std::string> representations;
    if (handle->hasRepresentation<VolumeDisk>())
        representations.push_back("Disk");
    if (handle->hasRepresentation<VolumeRAM>())
        representations.push_back("RAM");
    if (handle->hasRepresentation<VolumeOctreeBase>())
        representations.push_back("Octree");
    if (handle->hasRepresentation<VolumeGL>())
        representations.push_back("GL");
    descLines.push_back("MemSize: " + VolumeViewHelper::getVolumeMemorySize(handle) + " (" + strJoin(representations, ",") + ")");

    // format description lines
    QFontMetrics fontMetrics(QFontMetrics(treeItem->font(0)));
    int maxTextWidth = 175;
    QString descString;
    for (size_t i=0; i<descLines.size(); i++) {
        descString += fontMetrics.elidedText(QString::fromStdString(descLines.at(i)), Qt::ElideMiddle, maxTextWidth);
        if (i < descLines.size()-1)
            descString.append(QString(QChar::LineSeparator));
    }
    treeItem->setText(0, descString);

    // icon
    treeItem->setIcon(0, QIcon(VolumeViewHelper::generateBorderedPreview(handle, 63, 0)));
    treeItem->setSizeHint(0,QSize(65,65));

    // tooltip
    treeItem->setToolTip(0, QString::fromStdString("Port: " + port->getProcessor()->getGuiName() + "." + port->getGuiName()) + "\n"
        + QString::fromStdString("File: " + name) + "\n"
        + QString::fromStdString("Path: " + path) + "\n"
        + QString::fromStdString("Data Type: " + VolumeViewHelper::getVolumeType(handle)) + "\n"
        + QString::fromStdString("Dimensions: " + VolumeViewHelper::getVolumeDimension(handle)) + "\n"
        + QString::fromStdString("Voxel Size: " + VolumeViewHelper::getVolumeSpacing(handle)) + "\n"
        + QString::fromStdString("Memory Size: " + VolumeViewHelper::getVolumeMemorySize(handle)) + "\n"
        + QString::fromStdString("Representations: " + strJoin(representations, ", "))
    );

    return treeItem;
}

std::string VolumeViewer::calculateVolumeSizeString(const std::vector<const VolumeBase*>& handles) {

    uint64_t volumeSize = 0;
    uint64_t volumeSizeRam = 0;

    // use set to eliminate duplicates
    std::set<const VolumeBase*> handleSet;
    handleSet.insert(handles.begin(), handles.end());

    for (std::set<const VolumeBase*>::iterator it = handleSet.begin(); it != handleSet.end(); ++it) {
        volumeSize += VolumeViewHelper::getVolumeMemorySizeByte(*it);
        if ((*it)->hasRepresentation<VolumeRAM>())
            volumeSizeRam += VolumeViewHelper::getVolumeMemorySizeByte(*it);
    }
    int volumeCount = static_cast<int>(handleSet.size());

    std::stringstream out;
    /*if(volumeCount == 1) {
        out << volumeCount << " Volume (";
    }
    else if (volumeCount > 1) {
        out << volumeCount << " Volumes (";
    }*/

    out << "Total: " << formatMemorySize(volumeSize) << " | RAM: " << formatMemorySize(volumeSizeRam);

    return out.str();
}

void VolumeViewer::exportVolumes() {
    QList<QTreeWidgetItem *> allItems = volumeInfos_->selectedItems();
    QList<QTreeWidgetItem *>::iterator it = allItems.begin();
    while (it != allItems.end()) {
        int itemIndex = volumeInfos_->indexOfTopLevelItem((*it));
        tgtAssert(itemIndex < (int)volumeHandles_.size(), "invalid item index");
        tgtAssert(volumeHandles_.at(itemIndex), "volume null pointer");
        volumeIOHelper_.showFileSaveDialog(volumeHandles_.at(itemIndex));
        it++;
    }
}

void VolumeViewer::updateStateChanged(int /*state*/) {
    if (isVisible() && autoUpdateCheckbox_->isChecked()) {
        updateRequired_ = true;
        update();
    }
}

void VolumeViewer::volumeSaved(const VolumeBase* /*handle*/, const std::string& exportPath) {
    QMessageBox::information(this, "Volume Saved",
        QString::fromStdString("Saved volume to: " + exportPath));
}

} // namespace
