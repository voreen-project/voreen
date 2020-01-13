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

#include "segmentationplugin.h"

#include "voreen/core/properties/propertyvector.h"
#include "voreen/core/processors/volumeraycaster.h"
#include "voreen/core/utils/voreenpainter.h"

#include "voreen/core/network/networkevaluator.h"
#include "voreen/core/network/processornetwork.h"
#include "voreen/core/datastructures/callback/memberfunctioncallback.h"

#include <QCheckBox>
#include <QComboBox>
#include <QDialog>
#include <QFileDialog>
#include <QFrame>
#include <QHBoxLayout>
#include <QLabel>
#include <QMessageBox>
#include <QProcess>
#include <QPushButton>
#include <QSpinBox>
#include <QTextEdit>
#include <QVBoxLayout>
#include <QGroupBox>
#include <QMessageBox>
#include <QApplication>

namespace voreen {

using tgt::ivec2;

SegmentationPlugin::SegmentationPlugin(QWidget* parent, NetworkEvaluator* evaluator)
    : QWidget(parent),
      evaluator_(evaluator),
      intensityEditor_(0)
{
    setObjectName(tr("Segmentation"));
    setWindowIcon(QIcon(":/qt/icons/segmentation.png"));

    createWidgets();
    createConnections();
}

SegmentationPlugin::~SegmentationPlugin() {
}

bool SegmentationPlugin::usable(const std::vector<Processor*>& processors) {
    for (size_t i=0; i < processors.size(); i++) {
        if (dynamic_cast<SegmentationRaycaster*>(processors[i]))
            return true;
        if (dynamic_cast<RegionGrowingProcessor*>(processors[i]))
            return true;
    }
    return false;
}

void SegmentationPlugin::createWidgets() {
    QVBoxLayout* mainLayout = new QVBoxLayout();

    //
    // Visualization Box
    //
  /*  renderingBox_ = new QGroupBox(tr("Rendering"));

    QVBoxLayout* vboxLayout = new QVBoxLayout();
    QHBoxLayout* hboxLayout = new QHBoxLayout;

    // apply segmentation
    checkApplySegmentation_ = new QCheckBox(tr("Apply Segmentation"));
    hboxLayout->addWidget(checkApplySegmentation_);
    hboxLayout->addSpacing(20);

    // current segment
    spinCurrentSegment_ = new QSpinBox;
    spinCurrentSegment_->setRange(0, 1024);
    hboxLayout->addWidget(labelCurrentSegment_ = new QLabel(tr("Current Segment: ")));
    hboxLayout->addWidget(spinCurrentSegment_);
    hboxLayout->addSpacing(20);
    hboxLayout->addStretch();

    vboxLayout->addItem(hboxLayout);

    QFrame* line = new QFrame;
    line->setFrameShape(QFrame::HLine);
    line->setFrameShadow(QFrame::Sunken);
    vboxLayout->addWidget(line);

    // tf editor
    intensityEditor_ = new TransFunc1DKeysEditor(new TransFuncProperty("dummy","dummy"));
    intensityEditor_->createWidgets();
    intensityEditor_->createConnections();
    vboxLayout->addWidget(intensityEditor_);

    renderingBox_->setLayout(vboxLayout);
    mainLayout->addWidget(renderingBox_);

    //
    // Region Growing Box
    //
    regionGrowingBox_ = new QGroupBox(tr("Region Growing"));
    vboxLayout = new QVBoxLayout;
    hboxLayout = new QHBoxLayout;

    // mark seed button, output segment
    hboxLayout->addWidget(new QLabel(tr("Output Segment: ")));
    outputSegment_ = new QSpinBox;
    outputSegment_->setRange(0, 255);
    outputSegment_->setValue(1);
    hboxLayout->addWidget(outputSegment_);
    hboxLayout->addSpacing(10);
    markSeedButton_ = new QPushButton(tr(" Mark Seed "));
    markSeedButton_->setCheckable(true);
    hboxLayout->addWidget(markSeedButton_);
    hboxLayout->addSpacing(10);
    undoButton_ = new QPushButton(tr(" Undo Last "));
    hboxLayout->addWidget(undoButton_);
    hboxLayout->addSpacing(10);
    clearSegmentButton_ = new QPushButton(tr(" Clear Segment "));
    hboxLayout->addWidget(clearSegmentButton_);
    hboxLayout->addStretch();
    vboxLayout->addLayout(hboxLayout);
    line = new QFrame;
    line->setFrameShape(QFrame::HLine);
    line->setFrameShadow(QFrame::Sunken);
    vboxLayout->addWidget(line);

    hboxLayout = new QHBoxLayout;
    hboxLayout->addWidget(new QLabel("Strictness: "));
    spinStrictness_ = new QDoubleSpinBox();
    spinStrictness_->setValue(0.8);
    spinStrictness_->setRange(0.0, 65000.0);
    spinStrictness_->setSingleStep(0.05);
    hboxLayout->addWidget(spinStrictness_);
    hboxLayout->addSpacing(10);
    spinMaxSeedDist_ = new QSpinBox();
    spinMaxSeedDist_->setValue(0);
    spinMaxSeedDist_->setRange(0,999);
    hboxLayout->addWidget(new QLabel(tr("Max distance to seed:" )));
    hboxLayout->addWidget(spinMaxSeedDist_);
    hboxLayout->addSpacing(10);
    checkApplyThresholds_ = new QCheckBox(tr("Apply Thresholds"));
    checkApplyThresholds_->setChecked(true);
    hboxLayout->addWidget(checkApplyThresholds_);
    hboxLayout->addStretch();
    vboxLayout->addLayout(hboxLayout);

    hboxLayout = new QHBoxLayout;
    comboCostFunction_ = new QComboBox();
    comboCostFunction_->addItem(tr("Intensity"));
    comboCostFunction_->addItem(tr("Gradient Magnitude"));
    comboCostFunction_->addItem(tr("Weighted"));
    hboxLayout->addWidget(new QLabel(tr("Cost Function: ")));
    hboxLayout->addWidget(comboCostFunction_);
    hboxLayout->addSpacing(15);
    checkAdaptive_ = new QCheckBox(tr("Adaptive Growing Criteria"));
    hboxLayout->addWidget(checkAdaptive_);
    hboxLayout->addStretch();
    vboxLayout->addLayout(hboxLayout);


    // clear / save buttons
    line = new QFrame;
    line->setFrameShape(QFrame::HLine);
    line->setFrameShadow(QFrame::Sunken);
    vboxLayout->addWidget(line);
    hboxLayout = new QHBoxLayout;

    clearSegmentationButton_ = new QPushButton(tr(" Clear Segmentation "));
    hboxLayout->addWidget(clearSegmentationButton_);
    hboxLayout->addSpacing(5);
    saveSegmentationButton_ = new QPushButton(tr(" Save Segmentation "));
    hboxLayout->addWidget(saveSegmentationButton_);
    hboxLayout->addSpacing(5);
    loadSegmentationButton_ = new QPushButton(tr(" Load Segmentation "));
    hboxLayout->addWidget(loadSegmentationButton_);
    hboxLayout->addStretch();
    vboxLayout->addItem(hboxLayout);

    regionGrowingBox_->setLayout(vboxLayout);
    mainLayout->addWidget(regionGrowingBox_);

    // finish layout
    mainLayout->addStretch();
    setLayout(mainLayout);*/
}

void SegmentationPlugin::createConnections() {
    // connections
    connect(checkApplySegmentation_, SIGNAL(toggled(bool)), this, SLOT(toggleApplySegmentation(bool)));
    connect(spinCurrentSegment_, SIGNAL(valueChanged(int)), this, SLOT(setCurrentSegment(int)));

    connect(clearSegmentButton_, SIGNAL(clicked()), this, SLOT(clearSegment()));
    connect(undoButton_, SIGNAL(clicked()), this, SLOT(undoSegment()));
    connect(markSeedButton_, SIGNAL(clicked(bool)), this, SLOT(setSeed(bool)));

    connect(saveSegmentationButton_, SIGNAL(clicked()), this, SLOT(saveSegmentation()));
    connect(loadSegmentationButton_, SIGNAL(clicked()), this, SLOT(loadSegmentation()));
    connect(clearSegmentationButton_, SIGNAL(clicked()), this, SLOT(clearSegmentation()));

    connect(intensityEditor_, SIGNAL(transferFunctionChanged()), this, SLOT(repaintCanvas()));

}

void SegmentationPlugin::toggleApplySegmentation(bool useSegmentation) {

    if (!getSegmentationRaycaster())
        return;

    getSegmentationRaycaster()->getApplySegmentationProp().set(useSegmentation);

    repaintCanvas();
}

void SegmentationPlugin::setCurrentSegment(int /*segment*/) {
    updateIntensityEditor();
}

void SegmentationPlugin::undoSegment() {
    if (getRegionGrowingProcessor()) {
        getRegionGrowingProcessor()->undoLastGrowing();
        repaintCanvas();
    }
}

void SegmentationPlugin::clearSegment() {
    if (getRegionGrowingProcessor()) {
        getRegionGrowingProcessor()->clearSegment(outputSegment_->value());
        repaintCanvas();
    }
}

void SegmentationPlugin::setSeed(bool checked) {
    if (!getRegionGrowingProcessor()) {
        if (checked)
            markSeedButton_->setChecked(false);
        return;
    }
}

void SegmentationPlugin::mousePressEvent(tgt::MouseEvent* e) {
    e->ignore();

    if (e->button() != tgt::MouseEvent::MOUSE_BUTTON_LEFT ||
        !markSeedButton_->isChecked() ||
        !getRegionGrowingProcessor()) {
        return;
    }

    if (!evaluator_) {
        LERRORC("voreen.qt.SegmentationPlugin", "No network evaluator");
        return;
    }

    //if (!evaluator_->getTextureContainer()) {
        LERRORC("voreen.qt.SegmentationPlugin", "No texture container");
        return;
    //}

    e->accept();
    QApplication::setOverrideCursor(QCursor(Qt::WaitCursor));

    // propagate local settings to region growing processor
    propagateRegionGrowingParams();

    // propagate picking coords to RegionGrowing processor and start growing
    //ivec2 tcDims = evaluator_->getTextureContainer()->getSize();
    //getRegionGrowingProcessor()->startGrowing(ivec2(e->coord().x, tcDims.y - e->coord().y), outputSegment_->value());
    QApplication::restoreOverrideCursor();

    repaintCanvas();
}

void SegmentationPlugin::saveSegmentation() {
    if (!getRegionGrowingProcessor())
        return;

    QString fileName = QFileDialog::getSaveFileName(this, tr("Save Segmentation"),
                                                    QDir::currentPath() + "/../../data/segmentation.dat",
                                                    tr("Volumes (*.dat)"));
    if (!fileName.isEmpty())
        getRegionGrowingProcessor()->saveSegmentation(fileName.toStdString());
}

void SegmentationPlugin::loadSegmentation() {
    if (!getRegionGrowingProcessor())
        return;

    QString fileName = QFileDialog::getOpenFileName(this, tr("Load Segmentation"),
        QDir::currentPath() + "/../../data/segmentation.dat",
        tr("Volumes (*.dat)"));
    if (!fileName.isEmpty()) {
        getRegionGrowingProcessor()->loadSegmentation(fileName.toStdString());
        repaintCanvas();
    }
}


SegmentationRaycaster* SegmentationPlugin::getSegmentationRaycaster(bool suppressWarning) {

    if (!evaluator_) {
        LERRORC("voreen.qt.SegmentationPlugin", "No network evaluator.");
        return 0;
    }

    if (!evaluator_->getProcessorNetwork())
        return 0;

    SegmentationRaycaster* segproc = 0;
    const std::vector<Processor*>& procs = evaluator_->getProcessorNetwork()->getProcessors();
    for (size_t i=0; i < procs.size(); i++) {
        if ((segproc = dynamic_cast<SegmentationRaycaster*>(procs[i])) != 0)
            break;
    }

    if (!segproc && !suppressWarning)
        QMessageBox::warning(this, tr("No SegmentationRaycaster"), tr("No SegmentationRaycaster found in the current network."));

    return segproc;
}

RegionGrowingProcessor* SegmentationPlugin::getRegionGrowingProcessor(bool suppressWarning) {

    if (!evaluator_) {
        LERRORC("voreen.qt.SegmentationPlugin", "No network evaluator.");
        return 0;
    }

    if (!evaluator_->getProcessorNetwork())
        return 0;

    RegionGrowingProcessor* regionGrowing = 0;
    const std::vector<Processor*>& procs = evaluator_->getProcessorNetwork()->getProcessors();
    for (size_t i=0; i < procs.size(); i++) {
        if ((regionGrowing = dynamic_cast<RegionGrowingProcessor*>(procs[i])) != 0)
            break;
    }

    if (!regionGrowing && !suppressWarning)
        QMessageBox::warning(this, tr("No RegionGrowing processor"),
            tr("No RegionGrowing processor found in the current network."));

    return regionGrowing;
}

void SegmentationPlugin::transFuncChangedExternally() {
    intensityEditor_->updateFromProperty();
}

void SegmentationPlugin::clearSegmentation() {
    if (getRegionGrowingProcessor()) {
        getRegionGrowingProcessor()->clearSegmentation();
        repaintCanvas();
    }
}

void SegmentationPlugin::propagateRegionGrowingParams() {
    if (!getRegionGrowingProcessor())
        return;

    getRegionGrowingProcessor()->getStrictnessProp().set(spinStrictness_->value());
//    getRegionGrowingProcessor()->getCostFunctionProp().set(comboCostFunction_->currentIndex());
    getRegionGrowingProcessor()->getMaxSeedDistanceProp().set(spinMaxSeedDist_->value());
    getRegionGrowingProcessor()->getAdaptiveProp().set(checkAdaptive_->isChecked());
    if (checkApplyThresholds_->isChecked()) {
/*        TransFuncIntensity* tfIntensity = dynamic_cast<TransFuncIntensity*>(transFuncProp_.get());
        if (tfIntensity)
            getRegionGrowingProcessor()->getThresholdProp().set(tfIntensity->getThresholds());*/
    }
    getRegionGrowingProcessor()->getThresholdFillingProp().set(checkApplyThresholds_->isChecked());

}

void SegmentationPlugin::applySegmentationToggled() {
    if (!getSegmentationRaycaster())
        return;

    bool applySegmentation = getSegmentationRaycaster()->getApplySegmentationProp().get();

    checkApplySegmentation_->setChecked(applySegmentation);
    labelCurrentSegment_->setEnabled(applySegmentation);
    spinCurrentSegment_->setEnabled(applySegmentation);

    updateIntensityEditor();

}

void SegmentationPlugin::updateIntensityEditor() {
 /*   if (!getSegmentationRaycaster())
        return;

    bool useSegmentation = getSegmentationRaycaster()->getApplySegmentationProp().get();

    if (useSegmentation) {
        const PropertyVector& segmentationTransFuncs = getSegmentationRaycaster()->getSegmentationTransFuncs();
        int currentSegment = spinCurrentSegment_->value();
        if (currentSegment >= 0 && currentSegment < segmentationTransFuncs.size()) {
            TransFuncProperty* segTransFuncProp = segmentationTransFuncs.getProperty<TransFuncProperty*>(currentSegment);
            if (segTransFuncProp)
                intensityEditor_->setTransFuncProp(segTransFuncProp);
        }
    }
    else {
        intensityEditor_->setTransFuncProp(&getSegmentationRaycaster()->getTransFunc());
    }*/
}

void SegmentationPlugin::registerAsListener() {
    SegmentationRaycaster* segmentationRaycaster = getSegmentationRaycaster();
    if (!segmentationRaycaster)
        return;

    // register as change listener of segmentation raycaster's transfunc prop
    segmentationRaycaster->getTransFunc().onChange(
        MemberFunctionCallback<SegmentationPlugin>(this, &SegmentationPlugin::transFuncChangedExternally));

    // register as change listener of segmentation raycaster's segmentation transfunc props
/*    const PropertyVector& segTransFuncs = segmentationRaycaster->getSegmentationTransFuncs();
    for (int i=0; i<segTransFuncs.size(); ++i) {
        TransFuncProperty* tfProp = segTransFuncs.getProperty<TransFuncProperty*>(i);
        if (tfProp)
            tfProp->onChange(MemberFunctionCallback<SegmentationPlugin>(this,
                &SegmentationPlugin::transFuncChangedExternally));
    }

    // register as listener of the segmentation raycaster's apply segmentation prop
    segmentationRaycaster->getApplySegmentationProp().onChange(
        MemberFunctionCallback<SegmentationPlugin>(this, &SegmentationPlugin::applySegmentationToggled));*/
}

void SegmentationPlugin::repaintCanvas() {
    LWARNINGC("SegmentationPlugin", "SegmentationPlugin::repaintCanvas() called (unimplemented)");
}

} // namespace voreen
