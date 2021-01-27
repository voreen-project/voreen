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

#include "voreen/qt/widgets/property/transfunc/1d/1dkeys/transfunc1dkeyspropertyeditor.h"

#include "voreen/qt/widgets/property/transfunc/utils/doubleslider.h"
#include "voreen/qt/widgets/property/transfunc/1d/1dkeys/transfunc1dkeyspropertyeditorcanvas.h"

#include "voreen/core/datastructures/transfunc/1d/1dkeys/transfunc1dkeys.h"

#include <QGroupBox>
#include <QLayout>
#include <QToolButton>
#include <QSizePolicy>
#include <QCheckBox>

namespace voreen {

const std::string TransFunc1DKeysPropertyEditor::loggerCat_("voreen.qt.TransFunc1DKeysPropertyEditor");

//----------------------------------------------------------------------------------------------
//      Constructor and Qt Stuff
//----------------------------------------------------------------------------------------------
TransFunc1DKeysPropertyEditor::TransFunc1DKeysPropertyEditor(TransFunc1DKeysProperty* prop, QWidget* parent)
    : TransFunc1DPropertyEditor(prop, parent)
    //layout
    , mappingCanvas_(0) , thresholdSlider_(0)
    , showHistogramCB_(0), showTextureCB_(0)
    , makeRampButton_(0), invertButton_(0)
    , computeHistogram_(0)
    //member
    , transFunc1DKeysProp_(prop)
{
}

TransFunc1DKeysPropertyEditor::~TransFunc1DKeysPropertyEditor() {
}

void TransFunc1DKeysPropertyEditor::updateFromProperty() {
    TransFunc1DPropertyEditor::updateFromProperty();
    if(transFunc1DKeysProp_->get()) {
        //update slider
        thresholdSlider_->blockSignals(true);
        thresholdSlider_->setValues(transFunc1DKeysProp_->get()->getThreshold().x, transFunc1DKeysProp_->get()->getThreshold().y);
        thresholdSlider_->blockSignals(false);
        //update mapping canvas
        mappingCanvas_->updateFromProperty();

        //update computeHistogram
        computeHistogram_->blockSignals(true);
        computeHistogram_->setChecked(transFunc1DKeysProp_->getComputeHistogram());
        computeHistogram_->blockSignals(false);
    }
}

//----------------------------------------------------------------------------------------------
//      Layout
//----------------------------------------------------------------------------------------------
QWidget* TransFunc1DKeysPropertyEditor::layoutLeftComponents() {
    //main left and right widget
    QWidget* mainLeft = new QWidget();

    mappingCanvas_ = new TransFunc1DKeysPropertyEditorCanvas(0, transFunc1DKeysProp_);
    mappingCanvas_->setMinimumSize(300,100);
    mappingCanvas_->setSizePolicy(QSizePolicy::MinimumExpanding, QSizePolicy::MinimumExpanding);
    connect(mappingCanvas_, SIGNAL(toggleInteractionModeSignal(bool)), this, SLOT(toggleInteractionMode(bool)));
    connect(mappingCanvas_, SIGNAL(colorChangedSignal(const QColor&)), colorPicker_, SLOT(setHSByColorSlot(const QColor&)));
    connect(mappingCanvas_, SIGNAL(colorChangedSignal(const QColor&)), colorLumPicker_, SLOT(setHSVByColorSlot(const QColor&)));
    connect(colorLumPicker_, SIGNAL(newHSVSignal(int,int,int)), this, SLOT(markerColorChanged(int,int,int)));

    // threshold slider
    QHBoxLayout* hboxSlider = new QHBoxLayout();
    thresholdSlider_ = new DoubleSlider();
    thresholdSlider_->setSizePolicy(QSizePolicy::MinimumExpanding, QSizePolicy::Fixed);
    thresholdSlider_->setOffsets(12, 27);
    hboxSlider->addWidget(thresholdSlider_);
    connect(thresholdSlider_, SIGNAL(valuesChanged(float, float)), this, SLOT(sliderThresholdChanged()));
    connect(thresholdSlider_, SIGNAL(toggleInteractionMode(bool)), this, SLOT(toggleInteractionMode(bool)));

    //add special functions
        //toggle options
    computeHistogram_ = new QCheckBox("Auto-compute Histogram");
    computeHistogram_->setChecked(true);
    connect(computeHistogram_, SIGNAL(stateChanged(int)),this,SLOT(computeHistogramToggled(int)));

    showHistogramCB_ = new QCheckBox("Show Histogram");
    showHistogramCB_->setChecked(true);
    showTextureCB_ = new QCheckBox("Show Background");
    showTextureCB_->setChecked(true);
    connect(showHistogramCB_, SIGNAL(stateChanged(int)),this,SLOT(showHistogramToggled(int)));
    connect(showTextureCB_, SIGNAL(stateChanged(int)),this,SLOT(showTextureToggled(int)));
    QVBoxLayout* checkBoxLayout = new QVBoxLayout();
    checkBoxLayout->setAlignment(Qt::AlignCenter);
    checkBoxLayout->addWidget(computeHistogram_);
    checkBoxLayout->addWidget(showHistogramCB_);
    checkBoxLayout->addWidget(showTextureCB_);
    QGroupBox* viewBox = new QGroupBox("View Settings");
    viewBox->setLayout(checkBoxLayout);
    viewBox->setSizePolicy(QSizePolicy::Fixed,QSizePolicy::Minimum);
        //quick edit
    makeRampButton_ = new QToolButton();
    makeRampButton_->setToolButtonStyle(Qt::ToolButtonTextBesideIcon);
    makeRampButton_->setIcon(QIcon(":/qt/icons/ramp.png"));
    makeRampButton_->setText("Make Ramp");
    connect(makeRampButton_, SIGNAL(clicked()), this, SLOT(makeRampButtonClicked()));
    invertButton_ = new QToolButton();
    invertButton_->setToolButtonStyle(Qt::ToolButtonTextBesideIcon);
    invertButton_->setIcon(QIcon(":/qt/icons/arrow-leftright.png"));
    invertButton_->setText("Invert Map");
    connect(invertButton_, SIGNAL(clicked()), this, SLOT(invertMapButtonClicked()));
    QVBoxLayout* buttonLayout = new QVBoxLayout();
    buttonLayout->setAlignment(Qt::AlignCenter);
    buttonLayout->addWidget(makeRampButton_);
    buttonLayout->addWidget(invertButton_);
    QGroupBox* buttonBox = new QGroupBox("Quick Edit");
    buttonBox->setLayout(buttonLayout);
    buttonBox->setSizePolicy(QSizePolicy::Fixed,QSizePolicy::Minimum);

    QBoxLayout* specialFunctionLayout = new QHBoxLayout();
    specialFunctionLayout->addStretch();
    specialFunctionLayout->addWidget(viewBox);
    specialFunctionLayout->addSpacing(5);
    specialFunctionLayout->addWidget(buttonBox);
    specialFunctionLayout->addStretch();

    // put widgets in layout
    QVBoxLayout* vBox = new QVBoxLayout();
    vBox->setSpacing(1);
    vBox->addWidget(mappingCanvas_,1);
    vBox->addLayout(hboxSlider);
    vBox->addSpacing(10);
    vBox->addLayout(specialFunctionLayout);

    //set left layout
    mainLeft->setLayout(vBox);

    //updateThresholdFromProperty();
    return mainLeft;
}

//----------------------------------------------------------------------------------------------
//      layout slots
//----------------------------------------------------------------------------------------------
void TransFunc1DKeysPropertyEditor::sliderThresholdChanged() {
    tgtAssert(transFunc1DKeysProp_ && transFunc1DKeysProp_->get(), "no property or tf!");
    transFunc1DKeysProp_->get()->setThreshold(thresholdSlider_->getMinValue(), thresholdSlider_->getMaxValue());
    transFunc1DKeysProp_->invalidate();
}

void TransFunc1DKeysPropertyEditor::showHistogramToggled(int state) {
    mappingCanvas_->toggleHistogram(state);
}

void TransFunc1DKeysPropertyEditor::showTextureToggled(int state) {
    mappingCanvas_->toggleTexture(state);
}

void TransFunc1DKeysPropertyEditor::makeRampButtonClicked() {
    if (transFunc1DKeysProp_->get()) {
        transFunc1DKeysProp_->get()->makeRamp();
        transFunc1DKeysProp_->invalidate();
    }
}

void TransFunc1DKeysPropertyEditor::invertMapButtonClicked() {
    if (transFunc1DKeysProp_->get()) {
        transFunc1DKeysProp_->get()->invertKeys();
        transFunc1DKeysProp_->invalidate();
    }
}

void TransFunc1DKeysPropertyEditor::computeHistogramToggled(int /*state*/) {
    transFunc1DKeysProp_->setComputeHistogram(computeHistogram_->isChecked());
}

//----------------------------------------------------------------------------------------------
//      helper functions
//----------------------------------------------------------------------------------------------
void TransFunc1DKeysPropertyEditor::markerColorChanged(int h, int s, int v) {
    mappingCanvas_->changeCurrentColorSlot(QColor::fromHsv(h, s, v));
}



} // namespace voreen
