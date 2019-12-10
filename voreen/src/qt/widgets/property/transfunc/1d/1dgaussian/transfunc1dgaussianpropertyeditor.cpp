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

#include "voreen/qt/widgets/property/transfunc/1d/1dgaussian/transfunc1dgaussianpropertyeditor.h"

#include "voreen/qt/widgets/property/transfunc/utils/doubleslider.h"
#include "voreen/qt/widgets/property/transfunc/1d/1dgaussian/transfunc1dgaussianpropertyeditorcanvas.h"

#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/datastructures/volume/histogram.h"
#include "voreen/core/datastructures/transfunc/1d/1dgaussian/transfunc1dgaussian.h"

#include <QGroupBox>
#include <QLayout>
#include <QToolButton>
#include <QSizePolicy>
#include <QCheckBox>
#include <qtooltip.h>

namespace voreen {

const std::string TransFunc1DGaussianPropertyEditor::loggerCat_("voreen.qt.TransFunc1DGaussianPropertyEditor");

//----------------------------------------------------------------------------------------------
//      Constructor and Qt Stuff
//----------------------------------------------------------------------------------------------
TransFunc1DGaussianPropertyEditor::TransFunc1DGaussianPropertyEditor(TransFunc1DGaussianProperty* prop, QWidget* parent)
    : TransFunc1DPropertyEditor(prop, parent)
    //layout
    , mappingCanvas_(0) , thresholdSlider_(0)
    , showHistogramCB_(0), showTextureCB_(0)
    , autoTransFuncButton(0)
    , computeHistogram_(0)
    //member
    , transFunc1DGaussianProp_(prop)
{
}

TransFunc1DGaussianPropertyEditor::~TransFunc1DGaussianPropertyEditor() {
}

void TransFunc1DGaussianPropertyEditor::updateFromProperty() {
    TransFunc1DPropertyEditor::updateFromProperty();
    if(transFunc1DGaussianProp_->get()) {
        //update slider
        thresholdSlider_->blockSignals(true);
        thresholdSlider_->setValues(transFunc1DGaussianProp_->get()->getThreshold().x, transFunc1DGaussianProp_->get()->getThreshold().y);
        thresholdSlider_->blockSignals(false);
        //update mapping canvas
        mappingCanvas_->updateFromProperty();

        //update computeHistogram
        computeHistogram_->blockSignals(true);
        computeHistogram_->setChecked(transFunc1DGaussianProp_->getComputeHistogram());
        computeHistogram_->blockSignals(false);
    }
}

//----------------------------------------------------------------------------------------------
//      Layout
//----------------------------------------------------------------------------------------------
QWidget* TransFunc1DGaussianPropertyEditor::layoutLeftComponents() {
    //main left and right widget
    QWidget* mainLeft = new QWidget();

    mappingCanvas_ = new TransFunc1DGaussianPropertyEditorCanvas(0, transFunc1DGaussianProp_);
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
    autoTransFuncButton = new QToolButton();
    autoTransFuncButton->setToolButtonStyle(Qt::ToolButtonTextBesideIcon);
    autoTransFuncButton->setIcon(QIcon(":/qt/icons/tf_treshold.png"));
    autoTransFuncButton->setText("Auto Compute");
    connect(autoTransFuncButton, SIGNAL(clicked()), this, SLOT(autoTransFuncClicked()));
    QVBoxLayout* buttonLayout = new QVBoxLayout();
    buttonLayout->setAlignment(Qt::AlignCenter);
    buttonLayout->addWidget(autoTransFuncButton);
    QGroupBox* buttonBox = new QGroupBox("Quick Edit");
    buttonBox->setLayout(buttonLayout);
    buttonBox->setSizePolicy(QSizePolicy::Fixed,QSizePolicy::Minimum);

    QBoxLayout* specialFunctionLayout = new QHBoxLayout();
    specialFunctionLayout->addStretch();
    specialFunctionLayout->addWidget(viewBox);
    specialFunctionLayout->addSpacing(5);
    specialFunctionLayout->addWidget(buttonBox);
    specialFunctionLayout->addStretch();

    // put widgets in vertical layout on the left side of the window
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

QWidget* TransFunc1DGaussianPropertyEditor::createTooltip() {
    tooltipButton_ = new QToolButton();
    tooltipButton_->setText(tr(" ? "));
    connect(tooltipButton_, SIGNAL(clicked()), this, SLOT(showTooltip()));
    return tooltipButton_;
}


//----------------------------------------------------------------------------------------------
//      layout slots
//----------------------------------------------------------------------------------------------
void TransFunc1DGaussianPropertyEditor::sliderThresholdChanged() {
    tgtAssert(transFunc1DGaussianProp_ && transFunc1DGaussianProp_->get(), "no property or tf!");
    transFunc1DGaussianProp_->get()->setThreshold(thresholdSlider_->getMinValue(), thresholdSlider_->getMaxValue());
    transFunc1DGaussianProp_->invalidate();
}

void TransFunc1DGaussianPropertyEditor::showHistogramToggled(int state) {
    mappingCanvas_->toggleHistogram(state);
}

void TransFunc1DGaussianPropertyEditor::showTextureToggled(int state) {
    mappingCanvas_->toggleTexture(state);
}

void TransFunc1DGaussianPropertyEditor::autoTransFuncClicked() {
    if (transFunc1DGaussianProp_->get()) {
        if (const VolumeBase* volume = transFunc1DGaussianProp_->getVolume()) {
            if (volume->hasDerivedData<VolumeHistogramIntensity>()) {

                // automatic tresholding based on the volume histogram
                Histogram1D* histogram = &volume->getDerivedData<VolumeHistogramIntensity>()->getHistogram(transFunc1DGaussianProp_->getVolumeChannel());
                transFunc1DGaussianProp_->get()->autoTreshold(histogram);
                transFunc1DGaussianProp_->invalidate();
            }
            else {
                LERROR("volume histogram missing");
            }
        }

   }
}

void TransFunc1DGaussianPropertyEditor::computeHistogramToggled(int /*state*/) {
    transFunc1DGaussianProp_->setComputeHistogram(computeHistogram_->isChecked());
}

//----------------------------------------------------------------------------------------------
//      helper functions
//----------------------------------------------------------------------------------------------
void TransFunc1DGaussianPropertyEditor::markerColorChanged(int h, int s, int v) {
    mappingCanvas_->changeCurrentColorSlot(QColor::fromHsv(h, s, v));
}

void TransFunc1DGaussianPropertyEditor::showTooltip() {
    QToolTip::showText(QCursor::pos(), tr(
        "This transfer function is defined by multiple gauss functions.\n\n"
        "Insert new curves by clicking in a free space in the canvas.\n"
        "Dragging the curve's base markers changes its width and opacity base value.\n"
        "Curves are unicolored, i.e. use the same color for all three keys, by default.\n"
        "This behaviour can be changed in the context menu of a curve.\n"
        "The colors of superimposed curves are mixed additively.\n\n"
        "Hold [shift] to lock horizontal movement, [ctrl] to lock vertical movement.\n\n"
        "With the [Auto Compute] button a simple automatic classification can be performed.\n"
        "The currently visible histogram is used for this calculation. This means that the\n"
        "domain settings affect the outcome of the classification and can be used to\n"
        "exclude boundary intervals of the intensity values."
        ), tooltipButton_);
}


} // namespace voreen
