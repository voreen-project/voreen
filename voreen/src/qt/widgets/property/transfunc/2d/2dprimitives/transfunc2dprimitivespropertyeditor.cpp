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

#include "voreen/qt/widgets/property/transfunc/2d/2dprimitives/transfunc2dprimitivespropertyeditor.h"

#include "voreen/core/properties/transfunc/2d/2dprimitives/transfunc2dprimitivesproperty.h"
#include "voreen/core/datastructures/volume/histogram.h"

#include "voreen/qt/widgets/property/transfunc/2d/2dprimitives/transfunc2dprimitivespropertyeditorcanvas.h"


#include <QBoxLayout>
#include <QCheckBox>

namespace voreen {

const std::string TransFunc2DPrimitivesPropertyEditor::loggerCat_("TransFunc2DPrimitivesPropertyEditor");

using tgt::vec2;

TransFunc2DPrimitivesPropertyEditor::TransFunc2DPrimitivesPropertyEditor(TransFunc2DPrimitivesProperty* prop, QWidget* parent)
    : TransFunc2DPropertyEditor(prop, parent)
    , alphaPicker_(0), mappingCanvas_(0)
    , thresholdSliderX_(0), thresholdSliderY_(0)
    , showHistogramCB_(0), histogramBrightnessSlider_(0), showTextureCB_(0)
    , computeHistogram_(0)
    , addTriangleButton_(0), addQuadButton_(0),addBananaButton_(0)
    , primitiveFuzzinessSlider_(0), removePrimitiveButton_(0)
    , transFunc2DPrimitivesProp_(prop)
{
}

TransFunc2DPrimitivesPropertyEditor::~TransFunc2DPrimitivesPropertyEditor() {
}

void TransFunc2DPrimitivesPropertyEditor::updateFromProperty() {
    TransFunc2DPropertyEditor::updateFromProperty();
    if(transFunc2DPrimitivesProp_->get()) {
        //update slider
        thresholdSliderX_->blockSignals(true);
        thresholdSliderX_->setValues(transFunc2DPrimitivesProp_->get()->getThreshold(0).x, transFunc2DPrimitivesProp_->get()->getThreshold(0).y);
        thresholdSliderX_->blockSignals(false);
        thresholdSliderY_->blockSignals(true);
        thresholdSliderY_->setValues(1.f - transFunc2DPrimitivesProp_->get()->getThreshold(1).y, 1.f - transFunc2DPrimitivesProp_->get()->getThreshold(1).x);
        thresholdSliderY_->blockSignals(false);
        //update mapping canvas
        mappingCanvas_->updateFromProperty();

        //update computeHistogram
        computeHistogram_->blockSignals(true);
        computeHistogram_->setChecked(transFunc2DPrimitivesProp_->getComputeHistogram());
        computeHistogram_->blockSignals(false);
    }
}

//----------------------------------------------------------------------------------------------
//      update functions
//----------------------------------------------------------------------------------------------
void TransFunc2DPrimitivesPropertyEditor::updateVolumeDataBoundsFromProperty() {
    tgtAssert(transFunc2DProp_, "no Property");
    if(const VolumeBase* volume = transFunc2DProp_->getVolume()) {
        //use volumegradiant
        if(volume->hasDerivedData<VolumeHistogramIntensityGradient>()) {
            lowerVolumeBoundLabelX_->setText(QString::number(volume->getDerivedData<VolumeHistogramIntensityGradient>()->getMinValue(0,transFunc2DProp_->getVolumeChannel())));
            upperVolumeBoundLabelX_->setText(QString::number(volume->getDerivedData<VolumeHistogramIntensityGradient>()->getMaxValue(0,transFunc2DProp_->getVolumeChannel())));
            lowerVolumeBoundLabelY_->setText(QString::number(volume->getDerivedData<VolumeHistogramIntensityGradient>()->getMinValue(1,transFunc2DProp_->getVolumeChannel())));
            upperVolumeBoundLabelY_->setText(QString::number(volume->getDerivedData<VolumeHistogramIntensityGradient>()->getMaxValue(1,transFunc2DProp_->getVolumeChannel())));
        }
        else {
            if (transFunc2DProp_->getComputeHistogram())
                volume->getDerivedDataThreaded<VolumeHistogramIntensityGradient>();
            lowerVolumeBoundLabelX_->setText("0"); upperVolumeBoundLabelX_->setText("1");
            lowerVolumeBoundLabelY_->setText("0"); upperVolumeBoundLabelY_->setText("1");
            return;
        }
    } else {
        //set default, if no volume is present
        lowerVolumeBoundLabelX_->setText("0"); upperVolumeBoundLabelX_->setText("1");
        lowerVolumeBoundLabelY_->setText("0"); upperVolumeBoundLabelY_->setText("1");
        return;
    }
}

//----------------------------------------------------------------------------------------------
//      Layout
//----------------------------------------------------------------------------------------------
QGroupBox* TransFunc2DPrimitivesPropertyEditor::createColorPickerBox() {
    QGroupBox* box = TransFunc2DPropertyEditor::createColorPickerBox();
    alphaPicker_ = new AlphaPicker();
    alphaPicker_->setFixedWidth(20);
    alphaPicker_->setFixedHeight(100);
    connect(alphaPicker_, SIGNAL(toggleInteractionModeSignal(bool)), this, SLOT(toggleInteractionMode(bool)));
    connect(colorLumPicker_, SIGNAL(newHSVSignal(int,int,int)), alphaPicker_, SLOT(updateHSVSlot(int,int,int)));
    box->layout()->addWidget(alphaPicker_);
    return box;
}

QWidget* TransFunc2DPrimitivesPropertyEditor::layoutLeftComponents() {
    //main left and right widget
    QWidget* mainLeft = new QWidget();

    mappingCanvas_ = new TransFunc2DPrimitivesPropertyEditorCanvas(this, transFunc2DPrimitivesProp_);
    mappingCanvas_->setMinimumSize(300,100);
    mappingCanvas_->setSizePolicy(QSizePolicy::MinimumExpanding, QSizePolicy::MinimumExpanding);
    connect(mappingCanvas_, SIGNAL(toggleInteractionModeSignal(bool)), this, SLOT(toggleInteractionMode(bool)));
    connect(mappingCanvas_, SIGNAL(colorChangedSignal(const QColor&)), colorPicker_, SLOT(setHSByColorSlot(const QColor&)));
    connect(mappingCanvas_, SIGNAL(colorChangedSignal(const QColor&)), colorLumPicker_, SLOT(setHSVByColorSlot(const QColor&)));
    connect(mappingCanvas_, SIGNAL(colorChangedSignal(const QColor&)), alphaPicker_, SLOT(setHSVAByColorSlot(const QColor&)));
    connect(alphaPicker_, SIGNAL(newHSVASignal(int,int,int,int)), this, SLOT(markerColorChanged(int,int,int,int)));
    connect(this,SIGNAL(newPrimitiveAddedSignal(TransFuncPrimitive*)),mappingCanvas_,SLOT(newPrimitiveAddedSlot(TransFuncPrimitive*)));

    // threshold slider
    //QHBoxLayout* hboxSlider = new QHBoxLayout();
    thresholdSliderX_ = new DoubleSlider();
    thresholdSliderX_->setSizePolicy(QSizePolicy::MinimumExpanding, QSizePolicy::Fixed);
    thresholdSliderX_->setOffsets(12, 27);
    thresholdSliderY_ = new DoubleSlider(0,true);
    thresholdSliderY_->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::MinimumExpanding);
    thresholdSliderY_->setOffsets(27, 12);
    //hboxSlider->addWidget(thresholdSliderX_);
    connect(thresholdSliderX_, SIGNAL(valuesChanged(float, float)), this, SLOT(sliderThresholdChanged()));
    connect(thresholdSliderX_, SIGNAL(toggleInteractionMode(bool)), this, SLOT(toggleInteractionMode(bool)));
    connect(thresholdSliderY_, SIGNAL(valuesChanged(float, float)), this, SLOT(sliderThresholdChanged()));
    connect(thresholdSliderY_, SIGNAL(toggleInteractionMode(bool)), this, SLOT(toggleInteractionMode(bool)));

    //add special functions
        //toggle options
    computeHistogram_ = new QCheckBox("Auto-compute Histogram");
    computeHistogram_->setChecked(true);
    connect(computeHistogram_, SIGNAL(stateChanged(int)),this,SLOT(computeHistogramToggled(int)));

    showHistogramCB_ = new QCheckBox("Show Histogram");
    showHistogramCB_->setChecked(true);
    histogramBrightnessSlider_ = new QSlider(Qt::Horizontal);
    histogramBrightnessSlider_->setMinimum(10);
    histogramBrightnessSlider_->setMaximum(50);
    histogramBrightnessSlider_->setValue(20);
    showTextureCB_ = new QCheckBox("Show Background");
    showTextureCB_->setChecked(true);
    connect(showHistogramCB_, SIGNAL(stateChanged(int)),this,SLOT(showHistogramToggled(int)));
    connect(histogramBrightnessSlider_,SIGNAL(sliderMoved(int)),mappingCanvas_,SLOT(setHistogramBrightness(int)));
    connect(showTextureCB_, SIGNAL(stateChanged(int)),this,SLOT(showTextureToggled(int)));
    QVBoxLayout* checkBoxLayout = new QVBoxLayout();
    checkBoxLayout->setAlignment(Qt::AlignCenter);
    checkBoxLayout->addWidget(computeHistogram_);
    checkBoxLayout->addWidget(showHistogramCB_);
    QHBoxLayout* brightnessLayout = new QHBoxLayout();
    brightnessLayout->setContentsMargins(20,0,5,0);
    brightnessLayout->addWidget(histogramBrightnessSlider_);
    checkBoxLayout->addItem(brightnessLayout);
    checkBoxLayout->addWidget(showTextureCB_);
    QGroupBox* viewBox = new QGroupBox("View Settings");
    viewBox->setLayout(checkBoxLayout);
    viewBox->setSizePolicy(QSizePolicy::Fixed,QSizePolicy::Minimum);
        //quick edit
    addTriangleButton_ = new QToolButton();
    addTriangleButton_->setToolButtonStyle(Qt::ToolButtonTextBesideIcon);
    addTriangleButton_->setIcon(QIcon(":/qt/icons/triangle.png"));
    addTriangleButton_->setText("Add Triangle");
    addTriangleButton_->setSizePolicy(QSizePolicy::MinimumExpanding,QSizePolicy::Fixed);
    connect(addTriangleButton_, SIGNAL(clicked()), this, SLOT(addTriangleButtonClicked()));
    addQuadButton_ = new QToolButton();
    addQuadButton_->setToolButtonStyle(Qt::ToolButtonTextBesideIcon);
    addQuadButton_->setIcon(QIcon(":/qt/icons/quad.png"));
    addQuadButton_->setText("Add Quad");
    addQuadButton_->setSizePolicy(QSizePolicy::MinimumExpanding,QSizePolicy::Fixed);
    connect(addQuadButton_, SIGNAL(clicked()), this, SLOT(addQuadButtonClicked()));
    addBananaButton_ = new QToolButton();
    addBananaButton_->setToolButtonStyle(Qt::ToolButtonTextBesideIcon);
    addBananaButton_->setIcon(QIcon(":/qt/icons/banana.png"));
    addBananaButton_->setText("Add Banana");
    addBananaButton_->setSizePolicy(QSizePolicy::MinimumExpanding,QSizePolicy::Fixed);
    connect(addBananaButton_, SIGNAL(clicked()), this, SLOT(addBananaButtonClicked()));
    QVBoxLayout* buttonLayout = new QVBoxLayout();
    buttonLayout->setAlignment(Qt::AlignCenter);
    buttonLayout->addWidget(addTriangleButton_);
    buttonLayout->addWidget(addQuadButton_);
    buttonLayout->addWidget(addBananaButton_);
    QGroupBox* buttonBox = new QGroupBox("Add Primitives");
    buttonBox->setLayout(buttonLayout);
    buttonBox->setSizePolicy(QSizePolicy::Fixed,QSizePolicy::Minimum);
        //primitive edit
    QLabel* fuzzinessLabel = new QLabel("Border Fading:");
    primitiveFuzzinessSlider_ = new QSlider(Qt::Horizontal);
    primitiveFuzzinessSlider_->setMinimum(0);
    primitiveFuzzinessSlider_->setMaximum(100);
    primitiveFuzzinessSlider_->setValue(100);
    primitiveFuzzinessSlider_->setEnabled(false);
    connect(primitiveFuzzinessSlider_,SIGNAL(sliderMoved(int)),mappingCanvas_,SLOT(updateFuzzinessSlot(int)));
    connect(mappingCanvas_,SIGNAL(fuzzinessChangedSignal(int)),this,SLOT(fuzzinessChangedSlot(int)));
    removePrimitiveButton_ = new QToolButton();
    removePrimitiveButton_->setToolButtonStyle(Qt::ToolButtonTextBesideIcon);
    removePrimitiveButton_->setIcon(QIcon(":/qt/icons/eraser.png"));
    removePrimitiveButton_->setText("Remove Primitive");
    removePrimitiveButton_->setSizePolicy(QSizePolicy::MinimumExpanding,QSizePolicy::Fixed);
    connect(removePrimitiveButton_, SIGNAL(clicked()), mappingCanvas_, SLOT(deletePrimitiveSlot()));
    QVBoxLayout* primitiveLayout = new QVBoxLayout();
    primitiveLayout->setAlignment(Qt::AlignCenter);
    primitiveLayout->addWidget(fuzzinessLabel);
    QHBoxLayout* fuzzinessLayout = new QHBoxLayout();
    fuzzinessLayout->setContentsMargins(20,0,5,0);
    fuzzinessLayout->addWidget(primitiveFuzzinessSlider_);
    primitiveLayout->addItem(fuzzinessLayout);
    primitiveLayout->addWidget(removePrimitiveButton_);
    QGroupBox* primitiveBox = new QGroupBox("Primitive Settings");
    primitiveBox->setLayout(primitiveLayout);
    primitiveBox->setSizePolicy(QSizePolicy::Fixed,QSizePolicy::Minimum);

    QBoxLayout* specialFunctionLayout = new QHBoxLayout();
    specialFunctionLayout->addStretch();
    specialFunctionLayout->addWidget(viewBox);
    specialFunctionLayout->addSpacing(5);
    specialFunctionLayout->addWidget(buttonBox);
    specialFunctionLayout->addSpacing(5);
    specialFunctionLayout->addWidget(primitiveBox);
    specialFunctionLayout->addStretch();

    // put widgets in layout
    QGridLayout* gridLayout = new QGridLayout();
    gridLayout->setSpacing(1);
    gridLayout->addWidget(thresholdSliderY_,0,0);
    gridLayout->addWidget(mappingCanvas_,0,1);
    gridLayout->addWidget(thresholdSliderX_,1,1);
    gridLayout->addLayout(specialFunctionLayout,2,0,1,2);

    //set left layout
    mainLeft->setLayout(gridLayout);

    //updateThresholdFromProperty();
    return mainLeft;
}

//----------------------------------------------------------------------------------------------
//      layout slots
//----------------------------------------------------------------------------------------------
void TransFunc2DPrimitivesPropertyEditor::sliderThresholdChanged() {
    tgtAssert(transFunc2DPrimitivesProp_ && transFunc2DPrimitivesProp_->get(), "no property or tf!");
    //simply set threshold for intensity
    transFunc2DPrimitivesProp_->get()->setThreshold(thresholdSliderX_->getMinValue(), thresholdSliderX_->getMaxValue(),0);
    //for y take into account the inverse scale
    transFunc2DPrimitivesProp_->get()->setThreshold(1.f - thresholdSliderY_->getMaxValue(), 1.f - thresholdSliderY_->getMinValue(),1);
    transFunc2DPrimitivesProp_->invalidate();
}

void TransFunc2DPrimitivesPropertyEditor::showHistogramToggled(int state) {
    mappingCanvas_->toggleHistogram(state);
    histogramBrightnessSlider_->setEnabled(state);
}

void TransFunc2DPrimitivesPropertyEditor::showTextureToggled(int state) {
    mappingCanvas_->toggleTexture(state);
}

void TransFunc2DPrimitivesPropertyEditor::addTriangleButtonClicked() {
    if (transFunc2DPrimitivesProp_->get()) {
        TransFuncTriangle* triangley = new TransFuncTriangle(); //tgt::col4(128)
        triangley->getControlPoint(0).color_ = tgt::col4(255,255,0,255);
        triangley->getControlPoint(1).color_ = tgt::col4(255,255,0,255);
        triangley->getControlPoint(2).color_ = tgt::col4(255,255,0,255);
        triangley->setFuzziness(0.5f);
        transFunc2DPrimitivesProp_->get()->addPrimitive(triangley);
        emit newPrimitiveAddedSignal(triangley);
        transFunc2DPrimitivesProp_->invalidate();
    }
}

void TransFunc2DPrimitivesPropertyEditor::addQuadButtonClicked() {
    if (transFunc2DPrimitivesProp_->get()) {
        TransFuncQuad* quady = new TransFuncQuad(tgt::vec2(0.5f),0.2f,tgt::col4(128)); //tgt::col4(128)
        quady->getControlPoint(0).color_ = tgt::col4(255,255,0,255);
        quady->getControlPoint(1).color_ = tgt::col4(255,255,0,255);
        quady->getControlPoint(2).color_ = tgt::col4(255,255,0,255);
        quady->getControlPoint(3).color_ = tgt::col4(255,255,0,255);
        quady->setFuzziness(0.5f);
        transFunc2DPrimitivesProp_->get()->addPrimitive(quady);
        emit newPrimitiveAddedSignal(quady);
        transFunc2DPrimitivesProp_->invalidate();
    }
}

void TransFunc2DPrimitivesPropertyEditor::addBananaButtonClicked() {
    if (transFunc2DPrimitivesProp_->get()) {
        TransFuncBanana* bananay = new TransFuncBanana();
        bananay->getControlPoint(0).color_ = tgt::col4(255,255,0,255);
        bananay->getControlPoint(1).color_ = tgt::col4(255,255,0,255);
        bananay->getControlPoint(2).color_ = tgt::col4(255,255,0,255);
        bananay->getControlPoint(3).color_ = tgt::col4(255,255,0,255);
        bananay->setFuzziness(0.5f);
        transFunc2DPrimitivesProp_->get()->addPrimitive(bananay);
        emit newPrimitiveAddedSignal(bananay);
        transFunc2DPrimitivesProp_->invalidate();
    }
}

void TransFunc2DPrimitivesPropertyEditor::fuzzinessChangedSlot(int value) {
    if(value < 0)
        primitiveFuzzinessSlider_->setEnabled(false);
    else {
        primitiveFuzzinessSlider_->setEnabled(true);
        primitiveFuzzinessSlider_->blockSignals(true);
        primitiveFuzzinessSlider_->setValue(value);
        primitiveFuzzinessSlider_->blockSignals(false);
    }
}

void TransFunc2DPrimitivesPropertyEditor::computeHistogramToggled(int /*state*/) {
    transFunc2DPrimitivesProp_->setComputeHistogram(computeHistogram_->isChecked());
}

//----------------------------------------------------------------------------------------------
//      helper functions
//----------------------------------------------------------------------------------------------
void TransFunc2DPrimitivesPropertyEditor::markerColorChanged(int h, int s, int v, int a) {
    QColor tmp = QColor::fromHsv(h, s, v); tmp.setAlpha(a);
    mappingCanvas_->changeCurrentColorSlot(tmp);
}










//QLayout* TransFunc2DPrimitivesPropertyEditor::createMappingLayout() {
//    QLayout* layout = new QHBoxLayout();
//
//    transCanvas_ = new tgt::QtCanvas("", tgt::ivec2(1, 1), tgt::GLCanvas::RGBADD, 0, true);
//
//    painter_ = new TransFunc2DPrimitivesPropertyWidgetPainter(transCanvas_); //sets the painter in the canvas
//    transCanvas_->init();
//    if (transFuncGradient_)
//        painter_->setTransFunc(transFuncGradient_);
//
//    layout->addWidget(transCanvas_);
//    return layout;
//}
//
//QLayout* TransFunc2DPrimitivesPropertyEditor::createButtonLayout() {
//    QBoxLayout* buttonLayout;
//    if (orientation_ == Qt::Vertical)
//        buttonLayout = new QHBoxLayout();
//    else
//        buttonLayout = new QVBoxLayout();
//
//    clearButton_ = new QToolButton();
//    clearButton_->setIcon(QIcon(":/qt/icons/clear.png"));
//    clearButton_->setToolTip(tr("Reset transfer function to default"));
//
//    loadButton_ = new QToolButton();
//    loadButton_->setIcon(QIcon(":/qt/icons/open.png"));
//    loadButton_->setToolTip(tr("Load transfer function"));
//
//    saveButton_ = new QToolButton();
//    saveButton_->setIcon(QIcon(":/qt/icons/save.png"));
//    saveButton_->setToolTip(tr("Save transfer function"));
//
//    gridEnabledButton_ = new QToolButton();
//    gridEnabledButton_->setCheckable(true);
//    gridEnabledButton_->setChecked(false);
//    gridEnabledButton_->setIcon(QIcon(":/qt/icons/grid.png"));
//    gridEnabledButton_->setToolTip(tr("Show grid"));
//
//    histogramEnabledButton_ = new QToolButton();
//    histogramEnabledButton_->setCheckable(true);
//    histogramEnabledButton_->setChecked(false);
//    histogramEnabledButton_->setEnabled(false);
//    histogramEnabledButton_->setIcon(QIcon(":/qt/icons/histogram.png"));
//    histogramEnabledButton_->setToolTip(tr("Show data histogram"));
//
//    fitToDomainButton_ = new QToolButton();
//    fitToDomainButton_->setIcon(QIcon(":/qt/icons/histogram_fit.png"));
//    fitToDomainButton_->setToolTip(tr("Fit Domain to Data"));
//
//    buttonLayout->addWidget(clearButton_);
//    buttonLayout->addWidget(loadButton_);
//    buttonLayout->addWidget(saveButton_);
//    buttonLayout->addSpacing(4);
//    buttonLayout->addWidget(gridEnabledButton_);
//    buttonLayout->addWidget(histogramEnabledButton_);
//    buttonLayout->addWidget(fitToDomainButton_);
//
//    buttonLayout->addStretch();
//
//    return buttonLayout;
//}
//
//QLayout* TransFunc2DPrimitivesPropertyEditor::createPrimitivesButtonLayout() {
//    QBoxLayout* buttonLayout;
//    buttonLayout = new QHBoxLayout();
//
//    quadButton_ = new QToolButton();
//    quadButton_->setIcon(QIcon(":/qt/icons/quad.png"));
//    quadButton_->setToolTip(tr("Add a quad"));
//
//    bananaButton_ = new QToolButton();
//    bananaButton_->setIcon(QIcon(":/qt/icons/banana.png"));
//    bananaButton_->setToolTip(tr("Add a banana"));
//
//    deleteButton_ = new QToolButton();
//    deleteButton_->setIcon(QIcon(":/qt/icons/eraser.png"));
//    deleteButton_->setToolTip(tr("Delete selected primitive"));
//
//    colorButton_ = new QToolButton();
//    colorButton_->setIcon(QIcon(":/qt/icons/colorize.png"));
//    colorButton_->setToolTip(tr("Change the color of the selected primitive"));
//
//    buttonLayout->addWidget(quadButton_);
//    buttonLayout->addWidget(bananaButton_);
//    buttonLayout->addWidget(deleteButton_);
//    buttonLayout->addWidget(colorButton_);
//
//    buttonLayout->addStretch();
//
//    return buttonLayout;
//}
//
//QLayout* TransFunc2DPrimitivesPropertyEditor::createSliderLayout() {
//    QVBoxLayout* layout = new QVBoxLayout();
//
//    histogramBrightness_ = new QSlider(Qt::Horizontal);
//    histogramBrightness_->setEnabled(false);
//    histogramBrightness_->setMinimum(10);
//    histogramBrightness_->setMaximum(1000);
//    histogramBrightness_->setValue(200);
//
//    histogramLog_ = new QCheckBox(tr("Logarithmic Histogram"));
//    histogramLog_->setEnabled(false);
//    histogramLog_->setCheckState(Qt::Checked);
//
//    transparency_ = new QSlider(Qt::Horizontal);
//    transparency_->setEnabled(false);
//    transparency_->setMinimum(0);
//    transparency_->setMaximum(255);
//
//    fuzziness_ = new QSlider(Qt::Horizontal);
//    fuzziness_->setEnabled(false);
//    fuzziness_->setMinimum(0);
//    fuzziness_->setMaximum(100);
//
//    QHBoxLayout* intensityDomainLayout = new QHBoxLayout();
//    lowerIntensityDomainSpin_ = new QDoubleSpinBox();
//    upperIntensityDomainSpin_ = new QDoubleSpinBox();
//    intensityDomainLayout->addWidget(lowerIntensityDomainSpin_);
//    intensityDomainLayout->addWidget(new QLabel(tr("Intensity Domain")));
//    intensityDomainLayout->addWidget(upperIntensityDomainSpin_);
//
//    QHBoxLayout* gradientDomainLayout = new QHBoxLayout();
//    lowerGradientDomainSpin_ = new QDoubleSpinBox();
//    upperGradientDomainSpin_ = new QDoubleSpinBox();
//    gradientDomainLayout->addWidget(lowerGradientDomainSpin_);
//    gradientDomainLayout->addWidget(new QLabel(tr("Gradient Domain")));
//    gradientDomainLayout->addWidget(upperGradientDomainSpin_);
//
//    layout->addWidget(new QLabel(tr("Histogram Brightness:")));
//    layout->addWidget(histogramBrightness_);
//    layout->addWidget(histogramLog_);
//    layout->addWidget(new QLabel(tr("Transparency:")));
//    layout->addWidget(transparency_);
//    layout->addWidget(new QLabel(tr("Fuzziness:")));
//    layout->addWidget(fuzziness_);
//    layout->addLayout(intensityDomainLayout);
//    layout->addLayout(gradientDomainLayout);
//
//    layout->addStretch();
//
//    return layout;
//}
//
//void TransFunc2DPrimitivesPropertyEditor::createWidgets() {
//    QWidget* mapping = new QWidget();
//    QWidget* slider = new QWidget();
//
//    QLayout* mappingLayout = createMappingLayout();
//    QLayout* buttonLayout = createButtonLayout();
//    QLayout* primitivesButtonLayout = createPrimitivesButtonLayout();
//    QLayout* sliderLayout = createSliderLayout();
//
//    QSplitter* splitter = new QSplitter(orientation_);
//    QLayout* buttonSlider;
//    if (orientation_ == Qt::Vertical) {
//        buttonSlider = new QVBoxLayout();
//        buttonSlider->addItem(buttonLayout);
//        buttonSlider->addItem(primitivesButtonLayout);
//        buttonSlider->addItem(mappingLayout);
//        mapping->setLayout(buttonSlider);
//        slider->setLayout(sliderLayout);
//    }
//    else {
//        buttonSlider = new QHBoxLayout();
//        buttonSlider->addItem(buttonLayout);
//        buttonSlider->addItem(new QSpacerItem(5, 1));
//        buttonSlider->addItem(sliderLayout);
//
//        QVBoxLayout* vLayout = new QVBoxLayout();
//        vLayout->addItem(buttonSlider);
//        vLayout->addItem(primitivesButtonLayout);
//        vLayout->addStretch(10);
//
//        mapping->setLayout(mappingLayout);
//        slider->setLayout(vLayout);
//    }
//
//    splitter->setChildrenCollapsible(false);
//    splitter->addWidget(mapping);
//    splitter->addWidget(slider);
//
//    //mapping is more stretched then buttons and slider
//    splitter->setStretchFactor(0, 5);
//
//    QHBoxLayout* mainLayout = new QHBoxLayout();
//    mainLayout->addWidget(splitter);
//
//    setLayout(mainLayout);
//    updateFromProperty();
//    volumeChanged();
//}
//
//void TransFunc2DPrimitivesPropertyEditor::createConnections() {
//    // buttons
//    connect(loadButton_,  SIGNAL(clicked()), this, SLOT(loadTransferFunction()));
//    connect(saveButton_,  SIGNAL(clicked()), this, SLOT(saveTransferFunction()));
//    connect(clearButton_, SIGNAL(clicked()), painter_, SLOT(resetTransferFunction()));
//
//    connect(gridEnabledButton_,      SIGNAL(clicked()), this, SLOT(toggleShowGrid()));
//    connect(histogramEnabledButton_, SIGNAL(clicked()), this, SLOT(toggleShowHistogram()));
//    connect(fitToDomainButton_, SIGNAL(clicked()), this, SLOT(fitToDomain()));
//
//    connect(quadButton_,   SIGNAL(clicked()), painter_, SLOT(addQuadPrimitive()));
//    connect(bananaButton_, SIGNAL(clicked()), painter_, SLOT(addBananaPrimitive()));
//    connect(deleteButton_, SIGNAL(clicked()), painter_, SLOT(deletePrimitive()));
//    connect(colorButton_,  SIGNAL(clicked()), painter_, SLOT(colorizePrimitive()));
//
//    connect(histogramBrightness_, SIGNAL(sliderMoved(int)), painter_, SLOT(histogramBrightnessChanged(int)));
//    connect(histogramLog_, SIGNAL(stateChanged(int)), painter_, SLOT(toggleHistogramLogarithmic(int)));
//
//    // slider
//    connect(fuzziness_, SIGNAL(valueChanged(int)), painter_, SLOT(fuzzinessChanged(int)));
//    connect(transparency_, SIGNAL(valueChanged(int)), painter_, SLOT(transparencyChanged(int)));
//
//    connect(fuzziness_, SIGNAL(sliderPressed()), this, SLOT(startTracking()));
//    connect(transparency_, SIGNAL(sliderPressed()), this, SLOT(startTracking()));
//
//    connect(fuzziness_, SIGNAL(sliderReleased()), this, SLOT(stopTracking()));
//    connect(transparency_, SIGNAL(sliderReleased()), this, SLOT(stopTracking()));
//
//    connect(painter_, SIGNAL(setTransparencySlider(int)), this, SLOT(setTransparency(int)));
//    connect(painter_, SIGNAL(primitiveDeselected()), this, SLOT(primitiveDeselected()));
//    connect(painter_, SIGNAL(primitiveSelected()), this, SLOT(primitiveSelected()));
//    connect(painter_, SIGNAL(toggleInteractionMode(bool)), this, SLOT(toggleInteractionMode(bool)));
//    connect(painter_, SIGNAL(repaintSignal()), this, SLOT(repaintSignal()));
//
//    // domain
//    connect(lowerGradientDomainSpin_, SIGNAL(valueChanged(double)), this, SLOT(domainChanged()));
//    connect(upperGradientDomainSpin_, SIGNAL(valueChanged(double)), this, SLOT(domainChanged()));
//
//    connect(lowerIntensityDomainSpin_, SIGNAL(valueChanged(double)), this, SLOT(domainChanged()));
//    connect(upperIntensityDomainSpin_, SIGNAL(valueChanged(double)), this, SLOT(domainChanged()));
//}
//
//void TransFunc2DPrimitivesPropertyEditor::domainChanged() {
//    if(transFuncGradient_) {
//        transFuncGradient_->setDomain(vec2(lowerIntensityDomainSpin_->value(), upperIntensityDomainSpin_->value()), 0);
//        transFuncGradient_->setDomain(vec2(lowerGradientDomainSpin_->value(), upperGradientDomainSpin_->value()), 1);
//        //painter_->updateTF();
//        transCanvas_->update();
//    }
//}
//
//void TransFunc2DPrimitivesPropertyEditor::fitToDomain() {
//    if(transFuncGradient_) {
//        //painter_->fitToDomain();
//        property_->invalidate();
//    }
//}
//
//void TransFunc2DPrimitivesPropertyEditor::loadTransferFunction() {
//    if (!transFuncGradient_) {
//        LWARNING("No valid transfer function assigned");
//        return;
//    }
//
//    if(TransFuncIOHelperQt::loadTransferFunction(transFuncGradient_)) {
//        //painter_->updateTF();
//        transCanvas_->update();
//    }
//}
//
//void TransFunc2DPrimitivesPropertyEditor::saveTransferFunction() {
//    if (!transFuncGradient_) {
//        LWARNING("No valid transfer function assigned");
//        return;
//    }
//
//    TransFuncIOHelperQt::saveTransferFunction(transFuncGradient_);
//}
//
//void TransFunc2DPrimitivesPropertyEditor::toggleShowGrid() {
//    //painter_->toggleShowGrid(gridEnabledButton_->isChecked());
//    transCanvas_->update();
//}
//
//void TransFunc2DPrimitivesPropertyEditor::toggleShowHistogram() {
//    /*volume_ = property_->getVolume();
//
//    painter_->volumeChanged(volume_);
//    painter_->setHistogramVisible(histogramEnabledButton_->isChecked());
//    */histogramBrightness_->setEnabled(histogramEnabledButton_->isChecked());
//    histogramLog_->setEnabled(histogramEnabledButton_->isChecked());
//    transCanvas_->update();
//}
//
//void TransFunc2DPrimitivesPropertyEditor::primitiveSelected() {
//    /*const TransFuncPrimitive* p = painter_->getSelectedPrimitive();
//    fuzziness_->setValue(static_cast<int>(p->getFuzziness() * 100.f));
//    fuzziness_->setEnabled(true);
//    transparency_->setValue(p->getColor().a);
//    transparency_->setEnabled(true);*/
//}
//
//void TransFunc2DPrimitivesPropertyEditor::primitiveDeselected() {
//    fuzziness_->setValue(0);
//    fuzziness_->setEnabled(false);
//    transparency_->setValue(0);
//    transparency_->setEnabled(false);
//}
//
//void TransFunc2DPrimitivesPropertyEditor::setTransparency(int trans) {
//    transparency_->blockSignals(true);
//    transparency_->setValue(trans);
//    transparency_->blockSignals(false);
//}
//
//void TransFunc2DPrimitivesPropertyEditor::startTracking() {
//    toggleInteractionMode(true);
//}
//
//void TransFunc2DPrimitivesPropertyEditor::stopTracking() {
//    toggleInteractionMode(false);
//}
//
//void TransFunc2DPrimitivesPropertyEditor::updateFromProperty() {
//    const VolumeBase* newHandle = property_->getVolume();
//    /*if (newHandle != volume_) {
//        volume_ = newHandle;
//        volumeChanged();
//    }*/
//
//    transFuncGradient_ = dynamic_cast<TransFunc2DPrimitives*>(property_->get());
//    if (transFuncGradient_) {
//        transCanvas_->update();
//
//        // update domain spinboxes, disabling signals to prevent immediate change notification
//        lowerIntensityDomainSpin_->blockSignals(true);
//        lowerIntensityDomainSpin_->setValue(transFuncGradient_->getDomain(0).x);
//        lowerIntensityDomainSpin_->blockSignals(false);
//        upperIntensityDomainSpin_->blockSignals(true);
//        upperIntensityDomainSpin_->setValue(transFuncGradient_->getDomain(0).y);
//        upperIntensityDomainSpin_->blockSignals(false);
//        lowerGradientDomainSpin_->blockSignals(true);
//        lowerGradientDomainSpin_->setValue(transFuncGradient_->getDomain(1).x);
//        lowerGradientDomainSpin_->blockSignals(false);
//        upperGradientDomainSpin_->blockSignals(true);
//        upperGradientDomainSpin_->setValue(transFuncGradient_->getDomain(1).y);
//        upperGradientDomainSpin_->blockSignals(false);
//    }
//    //else
//        //resetEditor();
//}
//
//void TransFunc2DPrimitivesPropertyEditor::volumeChanged() {
//    // update control elements
//    histogramEnabledButton_->setEnabled(true);
//    histogramEnabledButton_->blockSignals(true);
//    histogramEnabledButton_->setChecked(false);
//    histogramEnabledButton_->blockSignals(false);
//
//    histogramBrightness_->setEnabled(false);
//    histogramLog_->setEnabled(false);
//    histogramBrightness_->blockSignals(true);
//    histogramBrightness_->setValue(100);
//    histogramBrightness_->blockSignals(false);
//
//    //if (volume_ && volume_->getRepresentation<VolumeRAM>()) {
//        //int bits = volume_->getRepresentation<VolumeRAM>()->getBitsAllocated() / volume_->getRepresentation<VolumeRAM>()->getNumChannels();
//        //maximumIntensity_ = static_cast<int>(pow(2.f, static_cast<float>(bits)))-1;
//    //}
//
//    // propagate volume to painter where the histogram is calculated
//    //painter_->volumeChanged(volume_);
//}
//
//void TransFunc2DPrimitivesPropertyEditor::resetEditor() {
//    // Need to make the GL context of this context current, as we use it for OpenGL calls
//    // with the transFuncGradient_ later on.
//    //transCanvas_->makeGLContextCurrent();
//    if (property_->get() != transFuncGradient_) {
//        LDEBUG("The pointers of property and transfer function do not match."
//                << "Creating new transfer function object.....");
//
//        transFuncGradient_ = new TransFunc2DPrimitives();
//        property_->set(transFuncGradient_);
//
//        painter_->setTransFunc(transFuncGradient_);
//    }
//
//    // reset transfer function to default, e.g. empty tf
//    //painter_->resetTransferFunction();
//
//    // update mapping widget, e.g. canvas
//    transCanvas_->update();
//
//    // cause repaint of volume rendering
//    /*property_->notifyChange();
//    emit transferFunctionChanged();*/
//    //transCanvas_->restorePreviuosGLContext();
//}
//
//
//void TransFunc2DPrimitivesPropertyEditor::repaintSignal() {
//    //property_->notifyChange();
//
//    //emit transferFunctionChanged();
//}

} // namespace voreen
