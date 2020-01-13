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

#include "voreen/qt/widgets/property/transfunc/2d/transfunc2dpropertyeditor.h"

#include "voreen/core/properties/transfunc/2d/transfunc2dproperty.h"

#include "voreen/core/datastructures/volume/volumeminmax.h"
#include "voreen/core/datastructures/volume/histogram.h"

#include <QMenu>
#include <QBoxLayout>

namespace {
        const int MAX_SPIN_BOX_DIGITS = 7;  ///< maximal digits of domain and threshold
}

namespace voreen {

    const std::string TransFunc2DPropertyEditor::loggerCat_("TransFunc2DPropertyEditor");

//----------------------------------------------------------------------------------------------
//      Constructor and Qt Stuff
//----------------------------------------------------------------------------------------------
TransFunc2DPropertyEditor::TransFunc2DPropertyEditor(TransFunc2DProperty* prop, QWidget* parent)
    : TransFuncPropertyEditorBase(prop, parent)
    , alphaButton_(0), alphaMenu_(0), gammaSpinX_(0), gammaSpinY_(0)
    , lowerVolumeBoundLabelX_(0), upperVolumeBoundLabelX_(0), lowerVolumeBoundLabelY_(0), upperVolumeBoundLabelY_(0)
    , lowerDomainSpinX_(0), upperDomainSpinX_(0), lowerDomainSpinY_(0), upperDomainSpinY_(0)
    , lowerThresholdSpinX_(0), upperThresholdSpinX_(0), lowerThresholdSpinY_(0), upperThresholdSpinY_(0)
    , cutoffDomainFittingX_(0), cutoffDomainFittingY_(0), domainFittingStrategy_(0), fitDomainToData_(0)
    , transFunc2DProp_(prop)
{
}

TransFunc2DPropertyEditor::~TransFunc2DPropertyEditor() {
}

void TransFunc2DPropertyEditor::updateFromProperty() {
    tgtAssert(transFunc2DProp_, "No property");
    tgtAssert(initialized_, "editor is not initialized!");
    if(transFunc2DProp_->get()) {
        //update color map settings
        updateAlphaButton(transFunc2DProp_->get()->getAlphaMode());
        gammaSpinX_->blockSignals(true);
        gammaSpinX_->setValue(transFunc2DProp_->get()->getGammaValue().x);
        gammaSpinX_->blockSignals(false);
        gammaSpinY_->blockSignals(true);
        gammaSpinY_->setValue(transFunc2DProp_->get()->getGammaValue().y);
        gammaSpinY_->blockSignals(false);

        //update domain and threshold settings
        if (transFunc2DProp_->getDomainFittingStrategy() < domainFittingStrategy_->count()) {
            domainFittingStrategy_->blockSignals(true);
            domainFittingStrategy_->setCurrentIndex(transFunc2DProp_->getDomainFittingStrategy());
            domainFittingStrategy_->blockSignals(false);
        }
        else {
            LERROR("Invalid domain fitting strategy: " << transFunc2DProp_->getDomainFittingStrategy());
        }
        updateVolumeDataBoundsFromProperty();
        updateDomainAndThresholdFromProperty();
    }
}

//----------------------------------------------------------------------------------------------
//      Layout
//----------------------------------------------------------------------------------------------
QGroupBox* TransFunc2DPropertyEditor::createColorMapSettingsBox() {
    //alpha
    QLabel* alphaLabel = new QLabel("Transparency:");
    alphaButton_ = new QToolButton();
    alphaButton_->setToolButtonStyle(Qt::ToolButtonTextBesideIcon);
    alphaButton_->setIcon(QIcon(":/qt/icons/alpha.png"));
    alphaButton_->setText("Alpha ");
    alphaButton_->setPopupMode(QToolButton::InstantPopup);
    alphaMenu_ = new QMenu();
    alphaMenu_->addAction(QIcon(":/qt/icons/alpha_trans.png"),"Transparent");
    alphaMenu_->addAction(QIcon(":/qt/icons/alpha_use.png"),"Use Alpha");
    alphaMenu_->addAction(QIcon(":/qt/icons/alpha_opaque.png"),"Opaque");
    alphaButton_->setMenu(alphaMenu_);
    connect(alphaMenu_, SIGNAL(triggered(QAction*)), this, SLOT(alphaButtonClicked(QAction*)));
    //gamma
    gammaSpinX_ = new QDoubleSpinBox();
    gammaSpinX_->setRange(0.1, 5.0);
    gammaSpinX_->setValue(1.0);
    gammaSpinX_->setSingleStep(0.1);
    gammaSpinX_->setKeyboardTracking(false);
    gammaSpinX_->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
    gammaSpinX_->setDecimals(2);
    gammaSpinX_->setFixedWidth(6*3+25);
    QLabel* gammaLabelX = new QLabel((" Gamma (" + getXGuiName() + "):").c_str());
    connect(gammaSpinX_, SIGNAL(valueChanged(double)), this, SLOT(gammaSpinChanged()));
    gammaSpinY_ = new QDoubleSpinBox();
    gammaSpinY_->setRange(0.1, 5.0);
    gammaSpinY_->setValue(1.0);
    gammaSpinY_->setSingleStep(0.1);
    gammaSpinY_->setKeyboardTracking(false);
    gammaSpinY_->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
    gammaSpinY_->setDecimals(2);
    gammaSpinY_->setFixedWidth(6*3+25);
    QLabel* gammaLabelY = new QLabel((" Gamma (" + getYGuiName() + "):").c_str());
    connect(gammaSpinY_, SIGNAL(valueChanged(double)), this, SLOT(gammaSpinChanged()));

    QBoxLayout* layoutH = new QHBoxLayout();
    layoutH->addWidget(alphaLabel);
    layoutH->addWidget(alphaButton_);
    //layout->addStretch();
    QGridLayout* layoutG = new QGridLayout();
    layoutG->addWidget(gammaLabelX,0,0); layoutG->addWidget(gammaSpinX_,0,1);
    layoutG->addWidget(gammaLabelY,1,0); layoutG->addWidget(gammaSpinY_,1,1);
    layoutH->addLayout(layoutG);

    QGroupBox* cmsBox = new QGroupBox("Color Map Settings");
    cmsBox->setLayout(layoutH);
    cmsBox->setSizePolicy(QSizePolicy::Minimum,QSizePolicy::Fixed);

    return cmsBox;
}

QGroupBox* TransFunc2DPropertyEditor::createDomainAndThresholdBox() {
    //data set domains X
    QHBoxLayout* dataSetDomainLayoutX = new QHBoxLayout();
    lowerVolumeBoundLabelX_ = new QLabel();
    lowerVolumeBoundLabelX_->setAlignment(Qt::AlignCenter);
    lowerVolumeBoundLabelX_->setFixedWidth(6*MAX_SPIN_BOX_DIGITS+25);
    upperVolumeBoundLabelX_ = new QLabel();
    upperVolumeBoundLabelX_->setAlignment(Qt::AlignCenter);
    upperVolumeBoundLabelX_->setFixedWidth(6*MAX_SPIN_BOX_DIGITS+25);
    QLabel* dataLabelX = new QLabel(("Data Set Bounds   (" + getXGuiName() + ")").c_str());
    dataSetDomainLayoutX->addWidget(lowerVolumeBoundLabelX_);
    dataSetDomainLayoutX->addStretch(1);
    dataSetDomainLayoutX->addWidget(dataLabelX);
    dataSetDomainLayoutX->addStretch(1);
    dataSetDomainLayoutX->addWidget(upperVolumeBoundLabelX_);
    //data set domains Y
    QHBoxLayout* dataSetDomainLayoutY = new QHBoxLayout();
    lowerVolumeBoundLabelY_ = new QLabel();
    lowerVolumeBoundLabelY_->setAlignment(Qt::AlignCenter);
    lowerVolumeBoundLabelY_->setFixedWidth(6*MAX_SPIN_BOX_DIGITS+25);
    upperVolumeBoundLabelY_ = new QLabel();
    upperVolumeBoundLabelY_->setAlignment(Qt::AlignCenter);
    upperVolumeBoundLabelY_->setFixedWidth(6*MAX_SPIN_BOX_DIGITS+25);
    QLabel* dataLabelY = new QLabel(("Data Set Bounds   (" + getYGuiName() + ")").c_str());
    dataSetDomainLayoutY->addWidget(lowerVolumeBoundLabelY_);
    dataSetDomainLayoutY->addStretch(1);
    dataSetDomainLayoutY->addWidget(dataLabelY);
    dataSetDomainLayoutY->addStretch(1);
    dataSetDomainLayoutY->addWidget(upperVolumeBoundLabelY_);

    //domains X
    QHBoxLayout* domainLayoutX = new QHBoxLayout();
    lowerDomainSpinX_ = new QDoubleSpinBox();
    upperDomainSpinX_ = new QDoubleSpinBox();
    upperDomainSpinX_->setRange(-9999999.0, 9999999.0);
    lowerDomainSpinX_->setRange(-9999999.0, 9999999.0);
    upperDomainSpinX_->setValue(1.0);
    lowerDomainSpinX_->setValue(0.0);
    upperDomainSpinX_->setKeyboardTracking(false);
    lowerDomainSpinX_->setKeyboardTracking(false);
    upperDomainSpinX_->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
    lowerDomainSpinX_->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
    upperDomainSpinX_->setDecimals(MAX_SPIN_BOX_DIGITS-1);
    lowerDomainSpinX_->setDecimals(MAX_SPIN_BOX_DIGITS-1);
    upperDomainSpinX_->setFixedWidth(6*MAX_SPIN_BOX_DIGITS+25);
    lowerDomainSpinX_->setFixedWidth(6*MAX_SPIN_BOX_DIGITS+25);
    QLabel* domainLabelX = new QLabel();
    domainLabelX->setText(("Color Map Domain (" + getXGuiName() + ")").c_str());
    domainLayoutX->addWidget(lowerDomainSpinX_);
    domainLayoutX->addStretch(1);
    domainLayoutX->addWidget(domainLabelX);
    domainLayoutX->addStretch(1);
    domainLayoutX->addWidget(upperDomainSpinX_);
    connect(lowerDomainSpinX_, SIGNAL(valueChanged(double)), this, SLOT(lowerDomainSpinChanged()));
    connect(upperDomainSpinX_, SIGNAL(valueChanged(double)), this, SLOT(upperDomainSpinChanged()));
    //domains Y
    QHBoxLayout* domainLayoutY = new QHBoxLayout();
    lowerDomainSpinY_ = new QDoubleSpinBox();
    upperDomainSpinY_ = new QDoubleSpinBox();
    upperDomainSpinY_->setRange(-9999999.0, 9999999.0);
    lowerDomainSpinY_->setRange(-9999999.0, 9999999.0);
    upperDomainSpinY_->setValue(1.0);
    lowerDomainSpinY_->setValue(0.0);
    upperDomainSpinY_->setKeyboardTracking(false);
    lowerDomainSpinY_->setKeyboardTracking(false);
    upperDomainSpinY_->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
    lowerDomainSpinY_->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
    upperDomainSpinY_->setDecimals(MAX_SPIN_BOX_DIGITS-1);
    lowerDomainSpinY_->setDecimals(MAX_SPIN_BOX_DIGITS-1);
    upperDomainSpinY_->setFixedWidth(6*MAX_SPIN_BOX_DIGITS+25);
    lowerDomainSpinY_->setFixedWidth(6*MAX_SPIN_BOX_DIGITS+25);
    QLabel* domainLabelY = new QLabel();
    domainLabelY->setText(("Color Map Domain (" + getYGuiName() + ")").c_str());
    domainLayoutY->addWidget(lowerDomainSpinY_);
    domainLayoutY->addStretch(1);
    domainLayoutY->addWidget(domainLabelY);
    domainLayoutY->addStretch(1);
    domainLayoutY->addWidget(upperDomainSpinY_);
    connect(lowerDomainSpinY_, SIGNAL(valueChanged(double)), this, SLOT(lowerDomainSpinChanged()));
    connect(upperDomainSpinY_, SIGNAL(valueChanged(double)), this, SLOT(upperDomainSpinChanged()));

    //autofit
    QVBoxLayout* autoFitLayout = new QVBoxLayout();
    QHBoxLayout* topAutoFitLayout = new QHBoxLayout();
    QGridLayout* leftAutoFitLayout = new QGridLayout();
    QHBoxLayout* bottomAutoFitLayout = new QHBoxLayout();
    QLabel* cutoffLabelX = new QLabel((getXGuiName() + " Cutoff (%):").c_str());
    cutoffDomainFittingX_ = new QSpinBox();
    cutoffDomainFittingX_->setRange(0,25);
    cutoffDomainFittingX_->setValue(0);
    cutoffDomainFittingX_->setKeyboardTracking(false);
    cutoffDomainFittingX_->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
    cutoffDomainFittingX_->setFixedWidth(37);
    cutoffDomainFittingX_->setEnabled(false);
    QLabel* cutoffLabelY = new QLabel((getYGuiName() + " Cutoff (%):").c_str());
    cutoffDomainFittingY_ = new QSpinBox();
    cutoffDomainFittingY_->setRange(0,25);
    cutoffDomainFittingY_->setValue(0);
    cutoffDomainFittingY_->setKeyboardTracking(false);
    cutoffDomainFittingY_->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
    cutoffDomainFittingY_->setFixedWidth(37);
    cutoffDomainFittingY_->setEnabled(false);
    fitDomainToData_ = new QToolButton();
    fitDomainToData_->setText(" Fit Domain To Data ");
    domainFittingStrategy_ = new QComboBox();
    domainFittingStrategy_->addItem("Never Fit Domain",     TransFuncPropertyBase::FIT_DOMAIN_NEVER);
    domainFittingStrategy_->addItem("Fit Initial Domain",   TransFuncPropertyBase::FIT_DOMAIN_INITIAL);
    domainFittingStrategy_->addItem("Always Fit Domain",    TransFuncPropertyBase::FIT_DOMAIN_ALWAYS);
    leftAutoFitLayout->addWidget(cutoffLabelX,0,0); leftAutoFitLayout->addWidget(cutoffDomainFittingX_,0,1);
    leftAutoFitLayout->addWidget(cutoffLabelY,1,0); leftAutoFitLayout->addWidget(cutoffDomainFittingY_,1,1);
    topAutoFitLayout->addWidget(domainFittingStrategy_);
    topAutoFitLayout->addLayout(leftAutoFitLayout);
    bottomAutoFitLayout->addWidget(fitDomainToData_);
    autoFitLayout->addLayout(topAutoFitLayout);
    autoFitLayout->addLayout(bottomAutoFitLayout);
    connect(domainFittingStrategy_, SIGNAL(currentIndexChanged(int)), this, SLOT(domainFittingStrategyChanged(int)));
    connect(fitDomainToData_, SIGNAL(clicked()), this, SLOT(fitDomainClicked()));

    //threshold X
    QHBoxLayout* thresholdLayoutX = new QHBoxLayout();
    lowerThresholdSpinX_ = new QDoubleSpinBox();
    upperThresholdSpinX_ = new QDoubleSpinBox();
    upperThresholdSpinX_->setRange(-9999999.0, 9999999.0);
    lowerThresholdSpinX_->setRange(-9999999.0, 9999999.0);
    upperThresholdSpinX_->setValue(1.0);
    lowerThresholdSpinX_->setValue(0.0);
    upperThresholdSpinX_->setKeyboardTracking(false);
    lowerThresholdSpinX_->setKeyboardTracking(false);
    upperThresholdSpinX_->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
    lowerThresholdSpinX_->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
    upperThresholdSpinX_->setDecimals(MAX_SPIN_BOX_DIGITS-1);
    lowerThresholdSpinX_->setDecimals(MAX_SPIN_BOX_DIGITS-1);
    upperThresholdSpinX_->setFixedWidth(6*MAX_SPIN_BOX_DIGITS+25);
    lowerThresholdSpinX_->setFixedWidth(6*MAX_SPIN_BOX_DIGITS+25);
    QLabel* thresholdLabelX = new QLabel();
    thresholdLabelX->setText(("Threshold Bounds (" + getXGuiName() + ")").c_str());
    thresholdLayoutX->addWidget(lowerThresholdSpinX_);
    thresholdLayoutX->addStretch(1);
    thresholdLayoutX->addWidget(thresholdLabelX);
    thresholdLayoutX->addStretch(1);
    thresholdLayoutX->addWidget(upperThresholdSpinX_);
    connect(lowerThresholdSpinX_, SIGNAL(valueChanged(double)), this, SLOT(lowerThresholdSpinChanged()));
    connect(upperThresholdSpinX_, SIGNAL(valueChanged(double)), this, SLOT(upperThresholdSpinChanged()));
    //threshold Y
    QHBoxLayout* thresholdLayoutY = new QHBoxLayout();
    lowerThresholdSpinY_ = new QDoubleSpinBox();
    upperThresholdSpinY_ = new QDoubleSpinBox();
    upperThresholdSpinY_->setRange(-9999999.0, 9999999.0);
    lowerThresholdSpinY_->setRange(-9999999.0, 9999999.0);
    upperThresholdSpinY_->setValue(1.0);
    lowerThresholdSpinY_->setValue(0.0);
    upperThresholdSpinY_->setKeyboardTracking(false);
    lowerThresholdSpinY_->setKeyboardTracking(false);
    upperThresholdSpinY_->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
    lowerThresholdSpinY_->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
    upperThresholdSpinY_->setDecimals(MAX_SPIN_BOX_DIGITS-1);
    lowerThresholdSpinY_->setDecimals(MAX_SPIN_BOX_DIGITS-1);
    upperThresholdSpinY_->setFixedWidth(6*MAX_SPIN_BOX_DIGITS+25);
    lowerThresholdSpinY_->setFixedWidth(6*MAX_SPIN_BOX_DIGITS+25);
    QLabel* thresholdLabelY = new QLabel();
    thresholdLabelY->setText(("Threshold Bounds (" + getYGuiName() + ")").c_str());
    thresholdLayoutY->addWidget(lowerThresholdSpinY_);
    thresholdLayoutY->addStretch(1);
    thresholdLayoutY->addWidget(thresholdLabelY);
    thresholdLayoutY->addStretch(1);
    thresholdLayoutY->addWidget(upperThresholdSpinY_);
    connect(lowerThresholdSpinY_, SIGNAL(valueChanged(double)), this, SLOT(lowerThresholdSpinChanged()));
    connect(upperThresholdSpinY_, SIGNAL(valueChanged(double)), this, SLOT(upperThresholdSpinChanged()));

    //line
    QFrame* line1 = new QFrame();
    line1->setFrameShape(QFrame::HLine);
    line1->setFrameShadow(QFrame::Sunken);
    QFrame* line2 = new QFrame();
    line2->setFrameShape(QFrame::HLine);
    line2->setFrameShadow(QFrame::Sunken);

    QVBoxLayout* vBox = new QVBoxLayout();
    vBox->addLayout(dataSetDomainLayoutX);
    vBox->addLayout(domainLayoutX);
    vBox->addLayout(dataSetDomainLayoutY);
    vBox->addLayout(domainLayoutY);
    vBox->addWidget(line1);
    vBox->addLayout(autoFitLayout);
    vBox->addWidget(line2);
    vBox->addLayout(thresholdLayoutX);
    vBox->addLayout(thresholdLayoutY);

    QGroupBox* domainBox = new QGroupBox("Domain Settings");
    domainBox->setLayout(vBox);
    domainBox->setSizePolicy(QSizePolicy::Minimum,QSizePolicy::Fixed);

    return domainBox;
}

//----------------------------------------------------------------------------------------------
//      layout slots
//----------------------------------------------------------------------------------------------
void TransFunc2DPropertyEditor::alphaButtonClicked(QAction* action) {
    if (transFunc2DProp_ && transFunc2DProp_->get()) {
        TransFuncBase::AlphaMode mode;
        if(action->text() == "Transparent") {
            mode = TransFuncBase::TF_ZERO_ALPHA;
        } else if(action->text() == "Opaque") {
            mode = TransFuncBase::TF_ONE_ALPHA;
        } else {
            mode = TransFuncBase::TF_USE_ALPHA;
        }
        transFunc2DProp_->get()->setAlphaMode(mode);
        updateAlphaButton(mode);
        transFunc2DProp_->invalidate();
    }
}

void TransFunc2DPropertyEditor::fitDomainClicked() {
    tgtAssert(transFunc2DProp_ && transFunc2DProp_->get(), "no property or tf!");
    transFunc2DProp_->applyDomainFromData();
}

void TransFunc2DPropertyEditor::domainFittingStrategyChanged(int index)  {
    tgtAssert(transFunc2DProp_ && transFunc2DProp_->get(), "no property or tf!");

    if (index < 0 || index >= domainFittingStrategy_->count()) {
        LERROR("domainFittingStrategyChanged(): invalid index " << index);
        return;
    }

    TransFuncPropertyBase::DomainAutoFittingStrategy strategy =
        (TransFuncPropertyBase::DomainAutoFittingStrategy)(domainFittingStrategy_->itemData(index).toUInt());
    transFunc2DProp_->setDomainFittingStrategy(strategy);
}

void TransFunc2DPropertyEditor::gammaSpinChanged() {
    tgtAssert(transFunc2DProp_ && transFunc2DProp_->get(), "no property or tf!");
    transFunc2DProp_->get()->setGammaValue(tgt::vec2(gammaSpinX_->value(),gammaSpinY_->value()));
    transFunc2DProp_->invalidate();
}

void TransFunc2DPropertyEditor::lowerDomainSpinChanged() {
    tgtAssert(transFunc2DProp_ && transFunc2DProp_->get(), "no property or tf!");

    //increment value of lower mapping spin when it equals value of upper mapping spin
    if (lowerDomainSpinX_->value() >= upperDomainSpinX_->value()) {
        upperDomainSpinX_->blockSignals(true);
        upperDomainSpinX_->setValue(lowerDomainSpinX_->value() + 0.01);
        upperDomainSpinX_->blockSignals(false);
    }
    if (lowerDomainSpinY_->value() >= upperDomainSpinY_->value()) {
        upperDomainSpinY_->blockSignals(true);
        upperDomainSpinY_->setValue(lowerDomainSpinY_->value() + 0.01);
        upperDomainSpinY_->blockSignals(false);
    }

    transFunc2DProp_->get()->setDomain(lowerDomainSpinX_->value(), upperDomainSpinX_->value(),0);
    transFunc2DProp_->get()->setDomain(lowerDomainSpinY_->value(), upperDomainSpinY_->value(),1);

    transFunc2DProp_->invalidate(); //updates the threshold
}

void TransFunc2DPropertyEditor::upperDomainSpinChanged() {
    tgtAssert(transFunc2DProp_ && transFunc2DProp_->get(), "no property or tf!");

    //increment value of upper mapping spin when it equals value of lower mapping spin
    if (upperDomainSpinX_->value() <= lowerDomainSpinX_->value()) {
        lowerDomainSpinX_->blockSignals(true);
        lowerDomainSpinX_->setValue(upperDomainSpinX_->value() - 0.01);
        lowerDomainSpinX_->blockSignals(false);
    }
    if (upperDomainSpinY_->value() <= lowerDomainSpinY_->value()) {
        lowerDomainSpinY_->blockSignals(true);
        lowerDomainSpinY_->setValue(upperDomainSpinY_->value() - 0.01);
        lowerDomainSpinY_->blockSignals(false);
    }

    transFunc2DProp_->get()->setDomain(lowerDomainSpinX_->value(), upperDomainSpinX_->value(),0);
    transFunc2DProp_->get()->setDomain(lowerDomainSpinY_->value(), upperDomainSpinY_->value(),1);

    transFunc2DProp_->invalidate(); //updates the threshold
}

void TransFunc2DPropertyEditor::lowerThresholdSpinChanged() {
    tgtAssert(transFunc2DProp_ && transFunc2DProp_->get(), "no property or tf!");

    //increment value of lower mapping spin when it equals value of upper mapping spin
    if (lowerThresholdSpinX_->value() > upperThresholdSpinX_->value()) {
        upperThresholdSpinX_->blockSignals(true);
        upperThresholdSpinX_->setValue(lowerThresholdSpinX_->value());
        upperThresholdSpinX_->blockSignals(false);
    }
    if (lowerThresholdSpinY_->value() > upperThresholdSpinY_->value()) {
        upperThresholdSpinY_->blockSignals(true);
        upperThresholdSpinY_->setValue(lowerThresholdSpinY_->value());
        upperThresholdSpinY_->blockSignals(false);
    }

    float lowerX = (lowerThresholdSpinX_->value()-lowerDomainSpinX_->value())/(upperDomainSpinX_->value()-lowerDomainSpinX_->value());
    float upperX = (upperThresholdSpinX_->value()-lowerDomainSpinX_->value())/(upperDomainSpinX_->value()-lowerDomainSpinX_->value());
    float lowerY = (lowerThresholdSpinY_->value()-lowerDomainSpinY_->value())/(upperDomainSpinY_->value()-lowerDomainSpinY_->value());
    float upperY = (upperThresholdSpinY_->value()-lowerDomainSpinY_->value())/(upperDomainSpinY_->value()-lowerDomainSpinY_->value());

    transFunc2DProp_->get()->setThreshold(lowerX, upperX,0);
    transFunc2DProp_->get()->setThreshold(lowerY, upperY,1);
    transFunc2DProp_->invalidate();
}

void TransFunc2DPropertyEditor::upperThresholdSpinChanged() {
    tgtAssert(transFunc2DProp_ && transFunc2DProp_->get(), "no property or tf!");

    //increment value of lower mapping spin when it equals value of upper mapping spin
    if (lowerThresholdSpinX_->value() > upperThresholdSpinX_->value()) {
        lowerThresholdSpinX_->blockSignals(true);
        lowerThresholdSpinX_->setValue(upperThresholdSpinX_->value());
        lowerThresholdSpinX_->blockSignals(false);
    }
    if (lowerThresholdSpinY_->value() > upperThresholdSpinY_->value()) {
        lowerThresholdSpinY_->blockSignals(true);
        lowerThresholdSpinY_->setValue(upperThresholdSpinY_->value());
        lowerThresholdSpinY_->blockSignals(false);
    }

    float lowerX = (lowerThresholdSpinX_->value()-lowerDomainSpinX_->value())/(upperDomainSpinX_->value()-lowerDomainSpinX_->value());
    float upperX = (upperThresholdSpinX_->value()-lowerDomainSpinX_->value())/(upperDomainSpinX_->value()-lowerDomainSpinX_->value());
    float lowerY = (lowerThresholdSpinY_->value()-lowerDomainSpinY_->value())/(upperDomainSpinY_->value()-lowerDomainSpinY_->value());
    float upperY = (upperThresholdSpinY_->value()-lowerDomainSpinY_->value())/(upperDomainSpinY_->value()-lowerDomainSpinY_->value());

    transFunc2DProp_->get()->setThreshold(lowerX, upperX,0);
    transFunc2DProp_->get()->setThreshold(lowerY, upperY,1);
    transFunc2DProp_->invalidate();
}


//----------------------------------------------------------------------------------------------
//      update functions
//----------------------------------------------------------------------------------------------
void TransFunc2DPropertyEditor::updateAlphaButton(TransFuncBase::AlphaMode mode) {
    alphaButton_->blockSignals(true);
    switch(mode) {
    case TransFuncBase::TF_ZERO_ALPHA:
        alphaButton_->setIcon(QIcon(":/qt/icons/alpha_trans.png"));
        alphaButton_->setText("   ");
    break;
    case TransFuncBase::TF_USE_ALPHA:
        alphaButton_->setIcon(QIcon(":/qt/icons/alpha_use.png"));
        alphaButton_->setText("   ");
    break;
    case TransFuncBase::TF_ONE_ALPHA:
        alphaButton_->setIcon(QIcon(":/qt/icons/alpha_opaque.png"));
        alphaButton_->setText("   ");
    break;
    }
    alphaButton_->blockSignals(false);
}

void TransFunc2DPropertyEditor::updateDomainAndThresholdFromProperty(){
   tgtAssert(transFunc2DProp_->get(), "no tf");

   //block signals
   lowerDomainSpinX_->blockSignals(true);    upperDomainSpinX_->blockSignals(true);
   lowerThresholdSpinX_->blockSignals(true); upperThresholdSpinX_->blockSignals(true);
   lowerDomainSpinY_->blockSignals(true);    upperDomainSpinY_->blockSignals(true);
   lowerThresholdSpinY_->blockSignals(true); upperThresholdSpinY_->blockSignals(true);

   //get values
   tgt::vec2 domainX    = transFunc2DProp_->get()->getDomain(0);
   tgt::vec2 domainY    = transFunc2DProp_->get()->getDomain(1);
   tgt::vec2 thresholdX = transFunc2DProp_->get()->getThreshold(0);
   tgt::vec2 thresholdY = transFunc2DProp_->get()->getThreshold(1);
   float domainRangeX = domainX.y - domainX.x;
   float domainRangeY = domainY.y - domainY.x;

   //set decimals X
   if(std::abs(domainX.x) < 1.f){
       lowerDomainSpinX_->setDecimals(MAX_SPIN_BOX_DIGITS-1);
       lowerThresholdSpinX_->setDecimals(MAX_SPIN_BOX_DIGITS-1);
   } else {
       lowerDomainSpinX_->setDecimals(MAX_SPIN_BOX_DIGITS-(static_cast<int>( log10( std::abs( domainX.x ) ) ) + 1));
       lowerThresholdSpinX_->setDecimals(MAX_SPIN_BOX_DIGITS-(static_cast<int>( log10( std::abs( domainX.x ) ) ) + 1));
   }
   if(std::abs(domainX.y) < 1.0) {
       upperDomainSpinX_->setDecimals(MAX_SPIN_BOX_DIGITS-1);
       upperThresholdSpinX_->setDecimals(MAX_SPIN_BOX_DIGITS-1);
   }
   else {
       upperDomainSpinX_->setDecimals(MAX_SPIN_BOX_DIGITS-(static_cast<int>( log10( std::abs( domainX.y ) ) ) + 1));
       upperThresholdSpinX_->setDecimals(MAX_SPIN_BOX_DIGITS-(static_cast<int>( log10( std::abs( domainX.y ) ) ) + 1));
   }
   //set decimals Y
   if(std::abs(domainY.x) < 1.f){
       lowerDomainSpinY_->setDecimals(MAX_SPIN_BOX_DIGITS-1);
       lowerThresholdSpinY_->setDecimals(MAX_SPIN_BOX_DIGITS-1);
   } else {
       lowerDomainSpinY_->setDecimals(MAX_SPIN_BOX_DIGITS-(static_cast<int>( log10( std::abs( domainY.x ) ) ) + 1));
       lowerThresholdSpinY_->setDecimals(MAX_SPIN_BOX_DIGITS-(static_cast<int>( log10( std::abs( domainY.x ) ) ) + 1));
   }
   if(std::abs(domainY.y) < 1.0) {
       upperDomainSpinY_->setDecimals(MAX_SPIN_BOX_DIGITS-1);
       upperThresholdSpinY_->setDecimals(MAX_SPIN_BOX_DIGITS-1);
   }
   else {
       upperDomainSpinY_->setDecimals(MAX_SPIN_BOX_DIGITS-(static_cast<int>( log10( std::abs( domainY.y ) ) ) + 1));
       upperThresholdSpinY_->setDecimals(MAX_SPIN_BOX_DIGITS-(static_cast<int>( log10( std::abs( domainY.y ) ) ) + 1));
   }

   //set stepsize
   lowerDomainSpinX_->setSingleStep(domainRangeX/1000.0);    upperDomainSpinX_->setSingleStep(domainRangeX/1000.0);
   lowerThresholdSpinX_->setSingleStep(domainRangeX/1000.0); upperThresholdSpinX_->setSingleStep(domainRangeX/1000.0);
   lowerDomainSpinY_->setSingleStep(domainRangeY/1000.0);    upperDomainSpinY_->setSingleStep(domainRangeY/1000.0);
   lowerThresholdSpinY_->setSingleStep(domainRangeY/1000.0); upperThresholdSpinY_->setSingleStep(domainRangeY/1000.0);

   //update domain slider
   lowerDomainSpinX_->setValue(domainX.x); upperDomainSpinX_->setValue(domainX.y);
   lowerDomainSpinY_->setValue(domainY.x); upperDomainSpinY_->setValue(domainY.y);

   //update threshold slider (relative to domain)
   lowerThresholdSpinX_->setRange(domainX.x,domainX.y);                 upperThresholdSpinX_->setRange(domainX.x,domainX.y);
   lowerThresholdSpinX_->setValue(domainX.x+domainRangeX*thresholdX.x); upperThresholdSpinX_->setValue(domainX.x+domainRangeX*thresholdX.y);
   lowerThresholdSpinY_->setRange(domainY.x,domainY.y);                 upperThresholdSpinY_->setRange(domainY.x,domainY.y);
   lowerThresholdSpinY_->setValue(domainY.x+domainRangeY*thresholdY.x); upperThresholdSpinY_->setValue(domainY.x+domainRangeY*thresholdY.y);

   lowerDomainSpinX_->blockSignals(false);    upperDomainSpinX_->blockSignals(false);
   lowerThresholdSpinX_->blockSignals(false); upperThresholdSpinX_->blockSignals(false);
   lowerDomainSpinY_->blockSignals(false);    upperDomainSpinY_->blockSignals(false);
   lowerThresholdSpinY_->blockSignals(false); upperThresholdSpinY_->blockSignals(false);
}


} // namespace voreen
