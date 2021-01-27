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

#include "voreen/qt/widgets/property/transfunc/1d/transfunc1dpropertyeditor.h"

#include "voreen/core/properties/transfunc/1d/transfunc1dproperty.h"

#include "voreen/core/datastructures/volume/volumeminmax.h"
#include "voreen/core/datastructures/volume/histogram.h"

#include <QMenu>
#include <QBoxLayout>

namespace {
        const int MAX_SPIN_BOX_DIGITS = 7;  ///< maximal digits of domain and threshold
}

namespace voreen {

    const std::string TransFunc1DPropertyEditor::loggerCat_("voreen.qt.TransFunc1DPropertyEditor");

//----------------------------------------------------------------------------------------------
//      Constructor and Qt Stuff
//----------------------------------------------------------------------------------------------
TransFunc1DPropertyEditor::TransFunc1DPropertyEditor(TransFunc1DProperty* prop, QWidget* parent)
    : TransFuncPropertyEditorBase(prop, parent)
    , lowerDomainSpin_(0), upperDomainSpin_(0), lowerThresholdSpin_(0), upperThresholdSpin_(0)
    , lowerVolumeBoundLabel_(0), upperVolumeBoundLabel_(0), cutoffDomainFitting_(0), domainFittingStrategy_(0)
    , fitDomainToData_(0)
    , transFunc1DProp_(prop)
{
}

TransFunc1DPropertyEditor::~TransFunc1DPropertyEditor() {
}

void TransFunc1DPropertyEditor::updateFromProperty() {
    tgtAssert(transFunc1DProp_, "No property");
    tgtAssert(initialized_, "editor is not initialized!");
    if(transFunc1DProp_->get()) {
        //update color map settings
        updateAlphaButton(transFunc1DProp_->get()->getAlphaMode());
        gammaSpin_->blockSignals(true);
        gammaSpin_->setValue(transFunc1DProp_->get()->getGammaValue());
        gammaSpin_->blockSignals(false);

        //update domain and threshold settings
        if (transFunc1DProp_->getDomainFittingStrategy() < domainFittingStrategy_->count()) {
            domainFittingStrategy_->blockSignals(true);
            domainFittingStrategy_->setCurrentIndex(transFunc1DProp_->getDomainFittingStrategy());
            domainFittingStrategy_->blockSignals(false);
        }
        else {
            LERROR("Invalid domain fitting strategy: " << transFunc1DProp_->getDomainFittingStrategy());
        }
        updateVolumeDataBoundsFromProperty();
        updateDomainAndThresholdFromProperty();
    }
}

//----------------------------------------------------------------------------------------------
//      Layout
//----------------------------------------------------------------------------------------------
QGroupBox* TransFunc1DPropertyEditor::createColorMapSettingsBox() {
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
    gammaSpin_ = new QDoubleSpinBox();
    gammaSpin_->setRange(0.1, 5.0);
    gammaSpin_->setValue(1.0);
    gammaSpin_->setSingleStep(0.1);
    gammaSpin_->setKeyboardTracking(false);
    gammaSpin_->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
    gammaSpin_->setDecimals(2);
    gammaSpin_->setFixedWidth(6*3+25);
    QLabel* gammaLabel = new QLabel("     Gamma Value:");
    connect(gammaSpin_, SIGNAL(valueChanged(double)), this, SLOT(gammaSpinChanged(double)));

    QBoxLayout* layout = new QHBoxLayout();
    layout->addWidget(alphaLabel);
    layout->addWidget(alphaButton_);
    //layout->addStretch();
    layout->addWidget(gammaLabel);
    layout->addWidget(gammaSpin_);

    QGroupBox* cmsBox = new QGroupBox("Color Map Settings");
    cmsBox->setLayout(layout);
    cmsBox->setSizePolicy(QSizePolicy::Minimum,QSizePolicy::Fixed);

    return cmsBox;
}

QGroupBox* TransFunc1DPropertyEditor::createDomainAndThresholdBox() {
    //data set domains
    QHBoxLayout* dataSetDomainLayout = new QHBoxLayout();
    lowerVolumeBoundLabel_ = new QLabel();
    lowerVolumeBoundLabel_->setAlignment(Qt::AlignCenter);
    lowerVolumeBoundLabel_->setFixedWidth(6*MAX_SPIN_BOX_DIGITS+25);
    upperVolumeBoundLabel_ = new QLabel();
    upperVolumeBoundLabel_->setAlignment(Qt::AlignCenter);
    upperVolumeBoundLabel_->setFixedWidth(6*MAX_SPIN_BOX_DIGITS+25);
    QLabel* dataLabel = new QLabel("Data Set Bounds");
    dataSetDomainLayout->addWidget(lowerVolumeBoundLabel_);
    dataSetDomainLayout->addStretch(1);
    dataSetDomainLayout->addWidget(dataLabel);
    dataSetDomainLayout->addStretch(1);
    dataSetDomainLayout->addWidget(upperVolumeBoundLabel_);
    //domains
    QHBoxLayout* domainLayout = new QHBoxLayout();
    lowerDomainSpin_ = new QDoubleSpinBox();
    upperDomainSpin_ = new QDoubleSpinBox();
    upperDomainSpin_->setRange(-9999999.0, 9999999.0);
    lowerDomainSpin_->setRange(-9999999.0, 9999999.0);
    upperDomainSpin_->setValue(1.0);
    lowerDomainSpin_->setValue(0.0);
    upperDomainSpin_->setKeyboardTracking(false);
    lowerDomainSpin_->setKeyboardTracking(false);
    upperDomainSpin_->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
    lowerDomainSpin_->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
    upperDomainSpin_->setDecimals(MAX_SPIN_BOX_DIGITS-1);
    lowerDomainSpin_->setDecimals(MAX_SPIN_BOX_DIGITS-1);
    upperDomainSpin_->setFixedWidth(6*MAX_SPIN_BOX_DIGITS+25);
    lowerDomainSpin_->setFixedWidth(6*MAX_SPIN_BOX_DIGITS+25);
    QLabel* domainLabel = new QLabel();
    domainLabel->setText("Color Map Domain");
    domainLayout->addWidget(lowerDomainSpin_);
    domainLayout->addStretch(1);
    domainLayout->addWidget(domainLabel);
    domainLayout->addStretch(1);
    domainLayout->addWidget(upperDomainSpin_);
    connect(lowerDomainSpin_, SIGNAL(valueChanged(double)), this, SLOT(lowerDomainSpinChanged(double)));
    connect(upperDomainSpin_, SIGNAL(valueChanged(double)), this, SLOT(upperDomainSpinChanged(double)));
    //autofit
    QVBoxLayout* autoFitLayout = new QVBoxLayout();
    QHBoxLayout* firstAutoFitLayout = new QHBoxLayout();
    QHBoxLayout* secondAutoFitLayout = new QHBoxLayout();
    QLabel* cutoffLabel = new QLabel("Fitting Cutoff (%):");
    cutoffDomainFitting_ = new QSpinBox();
    cutoffDomainFitting_->setRange(0,25);
    cutoffDomainFitting_->setValue(0);
    cutoffDomainFitting_->setKeyboardTracking(false);
    cutoffDomainFitting_->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
    cutoffDomainFitting_->setFixedWidth(37);
    cutoffDomainFitting_->setEnabled(false);
    fitDomainToData_ = new QToolButton();
    fitDomainToData_->setText(" Fit Domain To Data ");
    domainFittingStrategy_ = new QComboBox();
    domainFittingStrategy_->addItem("Never Fit Domain",     TransFuncPropertyBase::FIT_DOMAIN_NEVER);
    domainFittingStrategy_->addItem("Fit Initial Domain",   TransFuncPropertyBase::FIT_DOMAIN_INITIAL);
    domainFittingStrategy_->addItem("Always Fit Domain",    TransFuncPropertyBase::FIT_DOMAIN_ALWAYS);
    firstAutoFitLayout->addWidget(domainFittingStrategy_);
    firstAutoFitLayout->addWidget(cutoffLabel);
    firstAutoFitLayout->addWidget(cutoffDomainFitting_);
    secondAutoFitLayout->addWidget(fitDomainToData_);
    autoFitLayout->addLayout(firstAutoFitLayout);
    autoFitLayout->addLayout(secondAutoFitLayout);
    connect(domainFittingStrategy_, SIGNAL(currentIndexChanged(int)), this, SLOT(domainFittingStrategyChanged(int)));
    connect(fitDomainToData_, SIGNAL(clicked()), this, SLOT(fitDomainClicked()));
    //threshold
    QHBoxLayout* thresholdLayout = new QHBoxLayout();
    lowerThresholdSpin_ = new QDoubleSpinBox();
    upperThresholdSpin_ = new QDoubleSpinBox();
    upperThresholdSpin_->setRange(-9999999.0, 9999999.0);
    lowerThresholdSpin_->setRange(-9999999.0, 9999999.0);
    upperThresholdSpin_->setValue(1.0);
    lowerThresholdSpin_->setValue(0.0);
    upperThresholdSpin_->setKeyboardTracking(false);
    lowerThresholdSpin_->setKeyboardTracking(false);
    upperThresholdSpin_->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
    lowerThresholdSpin_->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
    upperThresholdSpin_->setDecimals(MAX_SPIN_BOX_DIGITS-1);
    lowerThresholdSpin_->setDecimals(MAX_SPIN_BOX_DIGITS-1);
    upperThresholdSpin_->setFixedWidth(6*MAX_SPIN_BOX_DIGITS+25);
    lowerThresholdSpin_->setFixedWidth(6*MAX_SPIN_BOX_DIGITS+25);
    QLabel* thresholdLabel = new QLabel();
    thresholdLabel->setText("Threshold Bounds");
    thresholdLayout->addWidget(lowerThresholdSpin_);
    thresholdLayout->addStretch(1);
    thresholdLayout->addWidget(thresholdLabel);
    thresholdLayout->addStretch(1);
    thresholdLayout->addWidget(upperThresholdSpin_);
    connect(lowerThresholdSpin_, SIGNAL(valueChanged(double)), this, SLOT(lowerThresholdSpinChanged(double)));
    connect(upperThresholdSpin_, SIGNAL(valueChanged(double)), this, SLOT(upperThresholdSpinChanged(double)));
    //line
    QFrame* line1 = new QFrame();
    line1->setFrameShape(QFrame::HLine);
    line1->setFrameShadow(QFrame::Sunken);
    QFrame* line2 = new QFrame();
    line2->setFrameShape(QFrame::HLine);
    line2->setFrameShadow(QFrame::Sunken);
    QVBoxLayout* vBox = new QVBoxLayout();
    vBox->addLayout(dataSetDomainLayout);
    vBox->addLayout(domainLayout);
    vBox->addWidget(line1);
    vBox->addLayout(autoFitLayout);
    vBox->addWidget(line2);
    vBox->addLayout(thresholdLayout);

    QGroupBox* domainBox = new QGroupBox("Domain Settings");
    domainBox->setLayout(vBox);
    domainBox->setSizePolicy(QSizePolicy::Minimum,QSizePolicy::Fixed);

    return domainBox;
}

//----------------------------------------------------------------------------------------------
//      layout slots
//----------------------------------------------------------------------------------------------
void TransFunc1DPropertyEditor::alphaButtonClicked(QAction* action) {
    if (transFunc1DProp_ && transFunc1DProp_->get()) {
        TransFuncBase::AlphaMode mode;
        if(action->text() == "Transparent") {
            mode = TransFuncBase::TF_ZERO_ALPHA;
        } else if(action->text() == "Opaque") {
            mode = TransFuncBase::TF_ONE_ALPHA;
        } else {
            mode = TransFuncBase::TF_USE_ALPHA;
        }
        transFunc1DProp_->get()->setAlphaMode(mode);
        updateAlphaButton(mode);
        transFunc1DProp_->invalidate();
    }
}

void TransFunc1DPropertyEditor::fitDomainClicked() {
    tgtAssert(transFunc1DProp_ && transFunc1DProp_->get(), "no property or tf!");
    transFunc1DProp_->applyDomainFromData();
}

void TransFunc1DPropertyEditor::domainFittingStrategyChanged(int index)  {
    tgtAssert(transFunc1DProp_ && transFunc1DProp_->get(), "no property or tf!");

    if (index < 0 || index >= domainFittingStrategy_->count()) {
        LERROR("domainFittingStrategyChanged(): invalid index " << index);
        return;
    }

    TransFuncPropertyBase::DomainAutoFittingStrategy strategy =
        (TransFuncPropertyBase::DomainAutoFittingStrategy)(domainFittingStrategy_->itemData(index).toUInt());
    transFunc1DProp_->setDomainFittingStrategy(strategy);
}

void TransFunc1DPropertyEditor::gammaSpinChanged(double gamma) {
    tgtAssert(transFunc1DProp_ && transFunc1DProp_->get(), "no property or tf!");
    transFunc1DProp_->get()->setGammaValue((float)gamma);
    transFunc1DProp_->invalidate();
}

void TransFunc1DPropertyEditor::lowerDomainSpinChanged(double value) {
    tgtAssert(transFunc1DProp_ && transFunc1DProp_->get(), "no property or tf!");

    //increment value of lower mapping spin when it equals value of upper mapping spin
    if (value >= upperDomainSpin_->value()) {
        upperDomainSpin_->blockSignals(true);
        upperDomainSpin_->setValue(value+0.01);
        upperDomainSpin_->blockSignals(false);
    }

    transFunc1DProp_->get()->setDomain(tgt::vec2(lowerDomainSpin_->value(), upperDomainSpin_->value()));
    transFunc1DProp_->invalidate(); //updates the threshold
}

void TransFunc1DPropertyEditor::upperDomainSpinChanged(double value) {
    tgtAssert(transFunc1DProp_ && transFunc1DProp_->get(), "no property or tf!");

    //increment value of upper mapping spin when it equals value of lower mapping spin
    if (value <= lowerDomainSpin_->value()) {
        lowerDomainSpin_->blockSignals(true);
        lowerDomainSpin_->setValue(value-0.01);
        lowerDomainSpin_->blockSignals(false);
    }

    transFunc1DProp_->get()->setDomain(tgt::vec2(lowerDomainSpin_->value(), upperDomainSpin_->value()));
    transFunc1DProp_->invalidate(); //updates the threshold
}

void TransFunc1DPropertyEditor::lowerThresholdSpinChanged(double value) {
    tgtAssert(transFunc1DProp_ && transFunc1DProp_->get(), "no property or tf!");

    //increment value of lower mapping spin when it equals value of upper mapping spin
    if (value > upperThresholdSpin_->value()) {
        upperThresholdSpin_->blockSignals(true);
        upperThresholdSpin_->setValue(value);
        upperThresholdSpin_->blockSignals(false);
    }

    float lower = (lowerThresholdSpin_->value()-lowerDomainSpin_->value())/(upperDomainSpin_->value()-lowerDomainSpin_->value());
    float upper = (upperThresholdSpin_->value()-lowerDomainSpin_->value())/(upperDomainSpin_->value()-lowerDomainSpin_->value());

    transFunc1DProp_->get()->setThreshold(lower, upper);
    transFunc1DProp_->invalidate();
}

void TransFunc1DPropertyEditor::upperThresholdSpinChanged(double value) {
    tgtAssert(transFunc1DProp_ && transFunc1DProp_->get(), "no property or tf!");

    //increment value of lower mapping spin when it equals value of upper mapping spin
    if (value < lowerThresholdSpin_->value()) {
        lowerThresholdSpin_->blockSignals(true);
        lowerThresholdSpin_->setValue(value);
        lowerThresholdSpin_->blockSignals(false);
    }

    float lower = (lowerThresholdSpin_->value()-lowerDomainSpin_->value())/(upperDomainSpin_->value()-lowerDomainSpin_->value());
    float upper = (upperThresholdSpin_->value()-lowerDomainSpin_->value())/(upperDomainSpin_->value()-lowerDomainSpin_->value());

    transFunc1DProp_->get()->setThreshold(lower, upper);
    transFunc1DProp_->invalidate();
}


//----------------------------------------------------------------------------------------------
//      update functions
//----------------------------------------------------------------------------------------------
void TransFunc1DPropertyEditor::updateAlphaButton(TransFuncBase::AlphaMode mode) {
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

void TransFunc1DPropertyEditor::updateVolumeDataBoundsFromProperty() {
    tgtAssert(transFunc1DProp_, "no Property");

    if(const VolumeBase* volume = transFunc1DProp_->getVolume()) {
        //use min max if present
        if(volume->hasDerivedData<VolumeMinMax>()) {
            lowerVolumeBoundLabel_->setText(QString::number(volume->getRealWorldMapping().normalizedToRealWorld(
                                            volume->getDerivedData<VolumeMinMax>()->getMinNormalized(transFunc1DProp_->getVolumeChannel()))));
            upperVolumeBoundLabel_->setText(QString::number(volume->getRealWorldMapping().normalizedToRealWorld(
                                            volume->getDerivedData<VolumeMinMax>()->getMaxNormalized(transFunc1DProp_->getVolumeChannel()))));
        }
        else {
            volume->getDerivedDataThreaded<VolumeMinMax>();
            lowerVolumeBoundLabel_->setText("0");
            upperVolumeBoundLabel_->setText("1");
            return;
        }
    } else {
        //set default, if no volume is present
        lowerVolumeBoundLabel_->setText("0");
        upperVolumeBoundLabel_->setText("1");
        return;
    }
}

void TransFunc1DPropertyEditor::updateDomainAndThresholdFromProperty(){
   tgtAssert(transFunc1DProp_->get(), "no tf");

   //block signals
   lowerDomainSpin_->blockSignals(true);
   upperDomainSpin_->blockSignals(true);
   lowerThresholdSpin_->blockSignals(true);
   upperThresholdSpin_->blockSignals(true);

   //get values
   tgt::vec2 domain    = transFunc1DProp_->get()->getDomain();
   tgt::vec2 threshold = transFunc1DProp_->get()->getThreshold();
   float domainRange = domain.y - domain.x;

   //set decimals
   if(std::abs(domain.x) < 1.f){
       lowerDomainSpin_->setDecimals(MAX_SPIN_BOX_DIGITS-1);
       lowerThresholdSpin_->setDecimals(MAX_SPIN_BOX_DIGITS-1);
   } else {
       lowerDomainSpin_->setDecimals(MAX_SPIN_BOX_DIGITS-(static_cast<int>( log10( std::abs( domain.x ) ) ) + 1));
       lowerThresholdSpin_->setDecimals(MAX_SPIN_BOX_DIGITS-(static_cast<int>( log10( std::abs( domain.x ) ) ) + 1));
   }
   if(std::abs(domain.y) < 1.0) {
       upperDomainSpin_->setDecimals(MAX_SPIN_BOX_DIGITS-1);
       upperThresholdSpin_->setDecimals(MAX_SPIN_BOX_DIGITS-1);
   }
   else {
       upperDomainSpin_->setDecimals(MAX_SPIN_BOX_DIGITS-(static_cast<int>( log10( std::abs( domain.y ) ) ) + 1));
       upperThresholdSpin_->setDecimals(MAX_SPIN_BOX_DIGITS-(static_cast<int>( log10( std::abs( domain.y ) ) ) + 1));
   }

   //set stepsize
   lowerDomainSpin_->setSingleStep(domainRange/1000.0);
   upperDomainSpin_->setSingleStep(domainRange/1000.0);
   lowerThresholdSpin_->setSingleStep(domainRange/1000.0);
   upperThresholdSpin_->setSingleStep(domainRange/1000.0);

   //update domain slider
   lowerDomainSpin_->setValue(domain.x);
   upperDomainSpin_->setValue(domain.y);

   //update threshold slider (relative to domain)
   lowerThresholdSpin_->setRange(domain.x,domain.y);
   upperThresholdSpin_->setRange(domain.x,domain.y);
   lowerThresholdSpin_->setValue(domain.x+domainRange*threshold.x);
   upperThresholdSpin_->setValue(domain.x+domainRange*threshold.y);

   lowerDomainSpin_->blockSignals(false);
   upperDomainSpin_->blockSignals(false);
   lowerThresholdSpin_->blockSignals(false);
   upperThresholdSpin_->blockSignals(false);
}


} // namespace voreen
