/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany.                        *
 * Visualization and Computer Graphics Group <http://viscg.uni-muenster.de>        *
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

#include "voreen/qt/widgets/property/transfunc/transfuncpropertyeditorbase.h"

#include "voreen/core/properties/transfunc/transfuncpropertybase.h"
#include "voreen/core/datastructures/transfunc/transfuncbase.h"

#include "voreen/qt/widgets/property/transfunc/utils/transfunciohelperqt.h"
#include "voreen/core/utils/voreenqualitymode.h"

#include <QBoxLayout>
#include <QSplitter>
#include <QLabel>

namespace voreen {

TransFuncPropertyEditorBase::TransFuncPropertyEditorBase(TransFuncPropertyBase* property, QWidget* parent)
    : QWidget(parent)
    , loadButton_(0) , saveButton_(0), clearButton_(0)
    , colorPicker_(0), colorLumPicker_(0)
    , baseProperty_(property)
    , initialized_(false)
{
    tgtAssert(property, "no property");
}

TransFuncPropertyEditorBase::~TransFuncPropertyEditorBase() {
}

void TransFuncPropertyEditorBase::initialize() {
    tgtAssert(!initialized_, "has been initialized already!");
    layoutComponents();
    initialized_ = true;
    updateFromProperty();
}

    //-----------------------------
    // Layout functions
    //-----------------------------
void TransFuncPropertyEditorBase::layoutComponents() {
    //create splitter for main layout
    QSplitter* splitter = new QSplitter(Qt::Horizontal);
    splitter->setStretchFactor(0, QSizePolicy::Expanding); // mapping should be stretched
    splitter->setStretchFactor(1, QSizePolicy::Fixed); // color should not be stretched
    splitter->setChildrenCollapsible(true);
    //first init right
    QWidget* right = layoutRightComponents();
    QWidget* left  = layoutLeftComponents();
    splitter->addWidget(left);
    splitter->addWidget(right);
    //create main layout
    QHBoxLayout* mainLayout = new QHBoxLayout();
    mainLayout->setMargin(4);
    mainLayout->addWidget(splitter);
    //set main layout
    setLayout(mainLayout);
}

QWidget* TransFuncPropertyEditorBase::layoutRightComponents() {
    //widget to return;
    QWidget* mainRight = new QWidget();

    QGroupBox* colorBox      = createColorPickerBox();
    QGroupBox* baseButtonBox = createBaseButtonBox();
    QWidget* toolTip = createTooltip();

    QGroupBox* settingsBox   = createColorMapSettingsBox();
    QGroupBox* domainBox     = createDomainAndThresholdBox();

    QGridLayout* gridLayout = new QGridLayout();
    gridLayout->addWidget(colorBox,0,0);
    gridLayout->addWidget(baseButtonBox,0,1);
    // add the tooltip widget if it exists
    if (toolTip) {
        gridLayout->addWidget(toolTip, 0, 2);
        gridLayout->setAlignment(toolTip, static_cast<Qt::Alignment>(Qt::AlignRight | Qt::AlignTop));
        // the tooltip column is smaller than the others
        gridLayout->setColumnStretch(0, 2);
        gridLayout->setColumnStretch(1, 2);
        gridLayout->setColumnStretch(2, 1);
    }
    gridLayout->addWidget(settingsBox,1,0,2,0);
    gridLayout->addWidget(domainBox,3,0,2,0);
    gridLayout->setSizeConstraint(QLayout::SetFixedSize);

    mainRight->setLayout(gridLayout);
    return mainRight;
}

QGroupBox* TransFuncPropertyEditorBase::createBaseButtonBox() {
    QBoxLayout* buttonLayout = new QVBoxLayout();

    clearButton_ = new QToolButton();
    clearButton_->setToolButtonStyle(Qt::ToolButtonTextBesideIcon);
    clearButton_->setText("Reset");
    clearButton_->setFixedWidth(70);
    clearButton_->setIcon(QIcon(":/qt/icons/clear.png"));
    connect(clearButton_, SIGNAL(clicked()), this, SLOT(clearButtonClicked()));

    loadButton_ = new QToolButton();
    loadButton_->setToolButtonStyle(Qt::ToolButtonTextBesideIcon);
    loadButton_->setText("Load");
    loadButton_->setFixedWidth(70);
    loadButton_->setIcon(QIcon(":/qt/icons/open.png"));
    connect(loadButton_, SIGNAL(clicked()), this, SLOT(loadButtonClicked()));

    saveButton_ = new QToolButton();
    saveButton_->setToolButtonStyle(Qt::ToolButtonTextBesideIcon);
    saveButton_->setText("Save");
    saveButton_->setFixedWidth(70);
    saveButton_->setIcon(QIcon(":/qt/icons/save.png"));
    connect(saveButton_, SIGNAL(clicked()), this, SLOT(saveButtonClicked()));

    //layout
    buttonLayout->setAlignment(Qt::AlignHCenter);
    buttonLayout->addWidget(clearButton_,Qt::AlignCenter);
    buttonLayout->addWidget(loadButton_,Qt::AlignCenter);
    buttonLayout->addWidget(saveButton_,Qt::AlignCenter);

    QGroupBox* baseBox = new QGroupBox("Load and Save");
    baseBox->setLayout(buttonLayout);
    baseBox->setSizePolicy(QSizePolicy::Minimum,QSizePolicy::Minimum);

    return baseBox;
}

QGroupBox* TransFuncPropertyEditorBase::createColorPickerBox() {
    // ColorPicker
    colorPicker_ = new ColorPicker();
    colorPicker_->setFrameStyle(QFrame::Panel | QFrame::Sunken);
    colorPicker_->setFixedWidth(100);
    colorPicker_->setFixedHeight(100);

    // ColorLuminacePicker
    colorLumPicker_ = new ColorLuminancePicker();
    colorLumPicker_->setFixedWidth(20);
    colorLumPicker_->setFixedHeight(100);

    QHBoxLayout* hBoxColor = new QHBoxLayout();
    hBoxColor->addWidget(colorPicker_);
    hBoxColor->addWidget(colorLumPicker_);

    QGroupBox* colorBox = new QGroupBox("Color Picking");
    colorBox->setLayout(hBoxColor);
    colorBox->setSizePolicy(QSizePolicy::Fixed,QSizePolicy::Fixed);

    connect(colorPicker_, SIGNAL(toggleInteractionModeSignal(bool)), this, SLOT(toggleInteractionMode(bool)));
    connect(colorLumPicker_, SIGNAL(toggleInteractionModeSignal(bool)), this, SLOT(toggleInteractionMode(bool)));
    connect(colorPicker_, SIGNAL(newHSSignal(int,int)), colorLumPicker_, SLOT(updateHSSlot(int,int)));

    return colorBox;
}

QWidget* TransFuncPropertyEditorBase::createTooltip() {
    // 0 will be returned per default. No tooltip will be added in this case.
    return 0;
}

    //-----------------------------
    // Layout slots
    //-----------------------------
void TransFuncPropertyEditorBase::clearButtonClicked() {
    baseProperty_->get()->reset();
    baseProperty_->applyDomainFromData();
    baseProperty_->invalidate();
}

void TransFuncPropertyEditorBase::saveButtonClicked() {
    tgtAssert(baseProperty_->get(), "No valid transfer function assigned");
    TransFuncIOHelperQt::saveTransferFunction(baseProperty_->get());
}

void TransFuncPropertyEditorBase::loadButtonClicked() {
    tgtAssert(baseProperty_->get(), "No valid transfer function assigned");

    if(TransFuncIOHelperQt::loadTransferFunction(baseProperty_->get())) {
        baseProperty_->invalidate();
    }
}

    //-----------------------------
    // General functions
    //-----------------------------
void TransFuncPropertyEditorBase::toggleInteractionMode(bool on) {
    tgtAssert(initialized_, "is not inizialized!");
    QualityMode.requestQualityMode((on ? VoreenQualityMode::RQ_INTERACTIVE : VoreenQualityMode::RQ_DEFAULT), this);
}

} // namespace voreen
