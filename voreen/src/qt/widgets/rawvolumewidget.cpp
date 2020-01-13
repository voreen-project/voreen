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

#include <QWidget>
#include <QComboBox>
#include <QDialog>
#include <QDoubleSpinBox>
#include <QGroupBox>
#include <QGridLayout>
#include <QLabel>
#include <QPushButton>
#include <QResizeEvent>
#include <QSpinBox>
#include <QVBoxLayout>

#include "voreen/core/voreenapplication.h"
#include "voreen/core/voreenmodule.h"
#include "voreen/core/properties/propertywidget.h"
#include "voreen/core/properties/matrixproperty.h"
#include "voreen/core/io/volumeserializerpopulator.h"

#include "voreen/qt/widgets/rawvolumewidget.h"
#include "voreen/qt/widgets/property/qpropertywidget.h"

namespace voreen {

RawVolumeWidget::RawVolumeWidget(QWidget* parent, const QString& filename, std::string& objectModel, std::string& format, int& numFrames,
                                 tgt::ivec3& dim, tgt::vec3& spacing, int& headerSkip, bool& bigEndian, tgt::mat4& trafoMat, int fixedZDim)
    : QDialog(parent)
    , objectModel_(objectModel)
    , format_(format)
    , numFrames_(numFrames)
    , dim_(dim)
    , spacing_(spacing)
    , trafoMat_(trafoMat)
    , headerSkip_(headerSkip)
    , bigEndian_(bigEndian)
{
    trafoMatProp_ = new FloatMat4Property("trafoMat", "trafoMat", trafoMat_);
    QVBoxLayout* mainLayout = new QVBoxLayout(this);
    QHBoxLayout* datatypeModelLayout = new QHBoxLayout();
    QVBoxLayout* datatypeBoxLayout = new QVBoxLayout();
    QVBoxLayout* objectModelBoxLayout = new QVBoxLayout();
    QVBoxLayout* numTimeFramesBoxLayout = new QVBoxLayout();
    QHBoxLayout* headerSkipBoxLayout = new QHBoxLayout();
    QHBoxLayout* endiannessBoxLayout_ = new QHBoxLayout();
    QHBoxLayout* headerSkipEndiannessBoxLayout = new QHBoxLayout();

    QGridLayout* dimensionBoxLayout = new QGridLayout();
    QGridLayout* spacingBoxLayout = new QGridLayout();

    headerSkipSpin_ = new QSpinBox();
    headerSkipSpin_->setValue(0);
    headerSkipSpin_->setMaximum(1<<20); // 1 MB

    endiannessCombo_ = new QComboBox();
    endiannessCombo_->addItem("Little Endian", false);
    endiannessCombo_->addItem("Big Endian", true);

    QGroupBox* datatypeBox = new QGroupBox(this);
    QGroupBox* objectModelBox = new QGroupBox(this);
    QGroupBox* numTimeFramesBox = new QGroupBox(this);
    QGroupBox* headerSkipBox = new QGroupBox(this);
    QGroupBox* endiannessBox_ = new QGroupBox(this);
    QGroupBox* dimensionBox = new QGroupBox(this);
    QGroupBox* spacingBox = new QGroupBox(this);

    QLabel* filenameLabel = new QLabel(filename);
    datatypeComboBox_ = new QComboBox();

    //datatypeComboBox_->addItems(generateVolumeTypeList());
    datatypeComboBox_->addItem("Unsigned Int 8-bit",    "UCHAR");
    datatypeComboBox_->addItem("Unsigned Int 12-bit",   "USHORT_12");
    datatypeComboBox_->addItem("Unsigned Int 16-bit",   "USHORT");
    datatypeComboBox_->addItem("Unsigned Int 32-bit",   "UINT");
    datatypeComboBox_->addItem("Integer 8-bit",         "CHAR");
    datatypeComboBox_->addItem("Integer 16-bit",        "SHORT");
    datatypeComboBox_->addItem("Integer 32-bit",        "INT");
    datatypeComboBox_->addItem("Float 32-bit",          "FLOAT");

    objectModelComboBox_ = new QComboBox();
    objectModelComboBox_->addItem("Intensity",              "I");
    objectModelComboBox_->addItem("RGB",                    "RGB");
    objectModelComboBox_->addItem("RGBA",                   "RGBA");
    objectModelComboBox_->addItem("Tensor Upper Diagonal",  "TENSOR_UP");
    objectModelComboBox_->addItem("Tensor Lower Diagonal",  "TENSOR_LOW");
    objectModelComboBox_->addItem("Tensor Diagonal Order",  "TENSOR_DIAG");
    objectModelComboBox_->addItem("Tensor Upper Diagonal Side-by-Side",  "TENSOR_FUSION_UP");
    objectModelComboBox_->addItem("Tensor Lower Diagonal Side-by-Side",  "TENSOR_FUSION_LOW");
    objectModelComboBox_->addItem("Tensor Diagonal Order Side-by-Side",  "TENSOR_FUSION_DIAG");

    numTimeFramesSpin_ = new QSpinBox();
    numTimeFramesSpin_->setMinimumWidth(60);
    numTimeFramesSpin_->setMinimum(1);
    numTimeFramesSpin_->setMaximum(512);


    objectModelBox->setTitle("Object Model");
    objectModelBox->setLayout(objectModelBoxLayout);
    objectModelBoxLayout->addWidget(objectModelComboBox_);

    datatypeBoxLayout->addWidget(datatypeComboBox_);
    datatypeBox->setLayout(datatypeBoxLayout);
    datatypeBox->setTitle("Data Type");

    numTimeFramesBox->setTitle("Number of Frames");
    numTimeFramesBox->setLayout(numTimeFramesBoxLayout);
    numTimeFramesBoxLayout->addWidget(numTimeFramesSpin_);

    headerSkipSpin_->setMinimumWidth(60);
    headerSkipBox->setTitle("Header Skip");
    headerSkipBoxLayout->addWidget(headerSkipSpin_);
    headerSkipBoxLayout->addWidget(new QLabel("Bytes"));
    headerSkipBoxLayout->addStretch();
    headerSkipBox->setLayout(headerSkipBoxLayout);

    endiannessCombo_->setMinimumWidth(80);
    endiannessBox_->setTitle("Byte Order");
    endiannessBoxLayout_->addWidget(endiannessCombo_);
    endiannessBox_->setLayout(endiannessBoxLayout_);

    QLabel* xText = new QLabel();
    xText->setText("x");
    QLabel* yText = new QLabel();
    yText->setText("y");
    QLabel* zText = new QLabel();
    zText->setText("z");
    QLabel* xText1 = new QLabel();
    xText1->setText("x");
    QLabel* yText1 = new QLabel();
    yText1->setText("y");
    QLabel* zText1 = new QLabel();
    zText1->setText("z");

    xDimension_ = new QSpinBox();
    yDimension_ = new QSpinBox();
    zDimension_ = new QSpinBox();
    xDimension_->setMinimum(1);
    yDimension_->setMinimum(1);
    zDimension_->setMinimum(1);
    xDimension_->setMaximum(8192);
    yDimension_->setMaximum(8192);
    zDimension_->setMaximum(8192);
    xDimension_->setValue(128);
    yDimension_->setValue(128);
    if (fixedZDim != -1) {
        zDimension_->setValue(fixedZDim);
        zDimension_->setEnabled(false);
    }
    else {
        zDimension_->setValue(128);
    }

    dimensionBoxLayout->addWidget(xText,0,0,1,1);
    dimensionBoxLayout->addWidget(yText,0,1,1,1);
    dimensionBoxLayout->addWidget(zText,0,2,1,1);
    dimensionBoxLayout->addWidget(xDimension_,1,0,1,1);
    dimensionBoxLayout->addWidget(yDimension_,1,1,1,1);
    dimensionBoxLayout->addWidget(zDimension_,1,2,1,1);
    dimensionBox->setLayout(dimensionBoxLayout);
    dimensionBox->setTitle("Dimensions");

    xSpacing_ = new QDoubleSpinBox();
    ySpacing_ = new QDoubleSpinBox();
    zSpacing_ = new QDoubleSpinBox();
    xSpacing_->setLocale(QLocale::English);
    ySpacing_->setLocale(QLocale::English);
    zSpacing_->setLocale(QLocale::English);
    xSpacing_->setDecimals(6);
    ySpacing_->setDecimals(6);
    zSpacing_->setDecimals(6);
    xSpacing_->setSingleStep(0.1);
    ySpacing_->setSingleStep(0.1);
    zSpacing_->setSingleStep(0.1);
    xSpacing_->setMinimum(0.000001);
    ySpacing_->setMinimum(0.000001);
    zSpacing_->setMinimum(0.000001);
    xSpacing_->setMaximum(100.0);
    ySpacing_->setMaximum(100.0);
    zSpacing_->setMaximum(100.0);
    xSpacing_->setMaximumWidth(150);
    ySpacing_->setMaximumWidth(150);
    zSpacing_->setMaximumWidth(150);
    xSpacing_->setValue(1.0);
    ySpacing_->setValue(1.0);
    zSpacing_->setValue(1.0);

    spacingBoxLayout->addWidget(xText1,0,0,1,1);
    spacingBoxLayout->addWidget(yText1,0,1,1,1);
    spacingBoxLayout->addWidget(zText1,0,2,1,1);
    spacingBoxLayout->addWidget(xSpacing_,1,0,1,1);
    spacingBoxLayout->addWidget(ySpacing_,1,1,1,1);
    spacingBoxLayout->addWidget(zSpacing_,1,2,1,1);

    spacingBox->setLayout(spacingBoxLayout);
    spacingBox->setTitle("Spacing");

    //trafoMat_ = new FloatMat4Property("TrafoMat", "Trafomat", tgt::mat4(1.0f));
    QHBoxLayout* matrixLayout = new QHBoxLayout();
    QLabel* trafoLabel = new QLabel("Transformation matrix", this);
    matrixLayout->addWidget(trafoLabel);
    trafoMatWidget_ = 0;
    const std::vector<VoreenModule*>& modules = VoreenApplication::app()->getModules();
    for (size_t m=0; m<modules.size(); m++) {
        PropertyWidget* propWidget = VoreenApplication::app()->createPropertyWidget(trafoMatProp_);
        if (propWidget) {
            trafoMatWidget_ = dynamic_cast<QPropertyWidget*>(propWidget);
            if (!trafoMatWidget_) {
                LERRORC("voreen.qt.RawVolumeWidget", "generated trafoMatProp_ widget is not of type QPropertyWidget");
                return;
            }
        }
    }
    tgtAssert(trafoMatWidget_, "trafoMatWidget_ is 0");
    connect(trafoMatWidget_, SIGNAL(valueModifiedByUser()), this, SLOT(updateValues()));
    matrixLayout->addWidget(trafoMatWidget_);
    mainLayout->addWidget(filenameLabel);

    datatypeModelLayout->addWidget(datatypeBox, 4);
    datatypeModelLayout->addWidget(objectModelBox, 3);
    datatypeModelLayout->addWidget(numTimeFramesBox, 2);
    mainLayout->addLayout(datatypeModelLayout);

    headerSkipEndiannessBoxLayout->addWidget(headerSkipBox, 4);
    headerSkipEndiannessBoxLayout->addWidget(endiannessBox_, 3);
    mainLayout->addLayout(headerSkipEndiannessBoxLayout);

    mainLayout->addWidget(dimensionBox);
    mainLayout->addWidget(spacingBox);

    mainLayout->addLayout(matrixLayout);

    connect(datatypeComboBox_, SIGNAL(currentIndexChanged(int)), this, SLOT(updateValues()));
    connect(objectModelComboBox_, SIGNAL(currentIndexChanged(int)), this, SLOT(updateValues()));
    connect(numTimeFramesSpin_, SIGNAL(valueChanged(int)), this, SLOT(updateValues()));
    connect(headerSkipSpin_, SIGNAL(valueChanged(int)), this, SLOT(updateValues()));
    connect(endiannessCombo_, SIGNAL(currentIndexChanged(int)), this, SLOT(updateValues()));

    connect(xDimension_, SIGNAL(valueChanged(int)), this, SLOT(updateValues()));
    connect(yDimension_, SIGNAL(valueChanged(int)), this, SLOT(updateValues()));
    connect(zDimension_, SIGNAL(valueChanged(int)), this, SLOT(updateValues()));

    connect(xSpacing_, SIGNAL(valueChanged(double)), this, SLOT(updateValues()));
    connect(ySpacing_, SIGNAL(valueChanged(double)), this, SLOT(updateValues()));
    connect(zSpacing_, SIGNAL(valueChanged(double)), this, SLOT(updateValues()));

    QPushButton* loadButton = new QPushButton();
    loadButton->setText("Load Volume");
    mainLayout->addWidget(loadButton);

    connect(loadButton, SIGNAL(clicked()), this, SLOT(accept()) );
    setWindowTitle(QString::fromStdString("Raw Volume Information"));

    setFixedSize(sizeHint());

    updateValues();
    setAttribute(Qt::WA_DeleteOnClose);
}

RawVolumeWidget::~RawVolumeWidget() {
    trafoMatWidget_->hide();
}

void RawVolumeWidget::resizeEvent(QResizeEvent*) {
}

void RawVolumeWidget::updateValues() {

    objectModel_ = objectModelComboBox_->itemData(objectModelComboBox_->currentIndex()).toString().toStdString();
    format_ = datatypeComboBox_->itemData(datatypeComboBox_->currentIndex()).toString().toStdString();
    numFrames_ = numTimeFramesSpin_->value();
    headerSkip_ = headerSkipSpin_->value();
    bigEndian_ = endiannessCombo_->itemData(endiannessCombo_->currentIndex()).toBool();

    dim_.x = xDimension_->value();
    dim_.y = yDimension_->value();
    dim_.z = zDimension_->value();
    spacing_.x = xSpacing_->value();
    spacing_.y = ySpacing_->value();
    spacing_.z = zSpacing_->value();
    trafoMat_ = trafoMatProp_->get();
}

} //namespace voreen
