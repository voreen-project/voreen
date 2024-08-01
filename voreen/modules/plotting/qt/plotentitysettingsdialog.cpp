/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2024 University of Muenster, Germany,                        *
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

#include "plotentitysettingsdialog.h"

#include "../properties/plotentitiesproperty.h"

#include "voreen/core/voreenapplication.h"
#include "voreen/qt/widgets/property/colorpropertywidget.h"

#include <QLabel>
#include <QGridLayout>
#include <QTabWidget>
#include <QComboBox>
#include <QCheckBox>
#include <QPushButton>
#include <QDoubleSpinBox>
#include <QDialogButtonBox>
#include <QDir>
#include <QFileDialog>

namespace voreen {

PlotEntitySettingsDialog::PlotEntitySettingsDialog(const PlotEntitySettings& entitySettings, const PlotData* plotData, const PlotEntitiesProperty* prop, QWidget* parent)
    : QDialog(parent)
    , entitySettings_(entitySettings)
    , plotData_(plotData)
    , prop_(prop)
{
    tgtAssert(plotData_, "PlotEntitySettingsDialog(): plotData is NULL");
    setObjectName("PlotEntitySettingsDialog");
    createWidgets();
    updateWidgets();
}

static const int MAXIMUM_WIDGET_WIDTH = 200;

void PlotEntitySettingsDialog::createWidgets() {

    buttonBox_ = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel);
    cbMainColumn_ = new QComboBox();
    cbMainColumn_->setMaximumWidth(MAXIMUM_WIDGET_WIDTH);
    fillComboBox(cbMainColumn_);
    cbMainColumn_->setCurrentIndex(cbMainColumn_->findData(entitySettings_.getMainColumnIndex()));
    cclrFirst_ = new ColorSelectorWidget("");
    cclrFirst_->setColor(toQColor(entitySettings_.getFirstColor()));
    switch (entitySettings_.getEntity()) {
        case PlotEntitySettings::LINE: {
            //create widgets
            QVBoxLayout* layout = new QVBoxLayout();
            tabWidget_ = new QTabWidget(this);
            QWidget* lineTab = new QWidget(this);
            QWidget* candleStickTab = new QWidget(this);
            QGridLayout* gridLayout = new QGridLayout(lineTab);
            QGridLayout* gridLayout2 = new QGridLayout(candleStickTab);
            cbCandleTopColumn_ = new QComboBox();
            cbCandleTopColumn_->setMaximumWidth(MAXIMUM_WIDGET_WIDTH);
            cbCandleBottomColumn_ = new QComboBox();
            cbCandleBottomColumn_->setMaximumWidth(MAXIMUM_WIDGET_WIDTH);
            cbStickTopColumn_ = new QComboBox();
            cbStickTopColumn_->setMaximumWidth(MAXIMUM_WIDGET_WIDTH);
            cbStickBottomColumn_ = new QComboBox();
            cbStickBottomColumn_->setMaximumWidth(MAXIMUM_WIDGET_WIDTH);
            cbOptionalColumn_ = new QComboBox();
            cbOptionalColumn_->setMaximumWidth(MAXIMUM_WIDGET_WIDTH);
            cbStyle_ = new QComboBox();
            cbStyle_->setMaximumWidth(MAXIMUM_WIDGET_WIDTH);
            cbStyle_->addItem("Continuous"); cbStyle_->addItem("Dotted");  cbStyle_->addItem("Dashed");
            cbStyle_->setCurrentIndex(static_cast<int>(entitySettings_.getLineStyle()));
            chbSpline_ = new QCheckBox();
            chbErrorbar_ = new QCheckBox();
            cclrSecond_ = new ColorSelectorWidget("");
            cclrSecond_->setColor(toQColor(entitySettings_.getSecondColor()));
            cclrThird_ = new ColorSelectorWidget("");
            cclrThird_->setColor(cclrFirst_->getColor());
            cclrFourth_ = new ColorSelectorWidget("");
            cclrFourth_->setColor(cclrFirst_->getColor());
            fillComboBox(cbCandleTopColumn_);
            fillComboBox(cbCandleBottomColumn_);
            fillComboBox(cbStickTopColumn_);
            fillComboBox(cbStickBottomColumn_);
            fillComboBox(cbOptionalColumn_, true);
            cbOptionalColumn_->setCurrentIndex(
                cbOptionalColumn_->findData(entitySettings_.getOptionalColumnIndex()));
            cbCandleTopColumn_->setCurrentIndex(
                cbCandleTopColumn_->findData(entitySettings_.getCandleTopColumnIndex()));
            cbCandleBottomColumn_->setCurrentIndex(
                cbCandleBottomColumn_->findData(entitySettings_.getCandleBottomColumnIndex()));
            cbStickTopColumn_->setCurrentIndex(
                cbStickTopColumn_->findData(entitySettings_.getStickTopColumnIndex()));
            cbStickBottomColumn_->setCurrentIndex(
                cbStickBottomColumn_->findData(entitySettings_.getStickBottomColumnIndex()));
            chbSpline_->setChecked(entitySettings_.getSplineFlag());
            chbErrorbar_->setChecked(entitySettings_.getErrorbarFlag());
            entitySettings_.getCandleStickFlag() ? tabWidget_->setCurrentIndex(1)
                                                 : tabWidget_->setCurrentIndex(0);
            //add widgets
            gridLayout->addWidget(new QLabel(tr("Line")), 0, 0);
            gridLayout->addWidget(cbMainColumn_, 0, 1, 1, 2);
            gridLayout->addWidget(new QLabel(tr("Color")), 1, 0);
            gridLayout->addWidget(cclrFirst_, 1, 1, 1, 2);
            gridLayout->addWidget(new QLabel(tr("Style")), 2, 0);
            gridLayout->addWidget(cbStyle_, 2, 1, 1, 2);
            gridLayout->addWidget(new QLabel(tr("Error")), 3, 0);
            gridLayout->addWidget(cbOptionalColumn_, 3, 1, 1, 2);
            gridLayout->addWidget(new QLabel(tr("Error Color")), 4, 0);
            gridLayout->addWidget(cclrSecond_, 4, 1, 1, 2);
            gridLayout->addWidget(new QLabel(tr("Errorbar")), 5, 0);
            gridLayout->addWidget(chbErrorbar_, 5, 1, 1, 2);
            gridLayout->addWidget(new QLabel(tr("Spline")), 6, 0);
            gridLayout->addWidget(chbSpline_, 6, 1, 1, 2);

            gridLayout2->addWidget(new QLabel(tr("Candle Top")), 0, 0);
            gridLayout2->addWidget(cbCandleTopColumn_, 0, 1, 1, 2);
            gridLayout2->addWidget(new QLabel(tr("Candle Bottom")), 1, 0);
            gridLayout2->addWidget(cbCandleBottomColumn_, 1, 1, 1, 2);
            gridLayout2->addWidget(new QLabel(tr("Stick Top")), 2, 0);
            gridLayout2->addWidget(cbStickTopColumn_, 2, 1, 1, 2);
            gridLayout2->addWidget(new QLabel(tr("Stick Bottom")), 3, 0);
            gridLayout2->addWidget(cbStickBottomColumn_, 3, 1, 1, 2);
            gridLayout2->addWidget(new QLabel(tr("Stick Color")), 4, 0);
            gridLayout2->addWidget(cclrThird_, 4, 1, 1, 2);
            gridLayout2->addWidget(new QLabel(tr("Candle Color")), 5, 0);
            gridLayout2->addWidget(cclrFourth_, 5, 1, 1, 2);

            tabWidget_->addTab(lineTab, tr("Line"));
            tabWidget_->addTab(candleStickTab, tr("Candlestick"));
            layout->addWidget(tabWidget_);
            layout->addWidget(buttonBox_);
            //create connection
            connect(cbCandleTopColumn_, SIGNAL(activated(int)), this, SLOT(cbCandleTopColumnChanged(int)));
            connect(cbCandleBottomColumn_, SIGNAL(activated(int)), this, SLOT(cbCandleBottomColumnChanged(int)));
            connect(cbStickTopColumn_, SIGNAL(activated(int)), this, SLOT(cbStickTopColumnChanged(int)));
            connect(cbOptionalColumn_, SIGNAL(activated(int)), this, SLOT(cbOptionalColumnChanged(int)));
            connect(cbStickBottomColumn_, SIGNAL(activated(int)), this, SLOT(cbStickBottomColumnChanged(int)));
            connect(cbStyle_, SIGNAL(activated(int)), this, SLOT(cbStyleChanged(int)));
            connect(chbSpline_, SIGNAL(stateChanged(int)), this, SLOT(chbSplineChanged()));
            connect(chbErrorbar_, SIGNAL(stateChanged(int)), this, SLOT(chbErrorbarChanged()));
            connect(cclrSecond_, SIGNAL(colorChangedSignal(QColor)), this, SLOT(secondColorChanged(QColor)));
            connect(cclrThird_, SIGNAL(colorChangedSignal(QColor)), this, SLOT(firstColorChanged(QColor)));
            connect(cclrFourth_, SIGNAL(colorChangedSignal(QColor)), this, SLOT(secondColorChanged(QColor)));
            connect(tabWidget_, SIGNAL(currentChanged(int)), this, SLOT(tabChanged(int)));
            setLayout(layout);
            setWindowTitle(tr("Line settings"));
            break;
            }
        case PlotEntitySettings::BAR: {
            //create widgets
            QGridLayout* gridLayout = new QGridLayout();
            //add widgets
            gridLayout->addWidget(new QLabel(tr("Bars")), 0, 0);
            gridLayout->addWidget(cbMainColumn_, 0, 1, 1, 2);
            gridLayout->addWidget(new QLabel(tr("Color")), 1, 0);
            gridLayout->addWidget(cclrFirst_, 1, 1, 1, 2);
            gridLayout->addWidget(buttonBox_, 2, 1, 1 , 2);
            setLayout(gridLayout);
            setWindowTitle(tr("Bar settings"));
            break;
            }
        case PlotEntitySettings::SURFACE: {
            //create widgets
            QGridLayout* gridLayout = new QGridLayout();
            cbColorMap_ = new QComboBox();
            std::vector<std::string> colorMapLabels = ColorMap::getColorMapLabels();
            for (size_t i = 0; i < colorMapLabels.size(); ++i)
                cbColorMap_->addItem(QString::fromStdString(colorMapLabels.at(i)));
            cbColorMap_->setVisible(false);
            cclrSecond_ = new ColorSelectorWidget("");
            cclrSecond_->setColor(toQColor(entitySettings_.getSecondColor()));
            chbWireOnly_ = new QCheckBox();
            chbHeightMap_ = new QCheckBox();
            cbOptionalColumn_ = new QComboBox();
            fillComboBox(cbOptionalColumn_, true);
            cbOptionalColumn_->setCurrentIndex(
                cbOptionalColumn_->findData(entitySettings_.getOptionalColumnIndex()));
            chbWireOnly_->setChecked(entitySettings_.getWireOnlyFlag());
            chbHeightMap_->setChecked(entitySettings_.getHeightmapFlag());
            //add widgets
            gridLayout->addWidget(new QLabel(tr("Surface")), 0, 0);
            gridLayout->addWidget(cbMainColumn_, 0, 1, 1, 2);
            gridLayout->addWidget(new QLabel(tr("Color")), 1, 0);
            gridLayout->addWidget(cclrFirst_, 1, 1, 1, 2);
            gridLayout->addWidget(cbColorMap_, 1, 1, 1, 2);
            gridLayout->addWidget(new QLabel(tr("Wire Color")), 2, 0);
            gridLayout->addWidget(cclrSecond_, 2, 1, 1, 2);
            gridLayout->addWidget(new QLabel(tr("Wire Only")), 3, 0);
            gridLayout->addWidget(chbWireOnly_, 3, 1, 1, 2);
            gridLayout->addWidget(new QLabel(tr("Height Map")), 4, 0);
            gridLayout->addWidget(chbHeightMap_, 4, 1, 1, 2);
            gridLayout->addWidget(new QLabel(tr("Color Data")), 5, 0);
            gridLayout->addWidget(cbOptionalColumn_, 5, 1, 1, 2);
            gridLayout->addWidget(buttonBox_, 6, 1, 1 , 2);
            //create connection
            connect(cbColorMap_, SIGNAL(activated(int)), this, SLOT(cbColorMapChanged(int)));
            connect(cclrSecond_, SIGNAL(colorChangedSignal(QColor)), this, SLOT(secondColorChanged(QColor)));
            connect(chbWireOnly_, SIGNAL(stateChanged(int)), this, SLOT(chbWireOnlyChanged()));
            connect(chbHeightMap_, SIGNAL(stateChanged(int)), this, SLOT(chbHeightMapChanged()));
            connect(cbOptionalColumn_, SIGNAL(activated(int)), this, SLOT(cbOptionalColumnChanged(int)));
            setLayout(gridLayout);
            break;
            }
        case PlotEntitySettings::SCATTER: {
            //create widgets
            QGridLayout* gridLayout = new QGridLayout();
            cbColorMap_ = new QComboBox();
            std::vector<std::string> colorMapLabels = ColorMap::getColorMapLabels();
            for (size_t i = 0; i < colorMapLabels.size(); ++i)
                cbColorMap_->addItem(QString::fromStdString(colorMapLabels.at(i)));
            cbColorMap_->setVisible(false);
            cbColorMap_->setCurrentIndex(cbColorMap_->findText(
                QString::fromStdString(entitySettings_.getColorMap().toString())));
            cbStyle_ = new QComboBox();
            cbStyle_->addItem("Point"); cbStyle_->addItem("Circle"); cbStyle_->addItem("Triangle");  cbStyle_->addItem("Quad");
            sbMinGlyphSize_ = new QDoubleSpinBox();
            sbMaxGlyphSize_ = new QDoubleSpinBox();
            sbMinGlyphSize_->setDecimals(2);
            sbMaxGlyphSize_->setDecimals(2);
            sbMinGlyphSize_->setVisible(false);
            sbMinGlyphSize_->setValue(static_cast<double>(entitySettings_.getMinGlyphSize()));
            sbMaxGlyphSize_->setValue(static_cast<double>(entitySettings_.getMaxGlyphSize()));
            cbOptionalColumn_ = new QComboBox();
            cbSecondOptionalColumn_ = new QComboBox();
            btOpenFileDialog_ = new QPushButton();
            fillComboBox(cbOptionalColumn_, true);
            fillComboBox(cbSecondOptionalColumn_, true);
            cbOptionalColumn_->setCurrentIndex(
                cbOptionalColumn_->findData(entitySettings_.getOptionalColumnIndex()));
            cbSecondOptionalColumn_->setCurrentIndex(
                cbSecondOptionalColumn_->findData(entitySettings_.getSecondOptionalColumnIndex()));
            cbStyle_->setCurrentIndex(static_cast<int>(entitySettings_.getGlyphStyle()));
            updateBtOpenFileDialogText(entitySettings_.getTexturePath());
            //add widgets
            gridLayout->addWidget(new QLabel(tr("Data")), 0, 0);
            gridLayout->addWidget(cbMainColumn_, 0, 1, 1, 2);
            gridLayout->addWidget(new QLabel(tr("Color")), 1, 0);
            gridLayout->addWidget(cclrFirst_, 1, 1, 1, 2);
            gridLayout->addWidget(cbColorMap_, 1, 1, 1, 2);
            gridLayout->addWidget(new QLabel(tr("Glyph")), 2, 0);
            gridLayout->addWidget(cbStyle_, 2, 1, 1, 2);
            gridLayout->addWidget(new QLabel(tr("Glyph Size")), 3, 0);
            gridLayout->addWidget(sbMinGlyphSize_, 3, 1, 1, 1);
            gridLayout->addWidget(sbMaxGlyphSize_, 3, 2, 1, 1);
            gridLayout->addWidget(new QLabel(tr("Color Data")), 4, 0);
            gridLayout->addWidget(cbOptionalColumn_, 4, 1, 1, 2);
            gridLayout->addWidget(new QLabel(tr("Size Data")), 5, 0);
            gridLayout->addWidget(cbSecondOptionalColumn_, 5, 1, 1, 2);
            gridLayout->addWidget(new QLabel(tr("Load Icon")), 6, 0);
            gridLayout->addWidget(btOpenFileDialog_, 6, 1, 1, 2);
            gridLayout->addWidget(buttonBox_, 7, 1, 1 , 2);
            //create connection
            connect(cbColorMap_, SIGNAL(activated(int)), this, SLOT(cbColorMapChanged(int)));
            connect(cbStyle_, SIGNAL(activated(int)), this, SLOT(cbStyleChanged(int)));
            connect(sbMinGlyphSize_, SIGNAL(valueChanged(double)), this, SLOT(sbMinGlyphSizeChanged(double)));
            connect(sbMaxGlyphSize_, SIGNAL(valueChanged(double)), this, SLOT(sbMaxGlyphSizeChanged(double)));
            connect(cbOptionalColumn_, SIGNAL(activated(int)), this, SLOT(cbOptionalColumnChanged(int)));
            connect(cbSecondOptionalColumn_, SIGNAL(activated(int)), this, SLOT(cbSecondOptionalColumnChanged(int)));
            connect(btOpenFileDialog_, SIGNAL(clicked(void)), this, SLOT(btOpenFileDialogClicked(void)));
            setLayout(gridLayout);
            setWindowTitle(tr("Scatter settings"));
            break;
            }
        default:
            break;
    }
    connect(cbMainColumn_, SIGNAL(activated(int)), this, SLOT(cbMainColumnChanged(int)));
    connect(cclrFirst_, SIGNAL(colorChangedSignal(QColor)), this, SLOT(firstColorChanged(QColor)));
    connect(buttonBox_, SIGNAL(accepted()), this, SLOT(accept()));
    connect(buttonBox_, SIGNAL(rejected()), this, SLOT(reject()));
}

void PlotEntitySettingsDialog::updateWidgets() {
    switch (entitySettings_.getEntity()) {
        case PlotEntitySettings::LINE:
            if (entitySettings_.getCandleStickFlag())
                tabWidget_->setCurrentIndex(1);
            break;
        case PlotEntitySettings::BAR:
            break;
        case PlotEntitySettings::SURFACE:
            cclrFirst_->setVisible((entitySettings_.getOptionalColumnIndex() == -1));
            cbColorMap_->setVisible((entitySettings_.getOptionalColumnIndex() != -1));
            break;
        case PlotEntitySettings::SCATTER:
            cclrFirst_->setVisible((entitySettings_.getOptionalColumnIndex() == -1));
            cbColorMap_->setVisible((entitySettings_.getOptionalColumnIndex() != -1));
            sbMinGlyphSize_->setVisible((entitySettings_.getSecondOptionalColumnIndex() != -1));
            break;
        default :
            break;
    }
}

void PlotEntitySettingsDialog::fillComboBox(QComboBox* cb, bool emptyEntry) {
    cb->clear();
    if (emptyEntry)
        cb->addItem(tr(""), QVariant(-1));
    const PlotData* data = prop_->getPlotData();
    for (int i = 0; i < data->getColumnCount(); ++i) {
        if (data->getColumnType(i) == PlotBase::NUMBER        //number
            && prop_->getXColumnIndex() != i                         //not current x axis
            && (entitySettings_.getEntity() == PlotEntitySettings::LINE ||
                entitySettings_.getEntity() == PlotEntitySettings::BAR ||
                prop_->getYColumnIndex() != i))                       //not current y axis
            cb->addItem(QString::fromStdString(data->getColumnLabel(i)), QVariant(i));
    }
    cb->setCurrentIndex(0);
}

tgt::Color PlotEntitySettingsDialog::toTgtColor(QColor color) {
    return tgt::Color(color.redF(),color.greenF(),color.blueF(), color.alphaF());
}

QColor PlotEntitySettingsDialog::toQColor(tgt::Color color) {
    return QColor(static_cast<int>(color.r * 255), static_cast<int>(color.g * 255),
                  static_cast<int>(color.b * 255), static_cast<int>(color.a * 255));
}

void PlotEntitySettingsDialog::firstColorChanged(QColor col) {
    if(col.isValid()) {
        cclrFirst_->setColor(col);
        //third ccclr only exist in line mode
        if (entitySettings_.getEntity() == PlotEntitySettings::LINE)
            cclrThird_->setColor(col);
        entitySettings_.setFirstColor(toTgtColor(col));
        entitySettings_.setUseTextureFlag(false);
    }
}

void PlotEntitySettingsDialog::secondColorChanged(QColor col) {
    if(col.isValid()) {
        cclrSecond_->setColor(col);
        if (entitySettings_.getEntity() == PlotEntitySettings::LINE)
            cclrFourth_->setColor(col);
        entitySettings_.setSecondColor(toTgtColor(col));
    }
}

void PlotEntitySettingsDialog::cbMainColumnChanged(int index) {
    int plotDataIndex = cbMainColumn_->itemData(index).toInt();

    entitySettings_.setMainColumnIndex(plotDataIndex);

    if(plotData_ && (plotDataIndex < plotData_->getColumnCount()) && plotData_->hasColumnColorHint(plotDataIndex)) {
        // we have a color hint for this column
        QColor col = toQColor(plotData_->getColumnColorHint(plotDataIndex));
        cclrFirst_->setColor(col);

        switch (entitySettings_.getEntity()) {
            case PlotEntitySettings::LINE: cclrFourth_->setColor(col);
            case PlotEntitySettings::SURFACE: cclrThird_->setColor(col);
            break;
            case PlotEntitySettings::BAR:
            case PlotEntitySettings::SCATTER:
            break;
            case PlotEntitySettings::NONE:
            ; //Continue on.
        }

        entitySettings_.setFirstColor(toTgtColor(col));
        entitySettings_.setSecondColor(toTgtColor(col));
    }
}

void PlotEntitySettingsDialog::cbCandleTopColumnChanged(int index) {
    entitySettings_.setCandleTopColumnIndex(cbCandleTopColumn_->itemData(index).toInt());
}

void PlotEntitySettingsDialog::cbCandleBottomColumnChanged(int index) {
    entitySettings_.setCandleBottomColumnIndex(cbCandleBottomColumn_->itemData(index).toInt());
}

void PlotEntitySettingsDialog::cbStickTopColumnChanged(int index) {
    entitySettings_.setStickTopColumnIndex(cbStickTopColumn_->itemData(index).toInt());
}

void PlotEntitySettingsDialog::cbStickBottomColumnChanged(int index) {
    entitySettings_.setStickBottomColumnIndex(cbStickBottomColumn_->itemData(index).toInt());
}

void PlotEntitySettingsDialog::cbOptionalColumnChanged(int index) {
    entitySettings_.setOptionalColumnIndex(cbOptionalColumn_->itemData(index).toInt());
    entitySettings_.setUseTextureFlag(false);
    updateWidgets();
}

void PlotEntitySettingsDialog::cbSecondOptionalColumnChanged(int index) {
    entitySettings_.setSecondOptionalColumnIndex(cbSecondOptionalColumn_->itemData(index).toInt());
    updateWidgets();
}

void PlotEntitySettingsDialog::cbStyleChanged(int index) {
    if (entitySettings_.getEntity() == PlotEntitySettings::LINE)
        entitySettings_.setLineStyle(static_cast<PlotEntitySettings::LineStyle>(index));
    else //if SCATTER
        entitySettings_.setGlyphStyle(static_cast<PlotEntitySettings::GlyphStyle>(index));
}

void PlotEntitySettingsDialog::chbSplineChanged() {
    entitySettings_.setSplineFlag(chbSpline_->isChecked());
}

void PlotEntitySettingsDialog::chbErrorbarChanged() {
    entitySettings_.setErrorbarFlag(chbErrorbar_->isChecked());
}

void PlotEntitySettingsDialog::tabChanged(int index) {
    index == 0 ? entitySettings_.setCandleStickFlag(false) : entitySettings_.setCandleStickFlag(true);
}

void PlotEntitySettingsDialog::chbWireOnlyChanged() {
    entitySettings_.setWireOnlyFlag(chbWireOnly_->isChecked());
}

void PlotEntitySettingsDialog::chbHeightMapChanged() {
    entitySettings_.setHeightmapFlag(chbHeightMap_->isChecked());
}

void PlotEntitySettingsDialog::cbColorMapChanged(int index) {
    entitySettings_.setColorMap(ColorMap::createColorMap(index));
    entitySettings_.setUseTextureFlag(false);
}

void PlotEntitySettingsDialog::sbMinGlyphSizeChanged(double value) {
    entitySettings_.setMinGlyphSize(static_cast<float>(value));
}

void PlotEntitySettingsDialog::sbMaxGlyphSizeChanged(double value) {
    entitySettings_.setMaxGlyphSize(static_cast<float>(value));
}

void PlotEntitySettingsDialog::btOpenFileDialogClicked() {
    QString directory;
    // use directory of current property value if any, default directory otherwise
    if (!entitySettings_.getTexturePath().empty())
        directory = QString::fromStdString(entitySettings_.getTexturePath());
    else
        directory = QString::fromStdString(VoreenApplication::app()->getUserDataPath("textures"));
    const QString fileFilter = tr("*.jpg *.png *.bmp;;All files (*)");
    QString filename = QFileDialog::getOpenFileName(QWidget::parentWidget(),
        tr("Select File"), directory, fileFilter);
    if (!filename.isEmpty()) {
        entitySettings_.setTexturePath(filename.toStdString());
        entitySettings_.setUseTextureFlag(true);
        updateBtOpenFileDialogText(entitySettings_.getTexturePath());
    }
}

void PlotEntitySettingsDialog::updateBtOpenFileDialogText(const std::string& text) {
    if (!text.empty()) {
        size_t index = text.find_last_of('/');
        if (index == text.npos)
            index = text.find_last_of('\\');
        std::string endFilename = text;
        if (index != text.npos)
            endFilename = text.substr(index + 1, text.length());
        btOpenFileDialog_->setText(QString::fromStdString(endFilename));
    }
    else {
        btOpenFileDialog_->setText(tr("Select File"));
    }
    btOpenFileDialog_->update();
}

const PlotEntitySettings& PlotEntitySettingsDialog::getEntitySettings() const {
    return entitySettings_;
}

} // namespace
