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

#include "plotentitiespropertywidget.h"

#include "../plotentitysettingsdialog.h"
#include "../../properties/plotentitiesproperty.h"
#include "../../datastructures/plotdata.h"

#include "voreen/qt/widgets/voreentoolwindow.h"

#include <QPushButton>
#include <QGridLayout>
#include <QLabel>
#include <QListWidget>
#include <QComboBox>

namespace voreen {

PlotEntitiesPropertyWidget::PlotEntitiesPropertyWidget(PlotEntitiesProperty* prop, QWidget* parent)
    : QPropertyWidget(prop, parent)
    , btAdd_(new QPushButton(tr("Add")))
    , btAddAll_(new QPushButton(tr("Add All")))
    , btDelete_(new QPushButton(tr("Delete")))
    , btApplyColormap_(new QPushButton(tr("Apply")))
    , lwEntities_(new QListWidget())
    , cbXAxis_(new QComboBox())
    , cbYAxis_(new QComboBox())
    , cbColormap_(new QComboBox())
    , cbLoadStrategy_(new QComboBox())
    , property_(prop)
{
    QWidget* mainWidget = new QWidget();

    QGridLayout* layout = new QGridLayout();
    layout->addWidget(new QLabel(tr("x Axis:")), 0, 0);
    layout->addWidget(cbXAxis_, 0, 1, 1, 2);

    //row counter
    int row = 1;

    if (property_->getEntities() == PlotEntitySettings::SURFACE || property_->getEntities() == PlotEntitySettings::SCATTER) {
        layout->addWidget(new QLabel(tr("y Axis:")), row, 0);
        layout->addWidget(cbYAxis_, row, 1, 1, 2);
        connect(cbYAxis_, SIGNAL(activated(int)), this, SLOT(updateYAxis(int)));
        ++row;
    }
    layout->addWidget(new QLabel(tr("Colormap:")), row, 0);
    layout->addWidget(cbColormap_, row, 1);
    std::vector<std::string> colorMapLabels = ColorMap::getColorMapLabels();
    for (size_t i = 0; i < colorMapLabels.size(); ++i)
        cbColormap_->addItem(QString::fromStdString(colorMapLabels.at(i)));

    layout->addWidget(btApplyColormap_, row, 2);
    ++row;
    layout->addWidget(lwEntities_, row, 0, 1, 3);
    ++row;
    layout->addWidget(btAdd_, row, 0, 1, 1);
    layout->addWidget(btAddAll_, row, 1, 1, 1);
    layout->addWidget(btDelete_, row, 2, 1, 1);
    ++row;
    layout->addWidget(new QLabel(tr("Loading:")), row, 0);
    layout->addWidget(cbLoadStrategy_, row, 1, 1, 2);
    cbLoadStrategy_->addItem(QString::fromStdString("No  Data"),static_cast<int>(PlotEntitiesProperty::LS_NON));
    cbLoadStrategy_->addItem(QString::fromStdString("New Data"),static_cast<int>(PlotEntitiesProperty::LS_NEW));
    cbLoadStrategy_->addItem(QString::fromStdString("All Data"),static_cast<int>(PlotEntitiesProperty::LS_ALL));

    connect(cbXAxis_, SIGNAL(activated(int)), this, SLOT(updateXAxis(int)));
    connect(cbColormap_, SIGNAL(activated(int)), this, SLOT(updateColormap(int)));
    connect(cbLoadStrategy_, SIGNAL(activated(int)), this, SLOT(updateLoadStrategy(int)));
    connect(btAdd_, SIGNAL(clicked()), this, SLOT(buttonAddClicked()));
    connect(btAddAll_, SIGNAL(clicked()), this, SLOT(buttonAddAllClicked()));
    connect(btDelete_, SIGNAL(clicked()), this, SLOT(buttonDeleteClicked()));
    connect(btApplyColormap_, SIGNAL(clicked()), this, SLOT(buttonApplyColormapClicked()));
    connect(lwEntities_, SIGNAL(itemDoubleClicked(QListWidgetItem*)), this, SLOT(listViewDoubleClicked(QListWidgetItem*)));

    mainWidget->setLayout(layout);
    addWidget(mainWidget);

    updateFromPropertySlot();
}

void PlotEntitiesPropertyWidget::updateFromPropertySlot() {
    reCreateWidgets();
    //first we clear all entries
    lwEntities_->clear();
    QVariant x = static_cast<QVariant>(property_->getLoadStrategy());
    cbLoadStrategy_->setCurrentIndex(cbLoadStrategy_->findData(static_cast<int>(property_->getLoadStrategy())));
    if (property_->dataValid()) {
        cbXAxis_->setCurrentIndex(cbXAxis_->findData(property_->getXColumnIndex()));
        cbColormap_->setCurrentIndex(cbColormap_->findText(QString::fromStdString(property_->getColorMap().toString())));
        const std::vector<PlotEntitySettings>& plotEntitySettings = property_->get();
        //fill the list widget
        for (size_t i = 0; i < plotEntitySettings.size(); ++i) {
            QListWidgetItem* item;
            if (plotEntitySettings.at(i).getEntity() == PlotEntitySettings::LINE && plotEntitySettings.at(i).getCandleStickFlag())
                item = new QListWidgetItem(QString::fromStdString(
                "Candle Stick: " + property_->getPlotData()->getColumnLabel(plotEntitySettings.at(i).getCandleTopColumnIndex())
                + ", " + property_->getPlotData()->getColumnLabel(plotEntitySettings.at(i).getCandleBottomColumnIndex())));
            else
                item = new QListWidgetItem(QString::fromStdString(
                property_->getPlotData()->getColumnLabel(plotEntitySettings.at(i).getMainColumnIndex())));
            item->setForeground(QBrush(toQColor(plotEntitySettings.at(i).getFirstColor())));
            lwEntities_->addItem(item);
        }
        if (property_->getEntities() == PlotEntitySettings::SURFACE || property_->getEntities() == PlotEntitySettings::SCATTER)
            cbYAxis_->setCurrentIndex(cbYAxis_->findData(property_->getYColumnIndex()));
    }
}

void PlotEntitiesPropertyWidget::reCreateWidgets() {
    //first we clear all entries
    cbXAxis_->clear();
    cbYAxis_->clear();
    if (property_->dataValid()) {
        if (property_->getEntities() == PlotEntitySettings::LINE || property_->getEntities() == PlotEntitySettings::BAR) {
            for (int i = 0; i < property_->getPlotData()->getColumnCount(); ++i) {
                //x axis can be label column
                if (property_->getPlotData()->getColumnType(i) == PlotBase::STRING)
                    cbXAxis_->addItem(QString::fromStdString(property_->getPlotData()->getColumnLabel(i)),QVariant(i));
                //x axis can be number column (in the line case it has be a key column)
                else if (property_->getPlotData()->getColumnType(i) == PlotBase::NUMBER
                    && (i < property_->getPlotData()->getKeyColumnCount() || property_->getEntities() == PlotEntitySettings::BAR))
                    cbXAxis_->addItem(QString::fromStdString(property_->getPlotData()->getColumnLabel(i)),QVariant(i));
            }
        }
        else {
            // the y axis (scatter) can be not selected
            if (property_->getEntities() == PlotEntitySettings::SCATTER)
                cbYAxis_->addItem(tr(""),QVariant(-1));
            // fill the x and y axis combo boxes, must be number column
            for (int i = 0; i < property_->getPlotData()->getColumnCount(); ++i) {
                if (property_->getPlotData()->getColumnType(i) == PlotBase::NUMBER) {
                    cbXAxis_->addItem(QString::fromStdString(property_->getPlotData()->getColumnLabel(i)),QVariant(i));
                    cbYAxis_->addItem(QString::fromStdString(property_->getPlotData()->getColumnLabel(i)),QVariant(i));
                }
            }
        }
    }
}

void PlotEntitiesPropertyWidget::buttonAddClicked() {
    if (property_->dataValid()) {
        PlotEntitySettingsDialog dialog(property_->createEntitySettings(), property_->getPlotData(), property_);
        if (dialog.exec() == QDialog::Accepted)
            property_->setPlotEntitySettings(dialog.getEntitySettings(), static_cast<int>(property_->get().size()));
    }
}

void PlotEntitiesPropertyWidget::buttonAddAllClicked() {
    if (property_->dataValid()) {
        std::vector<PlotEntitySettings> all = property_->createAllEntitySettings();
        for (size_t i = 0; i < all.size(); ++i)
            property_->setPlotEntitySettings(all.at(i), static_cast<int>(property_->get().size()));
    }
}

void PlotEntitiesPropertyWidget::buttonDeleteClicked() {
    if (property_->dataValid()) {
        int idx = lwEntities_->row(lwEntities_->currentItem());
        if (idx > -1 && idx < static_cast<int>(lwEntities_->count()))
            property_->deletePlotEntitySettings(idx);
    }
}

void PlotEntitiesPropertyWidget::buttonApplyColormapClicked() {
    if (disconnected_)
        return;
    property_->applyColormap();
}

void PlotEntitiesPropertyWidget::updateXAxis(int index) {
    if (disconnected_)
        return;
    property_->setXColumnIndex(cbXAxis_->itemData(index).toInt());
}

void PlotEntitiesPropertyWidget::updateYAxis(int index) {
    if (disconnected_)
        return;
    property_->setYColumnIndex(cbYAxis_->itemData(index).toInt());
}

void PlotEntitiesPropertyWidget::updateLoadStrategy(int index) {
    if (disconnected_)
        return;
    property_->setLoadStrategy(static_cast<PlotEntitiesProperty::loadStrategy>(cbLoadStrategy_->itemData(index).toInt()));
}

void PlotEntitiesPropertyWidget::updateColormap(int index) {
    if (disconnected_)
        return;
    property_->setColorMap(ColorMap::createColorMap(index));
}

void PlotEntitiesPropertyWidget::listViewDoubleClicked(QListWidgetItem*  item) {
    int row = lwEntities_->row(item);
    PlotEntitySettingsDialog dialog(property_->get().at(row), property_->getPlotData(), property_);
    if (dialog.exec() == QDialog::Accepted)
        property_->setPlotEntitySettings(dialog.getEntitySettings(), row);
}

tgt::Color PlotEntitiesPropertyWidget::toTgtColor(QColor color) {
    return tgt::Color(color.redF(),color.greenF(),color.blueF(), color.alphaF());
}

QColor PlotEntitiesPropertyWidget::toQColor(tgt::Color color) {
    return QColor(static_cast<int>(color.r * 255), static_cast<int>(color.g * 255),
                  static_cast<int>(color.b * 255), static_cast<int>(color.a * 255));
}

} // namespace
