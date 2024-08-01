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

#include "randomwalkeranalyzerwidget.h"

#include <QCheckBox>
#include <QVBoxLayout>
#include <QHBoxLayout>

#include <QPolygonF>
#include <QPointF>
#include <QVector>
#include <QBrush>
#include <QPen>
#include <QGraphicsEllipseItem>
#include <QGraphicsSceneMouseEvent>
#include <QGraphicsPixmapItem>
#include <QApplication>
#include <QMessageBox>
#include <QPushButton>


namespace voreen {

using tgt::vec3;
using tgt::vec2;
using tgt::ivec4;

const std::string RandomWalkerAnalyzerWidget::loggerCat_("voreen.qt.RandomWalkerAnalyzerWidget");

RandomWalkerAnalyzerWidget::RandomWalkerAnalyzerWidget(RandomWalkerAnalyzer* processor, QWidget* parent) :
    QProcessorWidget(processor, parent),
    blockProcessorUpdates_(false),
    curFocusRegion_(-1)
{

    QVBoxLayout* mainLayout = new QVBoxLayout();
    QHBoxLayout* hBoxLayout = new QHBoxLayout();

    // table widget
    hBoxLayout = new QHBoxLayout();
    hBoxLayout->addSpacing(2);
    hBoxLayout->addSpacing(5);

    tableUncertaintyRegions_ = new QTableWidget(0, 3);
    tableUncertaintyRegions_->setSizePolicy(QSizePolicy::MinimumExpanding, QSizePolicy::MinimumExpanding);
    QStringList headers;
    headers.push_back("Volume [ml]");
    headers.push_back("Position");
    headers.push_back("Dimensions");
//    headers.push_back("URB");

    tableUncertaintyRegions_->setHorizontalHeaderLabels(headers);
    tableUncertaintyRegions_->setColumnWidth(0, 80);
    tableUncertaintyRegions_->setColumnWidth(1, 90);
    tableUncertaintyRegions_->setColumnWidth(2, 90);
    //tableUncertaintyRegions_->setColumnWidth(3, 80);

    tableUncertaintyRegions_->setSelectionBehavior(QAbstractItemView::SelectRows);
    tableUncertaintyRegions_->setSelectionMode(QAbstractItemView::SingleSelection);

    mainLayout->addWidget(tableUncertaintyRegions_, 2);
    hBoxLayout->addSpacing(5);
    hBoxLayout->addStretch();
    mainLayout->addLayout(hBoxLayout);

    computeButton_= new QPushButton("Compute");
    hBoxLayout->addWidget(computeButton_);
    zoomInButton_ = new QPushButton("Zoom In");
    hBoxLayout->addWidget(zoomInButton_);
    resetZoomButton_ = new QPushButton("Zoom Out");
    hBoxLayout->addWidget(resetZoomButton_);

    setLayout(mainLayout);

    connect(tableUncertaintyRegions_, SIGNAL(itemSelectionChanged()), this, SLOT(tableSelectionChanged()));
    connect(tableUncertaintyRegions_, SIGNAL(clicked(QModelIndex)), this, SLOT(tableSelectionChanged()));
    connect(computeButton_, SIGNAL(clicked()), this, SLOT(computeRandomWalker()));
    connect(zoomInButton_, SIGNAL(clicked()), this, SLOT(zoomIn()));
    connect(resetZoomButton_, SIGNAL(clicked()), this, SLOT(resetZoom()));

    processorNameChanged();
    updateTableWidget();
}

RandomWalkerAnalyzerWidget::~RandomWalkerAnalyzerWidget() {
}

void RandomWalkerAnalyzerWidget::updateTableWidget() {

    RandomWalkerAnalyzer* rwAnalyzer = dynamic_cast<RandomWalkerAnalyzer*>(processor_);
    tgtAssert(rwAnalyzer, "No processor");

    float sizeMapping = 0.99f;

    const std::vector<RandomWalkerAnalyzer::UncertaintyRegion>& uncertaintyRegions = rwAnalyzer->getUncertaintyRegions();

    tableUncertaintyRegions_->blockSignals(true);
    tableUncertaintyRegions_->clearContents();
    tableUncertaintyRegions_->setRowCount(static_cast<int>(uncertaintyRegions.size()));
    for (int i=0; i<static_cast<int>(uncertaintyRegions.size()); i++) {

        const RandomWalkerAnalyzer::UncertaintyRegion& region = uncertaintyRegions.at(i);

        QTableWidgetItem* itemSize = new QTableWidgetItem();
        QString textSize;
        textSize.setNum(region.size_*sizeMapping);
        itemSize->setText(textSize);
        tableUncertaintyRegions_->setItem(i, 0, itemSize);

        QTableWidgetItem* itemPosition = new QTableWidgetItem();
        std::ostringstream stream;
        stream << tgt::iround(region.centerOfMass_);
        itemPosition->setText(QString::fromStdString(stream.str()));
        tableUncertaintyRegions_->setItem(i, 1, itemPosition);

        QTableWidgetItem* itemDim = new QTableWidgetItem();
        stream.str("");
        stream.clear();
        stream << (region.urb_ - region.llf_);
        itemDim->setText(QString::fromStdString(stream.str()));
        tableUncertaintyRegions_->setItem(i, 2, itemDim);
    }

    tableUncertaintyRegions_->blockSignals(false);
}


void RandomWalkerAnalyzerWidget::tableSelectionChanged() {

    RandomWalkerAnalyzer* rwAnalyzer = dynamic_cast<RandomWalkerAnalyzer*>(processor_);
    tgtAssert(rwAnalyzer, "No processor");

    const std::vector<RandomWalkerAnalyzer::UncertaintyRegion>& uncertaintyRegions =
        rwAnalyzer->getUncertaintyRegions();

    size_t row = tableUncertaintyRegions_->currentIndex().row();
    if (row < uncertaintyRegions.size()) {
        //tgt::ivec3 center = tgt::iround(uncertaintyRegions.at(row).centerOfMass_);
        //blockProcessorUpdates_ = true;
        rwAnalyzer->resetZoom();
        rwAnalyzer->setFocusRegion(static_cast<int>(row));
        curFocusRegion_ = static_cast<int>(row);
        qApp->processEvents();
        zoomIn();
        //rwAnalyzer->forceUpdate();
    }
}

void RandomWalkerAnalyzerWidget::updateFromProcessor() {

    if (!blockProcessorUpdates_) {
        QProcessorWidget::updateFromProcessor();
        updateTableWidget();

        if (qApp->overrideCursor())
            qApp->restoreOverrideCursor();
        tableUncertaintyRegions_->setEnabled(true);
    }
    else {
        blockProcessorUpdates_ = false;
    }
}

void RandomWalkerAnalyzerWidget::computeRandomWalker() {
    RandomWalkerAnalyzer* rwAnalyzer = dynamic_cast<RandomWalkerAnalyzer*>(processor_);
    tgtAssert(rwAnalyzer, "No processor");

    tableUncertaintyRegions_->setEnabled(false);
    qApp->setOverrideCursor(Qt::WaitCursor);
    //rwAnalyzer->resetZoom();
    rwAnalyzer->computeRandomWalker();
}

void RandomWalkerAnalyzerWidget::resetZoom() {
    RandomWalkerAnalyzer* rwAnalyzer = dynamic_cast<RandomWalkerAnalyzer*>(processor_);
    tgtAssert(rwAnalyzer, "No processor");
    rwAnalyzer->resetZoomAnimated(20);
}

void RandomWalkerAnalyzerWidget::zoomIn() {

    RandomWalkerAnalyzer* rwAnalyzer = dynamic_cast<RandomWalkerAnalyzer*>(processor_);
    tgtAssert(rwAnalyzer, "No processor");

    const std::vector<RandomWalkerAnalyzer::UncertaintyRegion>& uncertaintyRegions =
        rwAnalyzer->getUncertaintyRegions();

    if (curFocusRegion_ >= 0 && curFocusRegion_ < static_cast<int>(uncertaintyRegions.size())) {
        rwAnalyzer->zoomOnFocusRegion(curFocusRegion_, 20);
    }
}

void RandomWalkerAnalyzerWidget::processorNameChanged() {
    setWindowTitle(QString::fromStdString(processor_->getGuiName()) + QString(" - Uncertain Regions"));
}



} // namespace voreen

