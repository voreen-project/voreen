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

#include "plotdataprocessorwidget.h"

#include "../plotdatawidget.h"

#include "voreen/qt/widgets/voreentoolwindow.h"
#include "../../datastructures/plotdata.h"

#include <QPushButton>

namespace voreen {


PlotDataProcessorWidget::PlotDataProcessorWidget(QWidget* parent, Processor* proc, PlotDataProperty* data)
    : QProcessorWidget(proc, parent)
    , data_(data)
    , plotDataWidget_(0) {
    initialize();
}

PlotDataProcessorWidget::~PlotDataProcessorWidget() {
    delete plotDataWidget_;
}

void PlotDataProcessorWidget::initialize() {
    plotDataWidget_ = new PlotDataWidget(data_, parentWidget());
    plotDataWidget_->createWidgets();
}

void PlotDataProcessorWidget::updateFromProcessor() {
    if (plotDataWidget_)
        plotDataWidget_->updateFromPlotData();
}

} // namespace
