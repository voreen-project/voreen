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

#include "plotdatapropertywidget.h"

#include "../../properties/plotdataproperty.h"
#include "../plotdatawidget.h"

#include "voreen/qt/widgets/voreentoolwindow.h"

#include "../../datastructures/plotdata.h"

#include <QPushButton>

namespace voreen {


PlotDataPropertyWidget::PlotDataPropertyWidget(PlotDataProperty* prop, QWidget* parent)
    : QPropertyWidgetWithToolWindow(prop, parent, true, false)
    , editBt_(new QPushButton(tr("Show Data")))
    , property_(prop)
    , plotDataWidget_(0)
    , updateFromPropertyFlag_(true)
{
    createToolWindow(Qt::LeftDockWidgetArea);

    addWidget(editBt_);
    connect(editBt_, SIGNAL(clicked()), this, SLOT(toggleWidgetVisibility()));

    QFontInfo fontInfo(font());
    editBt_->setFont(QFont(fontInfo.family(), QPropertyWidget::fontSize_));
}

void PlotDataPropertyWidget::updateFromPropertySlot() {
    updateFromPropertyFlag_ = true;
}

void PlotDataPropertyWidget::toggleWidgetVisibility() {
    tgtAssert(toolWindow_, "No editor window");
    if (updateFromPropertyFlag_) {
        updateFromPropertyFlag_ = false;
        if (plotDataWidget_)
            plotDataWidget_->updateFromPlotData();
    }
    toolWindow_->setVisible(!toolWindow_->isVisible());
}

void PlotDataPropertyWidget::customizeToolWindow() {
    toolWindow_->setAllowedAreas(static_cast<QFlag>(Qt::RightDockWidgetArea | Qt::BottomDockWidgetArea));
    toolWindow_->setFloating(true);
}

QWidget* PlotDataPropertyWidget::createToolWindowWidget() {
    delete plotDataWidget_;

    plotDataWidget_ = new PlotDataWidget(property_, parentWidget());
    plotDataWidget_->createWidgets();

    return plotDataWidget_;
}

Property* PlotDataPropertyWidget::getProperty() {
    return property_;
}

} // namespace
