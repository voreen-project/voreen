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

#include "plotpredicatepropertywidget.h"

#include "../../datastructures/plotpredicate.h"
#include "../../properties/plotpredicateproperty.h"

#include <QListWidget>
#include <sstream>

namespace voreen {

PlotPredicatePropertyWidget::PlotPredicatePropertyWidget(PlotPredicateProperty* prop, QWidget* parent)
    : QPropertyWidget(prop, parent)
    , property_(prop)
    , lWidget_(new QListWidget)
{
    updateFromPropertySlot();
    QFontInfo fontInfo(font());
    lWidget_->setFont(QFont(fontInfo.family(), QPropertyWidget::fontSize_));
    addWidget(lWidget_);
    //connect(lWidget_, SIGNAL(currentIndexChanged(int)), this, SLOT(setProperty(int)));
    //connect(lWidget_, SIGNAL(currentIndexChanged(int)), this, SIGNAL(widgetChanged()));
}

void PlotPredicatePropertyWidget::updateFromPropertySlot() {
    lWidget_->blockSignals(true);
    lWidget_->clear();

    // build combo box from descriptions
    std::vector<std::pair<int, PlotPredicate*> >::const_iterator it;
    for (it = property_->get().begin(); it < property_->get().end(); ++it) {
        std::stringstream ss;
        ss << (*it).first << ": " << (*it).second->toString();
        lWidget_->addItem(QString::fromStdString(ss.str()));
    }

    lWidget_->blockSignals(false);
}


} // namespace
