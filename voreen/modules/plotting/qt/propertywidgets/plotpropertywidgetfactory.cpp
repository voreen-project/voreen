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

#include "plotpropertywidgetfactory.h"

#include "../../properties/colormapproperty.h"
#include "../../properties/plotdataproperty.h"
#include "../../properties/plotentitiesproperty.h"
#include "../../properties/plotpredicateproperty.h"
#include "../../properties/plotselectionproperty.h"

#include "colormappropertywidget.h"
#include "plotdatapropertywidget.h"
#include "plotentitiespropertywidget.h"
#include "plotpredicatepropertywidget.h"
#include "plotselectionpropertywidget.h"

namespace voreen {

PropertyWidget* PlotPropertyWidgetFactory::createAssociatedWidget(Property* prop) const {

    if (dynamic_cast<ColorMapProperty*>(prop))
        return new ColorMapPropertyWidget(static_cast<ColorMapProperty*>(prop), 0);

    if (dynamic_cast<PlotPredicateProperty*>(prop))
        return new PlotPredicatePropertyWidget(static_cast<PlotPredicateProperty*>(prop), 0);

    if (dynamic_cast<PlotDataProperty*>(prop))
        return new PlotDataPropertyWidget(static_cast<PlotDataProperty*>(prop), 0);

    if (dynamic_cast<PlotSelectionProperty*>(prop))
        return new PlotSelectionPropertyWidget(static_cast<PlotSelectionProperty*>(prop), 0);

    if (dynamic_cast<PlotEntitiesProperty*>(prop))
        return new PlotEntitiesPropertyWidget(static_cast<PlotEntitiesProperty*>(prop), 0);

    return 0;
}

} // namespace voreen
