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

#include "voreen/qt/widgets/property/processorpropertieswidget.h"

#include "voreen/core/processors/processor.h"
#include "voreen/core/ports/port.h"

namespace voreen {

ProcessorPropertiesWidget::ProcessorPropertiesWidget(Processor* processor, QWidget* parent,
                                                     bool expanded, bool userExpandable)
    : PropertyOwnerWidget(processor, parent, "",expanded, userExpandable, false)
{
    tgtAssert(processor, "null pointer passed");
}

std::vector<Property*>* ProcessorPropertiesWidget::createPropertyList() {
    Processor* proc = dynamic_cast<Processor*>(propertyOwner_);
    tgtAssert(proc, "no prcessor passed");
     std::vector<Property*>* propertyList =  new std::vector<Property*>(propertyOwner_->getProperties());
     for (size_t i=0; i<proc->getPorts().size(); i++) {
        const std::vector<Property*>& portProps =proc->getPorts().at(i)->getProperties();
        propertyList->insert(propertyList->end(), portProps.begin(), portProps.end());
    }
    return propertyList;
}

} // namespace
