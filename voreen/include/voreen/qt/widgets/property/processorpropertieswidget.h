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

#ifndef VRN_PROCESSORPROPERTIESWIDGET_H
#define VRN_PROCESSORPROPERTIESWIDGET_H

#include "propertyownerwidget.h"

namespace voreen {

class ExpandableHeaderButton;
class Processor;
class PropertyWidgetFactory;
class QPropertyWidget;
class GroupPropertyWidget;

/**
 * Widget containing each processor's property widgets.
 */
class VRN_QT_API ProcessorPropertiesWidget : public PropertyOwnerWidget {
    Q_OBJECT
public:

    ProcessorPropertiesWidget(Processor* processor, QWidget* parent = 0,
                              bool expanded = true, bool userExpandable = false);

protected:
    /** Get all properties of the prcessor and the ports */
    virtual std::vector<Property*>* createPropertyList();
};

} // namespace

#endif // VRN_PROCESSORPROPERTIESWIDGET_H
