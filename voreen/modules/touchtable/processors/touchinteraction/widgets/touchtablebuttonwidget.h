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

#ifndef VRN_TOUCHTABLEBUTTONWIDGET_H
#define VRN_TOUCHTABLEBUTTONWIDGET_H

#include "touchtablewidget.h"
#include "voreen/core/properties/buttonproperty.h"

namespace voreen {

/**
 * TouchTableWidget that has a button property which is triggered if a touch point associated with the widget is released.
 * This button property may be linked to other Voreen processors to trigger their functionality.
 */
class TouchTableButtonWidget : public TouchTableWidget {

public:

    TouchTableButtonWidget();

    virtual Processor* create() const;

    virtual std::string getClassName() const;

    virtual void pressed();

protected:

    virtual void setDescriptions();

    ButtonProperty buttonProp_;         ///< button property that is triggered by releasing a touch point associated with the widget

};

} // namespace

#endif // VRN_TOUCHTABLEBUTTONWIDGET_H
