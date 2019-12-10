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

#ifndef VRN_TOUCHTABLEBOOLWIDGET_H
#define VRN_TOUCHTABLEBOOLWIDGET_H

#include "touchtablewidget.h"
#include "voreen/core/properties/boolproperty.h"

namespace voreen {

/**
 * TouchScreenWidget that has a bool property which is triggered if a touch point associated with the widget is released.
 * This bool property may be linked to other processors to trigger functionality in Voreen.
 */
class TouchTableBoolWidget : public TouchTableWidget {
public:

    TouchTableBoolWidget();

    virtual Processor* create() const;
    virtual std::string getClassName() const;
    virtual void pressed();
protected:
    virtual void setDescriptions();

    //-----------
    // Callback
    //-----------
    /** Called from boolProperty changes.*/
    void updateWidget();

    BoolProperty boolProp_;         ///< bool property that is triggered by releasing a touchpoint associated with the widget and may be linked to other processors for functionality

};

} // namespace

#endif // VRN_TOUCHTALBEBOOLWIDGET_H
