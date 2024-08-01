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

#ifndef VRN_PLOTFUNCTIONSOURCE_H
#define VRN_PLOTFUNCTIONSOURCE_H

#include "../ports/plotport.h"
#include "voreen/core/processors/processor.h"

#include "voreen/core/properties/stringproperty.h"
#include "voreen/core/properties/buttonproperty.h"
#include "voreen/core/properties/optionproperty.h"
#include "voreen/core/properties/intproperty.h"

namespace voreen {

class PlotFunction;

class VRN_CORE_API PlotFunctionSource : public Processor {

public:
    PlotFunctionSource();
    ~PlotFunctionSource();
    virtual Processor* create() const;

    virtual std::string getCategory() const  { return "Data Source"; }
    virtual std::string getClassName() const { return "PlotFunctionSource"; }
    virtual CodeState getCodeState() const   { return CODE_STATE_TESTING; }

protected:
    virtual void setDescriptions() {
        setDescription("Creates a plot function from a string");
    }

    virtual void process();
    virtual void initialize();
    virtual void deinitialize();

private:
    void recalculate();
    void changeText();
    PlotFunction* readInputString();

    PlotPort outPort_;

    StringProperty functionInput_;
    StringProperty expressionNameInput_;
    IntOptionProperty expressionText_;

    StringProperty selfDescription_;
    IntProperty maxLength_;
    PlotFunction* pData_;

    ButtonProperty recalculate_;

    static const std::string loggerCat_;
};

} //namespace
#endif // VRN_PLOTFUNCTIONSOURCE_H

