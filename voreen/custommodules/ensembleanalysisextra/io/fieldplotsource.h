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

#ifndef VRN_FIELDPLOTSOURCE_H
#define VRN_FIELDPLOTSOURCE_H

#include "voreen/core/processors/processor.h"

#include "voreen/core/properties/filedialogproperty.h"
#include "voreen/core/properties/buttonproperty.h"
#include "voreen/core/properties/boolproperty.h"

#include "../ports/fieldplotdataport.h"

namespace voreen {

class VRN_CORE_API FieldPlotSource : public Processor {
public:
    FieldPlotSource();

    virtual Processor* create() const         { return new FieldPlotSource(); }
    virtual std::string getClassName() const  { return "FieldPlotSource"; }
    virtual std::string getCategory() const   { return "Input"; }
    virtual CodeState getCodeState() const    { return CODE_STATE_EXPERIMENTAL; }

protected:

    virtual void process();
    virtual void invalidate(int inv = 1);

private:

    void loadFieldPlot();

    FieldPlotDataPort outport_;

    FileDialogProperty filenameProp_;               ///< determines the name of the saved file
    ButtonProperty loadButton_;                     ///< triggers a load

    bool loadFieldPlot_;          ///< used to determine, if process should save or not

    static const std::string loggerCat_;
};

}   //namespace

#endif
