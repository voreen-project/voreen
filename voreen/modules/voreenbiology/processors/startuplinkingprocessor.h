/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany,                        *
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

#ifndef VRN_STARTUPLINKINGPROCESSOR_H
#define VRN_STARTUPLINKINGPROCESSOR_H

#include "voreen/core/processors/processor.h"
#include "voreen/core/datastructures/volume/volumeslicehelper.h"

#include "voreen/core/properties/optionproperty.h"
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/intproperty.h"

namespace voreen {


class StartUpLinkingProcessor : public Processor {
public:
    StartUpLinkingProcessor();
    virtual ~StartUpLinkingProcessor();
    virtual Processor* create() const { return new StartUpLinkingProcessor(); }

    virtual std::string getClassName() const    { return "StartUpLinkingProcessor";  }
    virtual std::string getCategory() const     { return "Utility";         }
    virtual CodeState getCodeState() const      { return CODE_STATE_STABLE; }
    virtual bool isUtility() const              { return true;              }

protected:
    virtual void setDescriptions() {
        setDescription("Used in the startup.vws workspace of VoreenBiology. Handles the property switching between different processors.");
    }
    virtual void process(){/*nothing to do here*/}
    virtual void initialize();

    void onGlobalAlignmentChange();
    void onGlobalNumberChange();

    void onXNumberChange();
    void onYNumberChange();
    void onZNumberChange();

private:
    //property gui output
    OptionProperty<SliceAlignment> globalSliceAlignment_;
    IntProperty globalSliceNumber_;
    //properties to be linked
        //compositors
    StringOptionProperty compositingMode_Y_X_;
    StringOptionProperty compositingMode_Z_YX_;
        //slice pos
            //x slice
    BoolProperty renderXSlice_;
    IntProperty xSliceIndexProp_;
            //y slice
    BoolProperty renderYSlice_;
    IntProperty ySliceIndexProp_;
            //z slice
    BoolProperty renderZSlice_;
    IntProperty zSliceIndexProp_;
};

}   // namespace

#endif
