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

#ifndef VRN_VOLUMEINPUTSWITCH_H
#define VRN_VOLUMEINPUTSWITCH_H

#include "voreen/core/processors/volumeprocessor.h"
#include "voreen/core/properties/intproperty.h"

namespace voreen {

class VRN_CORE_API VolumeInputSwitch : public VolumeProcessor {

public:
    VolumeInputSwitch();
    virtual Processor* create() const;

    virtual std::string getCategory() const   { return "Utility"; }
    virtual std::string getClassName() const  { return "VolumeInputSwitch";     }
    virtual CodeState getCodeState() const    { return CODE_STATE_TESTING;   }

    virtual bool isReady() const;

    virtual void process();

protected:

    virtual void setDescriptions() {
        setDescription("Routes the selected input volume to the outport (even if it does not contain any data).");
    }

private:

    VolumePort inport1_;
    VolumePort inport2_;
    VolumePort inport3_;
    VolumePort inport4_;

    VolumePort outport_;

    IntProperty inputSelect_;

    static const std::string loggerCat_; ///< category used in logging
};

}   //namespace

#endif // VRN_VOLUMEINPUTSWITCH_H
