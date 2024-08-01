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

#include "volumeinputswitch.h"

namespace voreen {

const std::string VolumeInputSwitch::loggerCat_("voreen.base.VolumeInputSwitch");

VolumeInputSwitch::VolumeInputSwitch()
    : VolumeProcessor()
    , inport1_(Port::INPORT, "volume.input1", "Volume Input 1")
    , inport2_(Port::INPORT, "volume.input2", "Volume Input 2")
    , inport3_(Port::INPORT, "volume.input3", "Volume Input 3")
    , inport4_(Port::INPORT, "volume.input4", "Volume Input 4")
    , outport_(Port::OUTPORT, "volume.output", "Volume Output")
    , inputSelect_("inputselect", "Select Input", 1, 1, 4)
{
    addPort(inport1_);
    addPort(inport2_);
    addPort(inport3_);
    addPort(inport4_);

    addPort(outport_);

    addProperty(inputSelect_);
}

Processor* VolumeInputSwitch::create() const {
    return new VolumeInputSwitch();
}

bool VolumeInputSwitch::isReady() const {
    if (!isInitialized())
        return false;
    else
        return true;
}
void VolumeInputSwitch::process() {
    switch (inputSelect_.get()) {
        case 1:
            {
                if (inport1_.getData())
                    outport_.setData(inport1_.getData(), false);
                else
                    outport_.setData(0);
                break;
            }
        case 2:
            {
                if (inport2_.getData())
                    outport_.setData(inport2_.getData(), false);
                else
                    outport_.setData(0);
                break;
            }
        case 3:
            {
                if (inport3_.getData())
                    outport_.setData(inport3_.getData(), false);
                else
                    outport_.setData(0);
                break;
            }
        case 4:
            {
                if (inport4_.getData())
                    outport_.setData(inport4_.getData(), false);
                else
                    outport_.setData(0);
                break;
            }
        default:
            outport_.setData(0);
    }
}

}   // namespace
