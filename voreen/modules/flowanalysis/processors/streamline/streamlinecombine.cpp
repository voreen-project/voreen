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

#include "streamlinecombine.h"

#include "modules/flowanalysis/datastructures/streamlinelist.h"

namespace voreen {

const std::string StreamlineCombine::loggerCat_("flowanalysis.StreamlineCombine");

StreamlineCombine::StreamlineCombine()
    : Processor()
    // ports
    , leftInport_(Port::INPORT, "leftInport", "Left Streamline Input", false)
    , rightInport_(Port::INPORT, "rightInport", "Right Streamline Input", false)
    , outport_(Port::OUTPORT, "outport", "Streamline Outport", false)
    // properties
    , combineProp_("combineProp", "Combine Mode:")
{
    addPort(leftInport_);
    addPort(rightInport_);
    addPort(outport_);

    addProperty(combineProp_);
        combineProp_.addOption("left","Use Only Left",SCB_LEFT);
        combineProp_.addOption("right","Use Only Right",SCB_RIGHT);
        combineProp_.addOption("combine","Combine Lists",SCB_COMBINE);
}

StreamlineCombine::~StreamlineCombine() {
}

bool StreamlineCombine::isReady() const {
    if(!leftInport_.isReady() && !rightInport_.isReady()) {
        setNotReadyErrorMessage("No input connected");
        return false;
    }

    return true;
}

void StreamlineCombine::process() {
    StreamlineListBase* output = nullptr;
    switch(combineProp_.getValue()) {
    case SCB_LEFT:
        output = leftInport_.hasData() ? leftInport_.getData()->clone() : nullptr;
        break;
    case SCB_RIGHT:
        output = rightInport_.hasData() ? rightInport_.getData()->clone() : nullptr;
        break;
    case SCB_COMBINE:
        if(leftInport_.hasData()) {
            output = leftInport_.getData()->clone();
            if(rightInport_.hasData()) {
                output->addStreamlineList(*rightInport_.getData()->clone());
            }
            break;
        }
        output = rightInport_.getData()->clone();
        break;
    default:
        LERROR("unknown combine option!");
        tgtAssert(false,"should not get here");
        break;
    }

    outport_.setData(output);
}

}   // namespace
