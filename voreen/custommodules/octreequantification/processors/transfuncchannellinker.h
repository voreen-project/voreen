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

#ifndef VRN_TRANSFUNCCHANNELLINKER_H
#define VRN_TRANSFUNCCHANNELLINKER_H

#include "voreen/core/processors/processor.h"

#include "voreen/core/properties/transfunc/1d/1dkeys/transfunc1dkeysproperty.h"
#include "voreen/core/properties/intproperty.h"

namespace voreen {

    /**
     * This processor can be used to select one of four linked input transfer functions which is then set to a linkable "output" transfer function.
     */
class VRN_CORE_API TransFuncChannelLinker : public Processor {
public:
    TransFuncChannelLinker();

    virtual Processor* create() const;

    virtual std::string getClassName() const    { return "TransFuncChannelLinker";      }
    virtual std::string getCategory() const     { return "Linker";  }
    virtual CodeState getCodeState() const      { return CODE_STATE_EXPERIMENTAL;  }

    virtual bool isEndProcessor() const { return true; }

protected:

    virtual void process();

    virtual void setDescriptions() {
        setDescription("This processor can be used to select one of four linked input transfer functions which is then set to a linkable \"output\" transfer function.");
    }

private:

    virtual void onChangeCallback();

    /*
     * input transfer functions
     */
    TransFunc1DKeysProperty channel1TransFuncProp_;     ///< first input TF
    TransFunc1DKeysProperty channel2TransFuncProp_;     ///< second input TF
    TransFunc1DKeysProperty channel3TransFuncProp_;     ///< third input TF
    TransFunc1DKeysProperty channel4TransFuncProp_;     ///< fourth input TF

    IntProperty channelProp_;                           ///< selected input TF

    TransFunc1DKeysProperty outputTransFuncProp_;       //< output TF

    static const std::string loggerCat_;        ///< miau
};

} // namespace

#endif 
