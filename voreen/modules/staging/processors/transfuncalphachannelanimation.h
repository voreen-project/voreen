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

#ifndef VRN_TRANSFUNCALPHACHANNELANIMATION_H
#define VRN_TRANSFUNCALPHACHANNELANIMATION_H

#include "voreen/core/processors/processor.h"

#include "voreen/core/properties/transfunc/1d/1dkeys/transfunc1dkeysproperty.h"
#include "voreen/core/properties/optionproperty.h"

namespace voreen {

    /**
     * Temporary workaround to be able to switch transfuncs on and off (i.e. set alpha to transparent) in the animation by manipulating option properties.
     */
class VRN_CORE_API TransFuncAlphaChannelAnimation : public Processor {
public:
    TransFuncAlphaChannelAnimation();

    virtual Processor* create() const;

    virtual std::string getClassName() const    { return "TransFuncAlphaChannelAnimation";      }
    virtual std::string getCategory() const     { return "Linker";  }
    virtual CodeState getCodeState() const      { return CODE_STATE_EXPERIMENTAL;  }

    virtual bool isReady() const;

protected:

    virtual void process() { 
        // do nothing 
    }

    virtual void setDescriptions() {
        setDescription("Use this processor to set the alpha mode of a linked transfer functions using the option property in the animation.");
    }

private:


    // callbacks for TF changes which sets the option properties
    void onTransFuncChange1();       
    void onTransFuncChange2();
    void onTransFuncChange3();
    void onTransFuncChange4();
    
    // callbacks for option property changes which sets the TF
    void onModeChange1();         
    void onModeChange2();
    void onModeChange3();
    void onModeChange4();

    // methods for converting between TF alpha mode and internal int representation 
    TransFuncBase::AlphaMode intToAlphaMode(int i);
    int alphaModeToInt(TransFuncBase::AlphaMode mode);

    IntOptionProperty channel1AlphaMode_; ///< use for animation of alpha mode 
    IntOptionProperty channel2AlphaMode_; ///< use for animation of alpha mode 
    IntOptionProperty channel3AlphaMode_; ///< use for animation of alpha mode 
    IntOptionProperty channel4AlphaMode_; ///< use for animation of alpha mode 

    TransFunc1DKeysProperty channel1TransFuncProp_;   ///< linked tf 
    TransFunc1DKeysProperty channel2TransFuncProp_;   ///< linked tf
    TransFunc1DKeysProperty channel3TransFuncProp_;   ///< linked tf
    TransFunc1DKeysProperty channel4TransFuncProp_;   ///< linked tf

    static const std::string loggerCat_;        ///< miau
};

} // namespace

#endif 
