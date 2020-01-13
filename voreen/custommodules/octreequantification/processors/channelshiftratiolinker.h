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

#ifndef VRN_CHANNELSHIFTRATIOLINKER_H
#define VRN_CHANNELSHIFTRATIOLINKER_H

#include "voreen/core/ports/volumeport.h"
#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/vectorproperty.h"
#include "voreen/core/properties/matrixproperty.h"

namespace voreen {

class VRN_CORE_API ChannelShiftRatioLinker : public Processor {
public:
    ChannelShiftRatioLinker();

    virtual Processor* create() const;

    virtual std::string getCategory() const   { return "Utility";              }
    virtual std::string getClassName() const  { return "ChannelShiftRatioLinker";   }
    virtual CodeState getCodeState() const    { return CODE_STATE_EXPERIMENTAL; }

protected:
    virtual void setDescriptions() {
        setDescription("Can be used to link the channel shift settings of an octree volume (i.e. the octree raycaster) to a downsampled single-channel version of the same data set. Also provides a transformation matrix for each channel which can be used with a VolumeTransform processor to apply the channel shift of the downsampled channel to a SingleVolumeRaycaster rendering. Caution: only transforms the input channel shift (input A) to channel shift B, not vice versa!");
    }

    virtual void process();

    /// sets the max values of the clipping regions 
    virtual void inputVolumesChanged();

    /// carries over changes from region A to region B
    virtual void channelShiftAChanged();

    static const std::string loggerCat_; ///< category used in logging

private:
    VolumePort inportA_;
    VolumePort inportB_;

    IntProperty channelSelect_;
    
    BoolProperty applyChannelShift_;

    FloatVec3Property channelShiftA1_;
    FloatVec3Property channelShiftA2_;
    FloatVec3Property channelShiftA3_;
    FloatVec3Property channelShiftA4_;

    FloatVec3Property channelShiftB_;
};

}   //namespace

#endif 
