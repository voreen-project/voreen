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

#ifndef VRN_VOLUMELISTMULTICHANNELADAPTER_H
#define VRN_VOLUMELISTMULTICHANNELADAPTER_H

#include "voreen/core/processors/processor.h"
#include "voreen/core/ports/genericport.h"
#include "voreen/core/ports/volumeport.h"

namespace voreen {

/**
 * Merges all separated channel volumes inside the input list to a single multi-channel volume in the output list.
 */
class VRN_CORE_API VolumeListMultiChannelAdapter : public Processor {
public:
    VolumeListMultiChannelAdapter();
    virtual ~VolumeListMultiChannelAdapter();
    virtual Processor* create() const;

    virtual std::string getClassName() const  { return "VolumeListMultiChannelAdapter"; }
    virtual std::string getCategory() const   { return "Input"; }
    virtual CodeState getCodeState() const    { return CODE_STATE_EXPERIMENTAL; }
    virtual bool isUtility() const { return true; }

protected:
    virtual void setDescriptions() {
        setDescription("Merges all separated channel volumes inside the input list to a single multi-channel volume in the output list.");
    }

    virtual void process();

private:

    void onChannelCountChanged();

    std::vector<std::unique_ptr<const VolumeBase>> volumes_;

    IntProperty numChannels_;
    StringOptionProperty layout_;
    BoolProperty invertChannel1_;
    BoolProperty invertChannel2_;
    BoolProperty invertChannel3_;
    BoolProperty invertChannel4_;

    VolumeListPort inport_;
    VolumeListPort outport_;
};

}   //namespace

#endif
