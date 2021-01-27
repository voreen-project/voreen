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

#ifndef VRN_FLOWMAPCREATOR_H
#define VRN_FLOWMAPCREATOR_H

#include "voreen/core/processors/processor.h"
#include "voreen/core/ports/genericport.h"
#include "voreen/core/ports/volumeport.h"

namespace voreen {

/**
 * Creates Flow Maps from separated x, y and z (volume) images.
 */
class VRN_CORE_API FlowMapCreator : public Processor {
public:
    FlowMapCreator();
    virtual ~FlowMapCreator();
    virtual Processor* create() const;

    virtual std::string getClassName() const  { return "FlowMapCreator";        }
    virtual std::string getCategory() const   { return "Volume Processing";     }
    virtual CodeState getCodeState() const    { return CODE_STATE_EXPERIMENTAL; }
    virtual bool isUtility() const            { return true; }

protected:
    virtual void setDescriptions() {
        setDescription("Creates Flow Maps from separated x, y and z (volume) images and provides additional options to "
                       "adjust certain properties. This is even supported for single-channel volumes. However, the "
                       "typical use-case is a series of x, y and z images converted into 3D vector fields to then be "
                       "used by other processors of the <br>FlowAnalysis</br> module.");
        numChannels_.setDescription("Number of channels of the output volume");
        layout_.setDescription("Layout in which channel volumes are stored in input list");
        mirrorX_.setDescription("Shall voxels be mirrored along x axis?");
        mirrorY_.setDescription("Shall voxels be mirrored along y axis?");
        mirrorZ_.setDescription("Shall voxels be mirrored along z axis?");
        swizzleChannel1_.setDescription("Map first channel to any other channel");
        swizzleChannel2_.setDescription("Map second channel to any other channel");
        swizzleChannel3_.setDescription("Map third channel to any other channel");
        swizzleChannel4_.setDescription("Map fourth channel to any other channel");
        negateChannel1_.setDescription("Shall the first channel (after swizzling) be negated?");
        negateChannel2_.setDescription("Shall the second channel (after swizzling) be negated?");
        negateChannel3_.setDescription("Shall the third channel (after swizzling) be negated?");
        negateChannel4_.setDescription("Shall the fourth channel (after swizzling) be negated?");
    }

    virtual void process();

private:

    void onChannelCountChanged();

    std::vector<std::unique_ptr<VolumeBase>> volumes_;

    IntProperty numChannels_;
    StringOptionProperty layout_;

    BoolProperty mirrorX_;
    BoolProperty mirrorY_;
    BoolProperty mirrorZ_;

    OptionProperty<size_t> swizzleChannel1_;
    OptionProperty<size_t> swizzleChannel2_;
    OptionProperty<size_t> swizzleChannel3_;
    OptionProperty<size_t> swizzleChannel4_;

    BoolProperty negateChannel1_;
    BoolProperty negateChannel2_;
    BoolProperty negateChannel3_;
    BoolProperty negateChannel4_;

    VolumeListPort inport_;
    VolumeListPort outport_;
};

}   //namespace

#endif
