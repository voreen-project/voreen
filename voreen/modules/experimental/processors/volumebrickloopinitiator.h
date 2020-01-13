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

#ifndef VRN_VOLUMEBRICKLOOPINITIATOR_H
#define VRN_VOLUMEBRICKLOOPINITIATOR_H

#include "voreen/core/ports/volumeport.h"
#include "voreen/core/ports/loopport.h"
#include "voreen/core/properties/optionproperty.h"
#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/vectorproperty.h"

namespace voreen {

/**
 * Create sub volumes to process through scene
 *
 * NOTICE: Presumes that largest allocation size is larger
 * on main memory than GPU memory.
 */
class VolumeBrickLoopInitiator : public Processor {
    /// Determines on which type of memory the volume should be stored.
    enum VolumeLocation {
        CPU = 0,
        GPU = 1
    };
    friend class OptionProperty<VolumeLocation>;

    /// Determines which property that determines how large the sub volumes are.
    enum AllocationDirective {
        MIN_NUM_BRICK = 0,
        DIV_BY_NUM = 1,
        SIZE_BRICK = 2
    };
    friend class OptionProperty<AllocationDirective>;

public:
    VolumeBrickLoopInitiator();
    virtual Processor* create() const;

    virtual std::string getCategory() const  { return "Utility"; }
    virtual std::string getClassName() const { return "VolumeBrickLoopInitiator"; }
    virtual CodeState getCodeState() const   { return Processor::CODE_STATE_EXPERIMENTAL; }
    virtual bool isUtility() const           { return true; }

    virtual bool isReady() const;

protected:
    virtual void setDescriptions() {
        setDescription("Bricks volume and sends bricks sequentially trough a volume port as new volumes.");
    }

    // calculate the maximum possible convenient allocation size
    void defineMaxAllocationSize();

    // return the next offset
    tgt::svec3 calculateNewOffset(int loopIteration, tgt::svec3 brickSize, tgt::svec3 dims, tgt::svec3 offset);

    // constrain brick size to volume dimensions
    tgt::svec3 calculateCorrectBrickSize(tgt::svec3 brickSize, tgt::svec3 dims, tgt::svec3 offset);

    // return the number of loops
    int calculateLoopIterations(tgt::svec3 brickSize, tgt::svec3 dims, bool hasData = true);

    // create volume if necessary (take spacing and transformation from incoming volume)
    void createNewVolume(Volume* volumeHandle, tgt::svec3 size, tgt::svec3 offset);

    virtual void process();
    virtual void initialize();
    virtual void deinitialize();

    virtual void showAndHideProperties();

    /// Property defining where the output volume is to be stored.
    OptionProperty<VolumeLocation> volOutput_;

    /// Property defining how bricks should be created.
    OptionProperty<AllocationDirective> allocDirective_;

    //IntProperty maxMemorySize_;
    IntProperty divideBy_;

    /// Power of two brick sizes
    OptionProperty<int> brickSize_;

    /// Volume Dimensions
    IntVec3Property dimension_;

    VolumePort inport_;
    VolumePort outport_;
    LoopPort loopInport_;

private:
    Volume* volHandle_;
    bool allFitsOnMainMem_;
    bool allFitsOnGpuMem_;
    tgt::svec3 currentBrickSize_;
    tgt::svec3 currentOffset_;
    static const std::string loggerCat_; ///< category used in logging
};

}

#endif //VRN_IMAGESEQUENCELOOPINITIATOR_H
