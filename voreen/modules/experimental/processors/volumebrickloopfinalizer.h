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

#ifndef VRN_VOLUMEBRICKLOOPFINALIZER_H
#define VRN_VOLUMEBRICKLOOPFINALIZER_H

#include "voreen/core/ports/volumeport.h"
#include "voreen/core/ports/loopport.h"
#include "voreen/core/properties/vectorproperty.h"

namespace voreen {

/**
 * Defines an brick-processing loop in combination with volumebrickloopInitiator.
 */
class VolumeBrickLoopFinalizer : public Processor {

public:
    VolumeBrickLoopFinalizer();
    ~VolumeBrickLoopFinalizer();
    virtual Processor* create() const;

    virtual std::string getCategory() const  { return "Utility"; }
    virtual std::string getClassName() const { return "VolumeBrickLoopFinalizer"; }
    virtual CodeState getCodeState() const   { return Processor::CODE_STATE_EXPERIMENTAL; }
    virtual bool isUtility() const           { return true; }

    virtual bool isReady() const;

protected:
    virtual void setDescriptions() {
        setDescription("Collects result from proccessed bricks and writes new volume to HDD.");
    }

    // determine the next offset
    void defineNewOffset();

    virtual void process();
    virtual void initialize();
    virtual void deinitialize();

    /// Volume Dimensions
    IntVec3Property dimension_;

    VolumePort inport_;
    VolumePort outport_;
    LoopPort loopOutport_;

    VolumeRAM* volume_;

private:
    tgt::svec3 currentOffset_;
    tgt::svec3 lastVolumeSize_;
    static const std::string loggerCat_; ///< category used in logging
};

}

#endif
