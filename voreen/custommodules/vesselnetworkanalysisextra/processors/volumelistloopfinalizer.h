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

#ifndef VRN_VOLUMELISTLOOPFINALIZER_H
#define VRN_VOLUMELISTLOOPFINALIZER_H

#include "voreen/core/processors/processor.h"
#include "voreen/core/ports/loopport.h"
#include "voreen/core/ports/volumeport.h"
//#include "voreen/core/ports/volumelistport.h"
#include "voreen/core/properties/progressproperty.h"

namespace voreen {

/**
 * Defines an volume-processing loop in combination with VolumeListLoopInitiator.
 * All volumes of the input list are processed by the loop.
 *
 * Note: This is very much experimental. Due to ownership issues volumes are copied from/into lists!
 */
class VRN_CORE_API VolumeListLoopFinalizer : public Processor {

public:
    VolumeListLoopFinalizer();
    virtual Processor* create() const;

    virtual std::string getCategory() const  { return "Utility"; }
    virtual std::string getClassName() const { return "VolumeListLoopFinalizer"; }
    virtual CodeState getCodeState() const   { return Processor::CODE_STATE_EXPERIMENTAL; }
    virtual bool isUtility() const           { return true; }

    virtual bool isReady() const;

protected:
    virtual void setDescriptions() {
        setDescription("Establishes a volume processing loop in combination with VolumeListLoopInitiator. All volumes of the input list are consecutively processed by the loop. Connect this processor's loop outport with the loop inport of an VolumeListLoopInitiator, in order to define the loop.");
    }

    virtual void process();

    VolumePort inport_;
    VolumeListPort outport_;
    LoopPort loopOutport_;
    ProgressProperty progress_;

    std::unique_ptr<VolumeList> volumeList_;
private:
    VolumeBase* getVolumeFromInport();

    static const std::string loggerCat_; ///< category used in logging

};

}

#endif //VRN_VOLUMELISTLOOPFINALIZER_H
