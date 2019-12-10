/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2019 University of Muenster, Germany,                        *
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

#ifndef VRN_RENDERLOOPINITIATOR_H
#define VRN_RENDERLOOPINITIATOR_H

#include "voreen/core/processors/renderprocessor.h"

#include "voreen/core/properties/intproperty.h"

namespace voreen {

/**
 * In combination with RenderLoopFinalizer, this class is used to define render loops.
 */
class VRN_CORE_API RenderLoopInitiator : public RenderProcessor {
public:
    RenderLoopInitiator();

    virtual std::string getCategory() const { return "Utility"; }
    virtual std::string getClassName() const { return "RenderLoopInitiator"; }
    virtual CodeState getCodeState() const { return CODE_STATE_STABLE; }
    virtual bool isUtility() const { return true; }

    virtual bool isReady() const;

protected:
    virtual void setDescriptions() {
        setDescription("In combination with RenderLoopFinalizer, this processor is used to establish render loops: The path between RenderLoopInitiator and RenderLoopFinalizer in the rendering network is executed as often as specified by this processor's \"Iterations\" property. At the end of each iteration, the temporary rendering result is passed back from the loop finalizer to the loop initiator where it is used as input for the next iteration.");
    }

    virtual void initialize();
    virtual void deinitialize();

    tgt::Shader* shader_;

    virtual Processor* create() const;
    virtual void process();

    RenderPort inport_;
    RenderPort outport_;
    RenderPort loopInport_;

    IntProperty numIterations_;

private:
    void updateIterationCount();

};

}

#endif //VRN_RENDERLOOPINITIATOR_H
