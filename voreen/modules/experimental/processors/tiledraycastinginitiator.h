/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany,                        *
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

#ifndef VRN_TILEDRAYCASTINGINITIATOR_H
#define VRN_TILEDRAYCASTINGINITIATOR_H

#include "voreen/core/processors/renderprocessor.h"

#include "voreen/core/ports/loopport.h"
#include "voreen/core/properties/intproperty.h"

namespace voreen {

class TiledRaycastingInitiator : public RenderProcessor {
public:
    TiledRaycastingInitiator();

    virtual std::string getCategory() const  { return "Utility"; }
    virtual std::string getClassName() const { return "TiledRaycastingInitiator"; }
    virtual CodeState getCodeState() const   { return CODE_STATE_EXPERIMENTAL; }
    virtual bool isUtility() const           { return true; }
    virtual bool usesExpensiveComputation() const { return true; }

    virtual bool isReady() const;

protected:
    virtual void setDescriptions() {
        setDescription("");
    }

    virtual void initialize();
    virtual void deinitialize();

    tgt::Shader* shader_;

    virtual Processor* create() const;
    virtual void process();

    RenderPort inportEntryPoints_;
    RenderPort inportExitPoints_;
    RenderPort outportEntryPoints_;
    RenderPort outportExitPoints_;

    LoopPort loopInport_;

    IntProperty subdivisions_;

    static const std::string loggerCat_;

private:
    void updateSubdivision();

};

}

#endif //VRN_TILEDRAYCASTINGINITIATOR_H
