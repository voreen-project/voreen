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

#ifndef VRN_TILEDRAYCASTINGFINALIZER_H
#define VRN_TILEDRAYCASTINGFINALIZER_H

#include "voreen/core/processors/renderprocessor.h"

#include "voreen/core/ports/loopport.h"
#include "voreen/core/properties/intproperty.h"

namespace voreen {

/**
 * In combination with RenderLoopInitiator, this class is used to define render loops.
 */
class TiledRaycastingFinalizer : public RenderProcessor {
public:
    TiledRaycastingFinalizer();
    virtual Processor* create() const;

    virtual std::string getCategory() const { return "Utility"; }
    virtual std::string getClassName() const { return "TiledRaycastingFinalizer"; }
    virtual CodeState getCodeState() const { return CODE_STATE_EXPERIMENTAL; }

    virtual bool isUtility() const { return true; }
    virtual bool isReady() const;

protected:
    virtual void setDescriptions() {
        setDescription("");
    }

    virtual void process();
    virtual void initialize();
    virtual void deinitialize();

    void updateInputSizeRequest();

    RenderPort inportTile_;
    RenderPort outportRendering_;
    LoopPort loopOutport_;

    IntProperty subdivisions_;

    tgt::Shader* shaderPrg_;

    static const std::string loggerCat_;
};

}

#endif //VRN_TILEDRAYCASTINGFINALIZER_H
