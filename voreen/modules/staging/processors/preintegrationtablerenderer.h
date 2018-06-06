/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany.                        *
 * Visualization and Computer Graphics Group <http://viscg.uni-muenster.de>        *
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

#ifndef VRN_PREINTEGRATIONTABLERENDERER_H
#define VRN_PREINTEGRATIONTABLERENDERER_H

#include "voreen/core/processors/imageprocessor.h"

#include "voreen/core/properties/transfunc/1d/1dkeys/transfunc1dkeysproperty.h"
#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/optionproperty.h"

#include "tgt/texturemanager.h"
#include "tgt/gpucapabilities.h"

namespace voreen {

class PreIntegrationTableRenderer : public ImageProcessor {
public:
    PreIntegrationTableRenderer();
    virtual Processor* create() const;

    virtual std::string getClassName() const  { return "PreIntegrationTableRenderer"; }
    virtual std::string getCategory() const   { return "Image Processing"; }
    virtual CodeState getCodeState() const    { return CODE_STATE_TESTING;  }

    virtual bool isReady() const;
protected:
    virtual void setDescriptions() {
        setDescription("Renders the pre-integration table of a 1D keys transfer function, mainly for testing and debugging purposes.");
    }

    virtual void process();

    virtual std::string generateHeader(const tgt::GpuCapabilities::GlVersion* version = 0);

private:
    RenderPort outport_;

    TransFunc1DKeysProperty transFunc_;
    StringOptionProperty piMode_;

    BoolProperty showGrid_;
    IntProperty gridDivisions_;
};

} // namespace

#endif
