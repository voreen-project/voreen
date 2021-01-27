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

#ifndef VRN_SCREENSPACEAMBIENTOCCLUSION_H
#define VRN_SCREENSPACEAMBIENTOCCLUSION_H

#include "voreen/core/processors/imageprocessorbypassable.h"


#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/cameraproperty.h"
#include "voreen/core/properties/floatproperty.h"
#include "voreen/core/properties/intproperty.h"

namespace voreen {

class ScreenSpaceAmbientOcclusion : public ImageProcessorBypassable {
public:
    ScreenSpaceAmbientOcclusion();
    virtual Processor* create() const;

    virtual std::string getCategory() const  { return "Image Processing";            }
    virtual std::string getClassName() const { return "ScreenSpaceAmbientOcclusion"; }
    virtual CodeState getCodeState() const   { return CODE_STATE_TESTING;            }

    virtual bool isReady() const;

protected:
    virtual void setDescriptions() {
        setDescription("Image processor which performs screen space ambient occlusion.");
    }

    void process();
    virtual std::string generateHeader(const tgt::GpuCapabilities::GlVersion* version);

    RenderPort inportOrig_;
    RenderPort outport_;

    CameraProperty camera_;
    IntProperty radius_;
    FloatProperty occFac_;
    FloatProperty colFac_;
};

} // namespace voreen

#endif //VRN_SCREENSPACEAMBIENTOCCLUSION_H
