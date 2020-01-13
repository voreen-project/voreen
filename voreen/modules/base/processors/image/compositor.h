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

#ifndef VRN_COMPOSITOR_H
#define VRN_COMPOSITOR_H

#include "voreen/core/processors/processor.h"
#include "voreen/core/processors/imageprocessor.h"

#include "voreen/core/properties/optionproperty.h"
#include "voreen/core/properties/floatproperty.h"
#include "voreen/core/properties/boolproperty.h"

namespace voreen {

/// Composites the the images at the two inports using the selected blending method.
class VRN_CORE_API Compositor : public ImageProcessor {
public:
    Compositor();
    ~Compositor();

    virtual std::string getCategory() const { return "Image Processing"; }
    virtual std::string getClassName() const { return "Compositor"; }
    virtual Processor::CodeState getCodeState() const { return CODE_STATE_STABLE; }
    virtual Processor* create() const {return new Compositor();}

    virtual bool isReady() const;

protected:
    virtual void setDescriptions() {
        setDescription("Composites two images with a selectable blending method.");
    }

    void process();
    virtual std::string generateHeader(const tgt::GpuCapabilities::GlVersion* version = 0);
    virtual void compile();

    void compositingModeChanged();

    virtual void initialize();
    virtual void deinitialize();

    StringOptionProperty compositingMode_;

    FloatProperty weightingFactor_;
    FloatProperty weightFirst_, weightSecond_;
    BoolProperty addDepth_;
    BoolProperty pipeThroughIfSecondNotReady_;

    RenderPort inport0_;
    RenderPort inport1_;
    RenderPort outport_;

    tgt::Shader* copyShader_;
};


} // namespace voreen

#endif //VRN_COMPOSITOR_H
