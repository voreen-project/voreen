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

#ifndef VRN_TOUCHPAINTER_H
#define VRN_TOUCHPAINTER_H

#include "voreen/core/processors/imageprocessor.h"

#include "voreen/core/properties/buttonproperty.h"
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/colorproperty.h"
#include "voreen/core/properties/floatproperty.h"
#include "voreen/core/properties/vectorproperty.h"

#include "voreen/core/ports/renderport.h"

#include "tgt/event/touchevent.h"
#include "tgt/event/touchpoint.h"

namespace voreen {

class TouchPainter : public ImageProcessor {

public:
    TouchPainter();

    virtual std::string getCategory() const { return "Image Processing"; }
    virtual std::string getClassName() const { return "TouchPainter"; }
    virtual Processor::CodeState getCodeState() const { return CODE_STATE_TESTING; }

    virtual Processor* create() const { return new TouchPainter(); }

    virtual void process();
    virtual bool isReady() const;
    virtual void initialize();
    virtual void deinitialize();

protected:
    virtual void setDescriptions() {
        setDescription("Paint some awesome stuff with your fingers!");
    }

    virtual void touchEvent(tgt::TouchEvent* e);
    virtual void mouseDoubleClickEvent(tgt::MouseEvent* e);
    virtual std::string generateHeader(const tgt::GpuCapabilities::GlVersion* version);

    void onClear();

    float getRandomFloat() const {
        return (float)std::rand() / RAND_MAX;
    }

    FloatProperty touchRadius_;
    ButtonProperty clearImage_;
    BoolProperty randomizeColors_;
    ColorProperty color_;
    FloatProperty opacity_;

    RenderPort inport_;
    RenderPort outport_;
    RenderPort privatePort_;

    std::deque<tgt::TouchPoint> toPaint_;
    tgt::Shader* copyShader_;
    std::map<unsigned int, tgt::vec3> colorMap_;

    static const std::string loggerCat_;
};

} // namespace voreen

#endif
