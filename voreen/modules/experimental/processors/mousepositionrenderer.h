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

#ifndef VRN_MOUSEPOSITIONRENDERER_H
#define VRN_MOUSEPOSITIONRENDERER_H

#include "voreen/core/processors/imageprocessor.h"

#include "voreen/core/properties/floatproperty.h"
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/colorproperty.h"
#include "voreen/core/properties/vectorproperty.h"
#include "voreen/core/properties/filedialogproperty.h"
#include "voreen/core/properties/cameraproperty.h"
#include "tgt/texturemanager.h"
#include "tgt/gpucapabilities.h"

#include <deque>

namespace voreen {

class MousePositionRenderer : public ImageProcessor {
public:
    MousePositionRenderer();
    virtual Processor* create() const;

    virtual std::string getClassName() const  { return "MousePositionRenderer";     }
    virtual std::string getCategory() const   { return "Image Processing"; }
    virtual CodeState getCodeState() const    { return CODE_STATE_EXPERIMENTAL;  }

protected:
    virtual void setDescriptions() {
        setDescription("");
    }

    virtual void process();

    virtual void initialize();
    virtual void deinitialize();
    virtual void onEvent(tgt::Event* e);

private:
    RenderPort inport_;
    RenderPort outport_;

    BoolProperty renderCursor_;
    BoolProperty renderTrail_;
    FloatProperty width_;
    FloatProperty opacity_;
    ColorProperty color_;
    FloatProperty cursorSize_;

    tgt::vec2 mpos_;
    tgt::vec2 mpos2_;
    bool isPressed_;
    bool isInside_;
    std::deque< std::pair<tgt::vec2, int> > history_;
    std::deque< std::pair<tgt::vec2, int> > history2_;
    tgt::Timer* timer_;
    tgt::EventHandler eventHandler_;    // A local eventhanlde which is added to the timer.

    tgt::Shader* copyShader_;
};

} // namespace

#endif
