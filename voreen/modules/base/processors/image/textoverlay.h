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

#ifndef VRN_TEXTOVERLAY_H
#define VRN_TEXTOVERLAY_H

#include "voreen/core/processors/imageprocessor.h"
#include "voreen/core/ports/textport.h"
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/optionproperty.h"
#include "voreen/core/properties/eventproperty.h"
#include "voreen/core/properties/fontproperty.h"
#include "voreen/core/properties/colorproperty.h"
#include "tgt/font.h"
#include "tgt/texturemanager.h"
#include "tgt/gpucapabilities.h"
#include "tgt/framebufferobject.h"

namespace voreen {

/**
 * Overlays text on top of the input image.
 */
class VRN_CORE_API TextOverlay : public ImageProcessor {
public:
    TextOverlay();
    ~TextOverlay();
    virtual Processor* create() const;

    virtual std::string getClassName() const { return "TextOverlay";      }
    virtual std::string getCategory() const  { return "Image Processing"; }
    virtual CodeState getCodeState() const   { return CODE_STATE_STABLE;  }

    /**
     * Returns true, if the processor has been initialized
     * and the RenderPors are connected.
     */
    virtual bool isReady() const;

protected:
    virtual void setDescriptions() {
        setDescription("Renders a text overlay on top of the input image.");
    }

    virtual void process();
    virtual void initialize();
    virtual void deinitialize();

    /**
     * Render text-overlay over scene
     */
    void renderOverlay();

    /**
     * Collect text to overlay from text-inports
     */
    std::string collectText(std::string key);

    int getNumberOfLines(std::string s);

    void mouseMove(tgt::MouseEvent* e);
    void mouseEnterExit(tgt::MouseEvent* e);

    RenderPort inport_;
    RenderPort outport_;
    RenderPort privatePort_;

    TextPort text0_;
    TextPort text1_;
    TextPort text2_;
    TextPort text3_;

    BoolProperty enabled_;
    StringOptionProperty layout0_;
    StringOptionProperty layout1_;
    StringOptionProperty layout2_;
    StringOptionProperty layout3_;

    EventProperty<TextOverlay> mouseMoveEventProp_;
    EventProperty<TextOverlay> mouseEnterExitEventProp_;

    tgt::ivec2 viewportSize_;
    tgt::ivec2 mousePos_;

    FontProperty fontProp_;

    OptionProperty<int> blendMode_;
        void blendModeOnChange();
    ColorProperty fontColor_;

    static const std::string loggerCat_; ///< category used in logging

private:
    bool renderFollowMouseText_;
};

} // namespace

#endif // VRN_TEXTOVERLAY_H
