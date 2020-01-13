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

#ifndef VRN_BUTTONOVERLAYPROCESSOR_H
#define VRN_BUTTONOVERLAYPROCESSOR_H

#include "voreen/core/voreenapplication.h"
#include "voreen/core/processors/volumeraycaster.h"
#include "voreen/core/properties/buttonproperty.h"
#include "voreen/core/properties/eventproperty.h"
#include "voreen/core/properties/stringproperty.h"
#include "voreen/core/interaction/interactionhandler.h"
#include "voreen/core/processors/imageprocessor.h"

#include "tgt/font.h"

namespace voreen {

/**
 * Performs a single pass multi volume raycasting.
 */
class ButtonOverlayProcessor : public ImageProcessor {

public:
    ButtonOverlayProcessor();

    virtual ~ButtonOverlayProcessor();

    virtual std::string getCategory() const { return "Image Processing"; }
    virtual std::string getClassName() const { return "ButtonOverlayProcessor"; }
    virtual Processor::CodeState getCodeState() const { return CODE_STATE_OBSOLETE; }
    virtual std::string getProcessorInfo() const;
    virtual Processor* create() const { return new ButtonOverlayProcessor(); }

    virtual void initialize();
    virtual void deinitialize();

    /**
     * Performs the raycasting.
     *
     * Initialize two texture units with the entry and exit params and renders
     * a screen aligned quad.
     */
    virtual void process();

    virtual bool isReady() const;

    void onMouseClick(tgt::MouseEvent* mouseEvent);

protected:
    virtual void setDescriptions() {
        setDescription("");
    }


    virtual std::string generateHeader(const tgt::GpuCapabilities::GlVersion* version = 0);

    void renderOverlayImage();
    void onOverlayChange();

    void renderFont(int i, tgt::ivec2& center);

private:

    RenderPort inport_;
    RenderPort outport_;
    RenderPort overlayPort_;
    RenderPort pickingPort_;

    IntProperty buttonRadiusX_;
    IntProperty buttonRadiusY_;
    IntProperty groupBorder_;
    IntProperty horzBorder_;
    IntProperty vertBorder_;
    IntProperty fontSize_;
    FloatProperty opacity_;        // between 0 and 1
    BoolProperty pushButtonMode_;        // between 0 and 1
    EventProperty<ButtonOverlayProcessor> buttonChecker_;

    tgt::Texture* buttonTexturePressed_;
    tgt::Texture* buttonTextureReleased_;
    std::vector<BoolProperty*> buttonProps_;
    std::vector<ButtonProperty*> pushButtonProps_;
    std::vector<BoolProperty*> isVisibleProps_;
    std::vector<StringProperty*> nameProps_;

    bool regenerateOverlay_;
    tgt::ivec2 oldSize_;

    tgt::Font* font_;
};


} // namespace voreen

#endif // VRN_BUTTONOVERLAYPROCESSOR_H
