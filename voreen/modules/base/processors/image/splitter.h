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

#ifndef VRN_SPLITTER_H
#define VRN_SPLITTER_H

#include "voreen/core/processors/renderprocessor.h"
#include "voreen/core/properties/eventproperty.h"
#include "voreen/core/properties/optionproperty.h"
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/floatproperty.h"
#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/colorproperty.h"

namespace voreen {

class VRN_CORE_API Splitter : public RenderProcessor {

    static const int HANDLE_GRAB_TOLERANCE;

public:
    Splitter();
    ~Splitter();
    virtual Processor* create() const;

    virtual std::string getClassName() const { return "Splitter";        }
    virtual std::string getCategory() const  { return "View";            }
    virtual CodeState getCodeState() const   { return CODE_STATE_STABLE; }

    virtual bool isReady() const;

    enum WindowConfiguration{
        FIRST,
        SECOND,
        HORIZONTAL,
        VERTICAL,
        HORIZONTAL_OVERLAY,
        VERTICAL_OVERLAY
    };

protected:
    virtual void setDescriptions() {
        setDescription("Composites to views next to each other horizontally or vertically in a configurable ratio.");
    }

    virtual void process();
    virtual void initialize();
    virtual void deinitialize();

    virtual void onEvent(tgt::Event* e);

    void updateSizes();
    //int getSplitIndex() const;
    void toggleMaximization(tgt::MouseEvent* me);

    void mouseEvent(tgt::MouseEvent* e);

    void distributeMouseEvent(int window, tgt::MouseEvent *newme);

    WindowConfiguration getWindowConfiguration() const;
    int getWindowForEvent(tgt::MouseEvent me, tgt::MouseEvent *translatedme);
    void renderSubview(RenderPort* image, tgt::vec2 minPx, tgt::vec2 maxPx, tgt::vec2 textureMin, tgt::vec2 textureMax);
    tgt::ivec2 getWindowViewport(WindowConfiguration windowConfiguration, int Window);

    BoolProperty showGrid_;
    ColorProperty gridColor_;
    FloatProperty lineWidth_;
    OptionProperty<bool> vertical_;
    BoolProperty overlay_;
    BoolProperty fixSplitPosition_;
    FloatProperty position_;
    IntProperty maximized_;
    BoolProperty maximizeOnDoubleClick_;
    EventProperty<Splitter> maximizeEventProp_;
    tgt::Shader* shader_;

    /// Inport whose rendering is mapped to the frame buffer.
    RenderPort outport_;

    RenderPort inport1_;
    RenderPort inport2_;
    int lastWindow_;
    int currentPort_;

    bool isDragging_;
};

} // namespace voreen

#endif
