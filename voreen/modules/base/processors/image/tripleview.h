/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2024 University of Muenster, Germany,                        *
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

#ifndef VRN_TRIPLEVIEW_H
#define VRN_TRIPLEVIEW_H

#include "voreen/core/processors/renderprocessor.h"
#include "voreen/core/properties/eventproperty.h"
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/optionproperty.h"
#include "voreen/core/properties/colorproperty.h"

namespace voreen {

class VRN_CORE_API TripleView : public RenderProcessor {
    enum WindowConfiguration {
        abc = 0,
        Abc = 1,
        Bac = 2,
        Cab = 3,
        // Internal for maximized views
        A = 4,
        B = 5,
        C = 6
    };

public:
    TripleView();
    ~TripleView();
    virtual Processor* create() const;

    virtual std::string getClassName() const { return "TripleView";       }
    virtual std::string getCategory() const  { return "View";             }
    virtual CodeState getCodeState() const   { return CODE_STATE_STABLE;  }

    virtual bool isReady() const;

protected:
    void renderPortQuad(RenderPort& rp, tgt::vec3 translate, tgt::vec3 scale);
    void renderLargeSmallSmall(RenderPort& large, RenderPort& small1, RenderPort& small2);

    virtual void setDescriptions() {
        setDescription("Combines three input images in a configurable layout.");
    }

    virtual void process();
    virtual void initialize();
    virtual void deinitialize();

    virtual void onEvent(tgt::Event* e);

    void distributeMouseEvent( int Window, tgt::MouseEvent *newme );

    void updateSizes();

    void toggleMaximization(tgt::MouseEvent* me);
    void mouseMove(tgt::MouseEvent* e);
    WindowConfiguration getWindowConfiguration() const;
    int getWindowForEvent(tgt::MouseEvent ev, tgt::MouseEvent* translatedMouseEvent);
    int getWindowForEventLargeSmallSmall(tgt::MouseEvent ev, int large, int small1, int small2, tgt::MouseEvent* translatedMouseEvent);

    tgt::ivec2 getWindowViewportLargeSmallSmall(bool isFirst);
    tgt::ivec2 getWindowViewport(WindowConfiguration configuration, int window);


    BoolProperty showGrid_;
    ColorProperty gridColor_;
    IntOptionProperty configuration_;
    EventProperty<TripleView> maximizeEventProp_;
    IntProperty maximized_;
    BoolProperty maximizeOnDoubleClick_;

    tgt::Shader* shader_;

    /// Inport whose rendering is mapped to the frame buffer.
    RenderPort outport_;

    RenderPort inport1_;
    RenderPort inport2_;
    RenderPort inport3_;

    int currentPort_;
    int lastWindow_;
    bool isDragging_;
};

} // namespace voreen

#endif
