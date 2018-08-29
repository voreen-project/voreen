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

#ifndef VRN_QUADVIEW_H
#define VRN_QUADVIEW_H

#include "voreen/core/processors/renderprocessor.h"
#include "voreen/core/properties/eventproperty.h"
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/colorproperty.h"

namespace voreen {

class VRN_CORE_API QuadView : public RenderProcessor {
public:
    QuadView();
    ~QuadView();

    virtual bool isReady() const;
    virtual void process();

    virtual std::string getCategory() const { return "View"; }
    virtual std::string getClassName() const { return "QuadView"; }
    virtual Processor::CodeState getCodeState() const { return CODE_STATE_STABLE; }
    virtual Processor* create() const;

    virtual void initialize();

    virtual void deinitialize();

    virtual void invalidate(int inv = INVALID_RESULT);

    virtual void onEvent(tgt::Event* e);
protected:
    virtual void setDescriptions() {
        setDescription("Allows to arrange four renderings in one multiview layout.");
    }

    void toggleMaximization(tgt::MouseEvent* me);
    void updateSizes();

    void portSizeReceiveChanged();

    void mouseMove(tgt::MouseEvent* e);
    void distributeMouseEvent(tgt::MouseEvent* me);
    void distributeTouchEvent(tgt::TouchEvent* te);
    void distributeKeyEvent(tgt::KeyEvent* ke);

    BoolProperty showGrid_;
    ColorProperty gridColor_;
    IntProperty maximized_;
    BoolProperty maximizeOnDoubleClick_;
    EventProperty<QuadView> maximizeEventProp_;
    EventProperty<QuadView> mouseMoveEventProp_;
    BoolProperty forwardDoubleClickEvent_;

    tgt::Shader* shader_;

    /// Inport whose rendering is mapped to the frame buffer.
    RenderPort outport_;

    RenderPort inport1_;
    RenderPort inport2_;
    RenderPort inport3_;
    RenderPort inport4_;

    int currentPort_;

    bool isDragging_;
};

} // namespace voreen

#endif
