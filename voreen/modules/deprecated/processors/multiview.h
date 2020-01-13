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

#ifndef VRN_MULTIVIEW_H
#define VRN_MULTIVIEW_H

#include "voreen/core/processors/renderprocessor.h"
#include "voreen/core/properties/eventproperty.h"
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/colorproperty.h"

namespace voreen {

class VRN_CORE_API MultiView : public RenderProcessor {
public:
    MultiView();
    ~MultiView();

    virtual bool isReady() const;

    virtual std::string getCategory() const { return "View"; }
    virtual std::string getClassName() const { return "MultiView"; }
    virtual Processor::CodeState getCodeState() const { return CODE_STATE_STABLE; }
    virtual Processor* create() const;

    virtual void invalidate(int inv = INVALID_RESULT);

    virtual void onEvent(tgt::Event* e);

    virtual void serialize(Serializer& s) const;
    virtual void deserialize(Deserializer& s);

protected:
    virtual void setDescriptions() {
        setDescription("Combines an arbitrary number of input images to a single view.");
    }

    virtual void beforeProcess();
    virtual void process();
    virtual void initialize();

    void toggleMaximization(tgt::MouseEvent* me);
    void updateSizes();
    void updateGridResolution();

    void mouseMove(tgt::MouseEvent* e);
    void renderQuad();

    void establishConnectionOrder();

    BoolProperty showGrid_;
    ColorProperty gridColor_;
    IntProperty maximized_;
    BoolProperty maximizeOnDoubleClick_;
    EventProperty<MultiView> maximizeEventProp_;
    EventProperty<MultiView> mouseMoveEventProp_;

    RenderPort outport_;

    /// Inports whose renderings are mapped to the frame buffer.
    RenderPort inport_;

    std::vector<RenderPort*> renderPorts_;
    tgt::ivec2 gridResolution_;

    int currentPort_;
    bool isDragging_;

    std::vector<std::string> portNameOrder_;
    std::vector<std::string> procNameOrder_;
};

} // namespace voreen

#endif
