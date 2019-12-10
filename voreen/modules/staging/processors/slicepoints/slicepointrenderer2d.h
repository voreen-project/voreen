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

#ifndef VRN_SLICEPOINTRENDERER2D_H
#define VRN_SLICEPOINTRENDERER2D_H

#include "voreen/core/processors/imageprocessor.h"
#include "modules/base/processors/render/sliceviewer.h"

#include "voreen/core/ports/renderport.h"

#include "voreen/core/properties/colorproperty.h"
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/shaderproperty.h"
#include "voreen/core/properties/vectorproperty.h"
#include "voreen/core/properties/matrixproperty.h"
#include "voreen/core/properties/optionproperty.h"
#include "voreen/core/properties/eventproperty.h"
#include "tgt/font.h"

namespace voreen {

    class TriangleMeshGeometryUInt16IndexedSimple;

class SlicePointRenderer2D : public ImageProcessor{
public:
    SlicePointRenderer2D();
    ~SlicePointRenderer2D();

    virtual std::string getCategory() const { return "Image Processing"; }
    virtual std::string getClassName() const { return "SlicePointRenderer2D"; }
    virtual Processor::CodeState getCodeState() const { return CODE_STATE_TESTING; }

    virtual Processor* create() const { return new SlicePointRenderer2D(); }

    virtual void process();

protected:
    virtual void setDescriptions() {
        setDescription("TODO.");
    }
    virtual std::string generateHeader(const tgt::GpuCapabilities::GlVersion* version = 0);

    void updateGeometry();
    virtual void deinitialize();

    /** Renders the selected slice points. */
    void renderPoints();
    void renderPointGeometryHelper(tgt::Shader* prog, const tgt::ivec3& pointPos, const tgt::vec4& screenPos, const tgt::vec4& col);
    void renderPointInfoHelper(const tgt::ivec3& pointPos, const tgt::vec4& screenPos, const tgt::vec4& col);
    /** Help function to blend the selected points over the input image. */
    void blendPointsOverInport();

    /** called, if the current mouse position of the linked slice view has been changed. */
    void linkedMousePositionInSliceOnChange();
private:
    //Ports
    RenderPort inport_;
    RenderPort outport_;
    RenderPort privatePort_;

    //Properties
        //general
    BoolProperty enable_;
    IntProperty pointRadius_;
    BoolProperty renderPointInfo_;
        //points
    BoolProperty renderPoint0_;
    ColorProperty pointColor0_;
    IntVec3Property pointPos0_;
    BoolProperty renderPoint1_;
    ColorProperty pointColor1_;
    IntVec3Property pointPos1_;
    BoolProperty renderPoint2_;
    ColorProperty pointColor2_;
    IntVec3Property pointPos2_;
    BoolProperty renderPoint3_;
    ColorProperty pointColor3_;
    IntVec3Property pointPos3_;
        //to be linked with sliceviewer
    IntVec3Property linkedMousePositionInSlice_;
    OptionProperty<SliceAlignment> linkedSliceAlignment_;
    IntProperty linkedSliceIndex_;
    FloatMat4Property linkedPickingMatrix_;
        //shader
    ShaderProperty shaderProp_;
        //member
    TriangleMeshGeometryUInt16IndexedSimple* equalGeometry_;
    TriangleMeshGeometryUInt16IndexedSimple* moreGeometry_;
    TriangleMeshGeometryUInt16IndexedSimple* lessGeometry_;
    tgt::Font* pointInfoFont_;
    float fontLineWidth_;

    EventProperty<SlicePointRenderer2D>* point0EnableKeyEvent_;
    EventProperty<SlicePointRenderer2D>* point1EnableKeyEvent_;
    EventProperty<SlicePointRenderer2D>* point2EnableKeyEvent_;
    EventProperty<SlicePointRenderer2D>* point3EnableKeyEvent_;
    EventProperty<SlicePointRenderer2D>* point0ActivateKeyEvent_;
    EventProperty<SlicePointRenderer2D>* point1ActivateKeyEvent_;
    EventProperty<SlicePointRenderer2D>* point2ActivateKeyEvent_;
    EventProperty<SlicePointRenderer2D>* point3ActivateKeyEvent_;
        //event callback functions
    void activatePoint0KeyEvent(tgt::KeyEvent* e);
    void activatePoint1KeyEvent(tgt::KeyEvent* e);
    void activatePoint2KeyEvent(tgt::KeyEvent* e);
    void activatePoint3KeyEvent(tgt::KeyEvent* e);
    void enablePoint0KeyEvent(tgt::KeyEvent* e);
    void enablePoint1KeyEvent(tgt::KeyEvent* e);
    void enablePoint2KeyEvent(tgt::KeyEvent* e);
    void enablePoint3KeyEvent(tgt::KeyEvent* e);
    void mainEventFunc(int point, BoolProperty* enableProp, bool enable = false);

    int activePoint_; ///< stores the currently active point. -1 if no point is active.

    static const std::string loggerCat_;
};

} // namespace voreen

#endif
