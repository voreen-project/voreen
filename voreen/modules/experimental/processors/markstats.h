/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2021 University of Muenster, Germany,                        *
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

#ifndef VRN_MARKSTATS_H
#define VRN_MARKSTATS_H

#include "voreen/core/processors/imageprocessor.h"
#include "voreen/core/datastructures/octree/volumeoctree.h"
#include "modules/base/processors/render/sliceviewer.h"

#include "voreen/core/ports/renderport.h"
#include "voreen/core/ports/volumeport.h"
#include "voreen/core/ports/textport.h"
#include "voreen/core/ports/geometryport.h"

#include "voreen/core/properties/colorproperty.h"
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/shaderproperty.h"
#include "voreen/core/properties/vectorproperty.h"
#include "voreen/core/properties/matrixproperty.h"
#include "voreen/core/properties/optionproperty.h"
#include "voreen/core/properties/eventproperty.h"
#include "voreen/core/properties/numeric/intervalproperty.h"
#include "voreen/core/properties/transfunc/1d/1dkeys/transfunc1dkeysproperty.h"
#include "tgt/font.h"

namespace voreen {

class MarkStats : public ImageProcessor{
public:
    MarkStats();
    ~MarkStats();

    virtual std::string getCategory() const { return "Image Processing"; }
    virtual std::string getClassName() const { return "MarkStats"; }
    virtual Processor::CodeState getCodeState() const { return CODE_STATE_EXPERIMENTAL; }

    virtual Processor* create() const { return new MarkStats(); }

    virtual bool isReady() const;

    virtual void process();
protected:
    virtual void setDescriptions();

    void setRenderState();
    void resetRenderState();
    void drawBrush();

    void createSliceBaseVectors(tgt::ivec3 &basex, tgt::ivec3 &basey);

    void mouseLocalization(tgt::MouseEvent* ev);
    void pickSegmentationId(tgt::MouseEvent* ev);

private:
    struct MarkPoint{
        tgt::ivec3 position_;
        int id_;
    };

    RenderPort inport_;

    RenderPort renderOutport_;
    TextPort textOutport_;
    GeometryPort geometryOutport_;

    IntProperty segmentationId_;
    FloatProperty brushAlpha_;

    IntVec3Property linkedMousePositionInSlice_;
    OptionProperty<SliceAlignment> linkedSliceAlignment_;
    IntProperty linkedSliceIndex_;
    FloatMat4Property linkedPickingMatrix_;

    EventProperty<MarkStats>* mouseEventPress_;
    //EventProperty<MarkStats>* mouseEventPick_; //Currently not used

    std::vector<MarkPoint> positions_;
    std::vector<tgt::vec4> colors_;

    static const std::string loggerCat_;
    bool shouldSetOutport_;
};

} // namespace voreen

#endif
