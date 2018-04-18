/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany.                        *
 * Visualization and Computer Graphics Group <http://viscg.uni-muenster.de>        *
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

#ifndef VRN_POIRENDERER_H
#define VRN_POIRENDERER_H

#include "poistorage.h"

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
#include "voreen/core/properties/fontproperty.h"
#include "voreen/core/properties/shaderproperty.h"
#include "voreen/core/properties/vectorproperty.h"
#include "voreen/core/properties/stringexpressionproperty.h"
#include "voreen/core/properties/matrixproperty.h"
#include "voreen/core/properties/optionproperty.h"
#include "voreen/core/properties/eventproperty.h"
#include "voreen/core/properties/numeric/intervalproperty.h"
#include "voreen/core/properties/transfunc/1d/1dkeys/transfunc1dkeysproperty.h"
#include "tgt/font.h"

#include <memory>

namespace voreen {

class POIRenderer2d : public ImageProcessor{
public:
    POIRenderer2d();
    ~POIRenderer2d();

    virtual std::string getCategory() const { return "Image Processing"; }
    virtual std::string getClassName() const { return "POIRenderer2d"; }
    virtual Processor::CodeState getCodeState() const { return CODE_STATE_EXPERIMENTAL; }

    virtual Processor* create() const { return new POIRenderer2d(); }

    virtual bool isReady() const;

    virtual void process();
    virtual void initialize() override;
protected:
    virtual void setDescriptions();
    
    void setRenderState();
    void resetRenderState();

    void createSliceBaseVectors(tgt::vec3 &basex, tgt::vec3 &basey);
    
    void mouseLocalization(tgt::MouseEvent* ev);

    POIPointID FindPOI(tgt::vec3 mousepos_voxel);

    void setPoint(tgt::MouseEvent* ev);
    
    void buildPointsPerSlice();
    void renderPOIs();
    void drawPOISlice(int slicenum); ///< needs a shader active

    tgt::mat4 voxelToScreenMatrix() const;

    void removePoint(tgt::KeyEvent* ev);
    float getRadiusVoxel();
    bool getRadiusPixel();
    void renderSelection( tgt::Bounds selection );
private:
    RenderPort inport_;
    VolumePort volumeInport_;
    RenderPort renderOutport_;
    GenericCoProcessorPort<POIStorage> cpPort_;
     
    FloatProperty radius_;

    BoolProperty showArrows_;
    IntProperty arrowRange_;

    ShaderProperty shader_;
    ShaderProperty arrowShader_;

    BoolProperty showIds_;
    FontProperty fontProp_;
    ColorProperty distanceMeasureColor_;
  
    OptionProperty<POIStorage::SelectionMode> selectionMode_;
    BoolProperty selectBackgroud_;

    IntVec3Property linkedMousePositionInSlice_;
    OptionProperty<SliceAlignment> linkedSliceAlignment_;
    
    IntProperty linkedSliceIndex_;
    FloatMat4Property linkedPickingMatrix_;

    std::unique_ptr<EventProperty<POIRenderer2d> > mouseEventPress_;
    std::unique_ptr<EventProperty<POIRenderer2d> > mouseEventPick_;
    std::unique_ptr<EventProperty<POIRenderer2d> > mouseEventSet_;
    std::unique_ptr<EventProperty<POIRenderer2d> > mouseEventLocation_;
    std::unique_ptr<EventProperty<POIRenderer2d> > removeEvent_;

    std::vector<std::vector<POIPoint>> pointsPerSlice_;

    static const std::string loggerCat_;
    bool shouldRebuildPointsPerSlice_;
    POIPointID draggedPoint_;
    bool isDragging_;
    bool isSelecting_;
    tgt::vec3 selectionStart_voxel_;
    tgt::vec3 selectionEnd_voxel_;
    tgt::vec3 dragOffset_;
};
}

#endif
