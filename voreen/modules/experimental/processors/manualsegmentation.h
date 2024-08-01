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

#ifndef VRN_MANUALSEGMENTATION_H
#define VRN_MANUALSEGMENTATION_H

#include "voreen/core/processors/imageprocessor.h"
#include "modules/base/processors/render/sliceviewer.h"

#include "voreen/core/ports/renderport.h"
#include "voreen/core/ports/volumeport.h"

#include "voreen/core/properties/colorproperty.h"
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/shaderproperty.h"
#include "voreen/core/properties/vectorproperty.h"
#include "voreen/core/properties/optionproperty.h"
#include "voreen/core/properties/eventproperty.h"
#include "tgt/timer.h"

#include "manualsegmentationstorage.h"

namespace voreen {

/**
 * The ManualSegmentation processor is used to manually segment a volume by painting over it.
 * It needs a sliceviewer to paint over it.
 *
 * A ManualSegmentationStorage is used to store the volume of the segmentation and need to be
 * connected as a coprocessor. This design allows to have multiple painters for different
 * slice alignments sharing one volume.
 */
class ManualSegmentation : public ImageProcessor{
public:
    ManualSegmentation();
    ~ManualSegmentation();

    virtual std::string getCategory() const { return "Image Processing"; }
    virtual std::string getClassName() const { return "ManualSegmentation"; }
    virtual Processor::CodeState getCodeState() const { return CODE_STATE_EXPERIMENTAL; }

    virtual Processor* create() const { return new ManualSegmentation(); }

    virtual void process();


    virtual void initialize() ;
protected:
    virtual void setDescriptions();

    // Renderstate managament
    void setRenderState();
    void resetRenderState();

    /**
     * Draws the Brush at the mouse Position
     */
    void drawBrush();
    /**
     * Draws the temporary line for the brush, before it is applied while painting.
     */

    void drawLines();

    /**
     * Creates two orthonormal base vectors of the current slice alignment
     */
    void createSliceBaseVectors(tgt::ivec3 &basex, tgt::ivec3 &basey);

    // Mouse events
    void mouseLocalization(tgt::MouseEvent* ev);
    void toggleEraser(tgt::KeyEvent* ev);
    void pickSegmentationId(tgt::MouseEvent* ev);
    void stopPainting(tgt::MouseEvent* ev);
    void blockPainting(tgt::KeyEvent* ev);

    void usingEraserChanged();
    void invalidateBrush();

    /**
     * Creates a brush with the currentently selected size and stores it in brush_
     * The size of it is not store explicitly
     */
    void createBrush();

    /**
     * This takes the currentId and traces the current lines with the current brush in the Volume
     * in the ManualSegmentationStorage
     */
    void applyBrush();

    /**
     * Traces a single line with the current brush between to points
     */
    void fillLine(tgt::vec3 start, tgt::vec3 end, uint16_t id);


    tgt::vec4 getColorForId(int id, bool usingEraser);

private:
    RenderPort inport_; ///< Inport for a sliceviewer

    RenderPort renderOutport_; ///< Inport overlayed with brush rendering
    GenericCoProcessorPort<ManualSegmentationStorage> storageCoprocessorPort_; ///< Port for the Storage of the Segmentation

    IntProperty segmentationId_; ///< Id used for painting
    IntProperty brushRadius_;  ///< Radius of the brush
    FloatProperty brushAlpha_; ///< alpha value for the brush rendering
    BoolProperty  delayedDraw_; ///< Brush only applied after each line and not while painting for performance reasons

    // Propertys that NEED to be linked to the sliceviewer
    IntVec3Property linkedMousePositionInSlice_;
    OptionProperty<SliceAlignment> linkedSliceAlignment_;
    IntProperty linkedSliceIndex_;
    FloatMat4Property linkedPickingMatrix_;
    BoolProperty usingEraser_;
    IntProperty doubleclickLatancy_;

    // Event Propertys for mouse events
    EventProperty<ManualSegmentation>* mouseEventPress_;
    EventProperty<ManualSegmentation>* keyEventToggleErase_;
    EventProperty<ManualSegmentation>* mouseEventPick_;
    EventProperty<ManualSegmentation>* mouseEventStopPainting_;
    EventProperty<ManualSegmentation>* keyEventBlockPainting_;

    bool isPainting_; ///< is the user painting with the mouse
    bool justFinishedPainting_; ///< painting is finished and will be applied in the process-method

    bool shouldCreateBrush_; ///< the brush in brush_ is invalid
    bool blockPainting_;


    struct Line{
        bool isErasing_;
        std::vector<tgt::vec3> line_;
    };
    /**
     * All the lines of the current painting operation.
     * If the mouse leave the screen and renters a new line is created
     */
    std::vector<Line > currentLines_;

    /**
     * The current brush in square 2d array.
     */
    std::vector<uint8_t> brush_;


    static const std::string loggerCat_;
};

} // namespace voreen

#endif

