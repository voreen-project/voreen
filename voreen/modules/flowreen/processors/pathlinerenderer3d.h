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

#ifndef VRN_PATHLINERENDERER3D_H

#include <string>
#include "voreen/core/processors/renderprocessor.h"
#include "voreen/core/ports/volumeport.h"
#include "voreen/core/ports/renderport.h"
#include "voreen/core/ports/genericcoprocessorport.h"
#include "flowreenprocessor.h"
#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/floatproperty.h"
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/optionproperty.h"
#include "voreen/core/properties/vectorproperty.h"
#include "voreen/core/properties/colorproperty.h"

namespace voreen {

class CameraInteractionHandler;
class Modality;
class FlowOrthogonalSliceRenderer;

/**
 * Processor for rendering pathlines from time-dependent flow data using
 * geometrical primitives like points, lines, tubes and arrows.
 */
class PathlineRenderer3D : public RenderProcessor, private FlowreenProcessor {
public:
    PathlineRenderer3D();
    virtual ~PathlineRenderer3D();

    virtual Processor* create() const { return new PathlineRenderer3D(); }
    virtual std::string getCategory() const { return "Flow Visualization"; }
    virtual std::string getClassName() const { return "PathlineRenderer3D"; }
    virtual Processor::CodeState getCodeState() const { return Processor::CODE_STATE_STABLE; }
    virtual void initialize();
    virtual void invalidate(int inv = INVALID_RESULT);
    virtual void process();

private:
    bool applyThresholds(std::vector<tgt::vec3>& pathline, const float& length) const;

    /**
     * Marks the data in the supplied VolumeSeries containing the contextual data by
     * checking whether the intesities fall within the range of the threshold.
     * Therefore Volumes of the same size consiting of bool are created in
     * intensityMasks_, indicating the test result.
     */
    size_t markContextData(const tgt::vec2& thresholds);

    /**
     * Verfies whether all input flow volumes are of the correct type (VolumeFlow3D)
     * and whether they are all of the same dimensions. The pointers to valid flow
     * volumes are stored in flows_ member which is replaced.
     *
     * @param   forceUpdate Determines whether to replace the existing flow volumes in
     *                      member flows_ or not.
     * @return  The number of valid flow volumes within flows_. If this member is empty,
     *          rendering cannot be performed and has to be canceled.
     */
    size_t checkFlows(const bool forceUpdate = false);

    void adjustTimestepProperty();
    void clearPathlines();
    void initPathlines(const size_t numPoints);
    void initPathlinesGrid(const size_t spacing);
    void initPathlinesSliceGrid(const size_t spacing);
    void onIntensityChange();
    void onLineStyleChange();
    void onSeedingStrategyChange();
    void onTimestepChange();
    void onThresholdingChange();

    void renderAsArrows(const tgt::vec4& lineColor);
    void renderAsLines(const tgt::vec4& lineColor);
    void renderAsLineSegments(const tgt::vec4& lineColor);
    void renderAsPoints(const tgt::vec4& lineColor);
    void renderAsTubes(const tgt::vec4& lineColor);

    bool setupShader();

private:
    virtual void setDescriptions() {
        setDescription("Processor for rendering pathlines from time-dependent flow data using geometrical primitives like points, lines, tubes and arrows. the PathlineRenderer3D can be used to visualize vector field data. It expects one flow volume and one scalar context volume. The flow volume is depicted by geometric primitives (e.g., arrows or lines). The processor provides different seeding strategies. Such that the geometric primitives can be arranged based on the grid or depend on the context volume. Thus, a slice-based seeding becomes possible.");
    }

    static const std::string loggerCat_;

    enum LineStyle {
        STYLE_POINTS,
        STYLE_LINES,
        STYLE_TUBES,
        STYLE_ARROWS,
        STYLE_SEGMENTS
    };
    friend class OptionProperty<LineStyle>;

    enum SeedingStrategy {
        SEED_RANDOM,
        SEED_GRID,
        SEED_SLICES_RANDOM,
        SEED_SLICES_GRID
    };
    friend class OptionProperty<SeedingStrategy>;

    enum Thresholding {
        THRESHOLDING_NONE,
        THRESHOLDING_LINELENGTH,
        THRESHOLDING_INTENSITY,
        THRESHOLDING_OR,
        THRESHOLDING_AND
    };
    friend class OptionProperty<Thresholding>;

    OptionProperty<SeedingStrategy>* seedingStrategyProp_;
    IntProperty numSeedpointsProp_;
    BoolProperty seedOnXYSliceProp_;
    BoolProperty seedOnXZSliceProp_;
    BoolProperty seedOnYZSliceProp_;
    IntProperty gridSpacingProp_;
    IntProperty objectSizeProp_;
    OptionProperty<LineStyle>* lineStyleProp_;
    OptionProperty<Thresholding>* thresholdingProp_;
    FloatVec2Property lineLengthProp_;
    FloatVec2Property intensityProp_;
    ColorProperty lineColorProp_;
    IntProperty segmentLengthProp_;
    FloatProperty integrationStepProp_;
    FloatProperty timestepProp_;          /** intended to be linked with a clock processor. */
    CameraProperty camProp_;
    CameraInteractionHandler* cameraHandler_;

    /** All flow volumes */
    std::vector<const Flow3D*> flows_;
    const VolumeList* contextCollection_;
    const VolumeList* flowCollection_;

    /** Bool Volumes holding "true" if the intensity at the location in the volume is
     * within the threshold range or "false" otherwise. */
    std::vector<VolumeBase*> intensityMasks_;

    tgt::Shader* shader_;

    LineStyle currentStyle_;
    SeedingStrategy seedingStrategy_;
    Thresholding thresholding_;
    size_t numPathlines_;
    std::vector<tgt::vec3>* pathlines_;
    size_t currentTimestep_;
    size_t previousTimestep_;
    tgt::ivec3 slicePositions_;

    GenericCoProcessorPort<FlowOrthogonalSliceRenderer> coInport_;
    RenderPort imgOutport_;
    VolumeListPort inportContext_;
    VolumeListPort inportFlows_;
};

}   // namespace

#endif  // VRN_PATHLINERENDERER3D_H

