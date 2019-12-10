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

#ifndef VRN_RANDOMWALKERANALYZER_H
#define VRN_RANDOMWALKERANALYZER_H

#include "voreen/core/processors/volumeprocessor.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"

#include "voreen/core/ports/geometryport.h"

#include "voreen/core/properties/optionproperty.h"
#include "voreen/core/properties/floatproperty.h"
#include "voreen/core/properties/vectorproperty.h"
#include "voreen/core/properties/cameraproperty.h"
#include "voreen/core/properties/buttonproperty.h"
#include "voreen/core/properties/intproperty.h"

#include <string>
#include "tgt/vector.h"


namespace voreen {

class Volume;

class VRN_CORE_API RandomWalkerAnalyzer : public VolumeProcessor {

public:
    struct UncertaintyRegion {
        int label_;
        int size_;
        tgt::vec3 centerOfMass_;
        tgt::ivec3 llf_;
        tgt::ivec3 urb_;
    };

    RandomWalkerAnalyzer();
    virtual ~RandomWalkerAnalyzer();
    virtual Processor* create() const;

    virtual std::string getCategory() const { return "Volume Processing"; }
    virtual std::string getClassName() const { return "RandomWalkerAnalyzer"; }
    virtual CodeState getCodeState() const { return CODE_STATE_EXPERIMENTAL; }
    virtual std::string getProcessorInfo() const;

    virtual void invalidate(int inv = INVALID_RESULT);
    virtual bool isReady() const;

    const std::vector<UncertaintyRegion>& getUncertaintyRegions() const;
    void setFocusRegion(int region);
    void setFocusPoint(const UncertaintyRegion& region);
    void zoomOnFocusRegion(int region, int numZooSteps = 40);

    void computeRandomWalker();

    void resetZoom();
    void resetZoomAnimated(int iterations);

    virtual void timerEvent(tgt::TimeEvent* e);
    void forceUpdate();

private:
    virtual void setDescriptions() {
        setDescription("Analyzes a random walker segmentation and detects regions of high uncertainty (experimental).");
    }

    virtual void process();
    virtual void initialize();
    virtual void deinitialize();

    void updateNumSlices();

    VolumeRAM_UInt8* probabilityThreshold(VolumeRAM_UInt8* volume, uint8_t min, uint8_t max);
    Volume* connectedComponents(Volume* volume, size_t& numLabels);
    void analyzeUncertaintyRegions(Volume* labelVolume, int numLabels);
    void stretchLabels(VolumeRAM_UInt16* labelVolume);
    void printLabelInformation() const;
    void updateFocusAreaPropVolume();

    VolumeRAM_UInt8* uncertaintyEdges(VolumeRAM_UInt8* volumeProb, VolumeRAM_UInt8* volumeFiltered);
    VolumeRAM_UInt8* erosion(VolumeRAM_UInt8* volume);
    VolumeRAM_UInt8* createDistanceMap(VolumeRAM_UInt8* volume);

    void setFocusAndZoomAnimated(tgt::vec2 focusXY, tgt::vec2 focusXZ, tgt::vec2 focusYZ, float zoom, int numIterations);
    void performAnimStep();

    VolumePort inportVolume_;
    VolumePort outportUncertaintyAreas_;
    VolumePort outportFocusArea_;
    VolumePort outportFocusAreaProbabilities_;
    GeometryPort outportFocusPoint_;

    FloatVec2Property uncertaintyRange_;

    IntOptionProperty connectivity_;
    IntProperty minUncertaintyAreaSize_;
    IntProperty maxUncertaintyAreas_;

    tgt::EventHandler* eventHandler_;
    tgt::Timer* timer_;
    IntProperty sliceXY_;
    IntProperty sliceXZ_;
    IntProperty sliceYZ_;
    FloatProperty zoomXY_;
    FloatProperty zoomXZ_;
    FloatProperty zoomYZ_;
    FloatVec2Property focusXY_;
    FloatVec2Property focusXZ_;
    FloatVec2Property focusYZ_;
    CameraProperty camera_;
    ButtonProperty computeButton_;

    std::vector<UncertaintyRegion> uncertaintyRegions_;
    tgt::vec3 focusPoint_;
    int focusRegion_;

    VolumeRAM_UInt8* probVolume_;
    Volume* probVolumeCC_;

    tgt::vec2 animFocusXYStart_;
    tgt::vec2 animFocusXZStop_;
    tgt::vec2 animFocusYZStart_;
    tgt::vec2 animFocusXYStop_;
    tgt::vec2 animFocusXZStart_;
    tgt::vec2 animFocusYZStop_;
    float animZoomStart_;
    float animZoomStop_;
    int numAnimSteps_;
    int animStep_;

    bool forceUpdate_;

    static const std::string loggerCat_; ///< category used in logging
};

} //namespace

#endif
