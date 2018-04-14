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

#ifndef VRN_MULTIVOLUMECROSSSECTIONANALYZER_H
#define VRN_MULTIVOLUMECROSSSECTIONANALYZER_H

#include "modules/plotting/ports/plotport.h"

#include "voreen/core/processors/volumeprocessor.h"
#include "voreen/core/ports/geometryport.h"
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/floatproperty.h"
#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/vectorproperty.h"

#include "modules/plotting/datastructures/plotdata.h"
#include "voreen/core/datastructures/geometry/meshlistgeometry.h"

#include <vector>

namespace voreen {

/**
 * Todo
 */
class MultiVolumeCrossSectionAnalyzer : public VolumeProcessor {
public:
    MultiVolumeCrossSectionAnalyzer();
    virtual ~MultiVolumeCrossSectionAnalyzer();
    virtual Processor* create() const;

    virtual std::string getClassName() const    { return "MultiVolumeCrossSectionAnalyzer"; }
    virtual std::string getCategory() const     { return "Utility"; }
    virtual CodeState getCodeState() const      { return CODE_STATE_EXPERIMENTAL; }

    virtual bool isReady() const;

protected:
    virtual void setDescriptions() {
        setDescription("todo");
    }

    virtual void beforeProcess();
    virtual void process();
    virtual void initialize();
    virtual void deinitialize();

private:
    void analyzeCrossSections();
    void putOutPlotData();

    void sampleSlice(const VolumeRAM* volume, const tgt::mat4& worldToVoxelMatrix,
        const tgt::vec3& normal, const MeshGeometry& clipFace, float spacing, int& occurrences) const;

    void computeMinMaxDistToPlane(const std::vector<const VolumeBase*>& volumes, const tgt::vec3& normal,
        float& minDist, float& maxDist) const;
    float computeWorldSpacing(const std::vector<const VolumeBase*>& volumes) const;
    int getCurrentSlabID(const std::vector<const VolumeBase*>& volumes, const tgt::vec3& normal,
        float slabPosition, float slabThickness) const;

    void thicknessChanged();

    void forceAnalysis();

    VolumePort volumeInport1_;
    VolumePort volumeInport2_;
    VolumePort volumeInport3_;
    VolumePort volumeInport4_;
    PlotPort plotOutport_;
    GeometryPort slicesOutport_;

    FloatVec3Property crossSectionNormal_;
    FloatProperty slabThickness_;
    IntProperty slicesPerSlab_;
    FloatProperty samplingRate_;
    FloatProperty slabPosition_;
    BoolProperty normalizeOutput_;

    PlotData plotData_;
    MeshListGeometry slicesGeometry_;

    std::vector< std::vector<int> > occurrences_;  ///< detection results (stores one entry per slab and volume)
    std::vector< int > maxOccurrences_;            ///< max occurrence per volume

    bool performAnalysis_;

    static const std::string loggerCat_;
};


} // namespace voreen

#endif // VRN_MULTIVOLUMECROSSSECTIONANALYZER_H
