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

#ifndef VRN_PATHSURFACECREATOR_H
#define VRN_PATHSURFACECREATOR_H

#include <list>
#include <memory>

#include "voreen/core/processors/asynccomputeprocessor.h"

#include "voreen/core/ports/geometryport.h"
#include "voreen/core/ports/volumeport.h"
#include "custommodules/sciviscontest2022/processors/surface/integrators.h"

#include "modules/flowanalysis/ports/streamlinelistport.h"
#include "voreen/core/datastructures/geometry/trianglemeshgeometry.h"
#include "voreen/core/properties/numeric/intervalproperty.h"
#include "voreen/core/datastructures/geometry/pointsegmentlistgeometry.h"

namespace voreen {

enum class PathSurfaceGeometryInputSpace {
    SPACE_VOXEL,
    SPACE_WORLD
};

static std::ostream& operator<<(std::ostream& stream, PathSurfaceGeometryInputSpace space) {
    if (space == PathSurfaceGeometryInputSpace::SPACE_VOXEL)
        stream << "Voxelspace";
    else
        stream << "Worldspace";
    return stream;
}

struct IntegralCurve;

using dense_output_type = tgt::vec3;

struct Interval {
    std::vector<float> s;               ///< sample positions (in parameter space)
    std::vector<dense_output_type> f;   ///< sample position
    std::vector<tgt::vec3> v;           ///< sample velocities
    bool done = false;
};

using QPredicate = std::function<bool(
    float, float, float, 
    const dense_output_type&, const dense_output_type&, const dense_output_type&)>;

struct PathSurfaceCreatorInput {
    tgt::mat4 worldToVoxelMatrix;                   ///< transforms world space to voxel space

    QPredicate refine;                              ///< predicate for the need of refinement
    QPredicate discontinuity;                       ///< predicate for discontinuities

    std::vector<Interval> intervals;                ///< Seed curves
    std::vector<float> time_parameters;             ///< sample times
    PortDataPointer<VolumeList> flowVolumes;        ///< volume list

    VolumeRAM::Filter filterMode;                   ///< Filtermode for spatial sampler
    float stepsize;                                 ///< stepsize for integration
    float velocityUnitConversion;                   ///< conversion for units
    float elapsedTimeBetweenVolumes;                ///< Elapsed time from one volume to another. Usually 1 second
};

struct PathSurfaceCreatorOutput {
    std::unique_ptr<Geometry> geometry;
    std::unique_ptr<Geometry> pointlist;
};

/**
 * @brief This class basically implements the pathsurface creation algorithm
 * described in 'Generation of Accurate Integral Surfaces in Time-Dependant Vector Fields' by
 * Christoph Garth, Hari Krishnan, Xavier Tricoche, et al. from 2008
 */
class VRN_CORE_API PathSurfaceCreator : public AsyncComputeProcessor<PathSurfaceCreatorInput, PathSurfaceCreatorOutput> {
private:
public:
    PathSurfaceCreator();

    Processor* create() const override { return new PathSurfaceCreator(); }
    bool isReady() const override;

    std::string getCategory() const override { return "Path surface Processing"; }
    std::string getClassName() const override { return "PathSurfaceCreator"; }
    Processor::CodeState getCodeState() const override { return CODE_STATE_EXPERIMENTAL; }
protected:
    void setDescriptions() override {
		setDescription("This processor is used to create path surfaces from a vec3 volume. The resulting pathlines can be visualized or modified "
					   "by other processors (especially the PathSurfaceRenderer). The input geometry must be a PointSegmentListGeometryVec3 and "
                       "describes seeding curve (one seeding curve per segment) in voxel or world space (as can be selected in the Processor).");
    }

    ComputeInput prepareComputeInput() override;
	ComputeOutput compute(ComputeInput input, ProgressReporter &progressReporter) const override;
	void processComputeOutput(ComputeOutput output) override;
private:
    static bool s_QRefineInterNodeDistance(const dense_output_type& f0, const dense_output_type& f1, float threshold);
    static bool s_QRefineInterSegmentAngle(const dense_output_type& f0, const dense_output_type& f1, const dense_output_type& f2, float threshold);
    static bool s_QRefineTriangleArea(const dense_output_type& f0, const dense_output_type& f1, const dense_output_type& f2, float threshold);

    static bool s_QDiscontinuity(
        float x0, float x1, float x2,
        const dense_output_type& f0, const dense_output_type& f1, const dense_output_type& f2,
        float threshold);
private:
    static const std::string loggerCat_;

    VolumeListPort volumeInport_;
	GeometryPort seedingInport_;
    GeometryPort geometryOutport_;
    GeometryPort pointlistOutport_;
	StreamlineListPort pathlineOutport_;

    BoolProperty refinementInterNodeDistance_;
    BoolProperty refinementInterSegmentAngle_;
    BoolProperty refinementTriangleArea_;
    FloatProperty interNodeDistanceThresholdProperty_;
    FloatProperty interSegmentAngleThresholdProperty_;
    FloatProperty triangleAreaThresholdProperty_;
    FloatProperty magnitudeUpperBoundProperty_;

    OptionProperty<VolumeRAM::Filter> filterProperty_;
    OptionProperty<PathSurfaceGeometryInputSpace> inputSpaceProperty_;
    FloatProperty elapsedTimeBetweenVolumes_;
    FloatOptionProperty velocityUnitConversionProperty_;

    IntIntervalProperty timeIntervalProperty_;
    IntProperty resolutionPerTimestepProperty_;
};

}

#endif
