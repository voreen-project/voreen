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

#include "pathsurfacecreator.h"
#include "integrators.h"

#include <Eigen/Dense>

namespace voreen {

const std::string PathSurfaceCreator::loggerCat_("sciviscontest2022.PathSurfaceCreator");

using Timeline = Streamline;

struct IntegralCurvePoint {
    tgt::vec3 p;
    float t;
};
struct IntegralCurve {
    bool left_edge;
    std::vector<IntegralCurvePoint> ps;
    float t_begin;
    float s;
};

using CurveSet = std::map<float, IntegralCurve>;
using InactiveCurveSet = std::map<float, std::pair<bool, IntegralCurve>>;


void advanceRibbon(
        const IntegralCurve& S0, 
        const IntegralCurve& S1, 
        TriangleMeshGeometryTexCoord* geometry,
        InactiveCurveSet& inactiveSet,
        std::pair<std::size_t, std::size_t> activeEdge);

PathSurfaceCreator::PathSurfaceCreator() :
      seedingInport_(Port::INPORT, "pathsurfacecreator.seedingInport", "Seeding Point Inport")
    , volumeInport_(Port::INPORT, "pathsurfacecreator.volumeInport", "Volume Inport")
    , geometryOutport_(Port::OUTPORT, "pathsurfacecreator.geometryOutport", "Geometry Outport")
    , pointlistOutport_(Port::OUTPORT, "pathsurfacecreator.pointlistOutport", "Pointlist Outport")
    , pathlineOutport_(Port::OUTPORT, "pathsurfacecreator.pathlineOutport", "Pathlines Outport")
    , interSegmentAngleThresholdProperty_("interSegmentAngleThreshold", "Angle threshold between timeline and pathline segments", .6f, 0.001f, 3.1415692f * 2)
    , interNodeDistanceThresholdProperty_("interNodeDistanceThreshold", "Distance threshold between timeline and pathline nodes", 10.0f, 1.f, 10000.0f)
    , magnitudeUpperBoundProperty_("magnitudeUpperBound", "Upper limit of the magnitude of the first and second derivative", 2000.0f, 1000.f, 100000.0f)
    , triangleAreaThresholdProperty_("triangleAreaThreshold", "Threshold for the triangle area of timelines and pathlines", 10.0f, 1.f, 10000.0f)
    , resolutionPerTimestepProperty_("numTimesteps", "Number of timesteps", 4, 2, 100)
    , timeIntervalProperty_("timeInterval", "Time Interval (for integration)", 0, 0, 0, 0)
    , filterProperty_("filterModeProp", "Filtering:", Processor::INVALID_RESULT, false, Property::LOD_DEVELOPMENT)
    , inputSpaceProperty_("inputSpace", "Input Space:", Processor::INVALID_RESULT, false, Property::LOD_DEFAULT)
    , refinementInterNodeDistance_("refineInterNodeDistance", "Refine using distance: ", true)
    , refinementInterSegmentAngle_("refineInterSegmentAngle", "Refine using angles: ", false)
    , refinementTriangleArea_("refineTriangleArea", "Refine using triangle area: ", false)
    , elapsedTimeBetweenVolumes_("elapsedTimeBetweenVolumes", "Elapsed time between volumes", 1.f, 0.0001f, std::numeric_limits<float>::max())
    , velocityUnitConversionProperty_("velocityUnitConversion", "Input Velocity Unit")
{
    addPort(seedingInport_);
    addPort(volumeInport_);
    addPort(geometryOutport_);
    addPort(pointlistOutport_);
    addPort(pathlineOutport_);

    addProperty(inputSpaceProperty_);
    
    addProperty(refinementInterNodeDistance_);
    addProperty(interNodeDistanceThresholdProperty_);
    addProperty(refinementInterSegmentAngle_);
    addProperty(interSegmentAngleThresholdProperty_);
    addProperty(refinementTriangleArea_);
    addProperty(triangleAreaThresholdProperty_);
    addProperty(magnitudeUpperBoundProperty_);

    addProperty(elapsedTimeBetweenVolumes_);
    addProperty(velocityUnitConversionProperty_);
    addProperty(resolutionPerTimestepProperty_);
    addProperty(timeIntervalProperty_);
    addProperty(filterProperty_);

    inputSpaceProperty_.addOption("voxel", "Voxelspace", PathSurfaceGeometryInputSpace::SPACE_VOXEL);
    inputSpaceProperty_.addOption("world", "Worldspace", PathSurfaceGeometryInputSpace::SPACE_WORLD);

    inputSpaceProperty_.setGroupID("Input");

    filterProperty_.addOption("linear", "Linear", VolumeRAM::LINEAR);
    filterProperty_.addOption("nearest", "Nearest", VolumeRAM::NEAREST);

    refinementInterNodeDistance_.setGroupID("Refinement");
    interNodeDistanceThresholdProperty_.setGroupID("Refinement");
    refinementInterSegmentAngle_.setGroupID("Refinement");
    interSegmentAngleThresholdProperty_.setGroupID("Refinement");
    refinementTriangleArea_.setGroupID("Refinement");
    triangleAreaThresholdProperty_.setGroupID("Refinement");

    magnitudeUpperBoundProperty_.setGroupID("Discontinuities");

    elapsedTimeBetweenVolumes_.setGroupID("Integration");
    resolutionPerTimestepProperty_.setGroupID("Integration");
    timeIntervalProperty_.setGroupID("Integration");
    filterProperty_.setGroupID("Integration");
    velocityUnitConversionProperty_.setGroupID("Integration");

    // Chose the values such that multiplying with real world values we get mm(/s)
    // which (for some reason) is the default voreen length unit.
    velocityUnitConversionProperty_.addOption("km/s", "km/s", 1000000.0f);
    velocityUnitConversionProperty_.addOption("m/s", "m/s", 1000.0f);
    velocityUnitConversionProperty_.addOption("dm/s", "dm/s", 100.0f); // Quite unusual.
    velocityUnitConversionProperty_.addOption("cm/s", "cm/s", 10.0f);
    velocityUnitConversionProperty_.addOption("mm/s", "mm/s", 1.0f);
    velocityUnitConversionProperty_.set("mm/s");

    ON_CHANGE_LAMBDA(volumeInport_, [&]() { 
        if (!volumeInport_.hasData())
            return;
        
        timeIntervalProperty_.setMaxValue(volumeInport_.getThreadSafeData()->size()); 
    });

    ON_CHANGE_LAMBDA(volumeInport_ , [&]() { if (isReady()) process(); });
    ON_CHANGE_LAMBDA(seedingInport_, [&]() { if (isReady()) process(); });
}

bool PathSurfaceCreator::isReady() const {
    auto ready = geometryOutport_.isConnected();
    ready &= seedingInport_.isReady();
    ready &= volumeInport_.isReady();
    return ready;
}

PathSurfaceCreatorInput PathSurfaceCreator::prepareComputeInput() {
    auto flowVolumes = volumeInport_.getThreadSafeData();

    const auto velocityUnitConversion = velocityUnitConversionProperty_.getValue();

    auto doInterNodeDistanceRefinement = refinementInterNodeDistance_.get();
    auto doInterSegmentAngleRefinement = refinementInterSegmentAngle_.get();
    auto doTriangleAreaRefinement = refinementTriangleArea_.get();
    auto interNodeDistanceThreshold = interNodeDistanceThresholdProperty_.get() * velocityUnitConversion;
    auto interSegmentAngleThreshold = interSegmentAngleThresholdProperty_.get();
    auto triangleAreaThreshold = triangleAreaThresholdProperty_.get() * velocityUnitConversion;

    auto magnitudeUpperBound = magnitudeUpperBoundProperty_.get() * velocityUnitConversion;

    if (!(doInterNodeDistanceRefinement || doInterSegmentAngleRefinement || doTriangleAreaRefinement)) {
        throw InvalidInputException("You need to enable at least one of the Refinement algorithms!", InvalidInputException::S_ERROR);
    }

    QPredicate q_refine = [=](
        float x0, float x1, float x2, 
        const dense_output_type& f0, const dense_output_type& f1, const dense_output_type& f2
    ) -> bool {
        return 
            (doInterNodeDistanceRefinement && s_QRefineInterNodeDistance(f0, f2, interNodeDistanceThreshold * interNodeDistanceThreshold)) ||
            (doTriangleAreaRefinement && s_QRefineTriangleArea(f0, f1, f2, triangleAreaThreshold * triangleAreaThreshold)) ||
            // do segment angle refinement last since it calculates two sqare roots
            (doInterSegmentAngleRefinement && s_QRefineInterSegmentAngle(f0, f1, f2, interSegmentAngleThreshold));
    };
    QPredicate q_discon = [=](
        float x0, float x1, float x2, 
        const dense_output_type& f0, const dense_output_type& f1, const dense_output_type& f2
    ) -> bool {
        return s_QDiscontinuity(x0, x1, x2, f0, f1, f2, magnitudeUpperBound * magnitudeUpperBound);
    };


    auto t0 = timeIntervalProperty_.get().x;
    auto t1 = timeIntervalProperty_.get().y;
    auto nT = resolutionPerTimestepProperty_.get() * (t1 - t0);

    auto dt = (t1 - t0) / static_cast<float>(nT - 1);
    std::vector<float> time_parameters;
    time_parameters.resize(nT);
    for (auto i = 0u; i < nT; ++i) time_parameters[i] = t0 + i * dt;

    const VolumeBase* referenceVolume = flowVolumes->first();

    auto geometryTransformation = seedingInport_.getData()->getTransformationMatrix();

    auto voxelToWorldMatrix = tgt::mat4::identity;

    if (inputSpaceProperty_.getValue() == PathSurfaceGeometryInputSpace::SPACE_VOXEL) {
        voxelToWorldMatrix = referenceVolume->getVoxelToWorldMatrix();
    }

    const tgt::mat4 worldToVoxelMatrix = referenceVolume->getWorldToVoxelMatrix();
    const tgt::Bounds roi = referenceVolume->getBoundingBox().getBoundingBox();
    const RealWorldMapping rwm = referenceVolume->getRealWorldMapping();

    auto baseVolume = std::floor(time_parameters[0]);
    VolumeRAMRepresentationLock current{ flowVolumes->at(baseVolume + 0) };
    VolumeRAMRepresentationLock next   { flowVolumes->at(baseVolume + 1) };
    if(!*current || !*next) {
        throw InvalidInputException("RAM representation not available", InvalidInputException::S_ERROR);
    }

    auto alpha = time_parameters[0] - baseVolume;
    SpatioTemporalSampler sampler(*current, *next, alpha, rwm, filterProperty_.getValue(), worldToVoxelMatrix);

    auto inGeometry = static_cast<const PointSegmentListGeometryVec3*>(seedingInport_.getData());
    std::vector<Interval> intervals;
    intervals.resize(inGeometry->getNumSegments());
    auto perSegment = 1.f / inGeometry->getNumSegments();
    for (auto i = 0; i < inGeometry->getNumSegments(); ++i) {
        const auto& segment = inGeometry->getSegment(i);
        const auto ds = perSegment / segment.size();
        for (auto j = 0u; j < segment.size(); ++j) {
            intervals[i].s.emplace_back(i * perSegment + j * ds);
            // Apply the geometry transformation then apply voxel to world
            // transformation because that is what the Sampler expects
            auto p = voxelToWorldMatrix * geometryTransformation * segment[j];
            intervals[i].f.emplace_back(p);
            intervals[i].v.emplace_back(sampler.sample(p));
        }        
    }
    
    return {
        worldToVoxelMatrix,
        std::move(q_refine),                        ///< Refinement algorithm
        std::move(q_discon),                        ///< Discontinuity detection algorithm
        std::move(intervals),                       ///< initial seeding curve intervals
        std::move(time_parameters),                 ///< time parameters 
        std::move(flowVolumes),                     ///< flow volumes
        filterProperty_.getValue(),                 ///< how to filter the volumes
        dt,                                         ///< timestep
        velocityUnitConversion,                     ///< velocity unit conversion
        elapsedTimeBetweenVolumes_.get()            ///< elapsed time between volumes
    };
}

void refine(std::vector<Interval>& intervals, const QPredicate& q_refine, const QPredicate& q_discon) {

    std::set<float> s_insert;
    auto i = 0u;
    // for every segment, basically
    while (true) {
        if (i >= intervals.size())
            break;

        tgtAssert(intervals[i].s.size() == intervals[i].f.size(), "Size of parameter array must match the size of the value array");

        if (intervals[i].done) {
            // we already processed this segment
            ++i;
            continue;
        }

        auto& s = intervals[i].s;
        auto& f = intervals[i].f;
        auto& v = intervals[i].v;

        s_insert.clear();
        for (auto j = 2u; j < s.size(); ++j) {
            // for every adjacent three points:
            // s_(i-1) = (x0, f0)
            // s_( i ) = (x1, f1)
            // s_(i+1) = (x2, f2)
            const auto x0 = s[j - 2];
            const auto x1 = s[j - 1];
            const auto x2 = s[j - 0];
            const auto f0 = f[j - 2];
            const auto f1 = f[j - 1];
            const auto f2 = f[j - 0];

            if (q_discon(x0, x1, x2, f0, f1, f2)) {
                // We approximate a discontinuity here, thus
                // we break execution up into two different passes
                const auto it0 = intervals.begin() + i;

                Interval i0;
                i0.s = std::vector<float>            { s.begin(), s.begin() + j };
                i0.f = std::vector<dense_output_type>{ f.begin(), f.begin() + j };
                i0.v = std::vector<tgt::vec3>        { v.begin(), v.begin() + j };

                Interval i1;
                i1.s = std::vector<float>            { s.begin() + (j + 1), s.end() };
                i1.f = std::vector<dense_output_type>{ f.begin() + (j + 1), f.end() };
                i1.v = std::vector<tgt::vec3>        { v.begin() + (j + 1), v.end() };
                
                intervals.erase(it0);
                // We ignore intervals that got too small
                if (i0.s.size() > 2)
                    intervals.emplace_back(std::move(i0));
                if (i1.s.size() > 2)
                    intervals.emplace_back(std::move(i1));
                // rerun with new intervals
                return refine(intervals, q_refine, q_discon);
            }

        }

        for (auto j = 2u; j < s.size(); ++j) {
            // for every adjacent three points:
            // s_(i-1) = (x0, f0)
            // s_( i ) = (x1, f1)
            // s_(i+1) = (x2, f2)
            const auto x0 = s[j - 2];
            const auto x1 = s[j - 1];
            const auto x2 = s[j - 0];
            const auto f0 = f[j - 2];
            const auto f1 = f[j - 1];
            const auto f2 = f[j - 0];
            if (q_refine(x0, x1, x2, f0, f1, f2)) {
                // the refinement threshold gods tell us to
                // evaluate at the midpoints as well
                auto v0 = (x0 + x1) / 2.f;
                auto v1 = (x1 + x2) / 2.f;

                s_insert.insert({ v0, v1 });
            }
        }

        if (s_insert.empty()) {
            // Advance to the next interval
            intervals[i++].done = true;
            continue;
        }
        
        for (auto s_ : s_insert) {
            tgt::vec3 f_;
            tgt::vec3 v_;
            for (auto k = 0u; k < s.size() - 1; ++k) {
                if (s[k] < s_ && s[k + 1] >= s_) {
                    auto t = (s_ - s[k]) / (s[k + 1] - s[k]);
                    f_ = t * f[k] + (1 - t) * f[k + 1];
                    v_ = t * v[k] + (1 - t) * v[k + 1];
                    break;
                }
            }

            // emplace s_, f_, v_ into the sorted arrays s, f, v (all sorted by the parameters in s).
            // do binary search to find the correct place for s_, f_ and v_
            auto upper_bound = s.size() - 1;
            auto lower_bound = 1u;
            while (true) {
                auto middle = (upper_bound + lower_bound) / 2;
                if (s[middle] <= s_) {
                    if (s.size() < middle + 1 || s[middle + 1] >= s_) {
                        // insert between middle and middle + 1
                        s.emplace(s.begin() + (middle + 1), s_);
                        f.emplace(f.begin() + (middle + 1), f_);
                        v.emplace(v.begin() + (middle + 1), v_);
                        break;
                    }
                    
                    lower_bound = middle;
                }
                else {
                    // s[middle] > s_
                    if (middle == 0 || s[middle - 1] <= s_) {
                        // insert between middle - 1 and middle
                        s.emplace(s.begin() + middle, s_);
                        f.emplace(f.begin() + middle, f_);
                        v.emplace(v.begin() + middle, v_);
                        break;
                    }

                    upper_bound = middle;
                }
            }
        }
        s_insert.clear();
    }
}

PathSurfaceCreatorOutput PathSurfaceCreator::compute(ComputeInput input, ProgressReporter &progressReporter) const {
    auto flowVolumes = std::move(input.flowVolumes);
    const VolumeBase* referenceVolume = flowVolumes->first();

    std::unique_ptr<VolumeRAMRepresentationLock> currentVolume;
    std::unique_ptr<VolumeRAMRepresentationLock> nextVolume;

    const tgt::mat4 worldToVoxelMatrix = referenceVolume->getWorldToVoxelMatrix();
    const tgt::Bounds roi = referenceVolume->getBoundingBox().getBoundingBox();
    const RealWorldMapping rwm = referenceVolume->getRealWorldMapping();

    auto inBounds = [roi](tgt::vec3 point) {
        return roi.containsPoint(point);
    };

    auto current_timeline = std::move(input.intervals);

    auto t_begin = input.time_parameters[0];
    auto t_end   = input.time_parameters[input.time_parameters.size() - 1];

    CurveSet activeSet;
    InactiveCurveSet inactiveSet;

    PathSurfaceCreatorOutput output;
    auto geometry = new TriangleMeshGeometryTexCoord();
    output.geometry = std::unique_ptr<Geometry>(geometry);
    refine(current_timeline, input.refine, input.discontinuity);

    for (const auto& interval : current_timeline) {
        for (auto i = 0u; i < interval.s.size(); ++i) {
            // advace active set for the first timeline
            activeSet.emplace(interval.s[i], IntegralCurve{ 
                i == 0u, 
                { 
                    IntegralCurvePoint{ interval.f[i], input.time_parameters[0] }
                }, 
                input.time_parameters[0],
                interval.s[i]
            });
        }
    }

    SubtaskProgressReporter pointGeneration(progressReporter, { 0.f, 0.8f });

    /**
     * @brief Algorithm:
     * Let I_0 = C [Seeding curve]
     * 
     * for i in 1...N:
     *   I_(i+1) = S(t_i, LI_i; t_(i+1))
     * 
     * where L is the linear interpolant of I
     */
    auto lastVolume = -1;
    bool finished;
    Streamline sl;
    for (auto k = 1u; k < input.time_parameters.size(); ++k) {
        auto t = input.time_parameters[k];
        if (std::floor(t) > lastVolume) {
            lastVolume = std::floor(t);

            if (flowVolumes->size() <= lastVolume + 1)
                break;

            if (nextVolume) {
                currentVolume = std::move(nextVolume);
            } else {
                // create currentVolume
                currentVolume.reset(new VolumeRAMRepresentationLock(flowVolumes->at(lastVolume)));
                if(!**currentVolume) {
                    throw InvalidInputException("RAM representation not available", InvalidInputException::S_ERROR);
                }
            }

            // create nextVolume
            nextVolume.reset(new VolumeRAMRepresentationLock(flowVolumes->at(lastVolume + 1)));
        }
        auto alpha = t - lastVolume;
        SpatioTemporalSampler sampler(**currentVolume, **nextVolume, alpha, rwm, input.filterMode, worldToVoxelMatrix);

        auto velocityModifier = input.velocityUnitConversion * input.elapsedTimeBetweenVolumes;

        // move current_timeline forward one timestep:
        for (auto& interval : current_timeline) {
            for (auto i = 0u; i < interval.s.size(); ++i) {
                const float epsilon = 1e-5f; // std::numeric_limits<float>::epsilon() is not enough.
                // Position.
                tgt::vec3& r =  interval.f[i];
                if (inBounds(r)) {
                    // Velocity.
                    tgt::vec3& v = interval.v[i];

                    // Execute 4th order Runge-Kutta step.
                    tgt::vec3 k1 = v * input.stepsize * velocityModifier;
                    tgt::vec3 k2 = sampler.sample(r + (k1 * 0.5f)) * input.stepsize * velocityModifier;
                    tgt::vec3 k3 = sampler.sample(r + (k2 * 0.5f)) * input.stepsize * velocityModifier;
                    tgt::vec3 k4 = sampler.sample(r + k3) * input.stepsize * velocityModifier;
                    tgt::vec3 dr = (k1 / 6.0f) + (k2 / 3.0f) + (k3 / 3.0f) + (k4 / 6.0f);

                    r += dr;
                    if (!inBounds(r))
                        v = tgt::vec3(0);
                    else
                        v = sampler.sample(r);
                }
            }
            interval.done = false;
        }

        refine(current_timeline, input.refine, input.discontinuity);
        for (const auto& interval : current_timeline) {
            for (auto i = 0u; i < interval.s.size(); ++i) {
                auto it = activeSet.find(interval.s[i]);
                if (it == activeSet.end()) {
                    // this is part of an initially inactive curve
                    auto it_ = inactiveSet.find(interval.s[i]);
                    if (it_ == inactiveSet.end()) {
                        inactiveSet.emplace(interval.s[i], std::pair<bool, IntegralCurve>{ true, IntegralCurve{
                            false,
                            { 
                                IntegralCurvePoint{ interval.f[i], t } 
                            }, 
                            t,
                            interval.s[i] 
                        }});
                    }
                    else {
                        it_->second.second.ps.emplace_back(IntegralCurvePoint{ interval.f[i], t});
                    }
                }
                else {
                    it->second.ps.emplace_back(IntegralCurvePoint{ interval.f[i], t});
                }
            }
        }
        pointGeneration.setProgress((t - t_begin) / (t_end - t_begin));
    }

    if (pointlistOutport_.isConnected()) {
        auto points = new PointSegmentListGeometryVec3();
        output.pointlist.reset(points);

        for (const auto& C : activeSet) {
            std::vector<tgt::vec3> segment;
            for (auto p : C.second.ps) {
                segment.emplace_back(p.p);
            }
            points->addSegment(segment);
        }
        for (const auto& C : inactiveSet) {
            std::vector<tgt::vec3> segment;
            for (auto p : C.second.second.ps) {
                segment.emplace_back(p.p);
            }
            points->addSegment(segment);
        }
    }

    tgtAssert(activeSet.size() > 1, "Too few seeding points!");

    SubtaskProgressReporter surfaceGeneration(progressReporter, {0.8f, 1.0f});

    if (activeSet.size() > 1) {
        // Create mesh geometry from active and inactive set:
        auto endPtr = std::prev(activeSet.end());
        auto i = 0.f;
        for (auto it = activeSet.begin(); it != endPtr; ++it) {
            const auto& S0 = it->second;
            const auto& S1 = std::next(it)->second;

            if (S1.left_edge)
                continue;

            advanceRibbon(S0, S1, geometry, inactiveSet, { 0, 0 });
            i += 1.f;
            surfaceGeneration.setProgress(i / activeSet.size());
        }
    }
    
    return output;
}

void advanceRibbon(
        const IntegralCurve& S0, 
        const IntegralCurve& S1, 
        TriangleMeshGeometryTexCoord* geometry,
        InactiveCurveSet& inactiveSet,
        std::pair<std::size_t, std::size_t> activeEdge) {
    
    tgtAssert(S0.s < S1.s, "Seeding parameter of first curve must be smaller");
    
    while (activeEdge.first < S0.ps.size() - 1 && activeEdge.second < S1.ps.size() - 1) {
        auto left_dg  = tgt::lengthSq(S0.ps[activeEdge.first + 1].p - S1.ps[activeEdge.second].p);
        auto right_dg = tgt::lengthSq(S1.ps[activeEdge.second + 1].p - S0.ps[activeEdge.first].p);

        auto v2 = left_dg >= right_dg
                ? VertexTexCoord{ S1.ps[activeEdge.second + 1].p, { S0.s, S1.ps[activeEdge.second + 1].t }, {} } 
                : VertexTexCoord{ S0.ps[activeEdge. first + 1].p, { S1.s, S0.ps[activeEdge.first + 1].t }, {} };

        bool found = false;
        auto previous = &S0;
        for (auto current = inactiveSet.begin(); current != inactiveSet.end(); ) {
            if (current->first >= S1.s)
                // parameter out of range
                break;

            bool deleted = false;
            if (current->second.first && current->first > S0.s) {
                if (current->second.second.t_begin <= S0.ps[activeEdge.first + 1].t) {
                    if (!found)
                        v2 = VertexTexCoord{ S0.ps[activeEdge.first + 1].p, { S0.s, S0.ps[activeEdge.first + 1].t }, {} };
                    
                    auto c0 = VertexTexCoord{ current->second.second.ps[0].p, { current->first, current->second.second.t_begin }, {} };
                    Triangle<VertexTexCoord> face;
                    face.v_[0] = VertexTexCoord{ S0.ps[activeEdge.first].p, { S0.s, S0.ps[activeEdge.first].t }, {} };
                    face.v_[1] = c0;
                    face.v_[2] = v2;
                    geometry->addTriangle(face);

                    v2 = c0;

                    advanceRibbon(
                        *previous,
                        current->second.second,
                        geometry,
                        inactiveSet,
                        { activeEdge.first + 1, 0 }
                    );
                    previous = &current->second.second;

                    found = true;
                    //deleted = true;
                    current->second.first = false;
                    //current = inactiveSet.erase(current);
                }
            }

            if (!deleted)
                current = std::next(current);
        }
        if (found) {
            Triangle<VertexTexCoord> face;
            face.v_[0] = VertexTexCoord{ S1.ps[activeEdge.second].p,     {S1.s, S1.ps[activeEdge.second].t }, {} };
            face.v_[1] = VertexTexCoord{ S1.ps[activeEdge.second + 1].p, {S1.s, S1.ps[activeEdge.second + 1].t }, {} };
            face.v_[2] = v2;
            geometry->addTriangle(face);
            face.v_[0] = VertexTexCoord{ S0.ps[activeEdge.first].p,  { S0.s, S0.ps[activeEdge.first].t }, {} };
            face.v_[1] = VertexTexCoord{ S1.ps[activeEdge.second].p, { S1.s,  S1.ps[activeEdge.second].t }, {} };
            face.v_[2] = v2;
            geometry->addTriangle(face);

            advanceRibbon(
                *previous,
                S1,
                geometry,
                inactiveSet,
                { 0, activeEdge.second + 1 }
            );

            return;
        }

        if (left_dg < right_dg) {
            Triangle<VertexTexCoord> face;
            face.v_[0] = VertexTexCoord{ S0.ps[activeEdge.first].p, { S0.s, S0.ps[activeEdge.first].t }, {} };
            face.v_[1] = VertexTexCoord{ S1.ps[activeEdge.second].p, { S1.s, S1.ps[activeEdge.second].t }, {} };
            face.v_[2] = VertexTexCoord{ S0.ps[activeEdge.first + 1].p, { S0.s, S0.ps[activeEdge.first + 1].t }, {} };
            geometry->addTriangle(face);
            face.v_[0] = VertexTexCoord{ S1.ps[activeEdge.second].p, { S1.s, S1.ps[activeEdge.second].t }, {} };
            face.v_[1] = VertexTexCoord{ S1.ps[activeEdge.second + 1].p, { S1.s, S1.ps[activeEdge.second + 1].t }, {} };
            face.v_[2] = VertexTexCoord{ S0.ps[activeEdge.first + 1].p, { S0.s, S0.ps[activeEdge.first + 1].t }, {} };
            geometry->addTriangle(face);
        }
        else {
            Triangle<VertexTexCoord> face;
            face.v_[0] = VertexTexCoord{ S0.ps[activeEdge.first].p, { S0.s, S0.ps[activeEdge.first].t }, {} };
            face.v_[1] = VertexTexCoord{ S1.ps[activeEdge.second].p, { S1.s, S1.ps[activeEdge.second].t }, {} };
            face.v_[2] = VertexTexCoord{ S1.ps[activeEdge.second + 1].p, { S1.s, S1.ps[activeEdge.second + 1].t }, {} };
            geometry->addTriangle(face);
            face.v_[0] = VertexTexCoord{ S0.ps[activeEdge.first].p, { S0.s, S0.ps[activeEdge.first].t }, {} };
            face.v_[1] = VertexTexCoord{ S1.ps[activeEdge.second + 1].p, { S1.s, S1.ps[activeEdge.second + 1].t }, {} };
            face.v_[2] = VertexTexCoord{ S0.ps[activeEdge.first + 1].p, { S0.s, S0.ps[activeEdge.first + 1].t }, {} };
            geometry->addTriangle(face);
        }

        activeEdge = { activeEdge.first + 1, activeEdge.second + 1 };
    }

}

void PathSurfaceCreator::processComputeOutput(ComputeOutput output) {
    geometryOutport_.setData(output.geometry.release(), true);
    if (pointlistOutport_.isConnected())
        pointlistOutport_.setData(output.pointlist.release(), true);
}

bool PathSurfaceCreator::s_QRefineInterNodeDistance(const dense_output_type& f0, const dense_output_type& f1, float threshold)
{
    return tgt::lengthSq(f1 - f0) >= threshold;
}

bool PathSurfaceCreator::s_QRefineInterSegmentAngle(const dense_output_type& f0, const dense_output_type& f1, const dense_output_type& f2, float threshold)
{
    auto v0 = f1 - f0;
    auto v1 = f2 - f0;
    auto angle = acos(tgt::dot(v0, v1) / (tgt::length(v0) * tgt::length(v1)));
    return angle >= threshold;
}

bool PathSurfaceCreator::s_QRefineTriangleArea(const dense_output_type& f0, const dense_output_type& f1, const dense_output_type& f2, float threshold)
{
    auto v0 = f1 - f0;
    auto v1 = f2 -f1;
    auto cross = tgt::cross(v0, v1);
    auto areaSq = 0.25f * tgt::lengthSq(cross);
    return areaSq >= threshold;
}

bool PathSurfaceCreator::s_QDiscontinuity(
    float x0, float x1, float x2,
    const dense_output_type& f0, const dense_output_type& f1, const dense_output_type& f2,
    float threshold) 
{
    auto df01 = (f1 - f0) / (x1 - x0);
    auto df12 = (f2 - f1) / (x2 - x1);
    auto ddf = (df12 - df01) / ((x2 + x1) / 2 - (x1 - x0) / 2);

    return 
        tgt::lengthSq(df01) >= threshold ||
        tgt::lengthSq(df12) >= threshold ||
        tgt::lengthSq(ddf) >= threshold;
}

}
