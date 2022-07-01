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

#include "flowprofilestacking.h"

#include "voreen/core/ports/conditions/portconditionvolumelist.h"

#include "voreen/core/datastructures/geometry/pointsegmentlistgeometry.h"

#include "modules/plotting/datastructures/plotcell.h"
#include "modules/plotting/datastructures/plotdata.h"
#include "modules/plotting/datastructures/plotrow.h"

#include "../../utils/utils.h"
#include "voreen/core/utils/statistics.h"

namespace voreen {

const std::string FlowProfileStacking::loggerCat_("voreen.flowsimulation.FlowProfileStacking");

FlowProfileStacking::FlowProfileStacking()
    : AsyncComputeProcessor()
    , volumeListPort_(Port::INPORT, "input.volumelist", "Volume List Port")
    , regionOfInterestPort_(Port::INPORT, "input.parameter", "Region Of Interest Port")
    , outport_(Port::OUTPORT,"output.plot", "Plot Port")
    , stackingMode_("stackingMode", "Stacking Mode")
    , selectedSegment_("selectedSegment", "Selected Segment", 0, 0, 0, Processor::VALID)
    , r_("r_", "#samples R", 16, 2, 256)
    , a_("a_", "#samples A", 32, 1, 256)
    , z_("z_", "#samples Z", 32, 1, 256)
{
    addPort(volumeListPort_);
    //volumeListPort_.addCondition(new PortConditionVolumeListEnsemble());
    ON_CHANGE(volumeListPort_, FlowProfileStacking, onInputDataChange);
    addPort(regionOfInterestPort_);
    ON_CHANGE(regionOfInterestPort_, FlowProfileStacking, onInputDataChange);
    addPort(outport_);

    addProperty(stackingMode_);
    stackingMode_.addOption("modeMean", "Stacking of means");
    stackingMode_.addOption("modeValue", "Stacking of values");
    addProperty(selectedSegment_);
    ON_CHANGE(selectedSegment_, FlowProfileStacking, onSelectedSegmentChange);
    addProperty(r_);
    r_.setTracking(false);
    addProperty(a_);
    a_.setTracking(false);
    addProperty(z_);
    z_.setTracking(false);
}

void FlowProfileStacking::beforeProcess() {
    AsyncComputeProcessor::beforeProcess();
}

FlowProfileStackingInput FlowProfileStacking::prepareComputeInput() {

    auto volumes = volumeListPort_.getThreadSafeData();
    if(!volumes || volumes->empty()) {
        throw InvalidInputException("No volumes", InvalidInputException::S_ERROR);
    }

    auto geometry = dynamic_cast<const PointSegmentListGeometryVec3*>(regionOfInterestPort_.getData());
    if(!geometry) {
        throw InvalidInputException("No geometry", InvalidInputException::S_ERROR);
    }

    auto segments = geometry->getData();
    auto geometryTransformationMatrix = geometry->getTransformationMatrix();

    auto resampleSegments = std::vector<std::vector<tgt::vec3>>();
    resampleSegments.reserve(segments.size());
    for(auto iter = segments.rbegin(); iter != segments.rend(); iter++) {
        auto segment = *iter;
        for(auto& pos : segment) {
            pos = geometryTransformationMatrix * pos;
        }
        resampleSegments.emplace_back(segment);
    }

    return FlowProfileStackingInput {
            std::move(volumes),
            std::move(resampleSegments),
            r_.get(),
            a_.get(),
            z_.get(),
            stackingMode_.get() == "modeValue"
    };
}

FlowProfileStackingOutput FlowProfileStacking::compute(FlowProfileStackingInput input, ProgressReporter& progressReporter) const {

    auto volumes = std::move(input.volumes);
    auto segments = std::move(input.segments);

    const int R = input.r_;
    const int A = input.a_;
    const int Z = input.z_;
    bool stackValues = input.values_;

    std::vector<std::unique_ptr<PlotData>> output;
    output.reserve(segments.size());
    std::vector<std::vector<std::vector<plot_t>>> rows(segments.size());
    size_t numColumns = stackValues ? A*Z : volumes->size();
    for(size_t j = 0; j < segments.size(); j++) {
        std::unique_ptr<PlotData> plotData(new PlotData(1, numColumns));
        plotData->setColumnLabel(0, "point");
        if (!stackValues) {
            for (size_t i = 0; i < volumes->size(); i++) {
                plotData->setColumnLabel(i + 1, volumes->at(i)->getModality().getName());
            }
        }
        else {
            for (size_t i = 0; i < numColumns; i++) {
                plotData->setColumnLabel(i + 1, std::to_string(i));
            }
        }
        output.emplace_back(std::move(plotData));

        rows[j].resize(R*2);
    }

    progressReporter.setProgress(0.0f);

    // Iterate volumes first to reduce disk io.
    size_t numVolumes = stackValues ? 1 : volumes->size();
    for(size_t i = 0; i < numVolumes; i++) {
        const VolumeBase* volume = volumes->at(i);
        RealWorldMapping rwm = volume->getRealWorldMapping();
        tgt::mat4 worldToVoxelMatrix = volume->getWorldToVoxelMatrix();
        VolumeRAMRepresentationLock lock(volume);

        for(size_t j = 0; j < segments.size(); j++) {
            auto& segment = segments.at(j);

            // Ali wants the first and last 20% of the domain to be removed.
            auto lower_z = static_cast<float>(lock->getDimensions().z * 0.2f);
            auto upper_z = static_cast<float>(lock->getDimensions().z * 0.8f);

            // Gather samples.
#ifdef VRN_MODULE_OPENMP
            #pragma omp parallel for
#endif
            for(int r=-R; r<R; r++) {
                Statistics stats;
                for(int a=0; a<A; a++) {
                    float radius = 1.0f;
                    tgt::vec3 center = segment.front();
                    if(segment.size() > 1) {
                        radius = tgt::distance(segment.front(), segment.back()) / 2.0f;
                        center = (segment.front() + segment.back()) / 2.0f;
                    }
                    float factor = radius / (R-1);

                    float angle = tgt::PIf*a/A;
                    tgt::vec3 p = center + (r * factor) * tgt::vec3(std::cos(angle), std::sin(angle), 0);
                    tgt::vec3 point = worldToVoxelMatrix * p;

                    for (size_t z = 0; z < Z; z++) {
                        point.z = z * (upper_z - lower_z) / Z + lower_z;

                        tgt::vec3 sample = tgt::vec3::zero;
                        for (size_t channel = 0; channel < lock->getNumChannels(); channel++) {
                            sample[channel] = rwm.normalizedToRealWorld(lock->getVoxelNormalizedLinear(point, channel));
                        }
                        float magnitude = tgt::length(sample);
                        if(stackValues)
                            rows[j][r+R].push_back(magnitude);
                        else
                            stats.addSample(magnitude);
                    }
                }
                if(!stackValues)
                    rows[j][r+R].push_back(stats.getMean());
            }
        }

        progressReporter.setProgress((i+1.0f) / volumes->size());
    }

    for(size_t j = 0; j < segments.size(); j++) {
        auto& segment = segments.at(j);
        for (size_t k = 0; k < R*2; k++) {
            std::vector<PlotCellValue> values;
            values.emplace_back(k);
            for(size_t i = 0; i < rows[j][k].size(); i++) {
                values.emplace_back(PlotCellValue(rows[j][k][i]));
            }
            output[j]->insert(values);
        }
    }

    progressReporter.setProgress(1.0f);

    return FlowProfileStackingOutput {
            std::move(output)
    };
}

void FlowProfileStacking::processComputeOutput(FlowProfileStackingOutput output) {
    plotData_ = std::move(output.output);
    onSelectedSegmentChange();
}

void FlowProfileStacking::onSelectedSegmentChange() {
    if(plotData_.empty()) {
        return;
    }

    if(selectedSegment_.get() >= plotData_.size()) {
        return;
    }

    auto* ptr = plotData_[selectedSegment_.get()].get();
    outport_.setData(ptr, false);
}

void FlowProfileStacking::onInputDataChange() {
    outport_.setData(nullptr);
    plotData_.clear();

    if(auto data = dynamic_cast<const PointSegmentListGeometryVec3*>(regionOfInterestPort_.getData())) {
        if(data->getNumSegments() == 0) {
            selectedSegment_.setReadOnlyFlag(true);
        }
        else {
            selectedSegment_.setReadOnlyFlag(false);
            selectedSegment_.setMaxValue(data->getNumSegments()-1);
        }
    }
}

}
