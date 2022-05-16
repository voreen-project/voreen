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

#include "regionofinterestanalysis.h"

#include "voreen/core/ports/conditions/portconditionvolumelist.h"

#include "voreen/core/datastructures/geometry/pointsegmentlistgeometry.h"

#include "modules/plotting/datastructures/plotcell.h"
#include "modules/plotting/datastructures/plotdata.h"
#include "modules/plotting/datastructures/plotrow.h"

#include "../../utils/utils.h"

namespace voreen {



static std::vector<tgt::vec3> resampleSegment(const std::vector<tgt::vec3>& segment) {

    std::vector<tgt::vec3> resampled;

    // Append the first element separately.
    resampled.push_back(segment.front());

    // In case we have more than two samples, we achieve the result by linear interpolation.
    if(segment.size() > 2) {

        std::vector<float> distances(segment.size());

        float totalLength = 0.0f;
        distances[0] = 0.0f;
        for(size_t i = 1; i < segment.size(); i++) {
            totalLength += tgt::distance(segment[i-1], segment[i]);
            distances[i] = totalLength;
        }

        const float segmentLength = totalLength / (segment.size() - 1);
        for (size_t i = 1, pos = 0; i < segment.size() - 1; i++) {

            while (distances[pos+1] < segmentLength * i) pos++;

            const float t = tgt::clamp((segmentLength * i - distances[pos]) / (distances[pos + 1] - distances[pos]), 0.0f, 1.0f);

            tgt::vec3 position = segment[pos] * (1.0f - t) + segment[pos + 1] * t;

            resampled.push_back(position);
        }
    }

    // Append the last element separately.
    resampled.push_back(segment.back());

    return resampled;
}



const std::string RegionOfInterestAnalysis::loggerCat_("voreen.flowsimulation.RegionOfInterestAnalysis");

RegionOfInterestAnalysis::RegionOfInterestAnalysis()
    : AsyncComputeProcessor()
    , volumeListPort_(Port::INPORT, "input.volumelist", "Volume List Port")
    , regionOfInterestPort_(Port::INPORT, "input.parameter", "Region Of Interest Port")
    , outport_(Port::OUTPORT,"output.plot", "Plot Port")
    , selectedSegment_("selectedSegment", "Selected Segment", 0, 0, 0, Processor::VALID)
{
    addPort(volumeListPort_);
    //volumeListPort_.addCondition(new PortConditionVolumeListEnsemble());
    ON_CHANGE(volumeListPort_, RegionOfInterestAnalysis, onInputDataChange);
    addPort(regionOfInterestPort_);
    ON_CHANGE(regionOfInterestPort_, RegionOfInterestAnalysis, onInputDataChange);
    addPort(outport_);

    addProperty(selectedSegment_);
    ON_CHANGE(selectedSegment_, RegionOfInterestAnalysis, onSelectedSegmentChange);
}

void RegionOfInterestAnalysis::beforeProcess() {
    AsyncComputeProcessor::beforeProcess();
}

RegionOfInterestAnalysisInput RegionOfInterestAnalysis::prepareComputeInput() {

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
        auto resampledSegment = resampleSegment(*iter);
        for(auto& pos : resampledSegment) {
            pos = geometryTransformationMatrix * pos;
        }
        resampleSegments.emplace_back(resampledSegment);
    }

    return RegionOfInterestAnalysisInput {
            std::move(volumes),
            std::move(resampleSegments),
    };
}

RegionOfInterestAnalysisOutput RegionOfInterestAnalysis::compute(RegionOfInterestAnalysisInput input, ProgressReporter& progressReporter) const {

    auto volumes = std::move(input.volumes);
    auto segments = std::move(input.segments);

    std::vector<std::unique_ptr<PlotData>> output;
    output.reserve(segments.size());
    std::vector<std::vector<std::vector<plot_t>>> rows(segments.size());
    for(size_t j = 0; j < segments.size(); j++) {
        std::unique_ptr<PlotData> plotData(new PlotData(1, volumes->size()));
        plotData->setColumnLabel(0, "point");

        for(size_t i = 0; i < volumes->size(); i++) {
            plotData->setColumnLabel(i+1, volumes->at(i)->getModality().getName());
        }
        output.emplace_back(std::move(plotData));

        auto& segment = segments.at(j);
        rows[j].resize(segment.size()); // Add id column.
    }

    progressReporter.setProgress(0.0f);

    // Iterate volumes first to reduce disk io.
    for(size_t i = 0; i < volumes->size(); i++) {
        const VolumeBase* volume = volumes->at(i);
        RealWorldMapping rwm = volume->getRealWorldMapping();
        tgt::mat4 worldToVoxelMatrix = volume->getWorldToVoxelMatrix();
        VolumeRAMRepresentationLock lock(volume);

        for(size_t j = 0; j < segments.size(); j++) {
            auto& segment = segments.at(j);

            // Gather samples.
            for(size_t k = 0; k < segment.size(); k++) {
                tgt::vec3 point = worldToVoxelMatrix * segment[k];

                tgt::vec3 sample = tgt::vec3::zero;
                for(size_t channel=0; channel < lock->getNumChannels(); channel++) {
                    sample[channel] = rwm.normalizedToRealWorld(lock->getVoxelNormalizedLinear(point, channel));
                }
                float magnitude = tgt::length(sample);
                rows[j][k].push_back(magnitude);
            }
        }

        progressReporter.setProgress((i+1.0f) / volumes->size());
    }

    for(size_t j = 0; j < segments.size(); j++) {
        auto& segment = segments.at(j);
        for (size_t k = 0; k < segment.size(); k++) {
            std::vector<PlotCellValue> values;
            values.emplace_back(k);
            for(size_t i = 0; i < volumes->size(); i++) {
                values.emplace_back(PlotCellValue(rows[j][k][i]));

            }
            output[j]->insert(values);
        }
    }

    progressReporter.setProgress(1.0f);

    return RegionOfInterestAnalysisOutput {
            std::move(output)
    };
}

void RegionOfInterestAnalysis::processComputeOutput(RegionOfInterestAnalysisOutput output) {
    plotData_ = std::move(output.output);
    onSelectedSegmentChange();
}

void RegionOfInterestAnalysis::onSelectedSegmentChange() {
    if(plotData_.empty()) {
        return;
    }

    if(selectedSegment_.get() >= plotData_.size()) {
        return;
    }

    auto* ptr = plotData_[selectedSegment_.get()].get();
    outport_.setData(ptr, false);
}

void RegionOfInterestAnalysis::onInputDataChange() {
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
