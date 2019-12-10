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

#include "multivolumecrosssectionanalyzer.h"

#include "voreen/core/datastructures/geometry/meshgeometry.h"
#include "voreen/core/datastructures/geometry/meshlistgeometry.h"
#include "voreen/core/datastructures/geometry/facegeometry.h"
#include "voreen/core/datastructures/callback/memberfunctioncallback.h"
#include "voreen/core/utils/stringutils.h"

#include "modules/plotting/datastructures/plotcell.h"

#include <sstream>

using tgt::vec3;

namespace voreen {

const std::string MultiVolumeCrossSectionAnalyzer::loggerCat_("voreen.experimental.MultiVolumeCrossSectionAnalyzer");

MultiVolumeCrossSectionAnalyzer::MultiVolumeCrossSectionAnalyzer()
    : VolumeProcessor()
    , volumeInport1_(Port::INPORT, "volume1", "volume1", false)
    , volumeInport2_(Port::INPORT, "volume2", "volume2", false)
    , volumeInport3_(Port::INPORT, "volume3", "volume3", false)
    , volumeInport4_(Port::INPORT, "volume4", "volume4", false)
    , plotOutport_(Port::OUTPORT, "crossSectionData")
    , slicesOutport_(Port::OUTPORT, "crossSections")
    , crossSectionNormal_("crossSectionNormal", "Cross Section Normal",
        tgt::vec3(1.f, 0.f, 0.f), tgt::vec3(-1.f), tgt::vec3(1.f))
    , slabThickness_("slabThickness", "Slab Thickness", 0.05f, 1e-4, 10.f)
    , slicesPerSlab_("slicesPerSlab", "Slices Per Slab", 1, 1, 100)
    , samplingRate_("samplingRate", "Slice Sampling Rate", 0.5f, 0.01f, 10.f)
    , slabPosition_("slabPosition", "Current Slab Position", 0.f, -10.f, 10.f)
    , normalizeOutput_("normalizeOutput", "Normalize Output", true)
    , plotData_(1, 1)
    , performAnalysis_(false)
{
    addPort(volumeInport1_);
    addPort(volumeInport2_);
    addPort(volumeInport3_);
    addPort(volumeInport4_);
    addPort(plotOutport_);
    addPort(slicesOutport_);

    crossSectionNormal_.onChange(MemberFunctionCallback<MultiVolumeCrossSectionAnalyzer>(this, &MultiVolumeCrossSectionAnalyzer::forceAnalysis));
    slabThickness_.onChange(MemberFunctionCallback<MultiVolumeCrossSectionAnalyzer>(this, &MultiVolumeCrossSectionAnalyzer::forceAnalysis));
    slicesPerSlab_.onChange(MemberFunctionCallback<MultiVolumeCrossSectionAnalyzer>(this, &MultiVolumeCrossSectionAnalyzer::forceAnalysis));
    samplingRate_.onChange(MemberFunctionCallback<MultiVolumeCrossSectionAnalyzer>(this, &MultiVolumeCrossSectionAnalyzer::forceAnalysis));

    slabThickness_.onChange(MemberFunctionCallback<MultiVolumeCrossSectionAnalyzer>(this, &MultiVolumeCrossSectionAnalyzer::thicknessChanged));

    addProperty(crossSectionNormal_);
    addProperty(slabThickness_);
    addProperty(slicesPerSlab_);
    addProperty(samplingRate_);
    addProperty(slabPosition_);
    addProperty(normalizeOutput_);
}

MultiVolumeCrossSectionAnalyzer::~MultiVolumeCrossSectionAnalyzer() {}

Processor* MultiVolumeCrossSectionAnalyzer::create() const {
    return new MultiVolumeCrossSectionAnalyzer();
}

void MultiVolumeCrossSectionAnalyzer::initialize() {
    VolumeProcessor::initialize();
}

void MultiVolumeCrossSectionAnalyzer::deinitialize() {
    occurrences_.clear();
    maxOccurrences_.clear();
    plotOutport_.setData(0);
    slicesOutport_.setData(0);

    VolumeProcessor::deinitialize();
}

bool MultiVolumeCrossSectionAnalyzer::isReady() const {

    // at least one inport and one outport needs to be connected
    bool ready = true;
    ready &= plotOutport_.isReady() || slicesOutport_.isReady();
    ready &= (volumeInport1_.isReady() || volumeInport2_.isReady() ||
              volumeInport3_.isReady() || volumeInport4_.isReady() );

    return ready;
}

void MultiVolumeCrossSectionAnalyzer::beforeProcess() {
    if (volumeInport1_.hasChanged() || volumeInport2_.hasChanged() ||
        volumeInport3_.hasChanged() || volumeInport4_.hasChanged() ) {
        performAnalysis_ = true;
    }
}

void MultiVolumeCrossSectionAnalyzer::process() {

    if (performAnalysis_) {
        analyzeCrossSections();
        performAnalysis_ = false;
    }

    putOutPlotData();
}

void MultiVolumeCrossSectionAnalyzer::analyzeCrossSections() {

    tgtAssert(isInitialized(), "not initialized");

    // retrieve volumes to be analyzed from inports
    std::vector<const VolumeBase*> volumes;
    if (volumeInport1_.isReady())
        volumes.push_back(volumeInport1_.getData());
    if (volumeInport2_.isReady())
        volumes.push_back(volumeInport2_.getData());
    if (volumeInport3_.isReady())
        volumes.push_back(volumeInport3_.getData());
    if (volumeInport4_.isReady())
        volumes.push_back(volumeInport4_.getData());

    // cross section normal
    tgt::vec3 normal = tgt::normalize(crossSectionNormal_.get());

    // extent of the common bounding box to be sliced in world coords
    tgt::vec3 llf, urb;
    float minDist, maxDist;
    computeMinMaxDistToPlane(volumes, normal, minDist, maxDist);

    // number of slabs (cross sections) to analyze and the number of slices per slab
    int numSlabs = tgt::iceil((maxDist - minDist) / slabThickness_.get());
    int numSlicesPerSlab = slicesPerSlab_.get();
    float sliceStepping = slabThickness_.get() / numSlicesPerSlab;

    // distance of sampling points within a slice (slice sampling rate)
    float innerSliceSpacing = computeWorldSpacing(volumes); //< minimum spacing of the volumes in world coords
    innerSliceSpacing /= samplingRate_.get();               //< modify "native" spacing by user-specified sampling rate
    tgtAssert(innerSliceSpacing > 0.f, "invalid innerSliceSpacing");

    // store volume's bounding boxes and worldToVoxelMatrices
    std::vector<MeshGeometry> boundingBoxes;
    std::vector<tgt::mat4> worldToVoxelMatrices;
    for (size_t i=0; i<volumes.size(); i++) {
        boundingBoxes.push_back(volumes.at(i)->getBoundingBox(true));
        worldToVoxelMatrices.push_back(volumes.at(i)->getWorldToVoxelMatrix());
    }

    // initialize result buffer
    occurrences_.clear();
    occurrences_ = std::vector< std::vector<int> >(numSlabs, std::vector<int>(volumes.size(), 0));
    maxOccurrences_.clear();
    maxOccurrences_ = std::vector<int>(volumes.size(), 0);

    // output geometry
    slicesGeometry_.clear();

    //
    // Iterate over cross sections and sample them
    //
    for (int slabID=0; slabID < numSlabs; slabID++) {
        // iterate over slab slices
        for (int sliceID=0; sliceID < numSlicesPerSlab; sliceID++) {
            float clipPlaneDist = minDist + slabID*slabThickness_.get() - slabThickness_.get()/2.f + (sliceID+0.5f)*sliceStepping;
            // iterate over volumes and sample them at the clip face
            for (size_t volID = 0; volID < volumes.size(); volID++) {
                const VolumeRAM* volume = volumes.at(volID)->getRepresentation<VolumeRAM>();
                MeshGeometry bb(boundingBoxes[volID]);
                MeshGeometry clipFace;
                bb.clip(tgt::plane(normal, clipPlaneDist), clipFace);
                if (!clipFace.empty()) {
                    slicesGeometry_.addMesh(MeshGeometry(clipFace));
                    sampleSlice(volume, worldToVoxelMatrices.at(volID), normal, clipFace, innerSliceSpacing, occurrences_[slabID][volID]);
                    maxOccurrences_[volID] = std::max(maxOccurrences_[volID], occurrences_[slabID][volID]);
                }
            }
        }
    }

    slicesOutport_.setData(&slicesGeometry_, false);
}


void MultiVolumeCrossSectionAnalyzer::putOutPlotData() {

    tgtAssert(isInitialized(), "not initialized");

    // retrieve volumes to be analyzed from inports
    std::vector<const VolumeBase*> volumes;
    if (volumeInport1_.isReady())
        volumes.push_back(volumeInport1_.getData());
    if (volumeInport2_.isReady())
        volumes.push_back(volumeInport2_.getData());
    if (volumeInport3_.isReady())
        volumes.push_back(volumeInport3_.getData());
    if (volumeInport4_.isReady())
        volumes.push_back(volumeInport4_.getData());

    tgtAssert(maxOccurrences_.size() == volumes.size(), "maxOccurrences vector size does not match number of volumes");

    bool normalize = normalizeOutput_.get();

    float maxValue;
    if (normalize)
        maxValue = 1.f;
    else {
        int imax = 0;
        for (size_t i=0; i<maxOccurrences_.size(); i++)
            imax = std::max(imax, maxOccurrences_[i]);
        maxValue = static_cast<float>(imax);
    }
    //maxValue *= 1.05f;

    // cross section normal
    tgt::vec3 normal = tgt::normalize(crossSectionNormal_.get());

    // construct plot data object from result buffer
    unsigned int currentSlab = getCurrentSlabID(volumes, normal, slabPosition_.get(), slabThickness_.get());
    plotData_.reset(1, static_cast<int>(volumes.size()) + 1);
    plotData_.setColumnLabel(0, "Slab");
    for (size_t volID = 0; volID < volumes.size(); volID++)
        plotData_.setColumnLabel(static_cast<int>(volID) + 1, "Volume " + itos(volID));
    plotData_.setColumnLabel(static_cast<int>(volumes.size()) + 1, "Slab Position");

    // generate one plot row per slab
    for (size_t slabID = 0; slabID < occurrences_.size(); slabID++) {

        tgtAssert(occurrences_[slabID].size() == volumes.size(), "invalid size of occurrences buffer");

        if (slabID != currentSlab) {
            // add row for the slab
            std::vector<PlotCellValue> cells = std::vector<PlotCellValue>();
            PlotCellValue cellValue(static_cast<plot_t>(slabID));
            cells.push_back(cellValue);
            for (size_t volID = 0; volID < volumes.size(); volID++) {
                float value = static_cast<float>(occurrences_[slabID][volID]);
                if (normalize)
                    value /= static_cast<float>(maxOccurrences_[volID]);
                cells.push_back(value);
            }
            cells.push_back(0.f);
            plotData_.insert(cells);
        }
        else {
            //
            // mark current slab by peak function
            //
            std::vector<PlotCellValue> cells = std::vector<PlotCellValue>();

            // rising edge
            if (slabID > 0) {
                cells.push_back(PlotCellValue(slabID - 0.5f - 1e-4f));
                for (size_t volID = 0; volID < volumes.size(); volID++) {
                    float value = static_cast<float>((occurrences_[slabID][volID] + occurrences_[slabID-1][volID]) / 2.f);
                    if (normalize)
                        value /= static_cast<float>(maxOccurrences_[volID]);
                    cells.push_back(value);
                }
                cells.push_back(0);
                plotData_.insert(cells);

                cells.clear();
                cells.push_back(PlotCellValue(slabID - 0.5f));
                for (size_t volID = 0; volID < volumes.size(); volID++) {
                    float value = static_cast<float>((occurrences_[slabID][volID] + occurrences_[slabID-1][volID]) / 2.f);
                    if (normalize)
                        value /= static_cast<float>(maxOccurrences_[volID]);
                    cells.push_back(value);
                }
                cells.push_back(maxValue);
                plotData_.insert(cells);
            }

            // center
            cells.clear();
            cells.push_back(PlotCellValue(static_cast<plot_t>(slabID)));
            for (size_t volID = 0; volID < volumes.size(); volID++) {
                float value = static_cast<float>(occurrences_[slabID][volID]);
                if (normalize)
                    value /= static_cast<float>(maxOccurrences_[volID]);
                cells.push_back(value);
            }
            cells.push_back(maxValue);
            plotData_.insert(cells);

            // falling edge
            if (slabID < occurrences_.size()-1) {
                cells.clear();
                cells.push_back(PlotCellValue(slabID + 0.5f));
                for (size_t volID = 0; volID < volumes.size(); volID++) {
                    float value = static_cast<float>((occurrences_[slabID][volID] + occurrences_[slabID+1][volID]) / 2.f);
                    if (normalize)
                        value /= static_cast<float>(maxOccurrences_[volID]);
                    cells.push_back(value);
                }
                cells.push_back(maxValue);
                plotData_.insert(cells);

                cells.clear();
                cells.push_back(PlotCellValue(slabID + 0.5f + 1e-4f));
                for (size_t volID = 0; volID < volumes.size(); volID++) {
                    float value = static_cast<float>((occurrences_[slabID][volID] + occurrences_[slabID+1][volID]) / 2.f);
                    if (normalize)
                        value /= static_cast<float>(maxOccurrences_[volID]);
                    cells.push_back(value);
                }
                cells.push_back(0.f);
                plotData_.insert(cells);
            }
        }
    }

    // output data
    plotOutport_.setData(&plotData_);
}

void MultiVolumeCrossSectionAnalyzer::computeMinMaxDistToPlane(const std::vector<const VolumeBase*>& volumes,
    const tgt::vec3& normal, float& minDist, float& maxDist) const {

    // min/max distances to the most outer planes intersecting the volumes' bounding boxes
    minDist = 1e6f;
    maxDist = -1e6f;
    for (size_t v=0; v<volumes.size(); v++) {
        tgtAssert(volumes.at(v), "no handle");
        const VolumeBase* volume = volumes.at(v);
        MeshGeometry boundingBox = MeshGeometry::createCube(volume->getLLF(), volume->getURB());
        boundingBox.transform(volume->getPhysicalToWorldMatrix());
        for (size_t faceID=0; faceID < boundingBox.getFaceCount(); faceID++) {
            for (size_t vertexID=0; vertexID < boundingBox.getFace(faceID).getVertexCount(); vertexID++) {
                tgt::vec3 vertex = boundingBox.getFace(faceID).getVertex(vertexID).getCoords();
                minDist = std::min(minDist, tgt::dot(vertex, normal));
                maxDist = std::max(maxDist, tgt::dot(vertex, normal));
            }
        }
    }

    // correct min/max distances so that their distances to the slab position are multiples of the slab thickness
    float distModMin = fmod(slabPosition_.get() - minDist, slabThickness_.get());
    minDist -= (slabThickness_.get() - distModMin);
    float distModMax = fmod(maxDist - slabPosition_.get(), slabThickness_.get());
    maxDist += (slabThickness_.get() - distModMax);
}

float MultiVolumeCrossSectionAnalyzer::computeWorldSpacing(const std::vector<const VolumeBase*>& volumes) const {
    float spacing = 1e6f;
    for (size_t v=0; v<volumes.size(); v++) {
        tgtAssert(volumes.at(v), "no handle");
        tgt::mat4 voxelToWorldMatrix = volumes.at(v)->getVoxelToWorldMatrix();
        tgt::vec3 spacingVec = tgt::abs(voxelToWorldMatrix*tgt::vec3(1.f) - voxelToWorldMatrix*tgt::vec3(0.f));
        spacing = std::min(spacing, tgt::min(spacingVec));
    }
    return spacing;
}

void MultiVolumeCrossSectionAnalyzer::sampleSlice(const VolumeRAM* volume, const tgt::mat4& worldToVoxelMatrix,
    const tgt::vec3& normal, const MeshGeometry& clipFace, float spacing, int& occurrence) const {

    occurrence = 0;

    tgtAssert(volume, "no volume");
    tgtAssert(!clipFace.empty(), "clip face empty");
    if (clipFace.getVertexCount() < 3) {
        LWARNING("Malformed clip face");
        return;
    }
    FaceGeometry face = clipFace.getFace(0);
    if (clipFace.getFaceCount() > 1) {
        LWARNING("Multiple clip faces. Considering the first one only.");
    }

    tgt::vec3 first = face.getVertex(0).getCoords();

    // determine baseline and up vector of clip face for sampling
    tgt::vec3 baseline = first - face.getVertex(1).getCoords();
    tgtAssert(baseline != tgt::vec3(0.f), "malformed clip face");
    tgt::vec3 upVector = tgt::cross(baseline, normal);
    baseline = tgt::normalize(baseline);
    upVector = tgt::normalize(upVector);

    // determine max distances of clip face vertices from first vertex
    // in both directions (baseline and updir)
    float maxBaseDist = 0.f;
    float maxUpDist = 0.f;
    for (size_t i=1; i<face.getVertexCount(); i++) {
        float baseDist = tgt::dot(face.getVertex(i).getCoords() - first, baseline);
        float upDist = tgt::dot(face.getVertex(i).getCoords() - first, upVector);
        if (std::abs(baseDist) > std::abs(maxBaseDist))
            maxBaseDist = baseDist;
        if (std::abs(upDist) > std::abs(maxUpDist))
            maxUpDist = upDist;
    }
    if (maxBaseDist < 0.f) {
        baseline = -baseline;
        maxBaseDist = -maxBaseDist;
    }
    if (maxUpDist < 0.f) {
        upVector = -upVector;
        maxUpDist = -maxUpDist;
    }

    // sample 2D bounding box of clip face
    for (float baseOffset = 0.f; baseOffset <= maxBaseDist; baseOffset += spacing) {
        for (float upOffset = 0.f; upOffset <= maxUpDist; upOffset += spacing) {
            tgt::vec3 vertex = first + baseOffset*baseline + upOffset*upVector;
            tgt::ivec3 voxel = tgt::iround(worldToVoxelMatrix*vertex);
            if (tgt::hand(tgt::greaterThanEqual(voxel, tgt::ivec3(0))) && tgt::hand(tgt::lessThan(voxel, tgt::ivec3(volume->getDimensions())))) {
                float intensity = volume->getVoxelNormalized(voxel);
                if (intensity > 0.f)
                    occurrence++;
            }
        }
    }

}

int MultiVolumeCrossSectionAnalyzer::getCurrentSlabID(const std::vector<const VolumeBase*>& volumes, const tgt::vec3& normal,
        float slabPosition, float slabThickness) const {
    float minDist, maxDist;
    computeMinMaxDistToPlane(volumes, normal, minDist, maxDist);
    return tgt::ifloor((slabPosition - minDist) / slabThickness);
}

void MultiVolumeCrossSectionAnalyzer::thicknessChanged() {
    slabPosition_.setStepping(slabThickness_.get());
}

void MultiVolumeCrossSectionAnalyzer::forceAnalysis()  {
    performAnalysis_ = true;
    invalidate();
}


} // namespace
