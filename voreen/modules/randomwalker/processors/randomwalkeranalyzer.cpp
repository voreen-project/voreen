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

#include "randomwalkeranalyzer.h"

#include "voreen/core/datastructures/volume/volumeram.h"
#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/processors/processorwidget.h"
#include "voreen/core/datastructures/geometry/pointlistgeometry.h"
#include "voreen/core/datastructures/volume/operators/volumeoperatorconvert.h"
#include "voreen/core/voreenapplication.h"

#include "tgt/vector.h"
#include "tgt/timer.h"

#ifdef VRN_MODULE_CONNEXE
#include "modules/connexe/ext/connexe/connexe.h"
#endif

inline size_t volumeCoordsToIndex(int x, int y, int z, const tgt::ivec3& dim) {
    return z*dim.y*dim.x + y*dim.x + x;
}

inline size_t volumeCoordsToIndex(const tgt::ivec3& coords, const tgt::ivec3& dim) {
    return coords.z*dim.y*dim.x + coords.y*dim.x + coords.x;
}

namespace voreen {

const std::string RandomWalkerAnalyzer::loggerCat_("voreen.RandomWalker.RandomWalkerAnalyzer");

RandomWalkerAnalyzer::RandomWalkerAnalyzer()
    : VolumeProcessor(),
    inportVolume_(Port::INPORT, "volume.input"),
    outportUncertaintyAreas_(Port::OUTPORT, "volume.analysis", "volume.analysis", false),
    outportFocusArea_(Port::OUTPORT, "volume.focusRegion", "volume.focusRegion", false),
    outportFocusAreaProbabilities_(Port::OUTPORT, "volume.focusRegionProbs", "volume.focusRegionProbs", false),
    outportFocusPoint_(Port::OUTPORT, "geometry.focusPoint"),
    uncertaintyRange_("unertaintyRange", "Uncertainty Range", tgt::vec2(0.2f, 0.8f)),
    connectivity_("connectivity", "Connectivity"),
    minUncertaintyAreaSize_("minUncertaintyAreaSize_", "Minimum Uncertainty Area Size", 1, 1, 1000),
    maxUncertaintyAreas_("maxUncertaintyAreas_", "Max Uncertainty Areas", 10000, 1, 10000),
    sliceXY_("sliceXY", "XY-Slice", 0, 0, 10000, Processor::VALID),
    sliceXZ_("sliceXZ", "XZ-Slice", 0, 0, 10000, Processor::VALID),
    sliceYZ_("sliceYZ", "YZ-Slice", 0, 0, 10000, Processor::VALID),
    zoomXY_("zoomXY", "Zoom XY", 1.f, 0.01, 1.f, Processor::VALID),
    zoomXZ_("zoomXZ", "Zoom XZ", 1.f, 0.01, 1.f, Processor::VALID),
    zoomYZ_("zoomYZ", "Zoom YZ", 1.f, 0.01, 1.f, Processor::VALID),
    focusXY_("focusXY", "XY-Focus", tgt::vec2(0.f), tgt::vec2(-10000.f), tgt::vec2(10000.f), Processor::VALID),
    focusXZ_("focusXZ", "XZ-Focus", tgt::vec2(0.f), tgt::vec2(-10000.f), tgt::vec2(10000.f), Processor::VALID),
    focusYZ_("focusYZ", "YZ-Focus", tgt::vec2(0.f), tgt::vec2(-10000.f), tgt::vec2(10000.f), Processor::VALID),
    camera_("camera", "Camera", tgt::Camera(), false, Processor::VALID),
    computeButton_("computeButton", "Compute", Processor::VALID),
    focusPoint_(-1.f),
    focusRegion_(-1),
    probVolume_(0),
    probVolumeCC_(0),
    forceUpdate_(true)
{
    eventHandler_ = new tgt::EventHandler();
    eventHandler_->addListenerToBack(this);
    timer_ = VoreenApplication::app()->createTimer(eventHandler_);

    addPort(inportVolume_);
    addPort(outportUncertaintyAreas_);
    addPort(outportFocusArea_);
    addPort(outportFocusAreaProbabilities_);
    addPort(outportFocusPoint_);

    connectivity_.addOption("6-neighborhood", "6 Neighborhood", 6);
    connectivity_.addOption("10-neighborhood", "10 Neighborhood", 10);
    connectivity_.addOption("18-neighborhood", "18 Neighborhood", 18);
    connectivity_.addOption("26-neighborhood", "26 Neighborhood", 26);
    connectivity_.select("18-neighborhood");

    addProperty(uncertaintyRange_);
    addProperty(connectivity_);
    addProperty(minUncertaintyAreaSize_);
    addProperty(maxUncertaintyAreas_);

    addProperty(sliceXY_);
    addProperty(sliceXZ_);
    addProperty(sliceYZ_);
    addProperty(zoomXY_);
    addProperty(zoomXZ_);
    addProperty(zoomYZ_);
    addProperty(focusXY_);
    addProperty(focusXZ_);
    addProperty(focusYZ_);

    addProperty(camera_);
    addProperty(computeButton_);

    sliceXY_.setReadOnlyFlag(true);
    sliceXZ_.setReadOnlyFlag(true);
    sliceYZ_.setReadOnlyFlag(true);
    camera_.setReadOnlyFlag(true);
    focusXY_.setReadOnlyFlag(true);
    focusXZ_.setReadOnlyFlag(true);
    focusYZ_.setReadOnlyFlag(true);
    computeButton_.setReadOnlyFlag(true);
}

RandomWalkerAnalyzer::~RandomWalkerAnalyzer() {
    delete timer_;
    delete eventHandler_;
}

std::string RandomWalkerAnalyzer::getProcessorInfo() const {
    return std::string("OpenCL-based computation of a Random Walker solution.");
}

Processor* RandomWalkerAnalyzer::create() const {
    return new RandomWalkerAnalyzer();
}

void RandomWalkerAnalyzer::initialize() {
    VolumeProcessor::initialize();
}

void RandomWalkerAnalyzer::deinitialize() {
    delete probVolume_;
    probVolume_ = 0;

    VolumeProcessor::deinitialize();
}

bool RandomWalkerAnalyzer::isReady() const {
    bool ready = false;
    ready |= outportUncertaintyAreas_.isConnected();
    ready |= outportFocusArea_.isConnected();
    ready |= outportFocusAreaProbabilities_.isConnected();
    ready &= inportVolume_.isReady();
    return ready;
}

const std::vector<RandomWalkerAnalyzer::UncertaintyRegion>& RandomWalkerAnalyzer::getUncertaintyRegions() const {
    return uncertaintyRegions_;
}

void RandomWalkerAnalyzer::setFocusRegion(int region) {

    focusRegion_ = region;
    setFocusPoint(uncertaintyRegions_.at(region));
    //invalidate(Processor::INVALID_RESULT);
    updateFocusAreaPropVolume();
}

void RandomWalkerAnalyzer::setFocusPoint(const UncertaintyRegion& region) {

    if (!inportVolume_.hasData())
        return;

    tgt::vec3 focus = region.centerOfMass_;

    const VolumeBase* volume = inportVolume_.getData();
    tgt::ivec3 volDim = volume->getDimensions();
    tgt::vec3 volCubeSize = volume->getCubeSize();

    // updates slices
    sliceXY_.set(0);
    sliceXY_.set(tgt::iround(focus.z));
    sliceXZ_.set(0);
    sliceXZ_.set(tgt::iround(focus.y));
    sliceYZ_.set(0);
    sliceYZ_.set(tgt::iround(focus.x));

    // update camera
    tgt::vec3 focusCorner;
    focusCorner.x = focus.x < volDim.x / 2 ? 0.f : 1.f;
    focusCorner.y = focus.y < volDim.y / 2 ? 0.f : 1.f;
    focusCorner.z = (focus.z / (volDim.z)); // + 0.01; //< volDim.z / 2 ? 0.f : 1.f;
    focusCorner -= 0.5f;
    focusCorner *= volCubeSize;
    focusCorner.xy() *= 2.0f;

    //tgt::vec3 focusPoint_ = (tgt::vec3(focus) / tgt::vec3(volDim-1)) - tgt::vec3(0.5f);
    //focusPoint_ *= volCubeSize;
    tgt::vec3 focusPoint_ = volume->getVoxelToPhysicalMatrix() * focus;

    camera_.setPosition(focusCorner);
    camera_.setFocus(focusPoint_);
    camera_.setUpVector(tgt::vec3(0.f, 0.f, 1.f));
    camera_.get().getViewMatrix();
    camera_.invalidate();

    //delete outportFocusPoint_.getData();
    PointListGeometryVec3* focusGeom = new PointListGeometryVec3();
    focusGeom->addPoint(focusPoint_);
    outportFocusPoint_.setData(focusGeom);
}


void RandomWalkerAnalyzer::computeRandomWalker() {
    computeButton_.clicked();
}

void RandomWalkerAnalyzer::process() {

    if (!inportVolume_.hasChanged() && !forceUpdate_) {
        if (getProcessorWidget())
            getProcessorWidget()->updateFromProcessor();
        return;
    }

    updateNumSlices();  // validate the currently set values and adjust them if necessary

    const VolumeRAM* inputVolume = inportVolume_.getData()->getRepresentation<VolumeRAM>();
    VolumeRAM_UInt8* volume = 0;
    if (!volume) {
        VolumeOperatorConvert voConvert;
        //VolumeOperatorConvert voConvert(inputVolume);
        //volume = new VolumeRAM_UInt8(inputVolume->getDimensions(), inputVolume->getSpacing(), inputVolume->getTransformation());
        Volume* vh = voConvert.apply<uint8_t>(inportVolume_.getData());
        volume = vh->getWritableRepresentation<VolumeRAM_UInt8>();
        vh->releaseAllRepresentations();
        delete vh;
        //LERROR("Volume of type VolumeRAM_UInt8 expected.");
        //return;
    }

    tgt::ivec2 rangeNormalized = tgt::iround(uncertaintyRange_.get() * 255.f);
    VolumeRAM_UInt8* volumeFiltered = probabilityThreshold(volume, rangeNormalized.x, rangeNormalized.y);

    if (volume != inputVolume)
        delete volume;
    volume = 0;

    if (!volumeFiltered) {
        LERROR("Filtering volume failed");
        return;
    }

    //VolumeOperatorErosion opErosion(3);
    //opErosion(volumeFiltered);

    VolumeRAM_UInt8* volumeEroded = erosion(volumeFiltered);
    if (!volumeEroded) {
        LERROR("Erosion failed");
        delete volumeFiltered;
        return;
    }
    delete volumeFiltered;
    volumeFiltered = volumeEroded;

    delete probVolume_;
    probVolume_ = volumeFiltered->clone();

    size_t numLabels;
    Volume temp(volumeFiltered, inportVolume_.getData());
    Volume* volumeCC = connectedComponents(&temp, numLabels);
    temp.releaseAllRepresentations();
    if (!volumeCC) {
        LERROR("Creating connected components volume failed");
        delete volumeFiltered;
        return;
    }

    //VolumeOperatorDilation opDilation(3);
    //opDilation(volumeCC);

    analyzeUncertaintyRegions(volumeCC, static_cast<int>(numLabels));
    stretchLabels(static_cast<VolumeRAM_UInt16*>(volumeCC->getWritableRepresentation<VolumeRAM>()));

    //printLabelInformation();

    outportUncertaintyAreas_.setData(volumeCC);
    probVolumeCC_ = volumeCC;

    delete volumeFiltered;

    updateFocusAreaPropVolume();

    /*VolumeRAM_UInt8* volumeDist = createDistanceMap(volumeFiltered);
        if (!volumeDist) {
            LERROR("Creating distance volume failed");
            delete volumeFiltered;
            return;
        }
        outportDistanceMap_.setData(new Volume(volumeDist), true);

        VolumeRAM_UInt8* volumeEdges = uncertaintyEdges(volume, volumeFiltered);
        if (!volumeEdges) {
            LERROR("Creating edge volume failed");
            delete volumeFiltered;
            return;
        }
        outportProbBoundaries_.setData(new Volume(volumeEdges), true); */

    if (getProcessorWidget())
        getProcessorWidget()->updateFromProcessor();

    forceUpdate_ = false;

}

void RandomWalkerAnalyzer::invalidate(int inv) {

    VolumeProcessor::invalidate(inv);
}

VolumeRAM_UInt8* RandomWalkerAnalyzer::probabilityThreshold(VolumeRAM_UInt8* volume, uint8_t min, uint8_t max) {

    tgtAssert(volume, "No volume");

    VolumeRAM_UInt8* result = 0;
    try {
        result = volume->clone();
    }
    catch (std::bad_alloc) {
        LERROR("bad allocation");
        return 0;
    }

    size_t numVoxels = result->getNumVoxels();
    for (size_t i=0; i<numVoxels; i++) {
        uint8_t value = result->voxel(i);
        if (value < min || value > max)
            result->voxel(i) = 0;
    }

    return result;
}

VolumeRAM_UInt8* RandomWalkerAnalyzer::createDistanceMap(VolumeRAM_UInt8* volume) {

    tgtAssert(volume, "No volume");

    VolumeRAM_UInt8* result = 0;
    try {
        result = volume->clone();
    }
    catch (std::bad_alloc) {
        LERROR("bad allocation");
        return 0;
    }

    tgt::ivec3 volDim = volume->getDimensions();
    tgt::ivec3 pos;
    for (pos.z=1; pos.z<volDim.z-1; pos.z++) {
        for (pos.y=1; pos.y<volDim.y-1; pos.y++) {
            for (pos.x=1; pos.x<volDim.x-1; pos.x++) {
                uint8_t value = volume->voxel(pos);
                if (value > 0) {
                    uint8_t v0, v1, v2, v3, v4, v5;
                    v0 = volume->voxel(pos + tgt::ivec3(1, 0, 0));
                    v1 = volume->voxel(pos + tgt::ivec3(0, 1, 0));
                    v2 = volume->voxel(pos + tgt::ivec3(0, 0, 1));

                    v3 = volume->voxel(pos + tgt::ivec3(-1, 0, 0));
                    v4 = volume->voxel(pos + tgt::ivec3(0, -1, 0));
                    v5 = volume->voxel(pos + tgt::ivec3(0, 0, -1));

                    float gradient = tgt::length(tgt::vec3(static_cast<float>(v3 - v0), static_cast<float>(v4 - v1), static_cast<float>(v5 - v2))) / 2.f;
                    gradient = std::max(gradient, 1.f) / 255.f;
                    float dev = std::abs((value/255.f) - 0.5f);
                    result->voxel(pos) = static_cast<uint8_t>(tgt::clamp(dev/gradient * 25.f, 0.f, 255.f));
                }
            }
        }
    }

    return result;
}

VolumeRAM_UInt8* RandomWalkerAnalyzer::erosion(VolumeRAM_UInt8* volume) {

    tgtAssert(volume, "No volume");

    VolumeRAM_UInt8* result = 0;
    try {
        result = volume->clone();
    }
    catch (std::bad_alloc) {
        LERROR("bad allocation");
        return 0;
    }
    tgt::ivec3 volDim = volume->getDimensions();
    tgt::ivec3 pos;
    for (pos.z=1; pos.z<volDim.z-1; pos.z++) {
        for (pos.y=1; pos.y<volDim.y-1; pos.y++) {
            for (pos.x=1; pos.x<volDim.x-1; pos.x++) {
                uint8_t value = volume->voxel(pos);
                if (value > 0) {
                    uint8_t v0, v1, v2, v3, v4, v5;
                    v0 = volume->voxel(pos + tgt::ivec3(1, 0, 0));
                    v1 = volume->voxel(pos + tgt::ivec3(0, 1, 0));
                    v2 = volume->voxel(pos + tgt::ivec3(0, 0, 1));

                    v3 = volume->voxel(pos + tgt::ivec3(-1, 0, 0));
                    v4 = volume->voxel(pos + tgt::ivec3(0, -1, 0));
                    v5 = volume->voxel(pos + tgt::ivec3(0, 0, -1));

                    if ((v0 == 0 && v3 == 0) || (v1 == 0 && v4 == 0) || (v2 == 0 && v5 == 0))
                        result->voxel(pos) = 0;
                }
            }
        }
    }

    return result;
}

VolumeRAM_UInt8* RandomWalkerAnalyzer::uncertaintyEdges(VolumeRAM_UInt8* volumeProb, VolumeRAM_UInt8* volumeFiltered) {

    tgtAssert(volumeProb, "No volume");

    VolumeRAM_UInt8* result = 0;
    try {
        result = volumeProb->clone();
    }
    catch (std::bad_alloc) {
        LERROR("bad allocation");
        return 0;
    }
    tgt::ivec3 volDim = volumeProb->getDimensions();
    tgt::ivec3 pos;
    for (pos.z=1; pos.z<volDim.z-1; pos.z++) {
        for (pos.y=1; pos.y<volDim.y-1; pos.y++) {
            for (pos.x=1; pos.x<volDim.x-1; pos.x++) {
                if (volumeFiltered->voxel(pos) > 0) {
                    uint8_t value = volumeProb->voxel(pos);
                    uint8_t v0, v1, v2, v3, v4, v5;
                    v0 = volumeFiltered->voxel(pos + tgt::ivec3(1, 0, 0));
                    v1 = volumeFiltered->voxel(pos + tgt::ivec3(0, 1, 0));
                    v2 = volumeFiltered->voxel(pos + tgt::ivec3(0, 0, 1));

                    v3 = volumeFiltered->voxel(pos + tgt::ivec3(-1, 0, 0));
                    v4 = volumeFiltered->voxel(pos + tgt::ivec3(0, -1, 0));
                    v5 = volumeFiltered->voxel(pos + tgt::ivec3(0, 0, -1));

                    if ( (value <= 127 && (v0 > 127 || v1 > 127 || v2 > 127 || v3 > 127 || v4 > 127 || v5 > 127)) ||
                        (value > 127 && (v0 <= 127 || v1 <= 127 || v2 <= 127 || v3 <= 127 || v4 <= 127 || v5 <= 127)) ) {

                    }
                    else {
                        result->voxel(pos) = 0;
                    }
                }
                else {
                    result->voxel(pos) = 0;
                }
            }
        }
    }

    return result;
}

Volume* RandomWalkerAnalyzer::connectedComponents(Volume* volumeHandle, size_t& numLabels) {

    tgtAssert(volumeHandle, "No volume");
    const VolumeRAM* volume = volumeHandle->getRepresentation<VolumeRAM>();

    // compute connected component labels
    uint16_t* labels = 0;
    try {
        labels = new uint16_t[volume->getNumVoxels()];
    }
    catch (std::bad_alloc) {
        LERROR("Failed to create label volume: bad allocation");
        return 0;
    }

#ifdef VRN_MODULE_CONNEXE
    clock_t start = clock();

    Connexe_verbose();
    numLabels = CountConnectedComponentsWithAllParams(const_cast<void*>(reinterpret_cast<const void*>(volume->getData())), UCHAR_TYPE,
        reinterpret_cast<void*>(labels), USHORT_TYPE, tgt::ivec3(volume->getDimensions()).elem, 1.0,
        connectivity_.getValue(), minUncertaintyAreaSize_.get(), maxUncertaintyAreas_.get(), 0);

    clock_t finish = clock();

    if (numLabels > 0) {
        float sec = (static_cast<float>(finish - start)/CLOCKS_PER_SEC);
        LINFO("CC-Analysis: " << sec << " sec");

        // sort components decreasingly by size
        RelabelConnectedComponentsBySize(labels, USHORT_TYPE, tgt::ivec3(volume->getDimensions()).elem, 0);
        LDEBUG("Num labels: " << numLabels);
    }
#else
    LWARNING("Voreen compiled without 'connexe' module.");
#endif
    VolumeRAM_UInt16* v = new VolumeRAM_UInt16(labels, volume->getDimensions());

    return new Volume(v, volumeHandle);
}

void RandomWalkerAnalyzer::analyzeUncertaintyRegions(Volume* labelHandle, int numLabels) {
    const VolumeRAM* labelVolume = labelHandle->getRepresentation<VolumeRAM>();
    tgtAssert(labelVolume, "No labels");

    const uint16_t* labels = reinterpret_cast<const uint16_t*>(labelVolume->getData());
    tgt::ivec3 volDim = labelVolume->getDimensions();

    uncertaintyRegions_.clear();
    std::vector<UncertaintyRegion> regions;
    for (int i=0; i<numLabels; i++) {
        UncertaintyRegion area;
        area.label_ = i+1;
        area.size_ = 0;
        area.centerOfMass_ = tgt::vec3(0.f);
        area.llf_ = volDim;
        area.urb_ = tgt::vec3(-1);
        regions.push_back(area);
    }

    // determine region sizes, bounding boxes and center of masses
    //for (int i=0; i<labelVolume->getNumVoxels(); i++) {
    tgt::ivec3 pos;
    for (pos.z=0; pos.z<volDim.z; pos.z++) {
        for (pos.y=0; pos.y<volDim.y; pos.y++) {
            for (pos.x=0; pos.x<volDim.x; pos.x++) {
                size_t index = volumeCoordsToIndex(pos, volDim);
                int label = labels[index];
                if (label > 0) {
                    tgtAssert(label <= numLabels, "Unexpected label id");
                    UncertaintyRegion& region = regions[label-1];
                    region.size_++;
                    region.centerOfMass_ += tgt::vec3(pos);
                    region.llf_.x = std::min(region.llf_.x, pos.x);
                    region.llf_.y = std::min(region.llf_.y, pos.y);
                    region.llf_.z = std::min(region.llf_.z, pos.z);
                    region.urb_.x = std::max(region.urb_.x, pos.x);
                    region.urb_.y = std::max(region.urb_.y, pos.y);
                    region.urb_.z = std::max(region.urb_.z, pos.z);
                }
            }
        }
    }
    // divide center of mass by size
    for (size_t i=0; i<regions.size(); i++) {
        tgtAssert(regions.at(i).size_ > 0, "Empty uncertainty region");
        regions.at(i).centerOfMass_ /= static_cast<float>(regions.at(i).size_);
    }

    // sanity checks
    for (size_t i=0; i<regions.size(); i++) {

        if (tgt::hor(tgt::greaterThan(regions[i].llf_, regions[i].urb_))) {
            LERROR("Region's LLF greater than its URB (Region with label:  " << regions[i].label_ << ")");
        }

        if (tgt::hor(tgt::lessThan(tgt::iceil(regions[i].centerOfMass_), regions[i].llf_)) ||
            tgt::hor(tgt::greaterThan(tgt::ifloor(regions[i].centerOfMass_), regions[i].urb_)) ) {
                LERROR("Center of mass outside bounding box detected (Region with label:  " << regions[i].label_ << ")");
        }

        if (tgt::hmul(regions[i].urb_ - regions[i].llf_ + tgt::ivec3(1)) < regions[i].size_) {
            LERROR("Region's volume is greater than volume of its bounding box (Region with label:  "
                << regions[i].label_ << ")");
        }
    }

    // revert order of regions (largest one first)
    for (size_t i=0; i<regions.size(); i++) {
        uncertaintyRegions_.push_back(regions.at(regions.size() - i - 1));
    }

}

void RandomWalkerAnalyzer::stretchLabels(VolumeRAM_UInt16* labelVolume) {

    tgtAssert(labelVolume, "No labels");

    if (uncertaintyRegions_.empty())
        return;

    uint16_t* labels =  reinterpret_cast<uint16_t*>(labelVolume->getData());

    float scale = static_cast<float>((1<<16) - 2) / uncertaintyRegions_.size();
    for (size_t i=0; i<uncertaintyRegions_.size(); i++) {
        uncertaintyRegions_[i].label_ = tgt::ifloor(uncertaintyRegions_.at(i).label_*scale);
    }
    for (size_t i=0; i<labelVolume->getNumVoxels(); i++) {
        if (labels[i] > 0)
            labels[i] = uncertaintyRegions_.at(uncertaintyRegions_.size() - labels[i]).label_;
    }
}

void RandomWalkerAnalyzer::printLabelInformation() const {

    std::ostringstream stream;
    for (size_t i=0; i<uncertaintyRegions_.size(); i++) {
        stream << "Label " << uncertaintyRegions_.at(i).label_ << ": "
               << uncertaintyRegions_.at(i).size_ << " , ";
    }
    LINFO(stream.str());
}

void RandomWalkerAnalyzer::updateNumSlices() {
    tgtAssert(inportVolume_.hasData() && inportVolume_.getData()->getRepresentation<VolumeRAM>(), "No volume");
    tgt::ivec3 volDim = inportVolume_.getData()->getRepresentation<VolumeRAM>()->getDimensions();

    // set number of slice for xy-plane (along z-axis)
    sliceXY_.setMaxValue(volDim.z);
    if (sliceXY_.get() >= volDim.z)
        sliceXY_.set(0);
    sliceXY_.updateWidgets();

    // set number of slices for zx-plane (along y-axis)
    sliceXZ_.setMaxValue(volDim.y);
    if (sliceXZ_.get() >= volDim.y)
        sliceXZ_.set(0);
    sliceXZ_.updateWidgets();

    // set number of slices for sagittal plane (along x-axis)
    sliceYZ_.setMaxValue(volDim.x);
    if (sliceYZ_.get() >= volDim.x)
        sliceYZ_.set(0);
    sliceYZ_.updateWidgets();
}

void RandomWalkerAnalyzer::resetZoom() {
    focusXY_.set(tgt::vec2(0.f));
    focusXZ_.set(tgt::vec2(0.f));
    focusYZ_.set(tgt::vec2(0.f));
    zoomXY_.set(1.f);
    zoomXZ_.set(1.f);
    zoomYZ_.set(1.f);

    if (inportVolume_.hasData()) {
        tgt::svec3 halfDim = inportVolume_.getData()->getRepresentation<VolumeRAM>()->getDimensions() / tgt::svec3(2);
        sliceXY_.set(static_cast<int>(halfDim.x));
        sliceXZ_.set(static_cast<int>(halfDim.y));
        sliceYZ_.set(static_cast<int>(halfDim.z));
    }
}

void RandomWalkerAnalyzer::setFocusAndZoomAnimated(tgt::vec2 focusXY, tgt::vec2 focusXZ, tgt::vec2 focusYZ, float zoom, int numIterations) {

    animFocusXYStop_ = focusXY;
    animFocusXZStop_ = focusXZ;
    animFocusYZStop_ = focusYZ;

    animFocusXYStart_ = focusXY_.get();
    animFocusXZStart_ = focusXZ_.get();
    animFocusYZStart_ = focusYZ_.get();

    animZoomStart_ = zoomXY_.get();
    animZoomStop_ = zoom;
    numAnimSteps_ = numIterations;
    animStep_ = 0;

    if (timer_) {
        timer_->start(10);
    }
}

void RandomWalkerAnalyzer::performAnimStep() {
    animStep_++;
    float blendFactor = 1.f;
    if (animStep_ >= numAnimSteps_) {
        timer_->stop();
    }
    else {
        blendFactor = (float)(animStep_) / numAnimSteps_;
    }

    tgt::vec2 curFocusXY = (1.f-blendFactor) * animFocusXYStart_ +  blendFactor * animFocusXYStop_;
    tgt::vec2 curFocusXZ = (1.f-blendFactor) * animFocusXZStart_ +  blendFactor * animFocusXZStop_;
    tgt::vec2 curFocusYZ = (1.f-blendFactor) * animFocusYZStart_ +  blendFactor * animFocusYZStop_;
    float curZoom = tgt::clamp((1.f-blendFactor) * animZoomStart_ + blendFactor * animZoomStop_,
        zoomXY_.getMinValue(), zoomXY_.getMaxValue());

    focusXY_.set(curFocusXY);
    focusXZ_.set(curFocusXZ);
    focusYZ_.set(curFocusYZ);
    zoomXY_.set(curZoom);
    zoomXZ_.set(curZoom);
    zoomYZ_.set(curZoom);
}

void RandomWalkerAnalyzer::timerEvent(tgt::TimeEvent* /*e*/) {
    performAnimStep();
}

void RandomWalkerAnalyzer::forceUpdate() {
    forceUpdate_ = true;
    invalidate(INVALID_RESULT);
}

void RandomWalkerAnalyzer::resetZoomAnimated(int iterations) {
    animFocusXYStop_ = tgt::vec2(0.f);
    animFocusXZStop_ = tgt::vec2(0.f);
    animFocusYZStop_ = tgt::vec2(0.f);

    animFocusXYStart_ = focusXY_.get();
    animFocusXZStart_ = focusXZ_.get();
    animFocusYZStart_ = focusYZ_.get();

    animZoomStart_ = zoomXY_.get();
    animZoomStop_ = 1.f;
    numAnimSteps_ = iterations;
    animStep_ = 0;

    if (timer_) {
        timer_->start(10);
    }
}

void RandomWalkerAnalyzer::zoomOnFocusRegion(int regionID, int numZooSteps) {

    if (!inportVolume_.hasData())
        return;

    UncertaintyRegion& region = uncertaintyRegions_.at(regionID);

    tgt::vec3 focus = region.centerOfMass_;

    const VolumeBase* volume = inportVolume_.getData();
    tgt::ivec3 volDim = volume->getDimensions();

    // updates slices
    sliceXY_.set(0);
    sliceXY_.set(tgt::iround(focus.z));
    sliceXZ_.set(0);
    sliceXZ_.set(tgt::iround(focus.y));
    sliceYZ_.set(0);
    sliceYZ_.set(tgt::iround(focus.x));

    // update focus points
    tgt::vec3 halfDim = (tgt::vec3(volDim) - 1.f) / 2.f;
    tgt::vec3 focusX, focusY, focusZ;
    tgt::vec2 focusXY = tgt::vec2(halfDim.x - focus.x, -1.f * (halfDim.y - focus.y));
    tgt::vec2 focusXZ = tgt::vec2(halfDim.x - focus.x, 1.f * (halfDim.z - focus.z));
    tgt::vec2 focusYZ = tgt::vec2(-1.f * (halfDim.y - focus.y), 1.f * (halfDim.z - focus.z));

    // update zoom factor
    float maxZoom = 0.35f;
    tgt::vec3 dim = tgt::vec3(volDim) - 1.f;
    tgt::vec3 bbDim = region.urb_ - region.llf_;
    tgt::vec2 zoomXY = tgt::vec2(bbDim.x / dim.x, bbDim.y / dim.y) * 1.5f;
    tgt::vec2 zoomXZ = tgt::vec2(bbDim.x / dim.x, bbDim.z / dim.z) * 1.5f;
    tgt::vec2 zoomYZ = tgt::vec2(bbDim.y / dim.y, bbDim.z / dim.z) * 1.5f;
    float zoom = std::max(tgt::max(zoomXY), tgt::max(zoomXZ));
    zoom = std::max(zoom, tgt::max(zoomYZ));
    zoom = std::max(zoom, maxZoom);

    setFocusAndZoomAnimated(focusXY, focusXZ, focusYZ, zoom, numZooSteps);
}

void RandomWalkerAnalyzer::updateFocusAreaPropVolume() {
    if (!inportVolume_.hasData() || !probVolume_ || !probVolumeCC_)
        return;

    VolumeRAM_UInt8* focusPropVolume = probVolume_->clone();
    const VolumeRAM_UInt8* probVol = static_cast<const VolumeRAM_UInt8*>(probVolumeCC_->getRepresentation<VolumeRAM>());

    if (focusRegion_ >= 0 && focusRegion_ < (int)uncertaintyRegions_.size()) {
        int label = uncertaintyRegions_.at(focusRegion_).label_;
        for (size_t i=0; i<focusPropVolume->getNumVoxels(); i++) {
            if (probVol->voxel(i) != label)
                focusPropVolume->voxel(i) = 0;
        }
    }
    else {
        focusPropVolume->clear();
    }
    outportFocusAreaProbabilities_.setData(new Volume(focusPropVolume, inportVolume_.getData()));
}


}   // namespace
