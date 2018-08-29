/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany,                        *
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

#include "nucleipositionquantification.h"

#include "voreen/core/datastructures/geometry/pointsegmentlistgeometry.h"
#include "voreen/core/datastructures/geometry/pointlistgeometry.h"

#include "voreen/core/datastructures/volume/operators/volumeoperatorinvert.h"
#include "../operators/volumeoperatordistancetransform.h"

#include "voreen/core/voreenapplication.h"

#include "tgt/stopwatch.h"

namespace voreen {

const std::string NucleiPositionQuantification::loggerCat_("voreen.nucleusdetection.NucleiPositionQuantification");

NucleiPositionQuantification::NucleiPositionQuantification()
    : Processor()
    , inportImage_(Port::INPORT, "volumehandle.input.image", "Volume Input Image")
    , inportSegmentation_(Port::INPORT, "volumehandle.input.binary", "Volume Input Region of Interest")
    , centroidInport_(Port::INPORT, "centroid.input", "Nuclei Centroids Input", false)
    , insideOutport_(Port::OUTPORT, "centroid_in.output", "Nuclei Centroids Output", true, voreen::Processor::VALID)
    , outsideOutport_(Port::OUTPORT, "centroid_out.output", "Nuclei Centroids Output", true, voreen::Processor::VALID)
    //, exportAsCSV_("exportCSV", "Export CSV File", true)
    , csvFile_("csvFile", "CSV File", "Save CSV file", VoreenApplication::app()->getUserDataPath(), "Commaseparated Files(*.csv);;Textfile (*.txt);;All TextFiles (*.csv *.txt)", FileDialogProperty::SAVE_FILE,Processor::VALID)
    , progressProperty_("progressProperty", "Progress")
    , compute_("compute", "Compute")
    , forceComputation_(false)
{
    addPort(inportImage_);
    addPort(inportSegmentation_);
    addPort(centroidInport_);

    addPort(insideOutport_);
    addPort(outsideOutport_);

    addProperty(csvFile_);
    //addProperty(exportAsCSV_);

    // progress and compute button
    addProperty(progressProperty_);
    addProgressBar(&progressProperty_);

    ON_CHANGE(compute_, NucleiPositionQuantification, performQuantification);
    addProperty(compute_);
}

Processor* NucleiPositionQuantification::create() const {
    return new NucleiPositionQuantification();
}

bool NucleiPositionQuantification::isReady() const {
    if (!isInitialized())
        return false;

    return (inportImage_.isReady() && inportSegmentation_.isReady() && centroidInport_.isReady());
}

void NucleiPositionQuantification::process() {
    // reset status
    setProgress(0.f);

    if (!forceComputation_) 
        return;
    else
        forceComputation_ = false;
    
    deactivateSettings();
    setProgress(0.f);

    insideOutport_.clear();
    outsideOutport_.clear();

    // try to cast the geometry into a PointSegmentListGeometry
    const PointSegmentListGeometry<tgt::vec3>* segmentListGeom = dynamic_cast<const PointSegmentListGeometry<tgt::vec3>* >(centroidInport_.getData());

    // message on invalid geometry
    if (!segmentListGeom) {
        LWARNING("Invalid geometry. PointSegmentListGeometry<vec3> expected.");
        reactivateSettings();
        return;
    }

    // get RAM representation of the segmentation volume
    const VolumeRAM* segmentationRAM = inportSegmentation_.getData()->getRepresentation<VolumeRAM>();
    if (!segmentationRAM) {
        LERROR("Cannot create RAM representation of the segmentation volume");
        reactivateSettings();
        return;
    }

    PointListGeometryVec3* insideList = new PointListGeometryVec3();
    PointListGeometryVec3* outsideList = new PointListGeometryVec3();

    tgt::Stopwatch stopwatch;
    LINFO("Starting position quantification...");
    stopwatch.start();

    // compute ratio for coordinate rescaling
    tgt::svec3 originalDimensions = inportImage_.getData()->getDimensions();
    tgt::svec3 segmentationDimensions = inportSegmentation_.getData()->getDimensions();
    tgt::vec3 ratio = tgt::vec3(segmentationDimensions) / tgt::vec3(originalDimensions);

    // here is the actually correct rescaling formula for geting the position in the segmentation volume:
    // (pos - (dim_orig - 1) / 2) * ratio + (dim_seg - 1) / 2 
    // define d_a := (dim_orig - 1) / 2, d_b := (dim_seg - 1) / 2
    tgt::vec3 d_a = tgt::vec3(originalDimensions - tgt::svec3::one) / tgt::vec3(2.f);
    tgt::vec3 d_b = tgt::vec3(segmentationDimensions - tgt::svec3::one) / tgt::vec3(2.f);

    // TODO: if the spacing is not uniform, perform a resampling of the segmentation to force it to be uniform

    // invert foreground volume to compute distance transform to foreground
    Volume* invertedSegmentation = VolumeOperatorInvert::APPLY_OP(inportSegmentation_.getData(), 0);
    if (!invertedSegmentation) {
        LERROR("Cannot apply VolumeOperatorInvert to segmentation!");
        reactivateSettings();
        return;
    }
    setProgress(0.15f);
    Volume* invertedDistanceTransform = VolumeOperatorEuclideanDistanceTransform::APPLY_OP(invertedSegmentation, 0);
    if (!invertedDistanceTransform) {
        LERROR("Cannot compute distance transform!");
        delete invertedSegmentation;
        setProgress(0.f);
        reactivateSettings();
        return;
    }
    setProgress(0.25f);
    const VolumeRAM* dtRAM = invertedDistanceTransform->getRepresentation<VolumeRAM>();
    if (!dtRAM) {
        LERROR("Cannot create RAM representstion of inverted distance transform");
        delete invertedSegmentation;
        delete invertedDistanceTransform;
        setProgress(0.f);
        reactivateSettings();
        return;
    }

    std::string filePath = tgt::FileSystem::cleanupPath(csvFile_.get());
    std::ofstream csvExportFile(filePath.c_str());
    if (csvExportFile.fail()) {
        LERROR("Cannot create CSV file");
        delete invertedSegmentation;
        delete invertedDistanceTransform;
        setProgress(0.f);
        reactivateSettings();
        return;
    }

    csvExportFile << "centroid (original volume),centroid (segmentation),outside,distance (voxels),distance (mm)\n";

    // Iterate over segments and inside the segments over the points
    const std::vector<std::vector<tgt::vec3> >& segmentList = segmentListGeom->getData();

    size_t currentSegment = 0;
    for (auto& segment : segmentList) {  // loop over segments
        for (auto& centroid : segment) { // loop over centroids within the segment   
            // rescale the coordinates with the same scheme as the VolumeOperatorResample
            tgt::vec3 rescaledPosition = (centroid - d_a) * ratio + d_b;

            // currently: work on floats that are linearly interpolated (TODO: improve)
            tgt::svec3 nearestIndex = tgt::svec3(tgt::clamp(tgt::iround(rescaledPosition), tgt::ivec3::zero, tgt::ivec3(originalDimensions - tgt::svec3::one)));
            tgt::vec3 linearIndex = tgt::clamp(rescaledPosition, tgt::vec3::zero, tgt::vec3(segmentationDimensions - tgt::svec3::one));

            csvExportFile << centroid << "," << linearIndex << ",";

            if (segmentationRAM->getVoxelNormalized(nearestIndex) == 0.f) {
                outsideList->addPoint(centroid);

                // compute distance... assumes uniform spacing!
                float distanceVoxels = dtRAM->getVoxelNormalizedLinear(linearIndex);
                float distanceMillimeters = distanceVoxels * inportSegmentation_.getData()->getSpacing().x;

                csvExportFile << "1" << "," << distanceVoxels << "," << distanceMillimeters << "\n";
            }
            else {
                insideList->addPoint(centroid);
                csvExportFile << "0" << "," << "0" << "," << "0" << "\n";
            }
        }
        currentSegment++;

        // update progress every 10th segment
        if (currentSegment % 10 == 0) {
            setProgress(std::min(0.99f, 0.25f + 0.75f * (static_cast<float>(currentSegment) / static_cast<float>(segmentList.size()))));
        }
    }

    stopwatch.stop();
    LINFO("Position quantification finished. Runtime: " << static_cast<double>(stopwatch.getRuntime()) / 1000.0 << " seconds");

    csvExportFile.close();

    delete invertedSegmentation;
    delete invertedDistanceTransform;
    
    insideList->setTransformationMatrix(segmentListGeom->getTransformationMatrix());
    outsideList->setTransformationMatrix(segmentListGeom->getTransformationMatrix());

    insideOutport_.setData(insideList, true);
    outsideOutport_.setData(outsideList, true);

    setProgress(1.f);
    reactivateSettings();
}

void NucleiPositionQuantification::deactivateSettings() {
    //exportAsCSV_.setReadOnlyFlag(true);
    csvFile_.setReadOnlyFlag(true);
    compute_.setReadOnlyFlag(true);
}

void NucleiPositionQuantification::reactivateSettings() {
    // re-activate all settings properties 
    //exportAsCSV_.setReadOnlyFlag(false);
    csvFile_.setReadOnlyFlag(false);
    compute_.setReadOnlyFlag(false);
}

void NucleiPositionQuantification::performQuantification() {
    if (isInitialized() && isReady()) {
        forceComputation_ = true;
        invalidate();
    }
}

}   // namespace
