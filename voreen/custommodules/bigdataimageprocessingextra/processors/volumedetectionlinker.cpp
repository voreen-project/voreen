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

#include "volumedetectionlinker.h"

namespace voreen {

const std::string VolumeDetectionLinker::loggerCat_("voreen.VolumeDetectionLinker");

VolumeDetectionLinker::VolumeDetectionLinker()
    : Processor()
    , volumePort_(Port::INPORT, "firstsegmentation", "First Segmentation Volume", true)
    , volumesDetected_("volumesDetected_", "Volumes Detected", false)
    , numVolumesRequiredForDetection_("numVolumesRequiredForDetection", "Number of Volumes required for detection", 1, 0, 1000)
    , numDetectionsRequired_("numDetectionsRequired", "Number of repeated volume detections required to trigger", 1, 1, 1000)
    , numDetections_("numDetections", "Number of detections", 0, 0, 1000, Processor::VALID)
    , reset_("reset", "Reset")
{
    addPort(volumePort_);
    ON_CHANGE(volumePort_, VolumeDetectionLinker, checkCondition);

    numDetections_.setReadOnlyFlag(true);

    addProperty(numVolumesRequiredForDetection_);
    ON_CHANGE_LAMBDA(numDetectionsRequired_, [&] () {
        numDetections_.set(0);
        volumesDetected_.set(false);
    });

    addProperty(numDetectionsRequired_);
    ON_CHANGE_LAMBDA(numDetectionsRequired_, [&] () {
        numDetections_.set(0);
        volumesDetected_.set(false);
    });
    addProperty(numDetections_);

    addProperty(volumesDetected_);
    volumesDetected_.setReadOnlyFlag(true);


    addProperty(reset_);
    ON_CHANGE_LAMBDA(reset_, [&] () {
        numDetections_.set(0);
        volumesDetected_.set(false);
    });
}

Processor* VolumeDetectionLinker::create() const {
    return new VolumeDetectionLinker();
}

void VolumeDetectionLinker::process() {
    // Nada
}

void VolumeDetectionLinker::checkCondition() {
    bool detection = volumePort_.getAllData().size() >= numVolumesRequiredForDetection_.get();
    if(detection) {
        numDetections_.set(numDetections_.get() + 1);
        if(numDetections_.get() >= numDetectionsRequired_.get()) {
            volumesDetected_.set(true);
        }
    }
}

} // namespace
