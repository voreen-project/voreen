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

#include "roidetector.h"
#include "voreen/core/io/progressbar.h"

#include "voreen/core/datastructures/volume/volumebase.h"
#include "voreen/core/datastructures/volume/volumeram.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"

namespace voreen {

const std::string RoiDetector::loggerCat_("voreen.nucleusdetection.RoiDetector");

RoiDetector::RoiDetector()
    : VolumeProcessor()
    , roiInport_(Port::INPORT, "roiinport", "ROI Input", false)
    , originalInport_(Port::INPORT, "originalinport", "Original Volumes Input", false)
    , computeButton_("save", "Compute")
    , startComputation_(false)
    , progressProperty_("progressProperty", "Progress")
{
    addPort(roiInport_);
    addPort(originalInport_);

    computeButton_.onChange(MemberFunctionCallback<RoiDetector>(this, &RoiDetector::computeDetection));

    addProperty(computeButton_);
    addProperty(progressProperty_);
    addProgressBar(&progressProperty_);
}

Processor* RoiDetector::create() const {
    return new RoiDetector();
}

void RoiDetector::process() {
    if (!startComputation_)
        return;

    startComputation_ = false;

    if (!roiInport_.isReady() || !roiInport_.getData() || roiInport_.getData()->empty()) {
        LERROR(" no ROI volume data ");
        return;
    }

    if (!originalInport_.isReady() || !originalInport_.getData() || originalInport_.getData()->empty()) {
        LERROR(" no ROI volume data ");
        return;
    }

    // deactivate settings
    computeButton_.setReadOnlyFlag(true);

    // loop over ROIs
    size_t numChecks = roiInport_.getData()->size() * originalInport_.getData()->size();
    for (size_t i = 0; i < roiInport_.getData()->size(); ++i) {

        LINFO("trying to detect ROI " << i);

        bool found = false;
        const VolumeRAM_UInt16* currentRoi = dynamic_cast<const VolumeRAM_UInt16*>(roiInport_.getData()->at(i)->getRepresentation<VolumeRAM>());
        if (!currentRoi) {
            LERROR("ROI volume " << i << " is not an uint16 volume - ignoring volume");
            continue;
        }
            
        // loop over original volumes
        for (size_t j = 0; j < originalInport_.getData()->size(); ++j) {

            LINFO("Searching original volume " << j );

            const VolumeRAM_UInt16* currentOriginal = dynamic_cast<const VolumeRAM_UInt16*>(originalInport_.getData()->at(j)->getRepresentation<VolumeRAM>());
            if (!currentOriginal) {
                LERROR("Original volume " << i << " is not an uint16 volume - ignoring volume");
                continue;
            }

            currentRoi = dynamic_cast<const VolumeRAM_UInt16*>(roiInport_.getData()->at(i)->getRepresentation<VolumeRAM>());

            if (!roiInport_.getData()->at(i)->hasRepresentation<VolumeRAM>() || !originalInport_.getData()->at(j)->hasRepresentation<VolumeRAM>()) {
                LERROR("VolumeRAM representation not available for both volumes at the same time!");
                continue;
            }

            // dimensions of original volume have to be greater equal
            if (!tgt::hand(tgt::greaterThanEqual(originalInport_.getData()->at(j)->getDimensions(), roiInport_.getData()->at(i)->getDimensions())))
                continue;

            tgt::svec3 lastStartVoxel = originalInport_.getData()->at(j)->getDimensions() - roiInport_.getData()->at(i)->getDimensions();
            tgt::svec3 startVoxel = tgt::svec3::zero;

            // now iterate over potential llf corners
            for (size_t startZ = 0; startZ <= lastStartVoxel.z; ++startZ) {
                if (found) break;
                LINFO("Checking startZ = " << startZ);
                for (size_t startY = 0; startY <= lastStartVoxel.y; ++startY) {
                    if (found) break;
                    for (size_t startX = 0; startX <= lastStartVoxel.x; ++startX) {
                        if (found) break;

                        // now we have a start voxel
                        startVoxel = tgt::svec3(startX, startY, startZ);

                        // now compare all voxels in the ROI volume to all of this ROI in the original volume until a voxel differs
                        bool differingVoxel = false;

                        for (size_t z = 0; z < roiInport_.getData()->at(i)->getDimensions().z; ++z) {
                            if (differingVoxel) break;
                            for (size_t y = 0; y < roiInport_.getData()->at(i)->getDimensions().y; ++y) {
                                if (differingVoxel) break;
                                for (size_t x = 0; x < roiInport_.getData()->at(i)->getDimensions().x; ++x) {
                                    if (differingVoxel) break;
                                    
                                    tgt::svec3 pos  = tgt::svec3(x,y,z);

                                    if (currentRoi->voxel(pos) != currentOriginal->voxel(startVoxel + pos))
                                        differingVoxel = true;
                                }
                            }
                        }

                        if (!differingVoxel)
                            found = true;
                    }
                }

            }
                
            if (found) {
                LINFO("Found ROI " << i << " in original volume " << j << " starting at voxel " << startVoxel);
                break;
            }

            // free memory
            originalInport_.getData()->at(j)->removeRepresentation<VolumeRAM>();
        }
        
        if (!found)
            LWARNING("Could not locate ROI volume " << i << " in list of original volumes");
        
        setProgress(std::min(static_cast<float>(i) / static_cast<float>(roiInport_.getData()->size()), 0.99f));
    }

    setProgress(1.f);

    computeButton_.setReadOnlyFlag(false);
}

void RoiDetector::computeDetection() {
    if (!isInitialized())
        return;

    startComputation_ = true;
    invalidate();
}

} // namespace
