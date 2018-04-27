/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2016 University of Muenster, Germany.                        *
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

#include "segmentationlistvalidation.h"

#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/datastructures/volume/operators/volumeoperatornumsignificant.h"
#include "voreen/core/ports/conditions/portconditionvolumetype.h"
#include "voreen/core/datastructures/callback/memberfunctioncallback.h"

#include "custommodules/bigdataimageprocessing/util/csvwriter.h"

#include <climits>

namespace voreen {

const std::string SegmentationListValidation::loggerCat_("voreen.base.SegmentationListValidation");

SegmentationListValidation::SegmentationListValidation()
    : Processor()
    , inportSegmentations_(Port::INPORT, "segmentationlist.in", "Segmentation Volume")
    , inportReference_(Port::INPORT, "segmentation.reference", "Reference Volume")
    , exportFilePath_("exportFilePath", "Export File Path", "Export File Path", "", "*.csv", FileDialogProperty::SAVE_FILE)
    , computeButton_("computeButton", "Compute")
    , forceUpdate_(false)
{
    inportReference_.addCondition(new PortConditionVolumeChannelCount(1));
    addPort(inportSegmentations_);
    addPort(inportReference_);

    addProperty(exportFilePath_);

    computeButton_.onClick(MemberFunctionCallback<SegmentationListValidation>(this, &SegmentationListValidation::forceUpdate));
    addProperty(computeButton_);

}

Processor* SegmentationListValidation::create() const {
    return new SegmentationListValidation();
}

void SegmentationListValidation::process() {
    if (!forceUpdate_)
        return;

    forceUpdate_ = false;

    if (!inportSegmentations_.getData() || !inportReference_.getData()) {
        resetOutput();
    }
    else {
        computeMetrics();
    }
}

void SegmentationListValidation::forceUpdate() {
    forceUpdate_ = true;
    process();
}

void SegmentationListValidation::computeMetrics() {
    const VolumeList* segmentationsHandle = inportSegmentations_.getData();
    const VolumeBase* refHandle = inportReference_.getData();
    tgtAssert(segmentationsHandle && refHandle, "No input data");

    try {
        CSVWriter<size_t, size_t, size_t, size_t, size_t, size_t, size_t, size_t, float, float, float, float, float> writer(exportFilePath_.get());
        writer.writeHeader("id", "numVoxels", "numSignificantSegmentation", "numSignificantRef", "truePositive", "falsePositive", "falseNegative", "trueNegative", "jaccardIndex", "diceIndex", "sensitivity", "specificity", "threshold");

        for(size_t i = 0; i < segmentationsHandle->size(); ++i) {
            const VolumeBase* segHandle = segmentationsHandle->at(i);
            tgtAssert(segHandle, "No input data");

            const VolumeRAM* segVolume = segHandle->getRepresentation<VolumeRAM>();
            const VolumeRAM* refVolume = refHandle->getRepresentation<VolumeRAM>();
            tgtAssert(segVolume && refVolume, "Missing RAM representations");

            // check input data
            if (segVolume->getDimensions() != refVolume->getDimensions()) {
                LERROR("Volumes must have the same dimensions");
                resetOutput();
                return;
            }

            // Hack: The threshold used is stored as offset in the RWM.
            float threshold = segHandle->getRealWorldMapping().getOffset();

            size_t numVoxels = refVolume->getNumVoxels();

            // compute measures
            size_t sizeSeg = VolumeOperatorNumSignificant::APPLY_OP(segHandle);
            size_t sizeRef = VolumeOperatorNumSignificant::APPLY_OP(refHandle);

            size_t truePositive = getNumCommonVoxels(segVolume, refVolume);
            size_t falsePositive = sizeSeg - truePositive;
            size_t falseNegative = sizeRef - truePositive;
            size_t trueNegative = numVoxels - truePositive - falsePositive - falseNegative;

            tgtAssert(truePositive <= sizeSeg && truePositive <= sizeRef, "Invalid result");
            tgtAssert(numVoxels == (truePositive + falsePositive + falseNegative + trueNegative), "Invalid result");
            tgtAssert(sizeSeg == (truePositive + falsePositive), "Invalid result");
            tgtAssert(sizeRef == (truePositive + falseNegative), "Invalid result");

            // compute derived indices
            float jaccard = truePositive / static_cast<float>(falsePositive + truePositive + falseNegative);
            float dice = 2.f*truePositive / static_cast<float>((falsePositive + truePositive) + (truePositive + falseNegative));
            float sensitivity = truePositive / static_cast<float>((truePositive + falseNegative));
            float specificity = trueNegative / static_cast<float>((trueNegative + falsePositive));


            writer.write(
                    i,
                    numVoxels,
                    sizeSeg,
                    sizeRef,
                    truePositive,
                    falsePositive,
                    falseNegative,
                    trueNegative,
                    jaccard,
                    dice,
                    sensitivity,
                    specificity,
                    threshold
                    );
        }
        LINFO("Writing Segmentation Stats: " << exportFilePath_.get());
    } catch(tgt::IOException& e) {
        LERROR("Error opening file " << exportFilePath_.get() << ": " << e.what());
    }
}

void SegmentationListValidation::resetOutput() {
}

size_t SegmentationListValidation::getNumCommonVoxels(const VolumeRAM* volumeA, const VolumeRAM* volumeB) {
    tgtAssert(volumeA && volumeB, "Null pointer passed");
    tgtAssert(volumeA->getDimensions() == volumeB->getDimensions(), "Volume dimensions differ");

    size_t result = 0;
    VRN_FOR_EACH_VOXEL(pos, tgt::ivec3(0), volumeA->getDimensions()) {
        if (volumeA->getVoxelNormalized(pos) > 0.f && volumeB->getVoxelNormalized(pos) > 0.f)
            result++;
    }
    return result;
}

} // namespace
