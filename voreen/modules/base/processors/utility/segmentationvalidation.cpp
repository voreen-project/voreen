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

#include "segmentationvalidation.h"

#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/datastructures/volume/operators/volumeoperatornumsignificant.h"
#include "voreen/core/ports/conditions/portconditionvolumetype.h"
#include "voreen/core/datastructures/callback/memberfunctioncallback.h"

#include "voreen/core/voreenapplication.h"

#include <climits>

namespace voreen {

const std::string SegmentationValidation::loggerCat_("voreen.base.SegmentationValidation");

SegmentationValidation::SegmentationValidation()
    : Processor()
    , inportSegmentation_(Port::INPORT, "segmentation.in", "Segmentation Volume")
    , inportReference_(Port::INPORT, "segmentation.reference", "Reference Volume")
    , computeButton_("computeButton", "Compute")
    , sizeSegmentation_("sizeSegmentation", "Size Segmentation |A|", 0, 0, INT_MAX)
    , sizeReference_("sizeReference", "Size Reference |B|", 0, 0, INT_MAX)
    , truePositive_("truePositive", "True Positive (TP)", 0, 0, INT_MAX)
    , falsePositive_("falsePositive", "False Positive (FP)", 0, 0, INT_MAX)
    , falseNegative_("falseNegative", "False Negative (FN)", 0, 0, INT_MAX)
    , trueNegative_("trueNegative", "True Negative (TN)", 0, 0, INT_MAX)
    , accuracy_("accuracy", "Accuracy")
    , accuracyForeground_("accuracyForeground", "Accuracy (Foreground)")
    , accuracyBackground_("accuracyBackground", "Accuracy (Background)")
    , jaccardIndex_("jaccardIndex", "Jaccard Index")
    , diceIndex_("diceIndex", "Dice Index")
    , sensitivity_("sensitivity", "Sensitivity")
    , specificity_("specificity", "Specificity")
    , exportFile_("exportfile", "Data Export (CSV)", "Save quantification results...", VoreenApplication::app()->getUserDataPath(), "Commaseparated Files(*.csv)", FileDialogProperty::SAVE_FILE, Processor::VALID)
    , exportButton_("exportbutton", "Export Data")
    , autoCompute_("autoCompute", "Auto Compute", false)
    , autoExport_("autoExport", "Auto Export CSV", false, Processor::INVALID_RESULT, Property::LOD_DEBUG)
    , forceUpdate_(false)
{
    inportSegmentation_.addCondition(new PortConditionVolumeChannelCount(1));
    inportReference_.addCondition(new PortConditionVolumeChannelCount(1));
    addPort(inportSegmentation_);
    addPort(inportReference_);

    accuracy_.setNumDecimals(4);
    accuracyForeground_.setNumDecimals(4);
    accuracyBackground_.setNumDecimals(4);
    jaccardIndex_.setNumDecimals(4);
    diceIndex_.setNumDecimals(4);
    sensitivity_.setNumDecimals(4);
    specificity_.setNumDecimals(4);

    computeButton_.onClick(MemberFunctionCallback<SegmentationValidation>(this, &SegmentationValidation::forceUpdate));
    sizeSegmentation_.setReadOnlyFlag(true);
    sizeReference_.setReadOnlyFlag(true);
    truePositive_.setReadOnlyFlag(true);
    falsePositive_.setReadOnlyFlag(true);
    falseNegative_.setReadOnlyFlag(true);
    trueNegative_.setReadOnlyFlag(true);
    accuracy_.setReadOnlyFlag(true);
    accuracyForeground_.setReadOnlyFlag(true);
    accuracyBackground_.setReadOnlyFlag(true);
    jaccardIndex_.setReadOnlyFlag(true);
    diceIndex_.setReadOnlyFlag(true);
    sensitivity_.setReadOnlyFlag(true);
    specificity_.setReadOnlyFlag(true);

    addProperty(computeButton_);
    addProperty(sizeSegmentation_);
    addProperty(sizeReference_);
    addProperty(truePositive_);
    addProperty(falsePositive_);
    addProperty(falseNegative_);
    addProperty(trueNegative_);
    addProperty(accuracy_);
    addProperty(accuracyForeground_);
    addProperty(accuracyBackground_);
    addProperty(jaccardIndex_);
    addProperty(diceIndex_);
    addProperty(sensitivity_);
    addProperty(specificity_);

    exportFile_.setGroupID("dataexport");
    exportButton_.setGroupID("dataexport");
    addProperty(exportFile_);
    addProperty(exportButton_);
    setPropertyGroupGuiName("dataexport", "Data Export");

    exportButton_.onClick(MemberFunctionCallback<SegmentationValidation>(this, &SegmentationValidation::exportData));

    addProperty(autoCompute_);
    addProperty(autoExport_);
}

Processor* SegmentationValidation::create() const {
    return new SegmentationValidation();
}

void SegmentationValidation::beforeProcess() {
    Processor::beforeProcess();
    if (autoCompute_.get())
        forceUpdate_ = true;
}

void SegmentationValidation::process() {
    if (!forceUpdate_)
        return;

    forceUpdate_ = false;

    if (!inportSegmentation_.getData() || !inportReference_.getData()) {
        resetOutput();
    }
    else {
        computeMetrics();
    }

    if (autoExport_.get())
        exportData();
}

void SegmentationValidation::forceUpdate() {
    forceUpdate_ = true;
    //process();
    invalidate();
}

void SegmentationValidation::computeMetrics() {
    const VolumeBase* segHandle = inportSegmentation_.getData();
    const VolumeBase* refHandle = inportReference_.getData();
    tgtAssert(segHandle && refHandle, "No input data");

    const VolumeRAM* segVolume = segHandle->getRepresentation<VolumeRAM>();
    const VolumeRAM* refVolume = refHandle->getRepresentation<VolumeRAM>();
    tgtAssert(segVolume && refVolume, "Missing volumes");

    // check input data
    if (segVolume->getDimensions() != refVolume->getDimensions()) {
        LWARNING("Volumes must have the same dimensions");
        resetOutput();
        return;
    }

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

    // compute accuracy
    float accuracy = static_cast<float>(truePositive + trueNegative) / static_cast<float>(numVoxels);
    float accuracyForeground = static_cast<float>(truePositive) / static_cast<float>(truePositive + falseNegative);
    float accuracyBackground = static_cast<float>(trueNegative) / static_cast<float>(trueNegative + falsePositive);

    // compute derived indices
    float jaccard = truePositive / static_cast<float>(falsePositive + truePositive + falseNegative);
    float dice = 2.f*truePositive / static_cast<float>((falsePositive + truePositive) + (truePositive + falseNegative));
    float sensitivity = truePositive / static_cast<float>((truePositive + falseNegative));
    float specificity = trueNegative / static_cast<float>((trueNegative + falsePositive));

    // update ranges
    sizeSegmentation_.setMaxValue(static_cast<int>(numVoxels));
    sizeReference_.setMaxValue(static_cast<int>(numVoxels));
    truePositive_.setMaxValue(static_cast<int>(numVoxels));
    falsePositive_.setMaxValue(static_cast<int>(numVoxels));
    falseNegative_.setMaxValue(static_cast<int>(numVoxels));
    trueNegative_.setMaxValue(static_cast<int>(numVoxels));

    // put out measures
    sizeSegmentation_.set(static_cast<int>(sizeSeg));
    sizeReference_.set(static_cast<int>(sizeRef));
    truePositive_.set(static_cast<int>(truePositive));
    falsePositive_.set(static_cast<int>(falsePositive));
    falseNegative_.set(static_cast<int>(falseNegative));
    trueNegative_.set(static_cast<int>(trueNegative));
    accuracy_.set(accuracy);
    accuracyForeground_.set(accuracyForeground);
    accuracyBackground_.set(accuracyBackground);
    jaccardIndex_.set(jaccard);
    diceIndex_.set(dice);
    sensitivity_.set(sensitivity);
    specificity_.set(specificity);
}

void SegmentationValidation::resetOutput() {
    sizeSegmentation_.set(0);
    sizeReference_.set(0);
    truePositive_.set(0);
    falsePositive_.set(0);
    falseNegative_.set(0);
    trueNegative_.set(0);
    accuracy_.set(0.f);
    accuracyForeground_.set(0.f);
    accuracyBackground_.set(0.f);
    jaccardIndex_.set(0.f);
    diceIndex_.set(0.f);
    sensitivity_.set(0.f);
    specificity_.set(0.f);
}

size_t SegmentationValidation::getNumCommonVoxels(const VolumeRAM* volumeA, const VolumeRAM* volumeB) {
    tgtAssert(volumeA && volumeB, "Null pointer passed");
    tgtAssert(volumeA->getDimensions() == volumeB->getDimensions(), "Volume dimensions differ");

    size_t result = 0;
    VRN_FOR_EACH_VOXEL(pos, tgt::svec3::zero, volumeA->getDimensions()) {
        if (volumeA->getVoxelNormalized(pos) > 0.f && volumeB->getVoxelNormalized(pos) > 0.f)
            result++;
    }
    return result;
}

void SegmentationValidation::exportData() {
    std::string exportPath = tgt::FileSystem::cleanupPath(exportFile_.get());
    std::ofstream exportFile(exportPath);
    if (!exportFile.good()) {
        LERROR("Cannot write to " << exportPath);
        exportFile.close();
        return;
    }

    // write header
    //exportFile << "size segmentation,size reference,true positives,false positives,false negatives,true negatives,jaccard index,dice index,sensitivity,specificity" << std::endl;
    exportFile <<   sizeSegmentation_.getGuiName()   << ",";
    exportFile <<   sizeReference_.getGuiName()      << ",";
    exportFile <<   truePositive_.getGuiName()       << ",";
    exportFile <<   falsePositive_.getGuiName()      << ",";
    exportFile <<   falseNegative_.getGuiName()      << ",";
    exportFile <<   trueNegative_.getGuiName()       << ",";
    exportFile <<   accuracy_.getGuiName()           << ",";
    exportFile <<   accuracyForeground_.getGuiName() << ",";
    exportFile <<   accuracyBackground_.getGuiName() << ",";
    exportFile <<   jaccardIndex_.getGuiName()       << ",";
    exportFile <<   diceIndex_.getGuiName()          << ",";
    exportFile <<   sensitivity_.getGuiName()        << ",";
    exportFile <<   specificity_.getGuiName()        << std::endl;

    // write values
    exportFile <<   sizeSegmentation_.get()   << ",";
    exportFile <<   sizeReference_.get()      << ",";
    exportFile <<   truePositive_.get()       << ",";
    exportFile <<   falsePositive_.get()      << ",";
    exportFile <<   falseNegative_.get()      << ",";
    exportFile <<   trueNegative_.get()       << ",";
    exportFile <<   accuracy_.get()           << ",";
    exportFile <<   accuracyForeground_.get() << ",";
    exportFile <<   accuracyBackground_.get() << ",";
    exportFile <<   jaccardIndex_.get()       << ",";
    exportFile <<   diceIndex_.get()          << ",";
    exportFile <<   sensitivity_.get()        << ",";
    exportFile <<   specificity_.get()        << std::endl;

    exportFile.close();

    LINFO("Finished data export");
}


} // namespace
