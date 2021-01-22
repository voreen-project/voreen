/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2020 University of Muenster, Germany,                        *
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

#include "volumecombine.h"

#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/datastructures/volume/operators/volumeoperatorresize.h"
#include "voreen/core/datastructures/geometry/meshlistgeometry.h"
#include "voreen/core/datastructures/callback/memberfunctioncallback.h"

using tgt::ivec3;
using tgt::svec3;
using tgt::vec3;

namespace voreen {

const std::string VolumeCombine::loggerCat_("voreen.VolumeCombine");

VolumeCombine::VolumeCombine()
    : CachingVolumeProcessor()
    , inportFirst_(Port::INPORT, "volume.first", "Volume1 Input")
    , inportSecond_(Port::INPORT, "volume.second", "Volume2 Input")
    , outport_(Port::OUTPORT, "outport", "Volume Output", true)
    , enableProcessing_("enabled", "Enable", true)
    , combineFunction_("combineFunction", "Combine Function")
    , factorC_("factorC", "Factor c", 0.5f, -2.f, 2.f)
    , factorD_("factorD", "Factor d", 0.5f, -2.f, 2.f)
    , filteringMode_("filteringMode", "Filtering")
    , referenceVolume_("referenceVolume", "Reference Volume")
    , pipeThroughIfSecondNotReady_("pipethrough", "Propagate First Volume if Second Not Ready", false)
{
    addPort(inportFirst_);
    addPort(inportSecond_);
    addPort(outport_);

    combineFunction_.addOption("max",              "Max (max(A,B))",             OP_MAX);
    combineFunction_.addOption("min",              "Min (min(A,B))",             OP_MIN);
    combineFunction_.addOption("add",              "Add (A+B)",                  OP_ADD);
    combineFunction_.addOption("a-minus-b",        "Subtract (A-B)",             OP_A_MINUS_B);
    combineFunction_.addOption("b-minus-a",        "Subtract (B-A)",             OP_B_MINUS_A);
    combineFunction_.addOption("mult",             "Multiply (A*B)",             OP_MULT);
    combineFunction_.addOption("avg",              "Avg ((A+B)/2)",              OP_AVG);
    combineFunction_.addOption("weightedSum",      "Weighted Sum (c*A+(1-c)*B)", OP_WEIGHTED_SUM);
    combineFunction_.addOption("weightedSum2p",    "Weighted Sum (c*A+d*B)",     OP_WEIGHTED_SUM_2P);
    combineFunction_.addOption("blend",            "Blend (A+B*(1-A))",          OP_BLEND);
    combineFunction_.addOption("mask-a-by-b",      "Mask A by B (B?A:0)",        OP_MASK_A_BY_B);
    combineFunction_.addOption("mask-b-by-a",      "Mask B by A (A?B:0)",        OP_MASK_B_BY_A);
    combineFunction_.addOption("priorityFirst",    "Priority First (A?A:B)",     OP_PRIORITY_FIRST);
    combineFunction_.addOption("prioritySecond",   "Priority Second (B?B:A)",    OP_PRIORITY_SECOND);
    combineFunction_.addOption("takeFirst",        "Take First (A)",             OP_TAKE_FIRST);
    combineFunction_.addOption("takeSecond",       "Take Second (B)",            OP_TAKE_SECOND);
    combineFunction_.select("max");

    filteringMode_.addOption("nearest", "Nearest");
    filteringMode_.addOption("linear",  "Linear ");
    filteringMode_.addOption("cubic",   "Cubic");
    filteringMode_.set("linear");

    referenceVolume_.addOption("first", "First");
    referenceVolume_.addOption("second", "Second");

    combineFunction_.onChange(MemberFunctionCallback<VolumeCombine>(this, &VolumeCombine::adjustPropertyVisibilities));

    factorC_.setTracking(false);
    factorD_.setTracking(false);

    addProperty(enableProcessing_);
    addProperty(combineFunction_);
    addProperty(factorC_);
    addProperty(factorD_);
    addProperty(filteringMode_);
    addProperty(referenceVolume_);
    addProperty(pipeThroughIfSecondNotReady_);

    adjustPropertyVisibilities();
}

VolumeCombine::~VolumeCombine() {}

Processor* VolumeCombine::create() const {
    return new VolumeCombine();
}

bool VolumeCombine::isReady() const {
    if(pipeThroughIfSecondNotReady_.get()) {
        return inportFirst_.isReady();
    } else if(combineFunction_.getValue() == OP_TAKE_FIRST) {
        return inportFirst_.isReady();
    } else if(combineFunction_.getValue() == OP_TAKE_SECOND) {
        return inportSecond_.isReady();
    } else {
        return inportFirst_.isReady() && inportSecond_.isReady();
    }
}

void VolumeCombine::process() {

    const VolumeBase* firstVolume = inportFirst_.getData();
    const VolumeBase* secondVolume = inportSecond_.getData();
    Volume* combinedVolume = 0;

    if (!enableProcessing_.get() || pipeThroughIfSecondNotReady_.get() && !inportSecond_.isReady()) {
        outport_.setData(inportFirst_.getData(), false);
        return;
    }

    tgtAssert(inportFirst_.getData() && inportFirst_.getData()->getRepresentation<VolumeRAM>(), "No input volume");
    tgtAssert(inportSecond_.getData() && inportSecond_.getData()->getRepresentation<VolumeRAM>(), "No input volume");

    if (firstVolume->getNumChannels() == secondVolume->getNumChannels()) {
        LINFO("Performing channel-wise combination.");
    }
    else {
        if (referenceVolume_.isSelected("first"))
            LINFO("Number of channels of input volumes do not match. Combining each channel of first volume with channel 0 of second volume");
        else
            LINFO("Number of channels of input volumes do not match. Combining each channel of second volume with channel 0 of first volume");
    }


    // optimized combination for volumes that share a common grid in world-space
    if (firstVolume->getDimensions() == secondVolume->getDimensions() &&
        firstVolume->getSpacing() == secondVolume->getSpacing()     &&
        firstVolume->getOffset() == secondVolume->getOffset()     &&
        firstVolume->getPhysicalToWorldMatrix() == secondVolume->getPhysicalToWorldMatrix()) {

        try {
            if (referenceVolume_.isSelected("first"))
                combinedVolume = firstVolume->clone();
            else
                combinedVolume = secondVolume->clone();
        }
        catch (const std::bad_alloc&) {
            LERROR("Failed to create combined volume with dimensions " << firstVolume->getDimensions()
                << " : bad allocation");
        }

        if (combinedVolume) {
            LINFO("Performing optimized combination on common grid with dimensions "
                << firstVolume->getDimensions() << "...");
            combineVolumesOnCommonGrid(combinedVolume, firstVolume, secondVolume, combineFunction_.getValue());
        }
    }
    // standard combination with resampling
    else {
        if (referenceVolume_.isSelected("first"))
            combinedVolume = createCombinedVolume(firstVolume, secondVolume);
        else
            combinedVolume = createCombinedVolume(secondVolume, firstVolume);

        if (combinedVolume) {
            LINFO("Creating combined volume with dimensions " << combinedVolume->getDimensions() << " using "
                << filteringMode_.get() << " filtering ...");
            combineVolumes(combinedVolume, firstVolume, secondVolume, combineFunction_.getValue());
        }
    }

    outport_.setData(combinedVolume);
}

inline bool withinRange(const tgt::vec3& pos, const tgt::vec3& llf, const tgt::vec3& urb) {
    return tgt::hand(tgt::greaterThanEqual(pos, llf)) &&
           tgt::hand(tgt::lessThanEqual(   pos, urb));
}

void VolumeCombine::combineVolumes(Volume* combinedVolume, const VolumeBase* firstVolume,
                                   const VolumeBase* secondVolume, CombineOperation operation) const {

    tgtAssert(combinedVolume && firstVolume && secondVolume, "Null pointer passed");

    // compute transformation from voxel coordinates of combined volume
    // to voxel coords of input volumes
    tgt::mat4 combinedToFirst = computeConversionMatrix(combinedVolume, firstVolume);
    tgt::mat4 combinedToSecond = computeConversionMatrix(combinedVolume, secondVolume);

    LDEBUG("Voxel-to-world (First): " << combinedToFirst);
    LDEBUG("Voxel-to-world (Second) " << combinedToSecond);

    // determine which volume is the other (non reference) one
    const VolumeBase* otherVolume = referenceVolume_.isSelected("first") ? secondVolume : firstVolume;

    //
    // voxel-wise combination
    //
    tgt::vec3 dimFirst(firstVolume->getDimensions() - svec3(1));
    tgt::vec3 dimSecond(secondVolume->getDimensions() - svec3(1));
    const float c = factorC_.get();
    const float d = factorD_.get();
    const bool nearestFiltering = filteringMode_.isSelected("nearest");
    const bool linearFiltering = filteringMode_.isSelected("linear");
    const bool cubicFiltering = filteringMode_.isSelected("cubic");

    const VolumeRAM* v1 = firstVolume->getRepresentation<VolumeRAM>();
    const VolumeRAM* v2 = secondVolume->getRepresentation<VolumeRAM>();
    VolumeRAM* vc = combinedVolume->getWritableRepresentation<VolumeRAM>();

    const tgt::svec3 combinedVolDim = combinedVolume->getDimensions();
    const size_t combinedVolChan = combinedVolume->getNumChannels();
    const size_t otherVolumeChan = otherVolume->getNumChannels();

    VRN_FOR_EACH_VOXEL_WITH_PROGRESS(pos, tgt::svec3::zero, combinedVolDim, const_cast<VolumeCombine*>(this)) {
        for (size_t referenceChannel = 0; referenceChannel < combinedVolChan; ++referenceChannel) {
            size_t otherChannel = (combinedVolChan == otherVolumeChan) ? referenceChannel : 0;

            // transform sampling pos to coordinate systems of input volumes
            tgt::vec3 posFirst = combinedToFirst*tgt::vec3(pos);
            tgt::vec3 posSecond = combinedToSecond*tgt::vec3(pos);

            // sample input volumes, if transformed voxel position lies inside respective volume
            float valFirst = 0.f;
            float valSecond = 0.f;
            if (nearestFiltering) {
                if (withinRange(posFirst, tgt::vec3::zero, dimFirst))
                    valFirst = v1->getVoxelNormalized(tgt::iround(posFirst), referenceChannel);
                if (withinRange(posSecond, tgt::vec3::zero, dimSecond))
                    valSecond = v2->getVoxelNormalized(tgt::iround(posSecond), otherChannel);
            }
            else if (linearFiltering) {
                if (withinRange(posFirst, tgt::vec3::zero, dimFirst))
                    valFirst = v1->getVoxelNormalizedLinear(posFirst, referenceChannel);
                if (withinRange(posSecond, tgt::vec3::zero, dimSecond))
                    valSecond = v2->getVoxelNormalizedLinear(posSecond, otherChannel);
            }
            else if (cubicFiltering) {
                if (withinRange(posFirst, tgt::vec3::zero, dimFirst))
                    valFirst = v1->getVoxelNormalizedCubic(posFirst, referenceChannel);
                if (withinRange(posSecond, tgt::vec3::zero, dimSecond))
                    valSecond = v2->getVoxelNormalizedCubic(posSecond, otherChannel);
            }
            else {
                LERROR("Unknown filter mode: " << filteringMode_.get());
                return;
            }

            // apply operation to sampled values (note: a switch-block within a volume traversal loop
            // should normally be avoided, however in this case the main operations are way more expensive)
            float result = 0.f;
            switch (operation) {
                case OP_MAX:
                    result = std::max(valFirst, valSecond);
                    break;
                case OP_MIN:
                    result = std::min(valFirst, valSecond);
                    break;
                case OP_ADD:
                    result = valFirst + valSecond;
                    break;
                case OP_A_MINUS_B:
                    result = valFirst - valSecond;
                    break;
                case OP_B_MINUS_A:
                    result = valSecond - valFirst;
                    break;
                case OP_MULT:
                    result = valFirst * valSecond;
                    break;
                case OP_AVG:
                    result = (valFirst + valSecond) / 2.f;
                    break;
                case OP_WEIGHTED_SUM:
                    result = c*valFirst + (1.f-c)*valSecond;
                    break;
                case OP_WEIGHTED_SUM_2P:
                    result = c*valFirst + d*valSecond;
                    break;
                case OP_BLEND:
                    result = valFirst + valSecond*(1.f - valFirst);
                    break;
                case OP_MASK_A_BY_B:
                    result = (valSecond > 0.f) ? valFirst : 0.f;
                    break;
                case OP_MASK_B_BY_A:
                    result = (valFirst > 0.f) ? valSecond : 0.f;
                    break;
                case OP_PRIORITY_FIRST:
                    result = (valFirst > 0.f) ? valFirst : valSecond;
                    break;
                case OP_PRIORITY_SECOND:
                    result = (valSecond > 0.f) ? valSecond : valFirst;
                    break;
                case OP_TAKE_FIRST:
                    result = valFirst;
                    break;
                case OP_TAKE_SECOND:
                    result = valSecond;
                    break;
                default:
                    LERROR("Unknown operation: " << combineFunction_.get());
                    return;
            }

            // assign clamped result to combined volume
            //vc->setVoxelNormalized(tgt::clamp(result, 0.f, 1.f), pos, referenceChannel); //FIXME: clamp prevents working with float volumes
            vc->setVoxelNormalized(result, pos, referenceChannel);
        }
    } // VRN_FOR_EACH_VOXEL_WITH_PROGRESS
}

void VolumeCombine::combineVolumesOnCommonGrid(Volume* combinedVolume, const VolumeBase* firstVolume,
                                               const VolumeBase* secondVolume, CombineOperation operation) const {
    tgtAssert(combinedVolume && firstVolume && secondVolume, "Null pointer passed");
    tgtAssert(combinedVolume->getDimensions() == firstVolume->getDimensions() &&
              combinedVolume->getDimensions() == secondVolume->getDimensions(), "Volume dimensions mismatch");

    // determine which volume is the other (non reference) one
    const VolumeBase* otherVolume = referenceVolume_.isSelected("first") ? secondVolume : firstVolume;

    const float c = factorC_.get();
    const float d = factorD_.get();
    const VolumeRAM* v1 = firstVolume->getRepresentation<VolumeRAM>();
    const VolumeRAM* v2 = secondVolume->getRepresentation<VolumeRAM>();
    VolumeRAM* vc = combinedVolume->getWritableRepresentation<VolumeRAM>();

    const tgt::svec3 combinedVolDim = combinedVolume->getDimensions();
    const size_t combinedVolChan = combinedVolume->getNumChannels();
    const size_t otherVolumeChan = otherVolume->getNumChannels();

    VRN_FOR_EACH_VOXEL_WITH_PROGRESS(pos, tgt::svec3::zero, combinedVolDim, const_cast<VolumeCombine*>(this)) {
        for (size_t referenceChannel = 0; referenceChannel < combinedVolChan; ++referenceChannel) {
            size_t otherChannel = (combinedVolChan == otherVolumeChan) ? referenceChannel : 0;
            float valFirst = v1->getVoxelNormalized(pos, referenceChannel);
            float valSecond = v2->getVoxelNormalized(pos, otherChannel);

            // apply operation to input voxel values
            float result = 0.f;
            switch (operation) {
                case OP_MAX:
                    result = std::max(valFirst, valSecond);
                    break;
                case OP_MIN:
                    result = std::min(valFirst, valSecond);
                    break;
                case OP_ADD:
                    result = valFirst + valSecond;
                    break;
                case OP_A_MINUS_B:
                    result = valFirst - valSecond;
                    break;
                case OP_B_MINUS_A:
                    result = valSecond - valFirst;
                    break;
                case OP_MULT:
                    result = valFirst * valSecond;
                    break;
                case OP_AVG:
                    result = (valFirst + valSecond) / 2.f;
                    break;
                case OP_WEIGHTED_SUM:
                    result = c*valFirst + (1.f-c)*valSecond;
                    break;
                case OP_WEIGHTED_SUM_2P:
                    result = c*valFirst + d*valSecond;
                    break;
                case OP_BLEND:
                    result = valFirst + valSecond*(1.f - valFirst);
                    break;
                case OP_MASK_A_BY_B:
                    result = (valSecond > 0.f) ? valFirst : 0.f;
                    break;
                case OP_MASK_B_BY_A:
                    result = (valFirst > 0.f) ? valSecond : 0.f;
                    break;
                case OP_PRIORITY_FIRST:
                    result = (valFirst > 0.f) ? valFirst : valSecond;
                    break;
                case OP_PRIORITY_SECOND:
                    result = (valSecond > 0.f) ? valSecond : valFirst;
                    break;
                case OP_TAKE_FIRST:
                    result = valFirst;
                    break;
                case OP_TAKE_SECOND:
                    result = valSecond;
                    break;
                default:
                    LERROR("Unknown operation: " << combineFunction_.get());
                    return;
            }

            // assign clamped result to combined volume
            //vc->setVoxelNormalized(tgt::clamp(result, 0.f, 1.f), pos, referenceChannel); //FIXME: clamp prevents working with float volumes
            vc->setVoxelNormalized(result, pos, referenceChannel);
        } // for referenceChannel
    } // VRN_FOR_EACH_VOXEL_WITH_PROGRESS
}

Volume* VolumeCombine::createCombinedVolume(const VolumeBase* refVolume, const VolumeBase* secondVolume) const {
    // compute untransformed bounding box of reference volume
    tgt::vec3 refLLF = refVolume->getLLF();  // all components negative
    tgt::vec3 refURB = refVolume->getURB();  // all components positive
    std::pair<vec3, vec3> refBB = std::pair<vec3, vec3>(refVolume->getLLF(), refVolume->getURB());

    // create bounding box for second volume and transform into reference coordinate system
    MeshGeometry secondBB = MeshGeometry::createCube(secondVolume->getLLF(), secondVolume->getURB());
    tgt::mat4 transformToRef = refVolume->getWorldToPhysicalMatrix();
    transformToRef *= secondVolume->getPhysicalToWorldMatrix();
    secondBB.transform(transformToRef);

    // determine combined bounding box of reference volume and transformed second volume
    tgt::vec3 combinedLLF = refLLF;
    tgt::vec3 combinedURB = refURB;
    for (MeshGeometry::const_iterator face = secondBB.begin(); face != secondBB.end(); ++face) {
        for (FaceGeometry::const_iterator vertex = face->begin(); vertex != face->end(); ++vertex) {
            combinedLLF = tgt::min(combinedLLF, transformToRef * vertex->getCoords());
            combinedURB = tgt::max(combinedURB, transformToRef * vertex->getCoords());
        }
    }
    std::pair<vec3, vec3> combinedBB = std::pair<vec3, vec3>(combinedLLF, combinedURB);
    LDEBUG("Combined BB: " << combinedBB.first << ", " << combinedBB.second);
    tgtAssert(tgt::hand(tgt::lessThanEqual(combinedBB.first, refBB.first)) &&
        tgt::hand(tgt::greaterThanEqual(combinedBB.second, refBB.second)), "Invalid combined bounding box");

    // derive resolution of combined volume from size ratios of combinedBB/refBB and the ref volume's
    // original resolution => spatial resolution of combined volume equals spatial res. of reference volume

    // the bounding box contains an additional voxel which we have to take into account
    tgt::vec3 llfVoxel = refVolume->getPhysicalToVoxelMatrix() * combinedBB.first + tgt::vec3(0.5f);
    tgt::vec3 urbVoxel = refVolume->getPhysicalToVoxelMatrix() * combinedBB.second - tgt::vec3(0.5f);

    tgtAssert(tgt::hand(tgt::greaterThan(vec3(urbVoxel - llfVoxel), vec3(0.f))), "invalid bounding box size");
    tgt::svec3 combinedDim = svec3(tgt::iceil(urbVoxel - llfVoxel)) + tgt::svec3::one;
    LDEBUG("Common Dim: " << combinedDim);
    tgtAssert(tgt::hand(tgt::greaterThanEqual(combinedDim, refVolume->getDimensions())), "Invalid combined volume dimensions");

    // create combined volume with proper dimensions and spacing
    Volume* combinedVolume = 0;
    try {
        combinedVolume = VolumeOperatorResize::APPLY_OP(refVolume, combinedDim);
    }
    catch (const std::bad_alloc&) {
        LERROR("Failed to create combined volume with dimensions " << combinedDim << " : bad allocation");
        delete combinedVolume;
        return 0;
    }

    // determine correct offset and spacing for combined volume
    combinedVolume->setOffset(refVolume->getVoxelToPhysicalMatrix() * llfVoxel);
    vec3 s = combinedBB.second - combinedBB.first;
    s /= vec3(combinedDim);
    combinedVolume->setSpacing(s);
    combinedVolume->setPhysicalToWorldMatrix(refVolume->getPhysicalToWorldMatrix());

    return combinedVolume;
}

void VolumeCombine::adjustPropertyVisibilities() {
    factorC_.setVisibleFlag(combineFunction_.isSelected("weightedSum") || combineFunction_.isSelected("weightedSum2p"));
    factorD_.setVisibleFlag(combineFunction_.isSelected("weightedSum2p"));
}

} // namespace
